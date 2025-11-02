#include <algorithm>
#include <minimesh/core/mohe/mesh_modifier_arap.hpp>
#include <minimesh/core/util/assert.hpp>

namespace minimesh
{
namespace mohecore
{


//
// Given two vertices, this function return the index of the half-edge going from v0 to v1.
// Returns -1 if no half-edge exists between the two vertices.
//
int
Mesh_modifier_arap::get_halfedge_between_vertices(const int v0, const int v1)
{
  // Get a ring iterator for v0
  Mesh_connectivity::Vertex_ring_iterator ring_iter = mesh().vertex_ring_at(v0);

  int answer = mesh().invalid_index;

  // Loop over all half-edges that end at v0.
  do
  {
    // Make sure that the half-edge does end at v0
    assert(ring_iter.half_edge().dest().index() == v0);

    // If the half-edge also starts and v1, then it's twin
    // goes from v0 to v1. This would be the half-edge that
    // we were looking for
    if(ring_iter.half_edge().origin().index() == v1)
    {
      answer = ring_iter.half_edge().twin().index();
    }
  } while(ring_iter.advance());

  if(answer != mesh().invalid_index)
  {
    assert(mesh().half_edge_at(answer).origin().index() == v0);
    assert(mesh().half_edge_at(answer).dest().index() == v1);
  }

  return answer;
}


bool
Mesh_modifier_arap::flip_edge(const int he_index)
{
  //
  // Take a reference to all involved entities
  //

  // HALF-EDGES
  Mesh_connectivity::Half_edge_iterator he0 = mesh().half_edge_at(he_index);
  Mesh_connectivity::Half_edge_iterator he1 = he0.twin();

  // meshes on the boundary are not flippable
  if(he0.face().is_equal(mesh().hole()) || he1.face().is_equal(mesh().hole()))
  {
    return false;
  }

  Mesh_connectivity::Half_edge_iterator he2 = he0.next();
  Mesh_connectivity::Half_edge_iterator he3 = he2.next();
  Mesh_connectivity::Half_edge_iterator he4 = he1.next();
  Mesh_connectivity::Half_edge_iterator he5 = he4.next();

  // VERTICES
  Mesh_connectivity::Vertex_iterator v0 = he1.origin();
  Mesh_connectivity::Vertex_iterator v1 = he0.origin();
  Mesh_connectivity::Vertex_iterator v2 = he3.origin();
  Mesh_connectivity::Vertex_iterator v3 = he5.origin();

  // FACES
  Mesh_connectivity::Face_iterator f0 = he0.face();
  Mesh_connectivity::Face_iterator f1 = he1.face();

  //
  // Now modify the connectivity
  //

  // HALF-EDGES
  he0.data().next = he3.index();
  he0.data().prev = he4.index();
  he0.data().origin = v3.index();
  //
  he1.data().next = he5.index();
  he1.data().prev = he2.index();
  he1.data().origin = v2.index();
  //
  he2.data().next = he1.index();
  he2.data().prev = he5.index();
  he2.data().face = f1.index();
  //
  he3.data().next = he4.index();
  he3.data().prev = he0.index();
  //
  he4.data().next = he0.index();
  he4.data().prev = he3.index();
  he4.data().face = f0.index();
  //
  he5.data().next = he2.index();
  he5.data().prev = he1.index();

  // VERTICES
  v0.data().half_edge = he2.index();
  v1.data().half_edge = he4.index();
  v2.data().half_edge = he1.index();
  v3.data().half_edge = he0.index();

  // FACES
  f0.data().half_edge = he0.index();
  f1.data().half_edge = he1.index();

  // operation successful
  return true;
} // All done


//
// ARAP Deformation Methods (stubs for now)
//

void
Mesh_modifier_arap::initialize()
{
  // Initialize the boolean array to track anchors
  _is_anchor.clear();
  _is_anchor.resize(mesh().n_total_vertices(), false);
}


bool
Mesh_modifier_arap::add_anchor(const int vertex_index)
{
  // Check if vertex exists and is active
  if(vertex_index < 0 || vertex_index >= static_cast<int>(_is_anchor.size()))
  {
    return false;
  }

  Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(vertex_index);
  if(!v.is_active())
  {
    return false;
  }

  // Check if already an anchor
  if(_is_anchor[vertex_index])
  {
    return false; // Already an anchor
  }

  // Mark as anchor
  _is_anchor[vertex_index] = true;
  return true;
}


bool
Mesh_modifier_arap::remove_anchor(const int vertex_index)
{
  // Check bounds
  if(vertex_index < 0 || vertex_index >= static_cast<int>(_is_anchor.size()))
  {
    return false;
  }

  // Check if it was an anchor
  if(!_is_anchor[vertex_index])
  {
    return false; // Not an anchor
  }

  // Remove anchor
  _is_anchor[vertex_index] = false;
  return true;
}


bool
Mesh_modifier_arap::is_anchor(const int vertex_index) const
{
  // Check bounds
  if(vertex_index < 0 || vertex_index >= static_cast<int>(_is_anchor.size()))
  {
    return false;
  }

  return _is_anchor[vertex_index];
}


void
Mesh_modifier_arap::clear_anchors()
{
  // Reset all anchors to false
  std::fill(_is_anchor.begin(), _is_anchor.end(), false);
}


std::vector<int>
Mesh_modifier_arap::get_anchors()
{
  std::vector<int> anchors;

  // Collect all indices where _is_anchor is true
  for(int i = 0; i < static_cast<int>(_is_anchor.size()); ++i)
  {
    if(_is_anchor[i])
    {
      // Also verify the vertex is still active
      Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(i);
      if(v.is_active())
      {
        anchors.push_back(i);
      }
    }
  }

  return anchors;
}

/* Core ARAP Deformation Solver (Private)
Uses the set anchors in original static position with one temporary anchor moved to a new position.
Assumes output matrix is already properly sized and initialized with positions.
All public deformation methods delegate to this.
*/
bool
Mesh_modifier_arap::_solve_arap(const int temp_anchor_index,
    const Eigen::Vector3d & pulled_position,
    Eigen::Matrix3Xd & output)
{
  // Check if vertex is active
  Mesh_connectivity::Vertex_iterator v_pulled = mesh().vertex_at(temp_anchor_index);
  if(!v_pulled.is_active())
  {
    return false;
  }

  printf("Deformed mesh with anchor vertex %d \n", temp_anchor_index);
  printf(
      " - New anchor position: (%f, %f, %f) \n", pulled_position.x(), pulled_position.y(), pulled_position.z());

  // TODO: Implement ARAP deformation algorithm
  // For now, this is a stub that doesn't modify positions (naive stub)
  // The actual ARAP solver will:
  // 1. Use the prefactorized system matrix (from initialize())
  // 2. Update constraint for moved anchor
  // 3. Solve for all vertex positions to minimize ARAP energy
  // 4. Store results in output matrix

  return true;
}


/* ARAP Deformation - In-place version (for GUI performance)
Uses the set anchors in original static position with one temporary anchor moved to a new position.
Modifies the passed set of deformed positions in-place (no allocation).
*/
bool
Mesh_modifier_arap::deform_with_temp_anchor(const int vertex_index,
    const Eigen::Vector3d & new_position,
    Eigen::Matrix3Xd & deformed_positions)
{
  // Direct pass-through to solver (GUI provides pre-allocated matrix)
  return _solve_arap(vertex_index, new_position, deformed_positions);
}


/* ARAP Deformation - Returning version (for testing & convenience)
Allocates a new matrix, initializes with current mesh positions, computes deformation, and returns it.
*/
Eigen::Matrix3Xd
Mesh_modifier_arap::compute_deformation(const int vertex_index,
    const Eigen::Vector3d & new_position)
{
  // Allocate new matrix
  Eigen::Matrix3Xd result(3, mesh().n_total_vertices());

  // Initialize with current mesh positions
  for(int i = 0; i < mesh().n_total_vertices(); ++i)
  {
    Mesh_connectivity::Vertex_iterator vi = mesh().vertex_at(i);
    if(vi.is_active())
    {
      result.col(i) = vi.xyz();
    }
  }

  // Solve
  bool success = _solve_arap(vertex_index, new_position, result);

  if(!success)
  {
    throw std::runtime_error("ARAP deformation failed for vertex " + std::to_string(vertex_index));
  }

  return result;
}


/* ARAP Deformation - Mesh-modifying version (for CLI & file output)
Computes deformation in a temporary matrix and applies directly to mesh vertex positions.
*/
bool
Mesh_modifier_arap::apply_deformation_to_mesh(const int vertex_index,
    const Eigen::Vector3d & new_position)
{

  Eigen::Matrix3Xd deformed_positions = compute_deformation(vertex_index, new_position);

  // Apply deformed positions back to mesh
  for(int i = 0; i < mesh().n_total_vertices(); ++i)
  {
    Mesh_connectivity::Vertex_iterator vi = mesh().vertex_at(i);
    if(vi.is_active())
    {
      vi.data().xyz = deformed_positions.col(i);
    }
  }

  return true;
}


} // end of mohecore
} // end of minimesh
