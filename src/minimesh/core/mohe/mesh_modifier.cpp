#include <cmath>
#include <minimesh/core/mohe/mesh_analysis.hpp>
#include <minimesh/core/mohe/mesh_modifier.hpp>
#include <minimesh/core/util/assert.hpp>
#include <queue>
#include <unordered_map>

namespace minimesh
{
namespace mohecore
{

namespace
{

struct LoopEvenCoeffs
{
  double center; // weight for the center vertex
  double neighbor; // weight for each 1-ring neighbor
};

// Compute and cache coefficients for a given valence n.
// Using formula: beta = (40 - (3 + 2*cos(2π/n))^2) / (64*n)
// then center = 1 - n*beta; neighbor = beta.
static const LoopEvenCoeffs &
coeffs_for_valence(int n)
{
  // Trivial fallback (also covers n<3)
  if(n < 3)
  {
    static const LoopEvenCoeffs kFallback{1.0, 0.0};
    return kFallback;
  }

  static std::vector<LoopEvenCoeffs> cache(4, {1.0, 0.0}); // 0..3 prefilled
  if((int)cache.size() <= n)
    cache.resize(n + 1, {0.0, 0.0});
  auto & c = cache[n];

  // If uninitialized, compute
  if(c.center == 0.0 && c.neighbor == 0.0)
  {
    const double ang = 2.0 * M_PI / double(n);
    const double t = 3.0 + 2.0 * std::cos(ang);
    const double beta = (40.0 - t * t) / (64.0 * double(n)); // no pow
    c.center = 1.0 - n * beta;
    c.neighbor = beta;
  }
  return c;
}

} // namespace


//
// Given two vertices, this function return the index of the half-edge going from v0 to v1.
// Returns -1 if no half-edge exists between the two vertices.
//
int
Mesh_modifier::get_halfedge_between_vertices(const int v0, const int v1)
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
Mesh_modifier::flip_edge(const int he_index)
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

int
Mesh_modifier::divide_edge(const int he_index, const double weight)
{
  // Validate input index
  if(he_index < 0 || he_index >= mesh().n_total_half_edges())
    return mesh().invalid_index;

  // Get the half-edge and its twin
  auto he = mesh().half_edge_at(he_index);
  if(!he.is_active())
    return mesh().invalid_index;

  auto twin = he.twin();
  if(!twin.is_active())
    return mesh().invalid_index;

  // Get the two vertices of the edge
  auto origin_vertex = he.origin();
  auto dest_vertex = twin.origin(); // Twin's origin is he's destination

  // Compute weighted position: weight * origin + (1-weight) * dest
  Eigen::Vector3d new_pos = weight * origin_vertex.xyz() + (1.0 - weight) * dest_vertex.xyz();

  // Create new vertex at weighted position
  auto new_vertex = mesh().add_vertex();
  new_vertex.data().xyz = new_pos;

  // Create two new half-edges to split the original edge
  auto new_he = mesh().add_half_edge();
  auto new_twin = mesh().add_half_edge();

  // Store original connectivity info before modification
  int he_next = he.data().next;
  int he_face = he.data().face;

  int twin_next = twin.data().next;
  int twin_face = twin.data().face;

  // Update original half-edge: now goes from origin to new vertex
  he.data().next = new_he.index();
  // prev stays the same
  // face stays the same
  // origin stays the same
  he.data().twin = new_twin.index();

  // Update new half-edge: goes from new vertex to destination
  new_he.data().next = he_next;
  new_he.data().prev = he.index();
  new_he.data().twin = twin.index();
  new_he.data().face = he_face;
  new_he.data().origin = new_vertex.index();

  // Update original twin: now goes from dest to new vertex
  twin.data().next = new_twin.index();
  // prev stays the same
  // face stays the same
  // origin stays the same (dest_vertex)
  twin.data().twin = new_he.index();

  // Update new twin: goes from new vertex to origin
  new_twin.data().next = twin_next;
  new_twin.data().prev = twin.index();
  new_twin.data().twin = he.index();
  new_twin.data().face = twin_face;
  new_twin.data().origin = new_vertex.index();

  // Update next/prev pointers of adjacent half-edges
  auto new_he_next_idx = new_he.data().next;
  if(new_he_next_idx != mesh().invalid_index)
  {
    new_he.next().data().prev = new_he.index();
  }
  auto new_twin_next_idx = new_twin.data().next;
  if(new_twin_next_idx != mesh().invalid_index)
  {
    mesh().half_edge_at(new_twin_next_idx).data().prev = new_twin.index();
  }

  // Update vertex half-edge pointers
  new_vertex.data().half_edge = new_he.index();

  // Origin and dest vertices retained the same half-edge pointers
  // Face retained the same half-edge pointers

  return new_vertex.index();
}

bool
Mesh_modifier::subdivide_loop()
{
  // 1. Validate triangular mesh
  if(!mohecore::analysis::is_triangular_mesh(mesh()))
    return false;

  // 2. Create a copy of the original mesh connectivity for reference
  Mesh_connectivity original_mesh;
  original_mesh.copy(mesh());

  // 3. Snapshot original faces and vertices, track visited edges and new vertices
  std::vector<int> original_active_face_ids;
  for(int f = 0; f < original_mesh.n_total_faces(); ++f)
  {
    auto face = original_mesh.face_at(f);
    if(face.is_active())
    {
      original_active_face_ids.push_back(f);
    }
  }

  // Store original active vertices for deformation step later
  std::vector<int> original_active_vertex_ids;
  for(int v = 0; v < original_mesh.n_total_vertices(); ++v)
  {
    auto vertex = original_mesh.vertex_at(v);
    if(vertex.is_active())
    {
      original_active_vertex_ids.push_back(v);
    }
  }

  // Create mapping from half-edge index to midpoint vertex index, invalid index means not visited
  std::vector<int> edge_to_midpoint(original_mesh.n_total_half_edges(), original_mesh.invalid_index);

  // Iterate over all faces in the original mesh
  for(int f_id : original_active_face_ids)
  {
    auto face = original_mesh.face_at(f_id);
    std::vector<int> edge_midpoint_vertices;
    edge_midpoint_vertices.reserve(3); // triangles

    // Iterate over the half-edges of the face to split edges and collect midpoints
    auto he_start = face.half_edge();
    auto he = he_start;
    do
    {
      int he_idx = he.index();
      int twin_idx = he.twin().index();

      // Check if edge has already been split
      int new_midpoint = edge_to_midpoint[he_idx];
      if(new_midpoint == original_mesh.invalid_index)
      {
        // Edge has not been split, create a new midpoint
        new_midpoint = divide_edge(he_idx, 0.5);
        if(new_midpoint == mesh().invalid_index)
        {
          return false; // Failed to divide edge
        }
        edge_to_midpoint[he_idx] = new_midpoint;
        edge_to_midpoint[twin_idx] = new_midpoint; // Twin edge shares the same midpoint
        // Update position of new vertex using Loop rules
        auto split_edge = original_mesh.half_edge_at(he_idx);
        auto twin_edge = split_edge.twin();
        const bool is_boundary_edge =
            split_edge.face().is_equal(original_mesh.hole()) || twin_edge.face().is_equal(original_mesh.hole());

        Eigen::Vector3d new_pos;
        if(is_boundary_edge)
        {
          // Boundary edge rule
          new_pos = 0.5 * (split_edge.origin().xyz() + split_edge.dest().xyz());
        }
        else
        {
          // Interior new vertex rule
          new_pos = 0.375 * (split_edge.origin().xyz() + split_edge.dest().xyz()) +
                    0.125 * (split_edge.next().dest().xyz() + twin_edge.next().dest().xyz());
        }
        mesh().vertex_at(new_midpoint).data().xyz = new_pos;
      }
      edge_midpoint_vertices.push_back(new_midpoint);
      he = he.next();
    } while(!he.is_equal(he_start));

    // Now create the four new faces from the original face and midpoints
    Mesh_connectivity::Face_iterator new_face = mesh().face_at(f_id); // get the face from the current mesh
    subdivide_face(new_face, edge_midpoint_vertices);
  }
  // Move the old vertices
  loop_update_old_vertices(original_mesh, original_active_vertex_ids);

  return true;
}

// At this stage triangle has divided edges, midpoint list available
// Create four new triangles and deactivate the original one
//
//
//           (4)
//           / \\
//          /   \\
//         /     \\
//       (5)------(3)
//       / \\central/\\
//      /   \\    /   \\
//     /     \\  /     \\
//   (0)-----(1)------(2)
//
bool
Mesh_modifier::subdivide_face(Mesh_connectivity::Face_iterator & face, const std::vector<int> & edge_midpoint_vertices)
{

  for(int v : edge_midpoint_vertices)
  {
    if(v < 0 || v >= mesh().n_total_vertices())
      return false;
  }

  // Get onto the new face with half edge
  auto he_start = face.half_edge();

  // Fetch first half edge that starts at vertex 0,
  auto he0 = he_start;
  while(he0.origin().index() != edge_midpoint_vertices[0])
  {
    he0 = he0.next();
    if(he0.is_equal(he_start))
      return false; // midpoint not found on face
  }


  // Create outer faces & the central face. We'll collect the central edges.
  auto central_face = mesh().add_face();
  std::vector<int> central_edges; // indices of central half-edges we create
  central_edges.reserve(edge_midpoint_vertices.size());

  auto edge_start_12 = he0; // outer triangle seed edge
  do
  {
    // the next outer triangle always 2 nexts away
    auto next_start = edge_start_12.next().next();
    auto outer_face = mesh().add_face();
    auto new_he_31 = mesh().add_half_edge(); // 3 -> 1 part of outer face

    // Close the outer triangle cycle: (1->2)->(2->3)->(3->1)
    auto next_edge_23 = edge_start_12.next();
    edge_start_12.data().prev = new_he_31.index();
    next_edge_23.data().next = new_he_31.index();
    new_he_31.data().next = edge_start_12.index();
    new_he_31.data().prev = next_edge_23.index();

    // Attach to outer face
    new_he_31.data().face = outer_face.index();
    edge_start_12.data().face = outer_face.index();
    next_edge_23.data().face = outer_face.index();
    outer_face.data().half_edge = new_he_31.index();

    // Twin the edge with the central face edge
    auto new_twin_13 = mesh().add_half_edge(); // 1 -> 3 part of central face
    new_he_31.data().twin = new_twin_13.index();
    new_twin_13.data().twin = new_he_31.index();
    new_twin_13.data().face = central_face.index();

    // Place origins correctly
    new_he_31.data().origin = next_edge_23.dest().index();
    new_twin_13.data().origin = edge_start_12.origin().index();

    // Collect central edges for cycle
    central_edges.push_back(new_twin_13.index());

    // Advance to next outer triangle
    edge_start_12 = next_start;
  } while(!(edge_start_12.index() == he0.index()));

  // We must have a ring of central edges now
  if(central_edges.size() != edge_midpoint_vertices.size())
    return false;

  // Set the face's reference edge ONCE
  central_face.data().half_edge = central_edges[0];

  // Close the central face cycle
  for(int i = 0; i < 3; ++i)
  {
    int he_idx = central_edges[i];
    int next_idx = central_edges[(i + 1) % 3];
    int prev_idx = central_edges[(i + 2) % 3];
    auto he = mesh().half_edge_at(he_idx);
    he.data().next = next_idx;
    he.data().prev = prev_idx;
    he.data().face = central_face.index();
  }

  // Deactivate the original face
  face.deactivate();

  return true;
}


void
Mesh_modifier::loop_update_old_vertices(Mesh_connectivity & original_mesh,
    const std::vector<int> & original_active_vertex_ids)
{
  // Find the umbrella of neighbors using the ring iterator
  for(int v_id : original_active_vertex_ids)
  {
    // Get ring iterator
    Mesh_connectivity::Vertex_ring_iterator v_ring = original_mesh.vertex_ring_at(v_id);

    // Check special boundary condition branch
    if(v_ring.reset_boundary())
    {
      int n1 = v_ring.half_edge().origin().index();
      v_ring.advance();
      int n2 = v_ring.half_edge().origin().index();

      const Eigen::Vector3d v = original_mesh.vertex_at(v_id).xyz();
      const Eigen::Vector3d nsum = original_mesh.vertex_at(n1).xyz() + original_mesh.vertex_at(n2).xyz();

      mesh().vertex_at(v_id).data().xyz = 0.75 * v + 0.125 * nsum;
      continue;
    }

    // Not a boundary vertex, proceed as normal
    // Grab all neighbors
    std::vector<int> neighbor_ids;
    neighbor_ids.reserve(12); // preallocate space
    do
    {
      // Record the origins since half-edges point to v_id
      neighbor_ids.push_back(v_ring.half_edge().origin().index());
    } while(v_ring.advance());

    // Compute beta
    int n = static_cast<int>(neighbor_ids.size());
    const LoopEvenCoeffs & coeffs = coeffs_for_valence(n);


    auto orig_v = original_mesh.vertex_at(v_id);
    // Compute weighted sum of neighbor positions
    Eigen::Vector3d neighbor_sum = Eigen::Vector3d::Zero();
    for(int nid : neighbor_ids)
    {
      neighbor_sum += original_mesh.vertex_at(nid).xyz();
    }
    // Apply Loop subdivision weights: new_pos = center_weight * old_pos + neighbor_weight * neighbor_sum
    Eigen::Vector3d new_pos = coeffs.center * orig_v.xyz() + coeffs.neighbor * neighbor_sum;
    // Update vertex position in current mesh
    mesh().vertex_at(v_id).data().xyz = new_pos;
  }
}

} // end of mohecore
} // end of minimesh
