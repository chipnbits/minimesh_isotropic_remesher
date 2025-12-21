// Implementation of Uniform Explicit Isotropic Remeshing
// Based on Botsch and Kobbelt (2004)
//
// Project Part 1: Baseline remeshing algorithm
//

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <minimesh/core/mohe/mesh_analysis.hpp>
#include <minimesh/core/mohe/remesher/remesher_isotropic.hpp>
#include <minimesh/core/util/assert.hpp>
#include <set>
#include <vector>

namespace minimesh
{
namespace mohecore
{

// ============================================================
// Main Entry Point
// ============================================================

void
Mesh_modifier_uniform_remeshing::remesh(double target_edge_length, int iterations)
{
  // TODO: Implement the main remeshing loop
  // For each iteration:
  //   1. split_long_edges(target_edge_length)
  //   2. collapse_short_edges(target_edge_length)
  //   3. flip_edges_to_optimize_valence()
  //   4. tangential_smoothing()
}

// ============================================================
// Core Remeshing Operations
// ============================================================

void
Mesh_modifier_uniform_remeshing::split_long_edges(double target_length, double alpha)
{
  double t = target_length * alpha;
  const int initial_he_count = mesh().n_total_half_edges();

  // Make a scratch list for half-edges to mark visted ones
  std::vector<bool> visited(initial_he_count, false);

  auto he = Mesh_connectivity::Half_edge_iterator();
  int twin_idx;
  double length;
  for(int i = 0; i < initial_he_count; ++i)
  {
    he = mesh().half_edge_at(i);
    twin_idx = he.twin().index();
    if(visited[i] || visited[twin_idx])
    {
      visited[i] = true;
      visited[twin_idx] = true;
      continue;
    }
    visited[i] = true;
    visited[twin_idx] = true;

    if(!he.is_active())
      continue;
    length = get_edge_length(he.index());
    printf("Found edge %d of length %f, against threshold %f\n", he.index(), length, t);
    if(length > t)
    {
      split_edge(he.index(), 0.5); // Split at midpoint
    }
  }
}

void
Mesh_modifier_uniform_remeshing::collapse_short_edges(double target_length, double beta)
{
  // TODO: Implement edge collapsing
  // Recommended beta = 4.0/5.0 or 0.8
  // For each edge: if length < beta * target_length AND is_legal_collapse(), collapse
}

void
Mesh_modifier_uniform_remeshing::flip_edges_to_optimize_valence()
{
  // TODO: Implement edge flipping
  // For each edge: if flipping reduces total valence deviation, flip it
}

void
Mesh_modifier_uniform_remeshing::tangential_smoothing(int smoothing_iters)
{
  // TODO: Implement tangential smoothing
  // For each vertex:
  //   1. Compute barycenter of one-ring neighbors
  //   2. Move vertex towards barycenter
  //   3. Project back to tangent plane (or skip projection for simple version)
}

// ============================================================
// Low-Level Mesh Modifiers
// ============================================================

int
Mesh_modifier_uniform_remeshing::split_edge(int he_index, double weight)
{
  // Validate input index
  if(he_index < 0 || he_index >= mesh().n_total_half_edges())
    return mesh().invalid_index;

  // Get the half-edge and its twin
  auto he_1 = mesh().half_edge_at(he_index);
  if(!he_1.is_active())
    return mesh().invalid_index;

  auto tw_1 = he_1.twin();
  if(!tw_1.is_active())
    return mesh().invalid_index;

  // Get all involved vertices, half-edges, and faces
  auto he_2 = he_1.next();
  ;
  auto tw_2 = tw_1.next();
  auto he_3 = he_1.prev();
  ;
  auto tw_3 = tw_1.prev();
  auto F1 = he_1.face();
  ;
  auto F2 = tw_1.face();

  // Get the two vertices of the edge
  auto v1 = he_1.origin();
  auto v2 = tw_1.origin(); // Twin's origin is he's destination
  auto v3 = he_2.dest(); // Opposite vertex in face 1
  auto v4 = tw_2.dest(); // Opposite vertex in face 2

  bool he_face_is_hole = F1.is_equal(mesh().hole());
  bool tw_face_is_hole = F2.is_equal(mesh().hole());

  // Compute weighted position: weight * origin + (1-weight) * dest
  Eigen::Vector3d new_pos = weight * v1.xyz() + (1.0 - weight) * v2.xyz();

  // Create new vertex at weighted position
  auto v5 = mesh().add_vertex();
  v5.data().xyz = new_pos;
  tw_1.data().origin = v5.index(); // Update twin to point to new vertex

  // Create two new half-edges to split the original edge (dont recycle to avoid resplitting)
  auto he_4 = mesh().add_half_edge(false);
  auto tw_4 = mesh().add_half_edge(false);
  he_4.data().twin = tw_4.index();
  he_4.data().origin = v5.index();
  v5.data().half_edge = he_4.index();

  tw_4.data().twin = he_4.index();
  tw_4.data().origin = v2.index();
  v2.data().half_edge = tw_4.index();

  if(he_face_is_hole)
  {
    he_1.data().next = he_4.index();
    he_4.data().prev = he_1.index();
    he_4.data().next = he_2.index();
    he_2.data().prev = he_4.index();
    he_4.data().face = F1.index();
  }
  else
  {
    // Add new face, new twinned edge splitting F1
    auto F3 = mesh().add_face();
    auto he_5 = mesh().add_half_edge(false);
    auto he_6 = mesh().add_half_edge(false);
    he_5.data().twin = he_6.index();
    he_5.data().origin = v3.index();
    v3.data().half_edge = he_5.index();
    he_6.data().twin = he_5.index();
    he_6.data().origin = v5.index();
    v5.data().half_edge = he_6.index();

    // Set up new F1 cycle and face
    // Next cycle
    he_1.data().next = he_6.index();
    he_6.data().next = he_3.index();
    he_3.data().next = he_1.index(); // already set

    // Prev cycle
    he_3.data().prev = he_6.index();
    he_6.data().prev = he_1.index();
    he_1.data().prev = he_3.index(); // already set

    // Face assignments
    he_1.data().face = F1.index();
    he_6.data().face = F1.index();
    he_3.data().face = F1.index();
    F1.data().half_edge = he_1.index();

    // Set up new F3 cycle and face
    // Next cycle
    he_4.data().next = he_2.index();
    he_2.data().next = he_5.index();
    he_5.data().next = he_4.index();

    // Prev cycle
    he_5.data().prev = he_2.index();
    he_2.data().prev = he_4.index();
    he_4.data().prev = he_5.index();

    // Face assignments
    he_4.data().face = F3.index();
    he_2.data().face = F3.index();
    he_5.data().face = F3.index();
    F3.data().half_edge = he_4.index();

    printf("Face %d: half-edge %d\n", F1.index(), F1.data().half_edge);
    printf("Face %d: half-edge %d\n", F3.index(), F3.data().half_edge);
  }

  if(tw_face_is_hole)
  {
    tw_1.data().prev = tw_4.index();
    tw_4.data().next = tw_1.index();
    tw_4.data().prev = tw_3.index();
    tw_3.data().next = tw_4.index();
    tw_4.data().face = F2.index();
  }
  else
  {
    // Add new face, new twinned edge splitting F2
    auto F4 = mesh().add_face();
    auto tw_5 = mesh().add_half_edge(false);
    auto tw_6 = mesh().add_half_edge(false);
    tw_5.data().twin = tw_6.index();
    tw_5.data().origin = v5.index();
    v5.data().half_edge = tw_5.index();
    tw_6.data().twin = tw_5.index();
    tw_6.data().origin = v4.index();
    v4.data().half_edge = tw_6.index();

    // Set up new F2 cycle and face
    // Next cycle
    // tw_1.data().next = tw_2.index(); already set
    tw_2.data().next = tw_6.index();
    tw_6.data().next = tw_1.index();

    // Prev cycle
    tw_1.data().prev = tw_6.index();
    tw_6.data().prev = tw_2.index();
    // tw_2.data().prev = tw_1.index(); already set

    // Face assignments
    tw_1.data().face = F2.index();
    tw_2.data().face = F2.index();
    tw_6.data().face = F2.index();
    F2.data().half_edge = tw_1.index();

    // Set up new F4 cycle and face
    // Next cycle
    tw_4.data().next = tw_5.index();
    tw_5.data().next = tw_3.index();
    tw_3.data().next = tw_4.index();

    // Prev cycle
    tw_5.data().prev = tw_4.index();
    tw_4.data().prev = tw_3.index();
    tw_3.data().prev = tw_5.index();

    // Face assignments
    tw_4.data().face = F4.index();
    tw_3.data().face = F4.index();
    tw_5.data().face = F4.index();
    F4.data().half_edge = tw_4.index();
    printf("Face %d: half-edge %d\n", F2.index(), F2.data().half_edge);
    printf("Face %d: half-edge %d\n", F4.index(), F4.data().half_edge);
  }

  return v5.index();
}

bool
Mesh_modifier_uniform_remeshing::collapse_edge(int he_index)
{
  int v1 = mesh().half_edge_at(he_index).origin().index();
  int v2 = mesh().half_edge_at(he_index).dest().index();
  
  if(!is_legal_collapse(v1, v2))
    return false;

  auto v1_iter = mesh().vertex_at(v1);
  auto v2_iter = mesh().vertex_at(v2);
  auto he_v1_v2 = mesh().half_edge_at(he_index);
  auto he_v2_v1 = he_v1_v2.twin();
  auto u1_iter = he_v1_v2.next().dest();
  auto u2_iter = he_v2_v1.prev().origin();

  auto face_top = he_v1_v2.face();
  auto face_bottom = he_v2_v1.face();

  auto face_f1 = he_v1_v2.next().twin().face();
  auto face_f2 = he_v2_v1.prev().twin().face();
  auto he_f1_associate = he_v1_v2.prev();
  auto he_f2_associate = he_v2_v1.next();

  // Get the half-edges that will need to be re-linked
  auto he_f1_bottom = he_v1_v2.next().twin().next();
  auto he_f2_top = he_v2_v1.prev().twin().prev();

  // outer half-edges of two adjacent v2 side faces
  // defined so that s1=s2 when collapsing edge of degenerate form
  auto he_s1 = he_v1_v2.next().twin().prev();
  auto he_s2 = he_v2_v1.prev().twin().next();

  Eigen::Vector3d midpoint = 0.5 * (v1_iter.xyz() + v2_iter.xyz());

  // Remap all half-edges originating from v2 to originate from v1
  relabel_vertex(v2, v1);

  // Deactivate faces
  face_top.deactivate();
  face_bottom.deactivate();
  // Deactivate 6 half-edges
  he_v1_v2.next().twin().deactivate();
  he_v1_v2.next().deactivate();
  he_v1_v2.deactivate();
  he_v2_v1.prev().twin().deactivate();
  he_v2_v1.prev().deactivate();
  he_v2_v1.deactivate();
  // Deactivate vertex v2
  v2_iter.deactivate();

  // Check if this is a degenerate case where face_f1 == face_f2
  bool is_degenerate = he_s1.is_equal(he_s2);

  if(is_degenerate)
  {
    // f1 and f2 the same and s1, s2 the same
    he_f1_associate.data().prev = he_s1.index();
    he_s1.data().next = he_f1_associate.index();

    he_f1_associate.data().next = he_f2_associate.index();
    he_f2_associate.data().prev = he_f1_associate.index();

    he_f2_associate.data().next = he_s2.index();
    he_s2.data().prev = he_f2_associate.index();

    he_f1_associate.data().face = face_f1.index();
    he_f2_associate.data().face = face_f1.index();
    face_f1.data().half_edge = he_f1_associate.index();

    // // Run assertions to ensure connectivity is correct
    // // Check f1 = f2
    // assert(face_f1.is_equal(face_f2));
    // // Check s1 = s2
    // assert(he_s1.is_equal(he_s2));
    // // Check that cycle of 3 half-edges is correct and closed
    // const int he_start_index = he_f1_associate.index();
    // auto he_iter = he_f1_associate;
    // do{
    //   // Print out half edge twin, origin, prev, next, face for debugging
    //   std::cout << "Half-edge " << he_iter.index() << ": "
    //             << "twin=" << he_iter.twin().index() << ", "
    //             << "origin=" << he_iter.origin().index() << ", "
    //             << "prev=" << he_iter.prev().index() << ", "
    //             << "next=" << he_iter.next().index() << ", "
    //             << "face=" << he_iter.face().index() << std::endl;
    //   assert(he_iter.face().is_equal(face_f1));
    //   he_iter = he_iter.next();
    // }
    // while(he_iter.index() != he_start_index);
  }
  else
  {
    // Top edges re-linking
    he_f1_associate.data().next = he_f1_bottom.index();
    he_f1_bottom.data().prev = he_f1_associate.index();
    he_f1_associate.data().prev = he_f1_bottom.next().index();
    he_f1_bottom.next().data().next = he_f1_associate.index();
    // Relink face f1
    he_f1_associate.data().face = face_f1.index();
    face_f1.data().half_edge = he_f1_associate.index();

    // Bottom edges re-linking
    he_f2_associate.data().next = he_f2_top.prev().index();
    he_f2_top.prev().data().prev = he_f2_associate.index();
    he_f2_associate.data().prev = he_f2_top.index();
    he_f2_top.data().next = he_f2_associate.index();
    // Relink face f2
    he_f2_associate.data().face = face_f2.index();
    face_f2.data().half_edge = he_f2_associate.index();
  }

  // Update vertices to point to valid half-edges
  v1_iter.data().half_edge = he_f1_associate.twin().index();
  u1_iter.data().half_edge = he_f1_associate.index();
  u2_iter.data().half_edge = he_f2_associate.twin().index();

  // Move vertex v1 to optimal position
  mesh().vertex_at(v1).data().xyz = midpoint;
  return true;
}

bool
Mesh_modifier_uniform_remeshing::flip_edge(int he_index)
{
  // TODO: Implement edge flip
  // Only works for edges shared by two triangles
  // Return true if successful
  return false;
}

// ============================================================
// Helpers & Topology Checks
// ============================================================

double
Mesh_modifier_uniform_remeshing::get_edge_length(int he_index) const
{
  Mesh_connectivity::Half_edge_iterator he = _m.half_edge_at(he_index);
  Mesh_connectivity::Vertex_iterator v0 = he.origin();
  Mesh_connectivity::Vertex_iterator v1 = he.dest();

  Eigen::Vector3d p0 = v0.xyz();
  Eigen::Vector3d p1 = v1.xyz();

  return (p1 - p0).norm();
}

int
Mesh_modifier_uniform_remeshing::get_vertex_valence(int v_index) const
{
  return analysis::vertex_valence(_m, v_index);
}

int
Mesh_modifier_uniform_remeshing::get_valence_deviation(int v_index) const
{
  int valence = get_vertex_valence(v_index);
  bool is_boundary = analysis::vertex_is_boundary(_m, v_index);

  int ideal_valence = is_boundary ? 4 : 6;
  return std::abs(valence - ideal_valence);
}

bool
Mesh_modifier_uniform_remeshing::is_legal_collapse(int v1, int v2)
{

  // Check for boundary conditions
  bool v1_boundary = analysis::vertex_is_boundary(_m, v1);
  bool v2_boundary = analysis::vertex_is_boundary(_m, v2);
  auto he = mesh().half_edge_at(get_halfedge_between_vertices(v1, v2));
  auto twin = he.twin();

  bool he_face_is_hole = he.face().is_equal(mesh().hole());
  bool tw_face_is_hole = twin.face().is_equal(mesh().hole());

  bool edge_is_boundary = (he_face_is_hole || tw_face_is_hole);

  // Boundary edge rule, do not collapse interior edge connecting two boundary vertices
  if(v1_boundary && v2_boundary && !edge_is_boundary)
  {
    return false;
  }

  // Gather neighbors (excluding each other)
  std::set<int> n1 = get_all_neighbors_from_vertex(v1);
  std::set<int> n2 = get_all_neighbors_from_vertex(v2);
  n1.erase(v2);
  n2.erase(v1);

  // Opposite vertices across edge (the two triangle tips)
  std::set<int> allowed;
  // If the main face is NOT a hole, u1 exists
  if(!he_face_is_hole)
  {
    int u1 = he.next().dest().index();
    allowed.insert(u1);
  }

  // If the twin face is NOT a hole, u2 exists
  if(!tw_face_is_hole)
  {
    int u2 = twin.next().dest().index();
    allowed.insert(u2);
  }

  // Common neighbors
  std::set<int> common;
  std::set_intersection(n1.begin(), n1.end(), n2.begin(), n2.end(), std::inserter(common, common.begin()));

  // Fail Condition A: the only common neighbors are u1 and u2 (link condition)
  std::set<int> common_minus_allowed;
  std::set_difference(common.begin(),
      common.end(),
      allowed.begin(),
      allowed.end(),
      std::inserter(common_minus_allowed, common_minus_allowed.begin()));
  if(!common_minus_allowed.empty())
    return false;

  // Condition B: there exists at least one neighbor NOT in common (tetrahedron check)
  std::set<int> symdiff;
  std::set_symmetric_difference(n1.begin(), n1.end(), n2.begin(), n2.end(), std::inserter(symdiff, symdiff.begin()));
  if(symdiff.empty())
    return false;

  return true;
}

void
Mesh_modifier_uniform_remeshing::relabel_vertex(int old_id, int new_id)
{
  auto he = mesh().vertex_at(old_id).half_edge();
  const int he_start_index = he.index();
  do
  {
    assert(he.origin().index() == old_id);
    he.data().origin = new_id;
    he = he.twin().next();
  } while(he.index() != he_start_index);
}

std::set<int>
Mesh_modifier_uniform_remeshing::get_all_neighbors_from_vertex(int vertex_id) const
{
  std::set<int> neighbors;
  auto ring_iter = _m.vertex_ring_at(vertex_id);
  do
  {
    // he always points to vertex, so origin is the neighbor
    neighbors.insert(ring_iter.half_edge().origin().index());
  } while(ring_iter.advance());
  return neighbors;
}

int
Mesh_modifier_uniform_remeshing::get_halfedge_between_vertices(const int v0, const int v1)
{
  // Traverse the one-ring of v0 to find the half-edge pointing to v1
  Mesh_connectivity::Vertex_ring_iterator ring = _m.vertex_ring_at(v0);
  do
  {
    Mesh_connectivity::Half_edge_iterator he = ring.half_edge();
    if(he.origin().index() == v1)
    {
      // Found half-edge from v1 to v0, return its twin (v0 to v1)
      return he.twin().index();
    }
  } while(ring.advance());

  return -1; // No such half-edge found
}

std::vector<int>
Mesh_modifier_uniform_remeshing::get_one_ring_neighbors(int v_index) const
{
  std::vector<int> neighbors;
  Mesh_connectivity::Vertex_ring_iterator ring = _m.vertex_ring_at(v_index);

  do
  {
    Mesh_connectivity::Half_edge_iterator he = ring.half_edge();
    neighbors.push_back(he.origin().index());
  } while(ring.advance());

  return neighbors;
}

Eigen::Vector3d
Mesh_modifier_uniform_remeshing::compute_barycenter(int v_index) const
{
  Eigen::Vector3d barycenter(0.0, 0.0, 0.0);
  int count = 0;

  Mesh_connectivity::Vertex_ring_iterator ring = _m.vertex_ring_at(v_index);
  do
  {
    Mesh_connectivity::Half_edge_iterator he = ring.half_edge();
    Mesh_connectivity::Vertex_iterator neighbor = he.origin();
    barycenter += neighbor.xyz();
    count++;
  } while(ring.advance());

  if(count > 0)
  {
    barycenter /= static_cast<double>(count);
  }

  return barycenter;
}

} // end of mohecore
} // end of minimesh
