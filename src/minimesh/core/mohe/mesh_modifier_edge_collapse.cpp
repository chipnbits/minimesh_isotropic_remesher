#include <cmath>
#include <iostream> // For debug output
#include <minimesh/core/mohe/mesh_modifier_edge_collapse.hpp>
#include <minimesh/core/util/assert.hpp>
#include <set>

namespace minimesh
{
namespace mohecore
{

//
// SymQuadric implementations
//

// Default constructor - initializes to zero matrix
SymQuadric::SymQuadric()
{
  setZero();
}

// Reset to zero matrix
void
SymQuadric::setZero()
{
  q.fill(0.0);
}

SymQuadric &
SymQuadric::operator+=(const SymQuadric & o)
{
  // Element-wise addition of the 10 unique entries in the matrices
  for(int i = 0; i < 10; ++i)
    q[i] += o.q[i];
  return *this;
}

SymQuadric
operator+(SymQuadric a, const SymQuadric & b)
{
  // Copy a, then add b to it and return new a+b
  a += b;
  return a;
}

void
SymQuadric::addPlane(const Eigen::Vector4d & p)
{
  const double nx = p[0], ny = p[1], nz = p[2], d = p[3];
  q[0] += nx * nx; // a00
  q[1] += nx * ny; // a01
  q[2] += nx * nz; // a02
  q[3] += nx * d; // a03
  q[4] += ny * ny; // a11
  q[5] += ny * nz; // a12
  q[6] += ny * d; // a13
  q[7] += nz * nz; // a22
  q[8] += nz * d; // a23
  q[9] += d * d; // a33
}

double
SymQuadric::evalMul_xt_Q_x(const Eigen::Vector3d & x) const
{
  // Unpack blocks
  const double a00 = q[0], a01 = q[1], a02 = q[2];
  const double a11 = q[4], a12 = q[5], a22 = q[7];
  const double b0 = q[3], b1 = q[6], b2 = q[8];
  const double c = q[9];

  // x^T Q x
  const double xtQx = x[0] * (a00 * x[0] + 2.0 * a01 * x[1] + 2.0 * a02 * x[2]) +
                      x[1] * (a11 * x[1] + 2.0 * a12 * x[2]) + x[2] * (a22 * x[2]);

  // 2 b^T x + c
  return xtQx + 2.0 * (b0 * x[0] + b1 * x[1] + b2 * x[2]) + c;
}

bool
SymQuadric::solveMinimizer(Eigen::Vector3d & x) const
{
  // A
  //   [ a00  a01  a02  ]  = [ q[0], q[1], q[2] ]
  //   [ a01  a11  a12  ]  = [ q[1], q[4], q[5] ]
  //   [ a02  a12  a22  ]  = [ q[2], q[5], q[7] ]
  Eigen::Matrix3d A;
  A << q[0], q[1], q[2], q[1], q[4], q[5], q[2], q[5], q[7];
  // b
  //  [ b0 ]  = [ q[3] ]
  //  [ b1 ]  = [ q[6] ]
  //  [ b2 ]  = [ q[8] ]
  Eigen::Vector3d b(q[3], q[6], q[8]);

  // Robust 3x3 solve and check for invertibility with LU decomposition
  Eigen::FullPivLU<Eigen::Matrix3d> lu(A);
  if(!lu.isInvertible() || lu.rcond() < 1e-15) // Tolerance for numerical stability
    return false;

  x = lu.solve(-b);
  return x.allFinite();
}

Eigen::Matrix4d
SymQuadric::toMatrix() const
{
  Eigen::Matrix4d M;
  M << q[0], q[1], q[2], q[3], q[1], q[4], q[5], q[6], q[2], q[5], q[7], q[8], q[3], q[6], q[8], q[9];
  return M;
}


//
// Given two vertices, this function return the index of the half-edge going from v0 to v1.
// Returns -1 if no half-edge exists between the two vertices.
//
int
Mesh_modifier_edge_collapse::get_halfedge_between_vertices(const int v0, const int v1)
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

//
// Initialize the mesh simplifier to the current state of the mesh
//
void
Mesh_modifier_edge_collapse::initialize()
{
  _vertex_quadrics.resize(mesh().n_total_vertices());
  _pair_versions.clear();
  while(!_pair_heap.empty())
    _pair_heap.pop();

  // First compute quadrics for all vertices
  initialize_quadrics();

  // Then build the initial set of valid pairs based on the quadrics
  initialize_valid_pairs();
}

//
// Initialize quadric error matrices for all vertices
// Traverses all faces and accumulates plane contributions to vertex quadrics
//
void
Mesh_modifier_edge_collapse::initialize_quadrics()
{
  // Reset all quadrics to zero
  for(auto & Q : _vertex_quadrics)
  {
    Q.setZero();
  }

  int total_faces = mesh().n_total_faces();
  for(int f_id = 0; f_id < total_faces; ++f_id)
  {
    // Get face iterator
    Mesh_connectivity::Face_iterator f = mesh().face_at(f_id);

    // Skip inactive faces
    if(!f.is_active())
      continue;

    // Get half-edge of the face
    Mesh_connectivity::Half_edge_iterator he = f.half_edge();

    // Collect vertices of the face
    std::vector<Mesh_connectivity::Vertex_iterator> face_vertices;
    do
    {
      if(he.origin().is_active())
      {
        face_vertices.push_back(he.origin());
      }
      he = he.next();
    } while(he.index() != f.half_edge().index());

    // assert all faces triangular or throw an error
    if(face_vertices.size() != 3)
    {
      throw std::runtime_error(
          "Non-triangular face encountered in initialize_quadrics(). Only triangular meshes are supported.");
    }

    // Compute face normal using cross product
    Eigen::Vector3d v0 = face_vertices[0].xyz();
    Eigen::Vector3d v1 = face_vertices[1].xyz();
    Eigen::Vector3d v2 = face_vertices[2].xyz();
    Eigen::Vector3d normal = (v1 - v0).cross(v2 - v0);
    normal.normalize();

    // nx * x + ny * y + nz * z + d = 0  =>  d = - (nx*x0 + ny*y0 + nz*z0)
    double d = -normal.dot(v0);
    // Full Kp vector for planar face
    Eigen::Vector4d plane_vector(normal[0], normal[1], normal[2], d);

    // Add plane contribution to each vertex's quadric
    for(auto & v : face_vertices)
    {
      int v_index = v.index();
      _vertex_quadrics[v_index].addPlane(plane_vector);
    }
  }
}

void
Mesh_modifier_edge_collapse::initialize_valid_pairs()
{

  // Track visited half-edges to avoid duplicates
  int total_half_edges = mesh().n_total_half_edges();
  std::vector<bool> visited(total_half_edges, false);

  // Traverse all half-edges to find unique edges
  for(int he_id = 0; he_id < total_half_edges; ++he_id)
  {
    // Skip if already visited
    if(visited[he_id])
      continue;

    // Get the half-edge iterator
    Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_id);

    // Skip inactive half-edges
    if(!he.is_active())
      continue;

    // Mark this half-edge and its twin as visited
    visited[he_id] = true;
    visited[he.twin().index()] = true;

    // Get the two vertices of this edge
    int v1 = he.origin().index();
    int v2 = he.dest().index();

    // Add to data structures (creates canonical pair internally)
    add_or_update_pair(v1, v2);
  }
}

//
// Add or update a pair in the valid pairs structure
// Fetches Q1 and Q2, computes x_opt and error, then packages into a versioned heap entry
//
void
Mesh_modifier_edge_collapse::add_or_update_pair(int v1, int v2)
{
  // Create canonical pair (always v1 < v2)
  VertexPair pair(v1, v2);

  // Compute combined quadric error and the optimal position
  SymQuadric Q = _vertex_quadrics[v1] + _vertex_quadrics[v2];
  double error = Q.evalMul_xt_Q_x(Eigen::Vector3d::Zero());
  Eigen::Vector3d x_opt = Eigen::Vector3d::Zero();
  if(!Q.solveMinimizer(x_opt))
  {
    // Fallback to midpoint if minimizer is not valid
    Eigen::Vector3d p1 = mesh().vertex_at(v1).xyz();
    Eigen::Vector3d p2 = mesh().vertex_at(v2).xyz();
    x_opt = 0.5 * (p1 + p2);
  }
  error = Q.evalMul_xt_Q_x(x_opt);

  // Increment versioning for heap invalidation of old entries
  _pair_versions[pair]++;

  // Add to heap with new version
  MergeCandidate entry{error, pair, x_opt, _pair_versions[pair]};
  _pair_heap.push(entry);
}

// Get the next best pair to collapse (minimum error)
bool
Mesh_modifier_edge_collapse::get_min_pair(MergeCandidate & top)
{
  // Pop stale entries until we find a valid one
  while(!_pair_heap.empty())
  {
    top = _pair_heap.top();
    _pair_heap.pop();

    // Check if this entry is still valid
    if(_pair_versions[top.pair] == top.version)
    {
      // Valid entry found
      return true;
    }
    // Otherwise, this is a stale entry - continue to next
  }

  // Heap is empty or all entries are stale
  return false;
}


//  
bool
Mesh_modifier_edge_collapse::collapse_edge(MergeCandidate & candidate)
{
  int v1 = candidate.pair.v1;
  int v2 = candidate.pair.v2;

  // Check if collapse is legal (topology-preserving)
  if(!is_legal_collapse(v1, v2))
  {
    return false; // Illegal collapse
  }

  auto v1_iter = mesh().vertex_at(v1);
  auto v2_iter = mesh().vertex_at(v2);
  auto he_v1_v2 = mesh().half_edge_at(get_halfedge_between_vertices(v1, v2));
  auto he_v2_v1 = he_v1_v2.twin();
  auto u1_iter =  he_v1_v2.next().dest();
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


  // Invalidate all pairs involving v2 from the heap
  auto pairs_to_invalidate = get_all_pairs_from_vertex(v2);
  for(const auto & p : pairs_to_invalidate)
  {
    invalidate_pair(p);
  }

  // Remap all half-edges originating from v2 to originate from v1
  relabel_vertex(v2, v1);

  // Debug: print half-edge info before deactivation
  printf("  DEBUG: Before deactivation:\n");
  printf("    he_f1_associate: idx=%d, prev=%d, next=%d, active=%d\n",
         he_f1_associate.index(), he_f1_associate.prev().index(), he_f1_associate.next().index(), he_f1_associate.is_active());
  printf("    he_f2_associate: idx=%d, prev=%d, next=%d, active=%d\n",
         he_f2_associate.index(), he_f2_associate.prev().index(), he_f2_associate.next().index(), he_f2_associate.is_active());
  printf("    he_f1_bottom: idx=%d, prev=%d, next=%d, active=%d\n",
         he_f1_bottom.index(), he_f1_bottom.prev().index(), he_f1_bottom.next().index(), he_f1_bottom.is_active());
  printf("    he_f2_top: idx=%d, prev=%d, next=%d, active=%d\n",
         he_f2_top.index(), he_f2_top.prev().index(), he_f2_top.next().index(), he_f2_top.is_active());

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

  printf("  DEBUG: After deactivation, deactivated half-edges: %d, %d, %d, %d, %d, %d\n",
         he_v1_v2.next().twin().index(), he_v1_v2.next().index(), he_v1_v2.index(),
         he_v2_v1.prev().twin().index(), he_v2_v1.prev().index(), he_v2_v1.index());
  // Deactivate vertex v2
  v2_iter.deactivate();

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

  // Update vertices to point to valid half-edges
  v1_iter.data().half_edge = he_f1_associate.twin().index();
  u1_iter.data().half_edge = he_f1_associate.index();
  u2_iter.data().half_edge = he_f2_associate.twin().index();

  // Move vertex v1 to optimal position
  mesh().vertex_at(v1).data().xyz = candidate.x_opt;
  // Update quadric of v1
  _vertex_quadrics[v1] += _vertex_quadrics[v2];

  // Recompute pairs involving v1
  auto new_pairs = get_all_pairs_from_vertex(v1);
  for(const auto & p : new_pairs)
  {
    add_or_update_pair(p.v1, p.v2);
  }

  return true;
}

// Relabel the origin on all the half edges originating from v2
void Mesh_modifier_edge_collapse::relabel_vertex(int old_id, int new_id)
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

//
// Get all pairs containing a specific vertex
//
std::set<Mesh_modifier_edge_collapse::VertexPair>
Mesh_modifier_edge_collapse::get_all_pairs_from_vertex(int vertex_id)
{
  // Use a ring iterator to gather all pairs (one pair per neighbor / two half edges)
  std::set<VertexPair> pairs;
  Mesh_connectivity::Vertex_ring_iterator ring_iter = mesh().vertex_ring_at(vertex_id);
  Mesh_connectivity::Half_edge_iterator he;

  do
  {
    he = ring_iter.half_edge();
    pairs.insert(VertexPair(he.origin().index(), he.dest().index()));
  } while(ring_iter.advance());
  return pairs;
}

std::set<int>
Mesh_modifier_edge_collapse::get_all_neighbors_from_vertex(int vertex_id)
{
  std::set<int> neighbors;
  auto ring_iter = mesh().vertex_ring_at(vertex_id);
  do
  {
    // he always points to vertex, so origin is the neighbor
    neighbors.insert(ring_iter.half_edge().origin().index());
  } while (ring_iter.advance());
  return neighbors;
}

// Peek top-N valid candidates (without changing queue contents),
// while removing any stale entries encountered.
std::vector<int>
Mesh_modifier_edge_collapse::get_top_n_candidates(int n)
{
    std::vector<MergeCandidate> result;
    result.reserve(n);

    // Temporarily hold valid entries we popped so we can restore them.
    std::vector<MergeCandidate> restore;
    restore.reserve(n);

    // Pop until we gather N valid entries or the heap is exhausted.
    while ((int)result.size() < n && !_pair_heap.empty())
    {
        MergeCandidate e = _pair_heap.top();
        _pair_heap.pop();

        // Stale? (version mismatch) -> drop it permanently (heap cleanup)
        if (_pair_versions[e.pair] != e.version)
            continue;

        // Optionally guard against deactivated vertices.
        // (Not strictly necessary if your versioning/invalidations are correct.)
        // auto v1_active = mesh().vertex_at(e.pair.v1).is_active();
        // auto v2_active = mesh().vertex_at(e.pair.v2).is_active();
        // if (!v1_active || !v2_active) continue;

        result.push_back(e);
        restore.push_back(e);
    }

    // Restore the valid entries so this truly behaves like "peek".
    for (const auto& e : restore)
        _pair_heap.push(e);

    // convert to half-edge indices    
    std::vector<int> half_edge_indices;
    half_edge_indices.reserve(result.size());
    for (const auto& e : result)
    {
        int he_index = get_halfedge_between_vertices(e.pair.v1, e.pair.v2);
        if (he_index != mesh().invalid_index)
            half_edge_indices.push_back(he_index);
    }

    return half_edge_indices;
}


void
Mesh_modifier_edge_collapse::invalidate_pair(VertexPair pair)
{
  VertexPair pair(v1, v2);
  _pair_versions[pair]++;
}

// Condition checks that edge (v1,v2) has exactly two shared neighbors on the two faces
// that it connects, and at least one unique neighbor on either side (tetrahedron catcher)
bool Mesh_modifier_edge_collapse::is_legal_collapse(int v1, int v2)
{
  // Gather neighbors (excluding each other)
  std::set<int> n1 = get_all_neighbors_from_vertex(v1);
  std::set<int> n2 = get_all_neighbors_from_vertex(v2);
  n1.erase(v2);
  n2.erase(v1);

  // Opposite vertices across edge (the two triangle tips)
  auto he = mesh().half_edge_at(get_halfedge_between_vertices(v1, v2));
  int u1 = he.next().dest().index();
  int u2 = he.twin().next().dest().index();

  // Common neighbors
  std::set<int> common;
  std::set_intersection(n1.begin(), n1.end(),
                        n2.begin(), n2.end(),
                        std::inserter(common, common.begin()));

  // Condition A: the only common neighbors allowed are u1 and u2
  std::set<int> common_minus_u;
  std::set<int> allowed = {u1, u2};
  std::set_difference(common.begin(), common.end(),
                      allowed.begin(), allowed.end(),
                      std::inserter(common_minus_u, common_minus_u.begin()));
  if (!common_minus_u.empty())
    return false;

  // Condition B: there exists at least one neighbor NOT in common
  // (i.e., the symmetric difference is non-empty)
  std::set<int> symdiff;
  std::set_symmetric_difference(n1.begin(), n1.end(),
                                n2.begin(), n2.end(),
                                std::inserter(symdiff, symdiff.begin()));
  if (symdiff.empty())
    return false;

  return true;
}


} // end of mohecore
} // end of minimesh
