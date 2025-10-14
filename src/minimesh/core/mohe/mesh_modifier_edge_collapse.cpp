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
// Get all pairs containing a specific vertex
//
std::set<Mesh_modifier_edge_collapse::VertexPair>
Mesh_modifier_edge_collapse::get_all_pairs_from_vertex(int vertex_id)
{
  // Use a ring iterator to find all incident half-edges
  std::set<VertexPair> pairs;
  Mesh_connectivity::Vertex_ring_iterator ring_iter = mesh().vertex_ring_at(vertex_id);
  do
  {
    auto he = ring_iter.half_edge();
    int v1 = he.origin().index();
    int v2 = he.dest().index();
    pairs.insert(VertexPair(v1, v2));
  } while(ring_iter.advance());
  return pairs;
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


//
// Testing helper: invalidate a pair by incrementing its version WITHOUT adding to heap
//
void
Mesh_modifier_edge_collapse::invalidate_pair(int v1, int v2)
{
  VertexPair pair(v1, v2);
  _pair_versions[pair]++;
}


} // end of mohecore
} // end of minimesh
