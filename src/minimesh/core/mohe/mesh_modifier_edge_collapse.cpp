#include <cmath>
#include <minimesh/core/mohe/mesh_modifier_edge_collapse.hpp>
#include <minimesh/core/util/assert.hpp>
#include <iostream> // For debug output

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
  if (!lu.isInvertible() || lu.rcond() < 1e-15) // Tolerance for numerical stability
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


bool
Mesh_modifier_edge_collapse::flip_edge(const int he_index)
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

  // TODO: Traverse all faces in the mesh
  // For each face:
  //   1. Compute face normal and plane equation [nx, ny, nz, d]
  //   2. Get all vertices of the face
  //   3. Add the plane contribution to each vertex's quadric
  //      using Q.addPlane(plane_vector, weight)

  // Stub implementation - will be filled in with face traversal
  int total_faces = mesh().n_total_faces();
  (void)total_faces; // Suppress unused warning for now
}


//
// Compute edge collapse metric for a given half-edge
// Uses the square root of the origin vertex ID as the metric
//
float
Mesh_modifier_edge_collapse::compute_edge_metric(Mesh_connectivity::Half_edge_iterator he)
{
  // Use square root of origin vertex index as metric
  int vertex_id = he.origin().index();
  return std::sqrt(static_cast<float>(vertex_id));
}


//
// Initialize the priority queue with all edges in the mesh
//
void
Mesh_modifier_edge_collapse::initialize_priority_queue()
{
  // Clear the priority queue
  while(!_edge_pq.empty())
  {
    _edge_pq.pop();
  }

  // Create a boolean vector to track visited half-edges
  int total_half_edges = mesh().n_total_half_edges();
  std::vector<bool> visited(total_half_edges, false);

  // Traverse all half-edges in the mesh
  for(int he_id = 0; he_id < total_half_edges; ++he_id)
  {
    // Skip if already visited
    if(visited[he_id])
    {
      continue;
    }

    // Get the half-edge iterator
    Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_id);

    // Skip inactive half-edges
    if(!he.is_active())
    {
      continue;
    }

    // Mark this half-edge and its twin as visited
    visited[he_id] = true;
    visited[he.twin().index()] = true;

    // Compute the metric for this edge
    float metric = compute_edge_metric(he);

    // Add to priority queue
    _edge_pq.push(std::make_pair(metric, he_id));
  }

  // Verify that PQ has exactly half the number of half-edges (one entry per full edge)
  int expected_edges = mesh().n_active_half_edges() / 2;
  int actual_edges = static_cast<int>(_edge_pq.size());
  assert(expected_edges == actual_edges && "Priority queue should contain exactly half the number of half-edges");
}


//
// Get the top N candidates from the priority queue
// Returns a vector of half-edge indices
//
std::vector<int>
Mesh_modifier_edge_collapse::get_top_n_candidates(int n)
{
  std::vector<int> candidates;

  // Make a copy of the priority queue so we don't modify the original
  auto pq_copy = _edge_pq;

  // Extract top N elements
  for(int i = 0; i < n && !pq_copy.empty(); ++i)
  {
    auto top = pq_copy.top();
    candidates.push_back(top.second); // second is the half-edge index
    pq_copy.pop();
  }

  return candidates;
}


} // end of mohecore
} // end of minimesh
