#pragma once

//
// mesh_modifier_edge_collapse.hpp
//
// Edge collapse mesh modifier.
//
// Author: Shayan Hoshyari
//


#include <array>
#include <limits>
#include <map>
#include <queue>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <minimesh/core/mohe/mesh_connectivity.hpp>

namespace minimesh
{
namespace mohecore
{

//
// Symmetric 4x4 quadric matrix for Quadric Error Metric (QEM)
// Stores only 10 unique entries due to symmetry
//
struct SymQuadric
{
  /*
  10 unique entries of symmetric 4x4 matrix in row-major "triangle" order
  q = [a00, a01, a02, a03, a11, a12, a13, a22, a23, a33]

  Represents the full matrix:
  [ a00  a01  a02  a03 ]
  [ a01  a11  a12  a13 ]
  [ a02  a12  a22  a23 ]
  [ a03  a13  a23  a33 ]
  */

  std::array<double, 10> q;

  SymQuadric();
  void setZero();

  // Add another quadric (Q += other)
  SymQuadric & operator+=(const SymQuadric & o);
  friend SymQuadric operator+(SymQuadric a, const SymQuadric & b);

  // Add a plane contribution outer product, where p = [nx, ny, nz, d] with unit normal vector
  void addPlane(const Eigen::Vector4d & p);

  // Evaluate the QEM energy at x (3D point, homogeneous [x y z 1])
  // cost(x) = x^T A x + 2 b^T x + c, with A 3x3 SPD-ish, b 3x1, c scalar
  double evalMul_xt_Q_x(const Eigen::Vector3d & x) const;

  // Solve A x = -b for the minimizer of the quadratic; returns true if A is PD enough.
  // Falls back to 'false' if singular/ill-conditioned; caller should try endpoints/midpoint.
  bool solveMinimizer(Eigen::Vector3d & x) const;

  // Optional: materialize full 4x4 (rarely needed at runtime)
  Eigen::Matrix4d toMatrix() const;
};


class Mesh_modifier_edge_collapse
{
public:
  // Constructor - allocates quadrics for all vertices
  Mesh_modifier_edge_collapse(Mesh_connectivity & mesh_in)
  : _m(mesh_in)
  {
    // Allocate quadrics for all vertices (including inactive ones for consistent indexing)
    _vertex_quadrics.resize(mesh_in.n_total_vertices());
  }

  // Get the underlying mesh
  Mesh_connectivity & mesh() { return _m; }
  const Mesh_connectivity & mesh() const { return _m; }

  //
  // Given two vertices, this function return the index of the half-edge going from v0 to v1.
  // Returns mesh::invalid_index if no half-edge exists between the two vertices.
  //
  int get_halfedge_between_vertices(const int v0, const int v1);

  //
  // Flip an edge in a mesh
  // Input: The mesh, and the index of a half-edge on the edge we wish to flip
  // Return true if operation successful, and false if operation not possible
  //
  // Assumption: mesh is all triangles
  //
  // NOTE: To see how this method works, take a look at edge-flip.svg
  //
  bool flip_edge(const int he_index);

  //
  // Initialize quadric error matrices for all vertices
  // Traverses all faces and accumulates plane contributions to vertex quadrics
  //
  void initialize_quadrics();

  // Retrieve the quadric for a given vertex
  SymQuadric & vertex_quadric(int v_id) { return _vertex_quadrics[v_id]; }
  const SymQuadric & vertex_quadric(int v_id) const { return _vertex_quadrics[v_id]; }

  // Initialize valid pairs for edge collapse
  void initialize_valid_pairs();

  //
  // Initialize the priority queue with edge metrics for all edges in the mesh
  //
  void initialize_priority_queue();

  //
  // Get the top N edge candidates from the priority queue
  // Returns a vector of half-edge indices
  //
  std::vector<int> get_top_n_candidates(int n);

private:
  // Priority queue entry: (metric, half_edge_index)
  // We use greater<> for min-heap (smallest metric has highest priority)
  using PQEntry = std::pair<float, int>;
  struct ComparePQEntry
  {
    bool operator()(const PQEntry & a, const PQEntry & b) const
    {
      return a.first > b.first; // Min-heap: smaller metric = higher priority
    }
  };

  // pointer to the mesh that we are working on.
  Mesh_connectivity & _m;

  // Priority queue for edge collapse candidates (min-heap based on metric)
  std::priority_queue<PQEntry, std::vector<PQEntry>, ComparePQEntry> _edge_pq;

  // Quadric error matrices for each vertex (indexed by vertex ID)
  // One quadric per vertex, stores sum of planes from incident faces
  std::vector<SymQuadric> _vertex_quadrics;

  //
  // Compute the collapse metric for a given edge
  // Returns a float value where lower values indicate better collapse candidates
  //
  float compute_edge_metric(Mesh_connectivity::Half_edge_iterator he);
};


} // end of mohecore
} // end of minimesh
