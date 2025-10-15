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
#include <set>
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
  // Vertex pair representation (public so it can be used in API)
  // Always stores v1 < v2 for canonical ordering
  struct VertexPair
  {
    int v1, v2;

    VertexPair(int a, int b)
    : v1(std::min(a, b))
    , v2(std::max(a, b))
    {
    }

    bool operator<(const VertexPair & other) const { return (v1 < other.v1) || (v1 == other.v1 && v2 < other.v2); }

    bool operator==(const VertexPair & other) const { return v1 == other.v1 && v2 == other.v2; }
  };

  // Heap entry for the priority queue
  struct MergeCandidate
  {
    double error;
    VertexPair pair;
    Eigen::Vector3d x_opt;
    int version; // For lazy deletion

    // Default constructor
    MergeCandidate()
    : error(0.0)
    , pair(0, 0)
    , x_opt(Eigen::Vector3d::Zero())
    , version(0)
    {
    }

    // Parameterized constructor
    MergeCandidate(double e, VertexPair p, Eigen::Vector3d x, int v)
    : error(e)
    , pair(p)
    , x_opt(x)
    , version(v)
    {
    }

    bool operator>(const MergeCandidate & other) const
    {
      return error > other.error; // For min-heap
    }
  };

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
  // Initialize the mesh simplifier to the current state of the mesh
  // This computes quadrics for all vertices and builds the initial set of valid pairs
  // Call this after construction or whenever you need to reset the simplifier
  //
  void initialize();

  // Retrieve the quadric for a given vertex
  SymQuadric & vertex_quadric(int v_id) { return _vertex_quadrics[v_id]; }
  const SymQuadric & vertex_quadric(int v_id) const { return _vertex_quadrics[v_id]; }

  // Add or update a pair in the valid pairs structure
  // If pair exists, updates its error and adds new entry to heap
  void add_or_update_pair(int v1, int v2);

  // Get the next best pair to collapse (minimum error)
  // Returns true if a valid pair exists, false if heap is empty, pops the entry from the heap
  // Automatically skips stale heap entries
  bool get_min_pair(MergeCandidate & top);

  // Collapse an edge between two vertices, merging v2 into v1
  // Updates connectivity and removes degenerate faces
  // Returns true if collapse was successful, false if illegal (topology-preserving)
  bool collapse_edge(MergeCandidate & candidate);

  // Invalidate a pair in the heap by incrementing the latest pair version WITHOUT adding to heap
  // This makes all heap entries for this pair stale
  void invalidate_pair(int v1, int v2);

  // Check if collapsing the edge between v1 and v2 is legal (topology-preserving)
  bool is_legal_collapse(int v1, int v2);

  // Get all pairs containing a specific vertex (returns copies for safe iteration)
  std::set<VertexPair> get_all_pairs_from_vertex(int vertex_id);

  // Get all neighboring vertices of a vertex
  std::set<int> get_all_neighbors_from_vertex(int vertex_id);

  //
  // Get the top N edge candidates from the priority queue
  // Returns a vector of half-edge indices
  //
  std::vector<int> get_top_n_candidates(int n);



private:
  // pointer to the mesh that we are working on.
  Mesh_connectivity & _m;

  // Quadric error matrices for each vertex (indexed by vertex ID)
  // One quadric per vertex, stores sum of planes from incident faces
  std::vector<SymQuadric> _vertex_quadrics;

  //
  // Initialize quadric error matrices for all vertices
  // Traverses all faces and accumulates plane contributions to vertex quadrics
  //
  void initialize_quadrics();

  // Initialize valid pairs for edge collapse
  void initialize_valid_pairs();

  //
  // Valid pairs data structure
  //
  // Min-heap ordered by error metric (supports lazy deletion via versions)
  std::priority_queue<MergeCandidate, std::vector<MergeCandidate>, std::greater<MergeCandidate>> _pair_heap;

  // Version tracking for lazy deletion in heap
  std::map<VertexPair, int> _pair_versions;
};


} // end of mohecore
} // end of minimesh
