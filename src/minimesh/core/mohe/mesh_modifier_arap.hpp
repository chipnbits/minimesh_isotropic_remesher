#pragma once

//
// mesh_modifier_edge_collapse.hpp
//
// Edge collapse mesh modifier.
//
// Author: Shayan Hoshyari
//


#include <string>
#include <unordered_set>
#include <vector>

#include <Eigen/Core>

#include <minimesh/core/mohe/mesh_connectivity.hpp>

namespace minimesh
{
namespace mohecore
{


class Mesh_modifier_arap
{
public:

  struct Neighbor {
    int index;      // Neighbor vertex index
    double weight;  // Cotangent weight w_ij
  };

  // Trivial constructor
  Mesh_modifier_arap(Mesh_connectivity & mesh_in)
  : _m(mesh_in)
  {
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
  // ARAP Deformation Methods (stubs for now)
  //

  // Initialize ARAP structures (call this after mesh is loaded)
  void initialize();

  // Build adjacency list with cotangent weights for all vertices
  // This computes w_ij = (cot(alpha) + cot(beta)) / 2 for each edge
  void build_cotangent_weights();

  // Add a vertex as an anchor point
  // Returns true if successful, false if vertex is already an anchor
  bool add_anchor(const int vertex_index);

  // Remove a vertex from anchors
  // Returns true if successful, false if vertex wasn't an anchor
  bool remove_anchor(const int vertex_index);

  // Check if a vertex is an anchor
  bool is_anchor(const int vertex_index) const;

  // Clear all anchors
  void clear_anchors();

  // Get all anchor vertex indices (returns vector for convenience)
  std::vector<int> get_anchors();

  // Perform ARAP deformation by moving an anchor (in-place version for GUI performance)
  // vertex_index: the anchor to move
  // new_position: the new position for the anchor
  // deformed_positions: output matrix (3 x n_vertices) with deformed positions
  // Returns true if deformation was successful
  // NOTE: This is a stub - actual ARAP implementation to be added later
  //       For now, just copies original positions and updates the anchor
  bool deform_with_temp_anchor(const int vertex_index,
      const Eigen::Vector3d & new_position,
      Eigen::Matrix3Xd & deformed_positions);

  // Compute ARAP deformation and return new positions (convenience version for testing)
  // vertex_index: the anchor to move
  // new_position: the new position for the anchor
  // Returns: Matrix (3 x n_vertices) with deformed positions
  // Throws: std::runtime_error if deformation fails
  Eigen::Matrix3Xd compute_deformation(const int vertex_index,
      const Eigen::Vector3d & new_position);

  // Apply ARAP deformation directly to mesh vertex positions (for CLI & file output)
  // vertex_index: the anchor to move
  // new_position: the new position for the anchor
  // Returns true if deformation was successful
  // NOTE: This modifies the actual mesh vertex positions in-place
  bool apply_deformation_to_mesh(const int vertex_index,
      const Eigen::Vector3d & new_position);


private:
  // pointer to the mesh that we are working on.
  Mesh_connectivity & _m;

  // Boolean array tracking which vertices are anchors (O(1) lookup)
  std::vector<bool> _is_anchor;

  // Adjacency list with cotangent weights
  // adj[i] holds (neighbor_index, cotangent_weight) for each neighbor of vertex i
  std::vector<std::vector<Neighbor>> _adj;

  Eigen::MatrixXd _vertices_rest;    // Rest positions of vertices

  // Core ARAP solver - all public deformation methods delegate to this
  // vertex_index: the temporary anchor to move
  // new_position: the new position for the temporary anchor
  // output: output matrix (3 x n_vertices) with deformed positions
  // Returns true if deformation was successful
  bool _solve_arap(const int vertex_index,
      const Eigen::Vector3d & new_position,
      Eigen::Matrix3Xd & output);
};


} // end of mohecore
} // end of minimesh
