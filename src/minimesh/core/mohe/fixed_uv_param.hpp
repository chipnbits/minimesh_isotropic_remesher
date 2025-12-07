#pragma once

//
// fixed_uv_param.hpp
//
// Fixed boundary UV parameterization.
//
// Author: Shayan Hoshyari
//


#include <string>
#include <Eigen/Dense>
#include <minimesh/core/mohe/mesh_connectivity.hpp>

namespace minimesh
{
namespace mohecore
{


class Fixed_boundary_uv_param
{
public:
  // Trivial constructor
  Fixed_boundary_uv_param(Mesh_connectivity & mesh_in)
  : _m(mesh_in), _is_computed(false)
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

  const std::vector<int> & get_boundary_loop() const { return _boundary_loop_indices; }

  //
  // Compute the UV parameterization
  // This populates the internal UV storage
  // Returns true if successful, false otherwise
  //
  bool compute_parameterization();

  //
  // Access the computed UV coordinates
  // Returns a reference to the UV vector (indexed by vertex index)
  // Each entry is a 2D point (u, v)
  //
  const std::vector<Eigen::Vector2d> & get_uv_coords() const { return _uv_coords; }

  //
  // Get UV coordinate for a specific vertex
  // Returns the (u, v) coordinate for the given vertex index
  //
  Eigen::Vector2d get_uv_at_vertex(int vertex_index) const;

  //
  // Check if parameterization has been computed
  //
  bool is_computed() const { return _is_computed; }


private:
  // pointer to the mesh that we are working on.
  Mesh_connectivity & _m;

  // Stores the boundary loop vertex indices
  std::vector<int> _boundary_loop_indices;

  // Stores interior vertex indices (non-boundary)
  std::vector<int> _interior_vertex_indices;

  // Stores the computed UV coordinates for each vertex
  // Index corresponds to vertex index in the mesh
  std::vector<Eigen::Vector2d> _uv_coords;

  // Flag to track if parameterization has been computed
  bool _is_computed;

  //
  // Helper: Initialize boundary and interior vertex lists
  //
  void initialize_vertex_lists();

  //
  // Helper: Map boundary vertices to a circle or square in UV space
  //
  void map_boundary_to_uv();

  //
  // Helper: Solve for interior vertex UV positions
  // (This is where you'd implement harmonic/Tutte/LSCM/etc.)
  //
  bool solve_interior_uvs();
};


} // end of mohecore
} // end of minimesh
