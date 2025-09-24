#pragma once

//
// mesh_modifier.hpp
//
// Some functionality for modifying meshes.
//
// Author: Shayan Hoshyari
//


#include <string>

#include <minimesh/core/mohe/mesh_connectivity.hpp>

namespace minimesh
{
namespace mohecore
{


class Mesh_modifier
{
public:
  // Trivial constructor
  Mesh_modifier(Mesh_connectivity & mesh_in)
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


  // Loop subdivision (placeholder implementation)
  bool subdivide_loop();

  // Subdivide an edge by adding a vertex at weighted position between endpoints
  // Returns the index of the newly created vertex, or invalid_index if operation failed
  int divide_edge(const int he_index, const double weight = 0.5);

  // Subdivide a into smaller faces using a set of ordered midpoint vertices that form a ring
  // Creates one new center face and n new faces around it
  // Returns true if operation successful, and false if operation failed
  bool subdivide_face(Mesh_connectivity::Face_iterator & face, const std::vector<int> & edge_midpoint_vertices);


private:
  // pointer to the mesh that we are working on.
  Mesh_connectivity & _m;

  // Helper functions for loop subdivision
  void loop_update_old_vertices(Mesh_connectivity & original_mesh,
      const std::vector<int> & original_active_vertex_ids);


};


} // end of mohecore
} // end of minimesh
