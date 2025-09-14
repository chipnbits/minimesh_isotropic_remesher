#pragma once

//
// mesh_analysis.hpp
//
// Functional utilities for analyzing mesh properties.
// These functions do not modify the mesh, only analyze its structure.
//

#include <minimesh/core/mohe/mesh_connectivity.hpp>

namespace minimesh
{
namespace mohecore
{
namespace analysis
{

// Check if mesh consists entirely of triangles
bool is_triangular_mesh(Mesh_connectivity& mesh);

// Check if mesh is manifold (each edge has exactly 2 adjacent faces, except boundaries)
bool is_manifold_mesh(const Mesh_connectivity& mesh);

// Check if mesh has boundary edges
bool has_boundary_edges(const Mesh_connectivity& mesh);

// Count boundary edges in the mesh
int count_boundary_edges(const Mesh_connectivity& mesh);

// Get vertex valence (number of edges connected to vertex)
int vertex_valence(Mesh_connectivity& mesh, int vertex_id);

// Check if vertex is on boundary
bool vertex_is_boundary(Mesh_connectivity& mesh, int vertex_id);

// Get average edge length in the mesh
double average_edge_length(const Mesh_connectivity& mesh);

// Get mesh bounding box diagonal length
double bounding_box_diagonal(const Mesh_connectivity& mesh);

// Check Euler characteristic (V - E + F = 2 - 2*genus for closed mesh)
bool check_euler_characteristic(const Mesh_connectivity& mesh, int expected_genus = 0);

// Count number of connected components in the mesh
int count_connected_components(Mesh_connectivity& mesh);

} // end of analysis
} // end of mohecore
} // end of minimesh