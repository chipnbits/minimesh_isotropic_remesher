# Minimesh Half-Edge Data Structure Reference

This document provides a concise reference for the half-edge mesh data structure used in the minimesh library for mesh traversal and modification operations. Other details helpful for context are found in `minimesh_extension_guide.md`. The cli can be run directly with the meshcli command that has been sourced into terminal. The GUI likewise can be run with the meshgui command and filepath argument.

## Core Data Structures

### Vertex_data
```cpp
class Vertex_data {
    int half_edge;           // Index of outgoing half-edge
    Eigen::Vector3d xyz;     // 3D position
    bool is_active;          // Active/deleted status
};
```

### Half_edge_data
```cpp
struct Half_edge_data {
    int next;     // Next half-edge in face loop
    int prev;     // Previous half-edge in face loop
    int twin;     // Twin (opposite) half-edge
    int face;     // Face this half-edge belongs to
    int origin;   // Vertex this half-edge originates from
    bool is_active;
};
```

### Face_data
```cpp
struct Face_data {
    int half_edge;    // One half-edge on this face
    bool is_active;
};
```

## Iterator Classes

### Vertex_iterator
**Access Methods:**
- `data()` - Get reference to vertex data for modification
- `index()` - Get vertex index
- `is_active()` - Check if vertex is active
- `deactivate()` - Delete the vertex
- `half_edge()` - Get outgoing half-edge iterator
- `xyz()` - Get position (const reference)
- `is_equal(other)` - Compare vertices

### Half_edge_iterator
**Navigation Methods:**
- `next()` - Next half-edge around face
- `prev()` - Previous half-edge around face
- `twin()` - Opposite half-edge on same edge
- `origin()` - Vertex where half-edge starts
- `dest()` - Vertex where half-edge ends
- `face()` - Face this half-edge bounds

**Access Methods:**
- `data()` - Get reference to half-edge data
- `index()` - Get half-edge index
- `is_active()` - Check if active
- `deactivate()` - Delete half-edge
- `is_equal(other)` - Compare half-edges

### Face_iterator
**Access Methods:**
- `data()` - Get reference to face data
- `index()` - Get face index
- `half_edge()` - Get one half-edge on face
- `n_vertices()` - Count vertices/edges on face
- `is_active()` - Check if active
- `deactivate()` - Delete face
- `is_equal(other)` - Compare faces

### Vertex_ring_iterator
**Traversal Methods:**
- `half_edge()` - Current half-edge pointing TO the vertex
- `advance()` - Move to next half-edge, returns false when done
- `reset_boundary()` - Reset to boundary half-edge if vertex on boundary

## Mesh_connectivity Methods

### Element Access
```cpp
Vertex_iterator vertex_at(int id);
Half_edge_iterator half_edge_at(int id);
Face_iterator face_at(int id);
Vertex_ring_iterator vertex_ring_at(int vertex_id);
Face_iterator hole();  // Dummy face representing holes/boundaries
```

### Counts
```cpp
int n_active_vertices();
int n_active_half_edges();
int n_active_faces();
int n_total_vertices();    // Including deleted
int n_total_half_edges();  // Including deleted
int n_total_faces();       // Including deleted
```

### Element Creation
```cpp
Vertex_iterator add_vertex(bool allow_recycling = true);
Half_edge_iterator add_half_edge(bool allow_recycling = true);
Face_iterator add_face(bool allow_recycling = true);
```

### Mesh Building
```cpp
void build_from_triangles(const std::vector<double>& xyz,
                         const std::vector<int>& triangle_verts);
void build_from_polygons(const std::vector<double>& xyz,
                        const std::vector<int>& polygon_verts,
                        const std::vector<int>& polygon_adj);
```

### Utilities
```cpp
bool check_sanity_slowly(bool verbose = true);
void copy(const Mesh_connectivity& other);
void swap(Mesh_connectivity& other);
void clear();
```

## Mesh_modifier Methods

### Edge Operations
```cpp
int get_halfedge_between_vertices(int v0, int v1);  // Returns half-edge from v0 to v1
bool flip_edge(int he_index);                       // Flip edge (triangular meshes only)
bool divide_edge(int he_index, double weight = 0.5); // Subdivide edge with weighted vertex placement
bool subdivide_loop();                               // Loop subdivision (placeholder)
```

## Mesh Analysis Utilities

**Location**: `mesh_analysis.hpp` in `src/minimesh/core/mohe/`

Contains functional utilities for analyzing mesh properties and making assertions about mesh structure. These functions are useful for debugging, testing, and validation throughout the codebase. Functions include topology checks, connectivity analysis, geometric measurements, and validation routines.

## Common Traversal Patterns

### Face Loop Traversal
```cpp
Half_edge_iterator start = face.half_edge();
Half_edge_iterator he = start;
do {
    // Process half-edge he
    he = he.next();
} while (!he.is_equal(start));
```

### Vertex Ring Traversal
```cpp
Vertex_ring_iterator ring = mesh.vertex_ring_at(vertex_id);
do {
    Half_edge_iterator he = ring.half_edge();  // Points TO the vertex
    Vertex_iterator neighbor = he.origin();    // Neighboring vertex
    // Process neighbor
} while (ring.advance());
```

### Edge Walking
```cpp
Half_edge_iterator he = /* some half-edge */;
Vertex_iterator start_vertex = he.origin();
Vertex_iterator end_vertex = he.dest();
Half_edge_iterator opposite = he.twin();
Face_iterator left_face = he.face();
Face_iterator right_face = opposite.face();
```

## Constants
```cpp
static constexpr int hole_index = -2;      // Index for boundary/hole
static constexpr int invalid_index = -3;   // Invalid/non-existent index
```

## Key Design Principles

1. **Iterator Pattern**: All mesh elements accessed via iterators, not direct objects
2. **Half-Edge Structure**: Each edge split into two half-edges for efficient traversal
3. **Connectivity Only**: Pure topological structure - geometry stored separately in vertex positions
4. **Recycling**: Deleted elements can be reused rather than creating new storage
5. **Boundary Handling**: Boundary edges have twin pointing to special "hole" face