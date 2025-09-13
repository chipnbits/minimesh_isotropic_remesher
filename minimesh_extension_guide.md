# Minimesh CLI and Extension Guide

## CLI Structure

The CLI in `src/minimesh/cli/main.cpp:75` is a simple executable that demonstrates basic mesh operations. It:

1. **Creates core objects**: `Mesh_connectivity`, `Mesh_io`, and `Mesh_modifier`
2. **Reads meshes**: Uses `mesh_io.read_obj_general()` for OBJ files
3. **Performs operations**: Currently demonstrates edge flipping with `mesh_modifier.flip_edge()`
4. **Writes output**: Saves to both VTK and OBJ formats

## Core Data Structures

**Mesh_connectivity** (`src/minimesh/core/mohe/mesh_connectivity.hpp:43`) - The main half-edge data structure:
- **Vertices** (`Vertex_data:57`): Store position (`xyz`) and outgoing half-edge
- **Half-edges** (`Half_edge_data:86`): Store `next`, `prev`, `twin`, `face`, and `origin`
- **Faces** (`Face_data:124`): Store reference to a half-edge on the face

**Iterators** for traversal:
- `Vertex_iterator` - Access vertex data and outgoing half-edges
- `Half_edge_iterator` - Navigate next/prev/twin/face relationships
- `Face_iterator` - Access face data and boundary half-edges
- `Vertex_ring_iterator` - Traverse all half-edges around a vertex

## Extending for Transforms and Traversals

### 1. **Create New Algorithms in `mesh_modifier.hpp`**
Add methods to the `Mesh_modifier` class following the pattern of `flip_edge()`:

```cpp
// Example: Add vertex smoothing
bool smooth_vertex(int vertex_id);

// Example: Add face subdivision
std::vector<int> subdivide_face(int face_id);
```

### 2. **Use Iterator Patterns for Traversals**
```cpp
// Traverse all faces
for(int i = 0; i < mesh.n_total_faces(); ++i) {
    auto face = mesh.face_at(i);
    if(face.is_active()) {
        // Process face
    }
}

// Ring traversal around vertex
auto ring = mesh.vertex_ring_at(vertex_id);
do {
    auto he = ring.half_edge();
    // Process half-edge pointing to vertex
} while(ring.advance());
```

### 3. **Mesh Transforms Implementation Pattern**
Follow the `flip_edge()` pattern in `mesh_modifier.cpp:46`:
- Get iterators to relevant entities
- Check validity/boundary conditions
- Modify connectivity via `.data()` access
- Update next/prev/twin/face/origin relationships

### 4. **Add to CLI**
Extend `main.cpp` with your new operations:
```cpp
// Add after line 114
modi.your_new_transform(parameters);
check_sanity_and_write_mesh();
```

### 5. **File I/O Extensions**
Use `Mesh_io` (`mesh_io.hpp:20`) for reading/writing:
- `read_obj_general()`, `read_off()` for input
- `write_obj()`, `write_vtk()` for output
- VTK format supports custom vertex/face data for visualization

## Key Implementation Notes

The half-edge structure provides O(1) navigation between adjacent mesh elements, making it ideal for implementing local mesh operations, smoothing algorithms, subdivision schemes, and topological modifications.

### Half-Edge Data Structure Benefits
- **Efficient traversal**: Navigate around vertices, faces, and edges in constant time
- **Manifold support**: Designed specifically for manifold meshes (with or without boundary)
- **Memory management**: Automatic recycling of deleted entities via internal stacks
- **Robustness**: Built-in sanity checking with `check_sanity_slowly()`

### Common Traversal Patterns
1. **Around a vertex**: Use `Vertex_ring_iterator` to visit all incident half-edges
2. **Around a face**: Use `Half_edge_iterator.next()` to walk face boundary
3. **Across mesh**: Iterate through all active entities using total count and `is_active()` checks
4. **Neighborhood queries**: Use `get_halfedge_between_vertices()` for connectivity queries

The codebase follows a clear separation between connectivity (`mesh_connectivity`), modification (`mesh_modifier`), and I/O (`mesh_io`), making it straightforward to extend with new geometric algorithms.