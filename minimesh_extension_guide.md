# Build and Compile
Latest version of CMake extension and some build configurations are probvided through the `CMakePresets.json` and the `launch.json` files.

There is a slower debug build and a faster optimized build. Since this is a development project, the debug build is the default. An automatic build and run with the cube.obj can be done using the CMake extension `Ctrl+Shift+B` from anywhere in the project. The default rebuild method is incremental and not a full clean rebuild.

New presets accessed through the bottom menu bar in VSCode can be added through the `CMakePresets.json` file.

### Command Line Interface
```bash
cd build-dbg
./bin/minimeshcli ../mesh/cube.obj
```

### Graphical User Interface
```bash
cd build-dbg
./bin/minimeshgui ../mesh/cube.obj
```

## Build Outputs

- **Executables**: `build-dbg/bin/` or `build-opt/bin/`
  - `minimeshcli` - Command line mesh processing
  - `minimeshgui` - Interactive GUI application

- **Libraries**:
  - `libminimeshcore.a` - Core mesh data structures
  - `libminimeshviz.a` - Visualization components

## File Structure

- `build-dbg/`, `build-opt/` - Build directories
- `src/minimesh/` - Source code
- `third-party/` - Dependencies
- `mesh/` - Example mesh files
- `.vscode/` - VSCode configuration

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

### 1. **Create New Mesh Modifier Classes**
The project uses specialized mesh modifier classes for different operations:

- `Mesh_modifier_loop_subdivision` (`mesh_modifier_loop_subdivision.hpp/cpp`) - Loop subdivision operations
- `Mesh_modifier_edge_collapse` (`mesh_modifier_edge_collapse.hpp/cpp`) - Edge collapse operations
- `Mesh_modifier_template` (`mesh_modifier_template.hpp/cpp`) - Template for new modifiers

To create a new modifier:
1. Copy `mesh_modifier_template.hpp/cpp` to a new name (e.g., `mesh_modifier_smoothing.hpp/cpp`)
2. Rename the class from `Mesh_modifier_template` to your new name
3. Add your new methods following the pattern of `flip_edge()`:

```cpp
// Example: Add vertex smoothing
bool smooth_vertex(int vertex_id);

// Example: Add face subdivision
std::vector<int> subdivide_face(int face_id);
```

4. Add the new files to `src/minimesh/core/CMakeLists.txt`

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
Follow the `flip_edge()` pattern in any of the mesh modifier implementation files (e.g., `mesh_modifier_template.cpp:46`):
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

The codebase follows a clear separation between connectivity (`mesh_connectivity`), modification (various `mesh_modifier_*` classes), and I/O (`mesh_io`), making it straightforward to extend with new geometric algorithms. Use `mesh_modifier_template` as a starting point for implementing new modifier classes.