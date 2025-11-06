# ARAP Mesh Deformation

**Author:** Simon Ghyselincks
**Student Number:** 12145124

## Overview

This project implements As-Rigid-As-Possible (ARAP) withs both GUI and CLI interfaces for mesh deformation.

## Running the Applications

The CLI and GUI can be run using the following aliased commands:

- **CLI**: `meshcli <filepath>`
  - Example: `meshcli ./mesh/camel_simple.obj`

- **GUI**: `meshgui <filepath>`
  - Example: `meshgui ./mesh/camel_simple.obj`

## GUI Usage

### Activating Deformation Mode

To begin deforming a mesh in the GUI:

1. **Load a mesh** via command line argument
2. **Enable Deform Mode** by clicking the `Deform` checkbox in the interface
3. **Select Static Anchors** using "Select Vertex" mode
   - Static anchors remain fixed during deformation
   - Multiple static anchors can be selected to constrain the deformation
   - Click on vertices in the viewport to add them as static anchors
   - Selected anchors are highlighted in the visualization as blue points
4. **Deselect Anchors** using the same "Select Vertex" mode on already-selected vertices
5. **Move the Dynamic Anchor** using "Move Vertex" mode
   - Only one dynamic anchor can be active at a time
   - The mesh deforms in real-time as you drag the dynamic anchor
   - ARAP optimization runs iteratively to minimize deformation energy

### Resetting Deformation

To reset all anchors and deformation state:
- **Deselect the `Deform` checkbox** - This exits ARAP deformation mode and clears all anchors
- **Reselect the checkbox** to start a fresh deformation session with the original mesh

## CLI Usage

### Basic Setup

The ARAP modifier is integrated into the minimesh framework and can be used programmatically in CLI applications:

```cpp
#include <minimesh/core/mohe/mesh_connectivity.hpp>
#include <minimesh/core/mohe/mesh_io.hpp>
#include <minimesh/core/mohe/mesh_modifier_arap.hpp>

// Load mesh
mohecore::Mesh_connectivity mesh;
mohecore::Mesh_io io(mesh);
io.read_obj("path/to/mesh.obj");

// Create and initialize ARAP modifier
mohecore::Mesh_modifier_arap modi_arap(mesh);
modi_arap.initialize();
```

### Managing Anchors

The ARAP modifier provides comprehensive anchor management through the following API:

```cpp
// Add static anchor points (vertices that remain fixed)
printf("Adding anchor points...\n");
for(int i = 0; i < 48 && i < mesh.n_active_vertices(); ++i)
{
  if(modi_arap.add_anchor(i))
  {
    auto v = mesh.vertex_at(i);
    Eigen::Vector3d pos = v.xyz();
    printf("  - Anchor %d at position (%.3f, %.3f, %.3f)\n",
           i, pos.x(), pos.y(), pos.z());
  }
}

// Check if a vertex is an anchor
if(modi_arap.is_anchor(vertex_id)) {
  printf("Vertex %d is anchored\n", vertex_id);
}

// Remove a specific anchor
modi_arap.remove_anchor(vertex_id);

// Get all static anchors
std::vector<int> anchors = modi_arap.get_static_anchors();
printf("Total anchors: %zu\n", anchors.size());

// Clear all anchors
modi_arap.clear_anchors();
```

### Applying Deformation

To deform the mesh, specify a vertex to move (the dynamic anchor) and its new position. The ARAP solver iteratively optimizes vertex positions and local rotations to achieve minimal deformation energy:

```cpp
// Select a vertex to move (dynamic anchor)
int test_vertex = 100;
Eigen::Vector3d pulled_pos(1.0, 2.0, 0.5);

// Apply deformation - iterates to 1% convergence
// Updates mesh vertex positions in-place
if(modi_arap.apply_deformation_to_mesh(test_vertex, pulled_pos))
{
  printf("Deformation applied successfully\n");

  // Save deformed mesh
  io.write_obj("deformed_mesh.obj");
}
else
{
  printf("Deformation failed\n");
}
```

### Alternative Deformation Methods

For more advanced use cases, the modifier provides additional deformation methods:

```cpp
// Compute deformation without modifying mesh (returns position matrix)
Eigen::Matrix3Xd deformed_positions =
  modi_arap.compute_deformation(vertex_index, new_position);

// In-place deformation with output buffer (for performance-critical GUI updates)
Eigen::Matrix3Xd output_buffer;
bool success = modi_arap.deform_with_temp_anchor(
  vertex_index, new_position, output_buffer);
```

## Implementation Details

### Core Algorithm

The ARAP implementation follows the algorithm from the original paper:

1. **Preprocessing**: Compute cotangents for all half-edges and then build edge weights and Laplacian matrices
2. **Constraint Modification**: Update Laplacian blocks to account for static anchor constraints. New laplacian system is formed as:
   ```
   [ Lff  Lfc ] [ Vf ] = [ Bf ]
   [  0    I  ] [ Vc ]   [ Vc ]
   ```
   where `Vf` are free vertices and `Vc` are constrained (anchored) vertices. The Lff and Lfc are fixed for a given set of anchors.
3. **Position Optimization**: The solver constructs the Bf term using the current rotations and constrained vertex positions.
4. **Rotation Update**: Each vertex's rotation matrix is updated to match best fit of its neighboring edges.

Efficiency: The GUI based implementation only restimates position and rotation once per frame, since deformations tend to be small and incremental. This provides better interactivity while using the previous solution as a warm start.

The CLI implementation converges to 1% change in energy using the equation from the ARAP paper. This may be over 40+ iterations for large deformations.

### Data Structures

- **Adjacency List**: Stores neighbor relationships with precomputed cotangent weights
- **Laplacian Blocks**: Sparse matrices `Lff` (free-free) and `Lfc` (free-constrained)
- **Rotation Matrices**: Per-vertex 3x3 rotation matrices updated each iteration
- **Anchor Management**: Efficient tracking of static anchors with versioning for incremental updates
- **Deformed Positions**: Internal storage for current deformed vertex positions to be reused for GUI interactivity
## API Reference

### Mesh_modifier_arap

**Initialization:**
- `Mesh_modifier_arap(Mesh_connectivity& mesh)` - Constructor
- `void initialize()` - Build ARAP data structures (cotangent weights, Laplacian)

**Anchor Management:**
- `bool add_anchor(int vertex_index)` - Add static anchor, returns false if already anchored
- `bool remove_anchor(int vertex_index)` - Remove static anchor
- `bool is_anchor(int vertex_index) const` - Check anchor status
- `void clear_anchors()` - Remove all static anchors
- `std::vector<int> get_static_anchors()` - Get list of all static anchor indices

**Deformation:**
- `bool apply_deformation_to_mesh(int vertex_index, const Eigen::Vector3d& new_position)` - Apply deformation to mesh vertices (modifies mesh in-place)
- `Eigen::Matrix3Xd compute_deformation(int vertex_index, const Eigen::Vector3d& new_position)` - Compute and return deformed positions without modifying mesh
- `bool deform_with_temp_anchor(int vertex_index, const Eigen::Vector3d& new_position, Eigen::Matrix3Xd& output)` - Efficient in-place deformation for GUI
