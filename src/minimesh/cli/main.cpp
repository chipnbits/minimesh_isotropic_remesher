//
// This is a bare executable that gets linked to the core library.
// You can use it to test your code, or start learning about the library.
//
// Here I have put an executable which reads a predefined mesh, flips a specific edge in that
// mesh, and then writes it back to .vtk and .obj formats. Feel free to play with the example, change
// it, or move it to a completely different file.
//

#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <unordered_map>

#include <minimesh/core/util/assert.hpp>
#include <minimesh/core/util/macros.hpp>

#include <minimesh/core/mohe/mesh_analysis.hpp>
#include <minimesh/core/mohe/mesh_connectivity.hpp>
#include <minimesh/core/mohe/mesh_io.hpp>
#include <minimesh/core/mohe/mesh_modifier_arap.hpp>
#include <minimesh/core/mohe/mesh_modifier_edge_collapse.hpp>
#include <minimesh/core/mohe/mesh_modifier_loop_subdivision.hpp>

using namespace minimesh;

// ===================
// EXAMPLE UTILITY FUNCTIONS
// ===================

namespace // Mark the start of anonymous namespace (cant be called from outside)
{

// Writes OBJ (and optional VTK) with a consistent naming scheme and sanity check.
// If `counter` is provided, it's incremented and used to suffix the file name.
// Otherwise a per-label static counter is used.
// Returns the full OBJ path written.
std::string
write_mesh_checked(mohecore::Mesh_connectivity & mesh,
    mohecore::Mesh_io & io,
    const std::string & label,
    int * counter = nullptr,
    bool also_vtk = false)
{
  force_assert(mesh.check_sanity_slowly());

  std::string base = "exports/" + label;

  if(counter)
  {
    if(*counter > 0)
      base += "_" + std::to_string(*counter);
    ++(*counter);
  }
  else
  {
    static std::unordered_map<std::string, int> per_label_counts;
    int & c = per_label_counts[label];
    if(c > 0)
      base += "_" + std::to_string(c);
    ++c;
  }

  const std::string obj_path = base + ".obj";
  printf("writing %s\n", obj_path.c_str());
  io.write_obj(obj_path);

  if(also_vtk)
  {
    const std::string vtk_path = base + ".vtk";
    printf("writing %s\n", vtk_path.c_str());
    io.write_vtk(vtk_path);
  }

  return obj_path;
}

///
/// Test ARAP deformation functionality
///
void
test_arap_deformation(mohecore::Mesh_connectivity & mesh, mohecore::Mesh_io & io)
{
  printf("\n=== ARAP DEFORMATION TEST === \n");

  // Initialize ARAP modifier
  mohecore::Mesh_modifier_arap modi_arap(mesh);
  modi_arap.initialize();
  printf("ARAP modifier initialized.\n");

  // Add some anchor points (first few vertices)
  printf("\nAdding anchor points...\n");
  for(int i = 0; i < 5 && i < mesh.n_active_vertices(); ++i)
  {
    if(modi_arap.add_anchor(i))
    {
      auto v = mesh.vertex_at(i);
      Eigen::Vector3d pos = v.xyz();
      printf("  - Anchor %d at position (%.3f, %.3f, %.3f)\n", i, pos.x(), pos.y(), pos.z());
    }
  }

  std::vector<int> anchors = modi_arap.get_static_anchors();
  printf("Total anchors set: %d\n", static_cast<int>(anchors.size()));

  // Test 1: compute_deformation() - Returns new matrix without modifying mesh
  printf("\n--- Test 1: compute_deformation() ---\n");
  int test_vertex = 10; // Pick a vertex to pull
  if(test_vertex < mesh.n_active_vertices())
  {
    auto v = mesh.vertex_at(test_vertex);
    Eigen::Vector3d original_pos = v.xyz();
    Eigen::Vector3d pulled_pos = original_pos + Eigen::Vector3d(0.5, 0.5, 0.0);

    printf("Testing compute_deformation() on vertex %d\n", test_vertex);
    printf("  Original position: (%.3f, %.3f, %.3f)\n", original_pos.x(), original_pos.y(), original_pos.z());
    printf("  Target position: (%.3f, %.3f, %.3f)\n", pulled_pos.x(), pulled_pos.y(), pulled_pos.z());

    try
    {
      Eigen::Matrix3Xd deformed = modi_arap.compute_deformation(test_vertex, pulled_pos);
      printf("  ✓ Deformation computed successfully (mesh unchanged)\n");
      printf("  Deformed matrix size: 3 x %ld\n", deformed.cols());
    }
    catch(const std::exception & e)
    {
      printf("  ✗ Deformation failed: %s\n", e.what());
    }
  }

  // Test 2: apply_deformation_to_mesh() - Modifies mesh directly
  printf("\n--- Test 2: apply_deformation_to_mesh() ---\n");
  test_vertex = 10;
  if(test_vertex < mesh.n_active_vertices())
  {
    auto v = mesh.vertex_at(test_vertex);
    Eigen::Vector3d original_pos = v.xyz();
    Eigen::Vector3d pulled_pos = original_pos + Eigen::Vector3d(0.5, 0.5, 0.0);

    printf("Testing apply_deformation_to_mesh() on vertex %d\n", test_vertex);
    printf("  Original position: (%.3f, %.3f, %.3f)\n", original_pos.x(), original_pos.y(), original_pos.z());
    printf("  Target position: (%.3f, %.3f, %.3f)\n", pulled_pos.x(), pulled_pos.y(), pulled_pos.z());

    if(modi_arap.apply_deformation_to_mesh(test_vertex, pulled_pos))
    {
      printf("  ✓ Deformation applied to mesh successfully\n");

      // Write out the deformed mesh
      printf("\nWriting deformed mesh...\n");
      write_mesh_checked(mesh, io, "arap_deformed", nullptr, false);
      printf("Deformed mesh written to exports/arap_deformed.{obj,vtk}\n");
    }
    else
    {
      printf("  ✗ Deformation failed\n");
    }
  }

  printf("\n=== ARAP TEST COMPLETE === \n");
}

//
// Create an example mesh file that we can read later.
//
// void
// write_example_mesh()
// {
//   FILE *fl = fopen("exports/example_mesh.obj", "w");

//   //
//   //
//   //  (6) ----- (7) --- (8)
//   //   |\      / \      /|
//   //   | \    /   \    / |
//   //   |  \  /     \  /  |
//   //   |   (4)------(5)  |
//   //   |   / \      / \  |
//   //   |  /   \    /   \ |
//   //   | /     \  /     \|
//   //   (1)----(2)------(3)
//   //

//   // Write the vertex coordinates
//   fprintf(fl, "v 0 0 0 \n");
//   fprintf(fl, "v 1 0 0 \n");
//   fprintf(fl, "v 2 0 0 \n");
//   fprintf(fl, "v 0.5 0.5 0 \n");
//   fprintf(fl, "v 1.5 0.5 0 \n");
//   fprintf(fl, "v 0 1 0 \n");
//   fprintf(fl, "v 1 1 0 \n");
//   fprintf(fl, "v 2 1 0 \n");

//   // Write the faces (vertices are index1-based in .obj format)
//   fprintf(fl, "\n");
//   fprintf(fl, "f 1 2 4 \n");
//   fprintf(fl, "f 2 3 5 \n");
//   fprintf(fl, "f 1 4 6 \n");
//   fprintf(fl, "f 4 2 5 \n");
//   fprintf(fl, "f 5 3 8 \n");
//   fprintf(fl, "f 6 4 7 \n");
//   fprintf(fl, "f 7 4 5 \n");
//   fprintf(fl, "f 7 5 8 \n");

//   fclose(fl);
// }

} // end of anonymus namespace

int
main(int argc, char ** argv)
{
  // Create a mesh_connectivity and a mesh reader
  mohecore::Mesh_connectivity mesh;
  mohecore::Mesh_io io(mesh);

  printf("=== MESH DEFORMATION TESTING === \n");

  // Check if filename provided as argument
  std::string filename;
  if(argc > 1)
  {
    filename = argv[1];
    printf("Processing file: %s\n", filename.c_str());
  }
  else
  {
    // Default behavior - use camel_simple.obj for testing
    filename = "./mesh/camel_simple.obj";
    printf("No filename provided, using default: %s\n", filename.c_str());
  }

  // Read the specified mesh file
  printf("reading %s \n", filename.c_str());
  io.read_obj_general(filename);

  printf("Total vertices: %d \n", mesh.n_active_vertices());
  printf("Total faces: %d \n", mesh.n_active_faces());
  printf("Total half-edges: %d \n", mesh.n_active_half_edges());

  // Run full ARAP deformation test
  test_arap_deformation(mesh, io);
  return 0;
} // end of main()
