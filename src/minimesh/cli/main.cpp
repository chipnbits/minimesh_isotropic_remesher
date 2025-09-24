//
// This is a bare executable that gets linked to the core library.
// You can use it to test your code, or start learning about the library.
//
// Here I have put an executable which reads a predefined mesh, flips a specific edge in that 
// mesh, and then writes it back to .vtk and .obj formats. Feel free to play with the example, change
// it, or move it to a completely different file.
//

#include <cstdio>
#include <string>
#include <chrono>

#include <minimesh/core/util/assert.hpp>
#include <minimesh/core/util/macros.hpp>

#include <minimesh/core/mohe/mesh_connectivity.hpp>
#include <minimesh/core/mohe/mesh_io.hpp>
#include <minimesh/core/mohe/mesh_modifier.hpp>
#include <minimesh/core/mohe/mesh_analysis.hpp>

using namespace minimesh;

// ===================
// EXAMPLE UTILITY FUNCTIONS
// ===================

namespace  // Mark the start of anonymous namespace (cant be called from outside)
{

//
// Create an example mesh file that we can read later.
//
void
write_example_mesh()
{
  FILE *fl = fopen("exports/example_mesh.obj", "w");

  //
  //
  //  (6) ----- (7) --- (8) 
  //   |\      / \      /|
  //   | \    /   \    / |
  //   |  \  /     \  /  |
  //   |   (4)------(5)  |
  //   |   / \      / \  |
  //   |  /   \    /   \ |
  //   | /     \  /     \|
  //   (1)----(2)------(3)
  //

  // Write the vertex coordinates
  fprintf(fl, "v 0 0 0 \n");
  fprintf(fl, "v 1 0 0 \n");
  fprintf(fl, "v 2 0 0 \n");
  fprintf(fl, "v 0.5 0.5 0 \n");
  fprintf(fl, "v 1.5 0.5 0 \n");
  fprintf(fl, "v 0 1 0 \n");
  fprintf(fl, "v 1 1 0 \n");
  fprintf(fl, "v 2 1 0 \n");

  // Write the faces (vertices are index1-based in .obj format)
  fprintf(fl, "\n");
  fprintf(fl, "f 1 2 4 \n");
  fprintf(fl, "f 2 3 5 \n");
  fprintf(fl, "f 1 4 6 \n");
  fprintf(fl, "f 4 2 5 \n");
  fprintf(fl, "f 5 3 8 \n");
  fprintf(fl, "f 6 4 7 \n");
  fprintf(fl, "f 7 4 5 \n");
  fprintf(fl, "f 7 5 8 \n");

  fclose(fl);
}

} // end of anonymus namespace

int main(int argc, char **argv)
{
  // Create a mesh_connectivity and a mesh reader
  mohecore::Mesh_connectivity mesh;
  mohecore::Mesh_io io(mesh);
  mohecore::Mesh_modifier modi(mesh);

  printf("=== MESH EDITTING EXAMPLE === \n");

  // Check if filename provided as argument
  std::string filename;
  if (argc > 1) {
    filename = argv[1];
    printf("Processing file: %s\n", filename.c_str());
  } else {
    // Default behavior - write and use example mesh
    printf("writing example_mesh.obj \n");
    write_example_mesh();
    filename = "./mesh/tetra.obj";
    printf("No filename provided, using default: %s\n", filename.c_str());
  }

  // Read the specified mesh file
  printf("reading %s \n", filename.c_str());
  io.read_obj_general(filename);

  // A lambda for checking the mesh sanity and writing it
  int write_count = 0;
  auto check_sanity_and_write_mesh = [&io, &mesh, &write_count](const std::string &label = "mesh") {
    force_assert(mesh.check_sanity_slowly());

    std::string base = "exports/" + label;
    if (write_count > 0) base += "_" + std::to_string(write_count);

    printf("writing %s.vtk and %s.obj\n", base.c_str(), base.c_str());
    io.write_obj(base + ".obj");
    io.write_vtk(base + ".vtk");

    ++write_count;

    // return the filepath
    return base + ".obj";
  };

  // Now check that the mesh is sane and write it  in both 
  // .vtk and .obj formats
  check_sanity_and_write_mesh();

  // Test edge division functionality
  printf("Total vertices: %d \n", mesh.n_active_vertices());
  printf("Total faces: %d \n", mesh.n_active_faces());
  printf("Total half-edges: %d \n", mesh.n_active_half_edges());

  printf("dividing with subdivision\n");

  auto start = std::chrono::high_resolution_clock::now();
  modi.subdivide_loop();
  modi.subdivide_loop();
  modi.subdivide_loop();
  modi.subdivide_loop();
  modi.subdivide_loop();
  modi.subdivide_loop();
  modi.subdivide_loop();
  modi.subdivide_loop();
  auto end   = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = end - start;
  printf("Subdivision took %.6f seconds\n", elapsed.count());

  printf("Total vertices: %d \n", mesh.n_active_vertices());
  printf("Total faces: %d \n", mesh.n_active_faces());
  printf("Total half-edges: %d \n", mesh.n_active_half_edges());

  check_sanity_and_write_mesh("subdivided_mesh");
  return 0;
} // end of main()
