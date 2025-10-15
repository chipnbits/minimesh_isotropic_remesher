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
#include <cmath>

#include <minimesh/core/util/assert.hpp>
#include <minimesh/core/util/macros.hpp>

#include <minimesh/core/mohe/mesh_connectivity.hpp>
#include <minimesh/core/mohe/mesh_io.hpp>
#include <minimesh/core/mohe/mesh_modifier_loop_subdivision.hpp>
#include <minimesh/core/mohe/mesh_modifier_edge_collapse.hpp>
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

int main(int argc, char **argv)
{
  // Create a mesh_connectivity and a mesh reader
  mohecore::Mesh_connectivity mesh;
  mohecore::Mesh_io io(mesh);

  printf("=== MESH EDGE COLLAPSE PRIORITY QUEUE TEST === \n");

  // Check if filename provided as argument
  std::string filename;
  if (argc > 1) {
    filename = argv[1];
    printf("Processing file: %s\n", filename.c_str());
  } else {
    // Default behavior - use camel_simple.obj for testing
    filename = "./mesh/camel_simple.obj";
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

  printf("Total vertices: %d \n", mesh.n_active_vertices());
  printf("Total faces: %d \n", mesh.n_active_faces());
  printf("Total half-edges: %d \n", mesh.n_active_half_edges());

  // Initialize the mesh simplifier (computes quadrics and builds valid pairs)
  printf("\n=== Initializing Mesh Simplifier === \n");
  mohecore::Mesh_modifier_edge_collapse modi(mesh);
  modi.initialize();
  printf("Mesh simplifier initialized (quadrics and valid pairs computed).\n");

  // Count valid pairs by popping all entries
  int num_valid_pairs = 0;
  int num_edges = mesh.n_active_half_edges() / 2;
  mohecore::Mesh_modifier_edge_collapse::MergeCandidate candidate{0.0, {0, 0}, Eigen::Vector3d::Zero(), 0};

  while (modi.get_min_pair(candidate)) {
    num_valid_pairs++;
  }

  printf("Valid pairs found: %d\n", num_valid_pairs);
  printf("Expected (number of edges): %d\n", num_edges);

  if (num_valid_pairs == num_edges) {
    printf("✓ SUCCESS: Valid pairs count matches edge count!\n");
  } else {
    printf("✗ MISMATCH: Valid pairs count does not match edge count!\n");
  }

  // Print first few candidates
  printf("\n=== First 5 Edge Collapse Candidates (by QEM error) === \n");
  modi.initialize(); // Re-initialize since we popped everything

  for (int i = 0; i < 10 && modi.get_min_pair(candidate); ++i) {
    int v1 = candidate.pair.v1;
    int v2 = candidate.pair.v2;
    double error = candidate.error;
    printf("Candidate %d: Edge (%d, %d), error=%.9f, x_opt=(%.3f, %.3f, %.3f)\n",
           i + 1, v1, v2, error,
           candidate.x_opt[0], candidate.x_opt[1], candidate.x_opt[2]);
  }

  // use modi to collapse 10 edges or until no valid pairs remain
  printf("\n=== Collapsing 10 Edges === \n");
  int collapses = 0;
  while (collapses < 93 && modi.get_min_pair(candidate)) {
    printf("Collapsing edge (%d, %d) with error %.9f [collapse #%d]\n",
           candidate.pair.v1, candidate.pair.v2, candidate.error, collapses + 1);
    if (modi.collapse_edge(candidate)) {
      collapses++;
    }
  }

  printf("\n");
  check_sanity_and_write_mesh("cli_mesh");
  return 0;
} // end of main()
