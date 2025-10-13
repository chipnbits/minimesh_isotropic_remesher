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
  mohecore::Mesh_modifier_edge_collapse modi(mesh);

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

  // Initialize quadrics for all vertices
  printf("\n=== Initializing Quadrics === \n");
  modi.initialize_quadrics();
  printf("Quadrics initialized for all vertices.\n");

  // Initialize the priority queue
  printf("\n=== Initializing Priority Queue === \n");
  modi.initialize_priority_queue();
  printf("Priority queue initialized.\n");
  printf("Verification: PQ has %d edges (should be half of %d half-edges)\n",
         static_cast<int>(modi.get_top_n_candidates(mesh.n_active_half_edges()).size()),
         mesh.n_active_half_edges());

  // Get and print the top 10 candidates
  printf("\n=== Top 10 Edge Collapse Candidates (sqrt of origin vertex ID) === \n");
  std::vector<int> top_candidates = modi.get_top_n_candidates(10);

  for (size_t i = 0; i < top_candidates.size(); ++i) {
    int he_idx = top_candidates[i];
    auto he = mesh.half_edge_at(he_idx);
    int v_origin = he.origin().index();
    float metric = std::sqrt(static_cast<float>(v_origin));
    printf("Candidate %zu: Half-edge %d (vertex %d -> %d), metric=%.3f\n",
           i + 1, he_idx, v_origin, he.dest().index(), metric);
  }

  printf("\n");
  check_sanity_and_write_mesh("cli_mesh");
  return 0;
} // end of main()
