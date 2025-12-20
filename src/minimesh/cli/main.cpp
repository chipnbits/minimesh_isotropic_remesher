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
#include <minimesh/core/mohe/fixed_uv_param.hpp>
#include <minimesh/core/mohe/lscm_uv_param.hpp>

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

  // Check if filename provided as argument
  std::string filename;
  std::string algorithm = "fixed";
  if(argc > 1)
  {
    filename = argv[1];
    printf("Processing file: %s\n", filename.c_str());
  }
  else
  {
    // Default behavior - use camel_simple.obj for testing
    filename = "./hw4_mesh/cat.obj";
    printf("No filename provided, using default: %s\n", filename.c_str());
  }
  if (argc > 2)
  {
    algorithm = argv[2];
    printf("Using algorithm: %s\n", algorithm.c_str());
  }

  // Read the specified mesh file
  printf("reading %s \n", filename.c_str());
  io.read_obj_general(filename);

  printf("Total vertices: %d \n", mesh.n_active_vertices());
  printf("Total faces: %d \n", mesh.n_active_faces());
  printf("Total half-edges: %d \n", mesh.n_active_half_edges());

  if (algorithm == "fixed")
  {
    // Create a UV parameterization object
    mohecore::Fixed_boundary_uv_param uv_param(mesh);

    uv_param.compute_parameterization();

    // overwrite coords into the 2D plane
    for(int v = 0; v < mesh.n_total_vertices(); ++v)
    {
      mohecore::Mesh_connectivity::Vertex_iterator vertex = mesh.vertex_at(v);
      if(vertex.is_active())
      {
        Eigen::Vector2d uv = uv_param.get_uv_at_vertex(vertex.index());
        vertex.data().xyz = Eigen::Vector3d(uv[0], uv[1], 0.0);
      }
    }
  }
  else if (algorithm == "lscm")
  {
    // Create a UV parameterization object
    mohecore::LSCM_uv_param uv_param(mesh, mohecore::LSCM_uv_param::PinningStrategy::MAX_DISTANCE);
    uv_param.compute_parameterization();

    for (int v = 0; v < mesh.n_total_vertices(); ++v)
    {
      mohecore::Mesh_connectivity::Vertex_iterator vertex = mesh.vertex_at(v);
      if(vertex.is_active())
      {
        Eigen::Vector2d uv = uv_param.get_uv_at_vertex(vertex.index());
        vertex.data().xyz = Eigen::Vector3d(uv[0], uv[1], 0.0);
      }
    }
  }
  else
  {
    printf("Unknown algorithm: %s\n", algorithm.c_str());
    return -1;
  }
  // Reuse the same filename but to export the result (take only filename without path or .obj)
  std::string mesh_out_path = filename.substr(filename.find_last_of("/\\") + 1);
  mesh_out_path = mesh_out_path.substr(0, mesh_out_path.find_last_of('.'));
  mesh_out_path = mesh_out_path + "_cli";
  write_mesh_checked(mesh, io, mesh_out_path, nullptr, false);

  return 0;
} // end of main()

