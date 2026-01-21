// From standard library
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>

// eigen
#include <Eigen/Core>
#include <Eigen/Geometry>

// core
#include <minimesh/core/mohe/mesh_connectivity.hpp>
#include <minimesh/core/mohe/mesh_io.hpp>
#include <minimesh/core/mohe/mesh_modifier_arap.hpp>
#include <minimesh/core/mohe/mesh_modifier_edge_collapse.hpp>
#include <minimesh/core/mohe/mesh_modifier_loop_subdivision.hpp>
#include <minimesh/core/mohe/remesher/remesher_isotropic.hpp>
#include <minimesh/core/mohe/fixed_uv_param.hpp>
#include <minimesh/core/mohe/lscm_uv_param.hpp>
#include <minimesh/core/util/assert.hpp>
#include <minimesh/core/util/foldertools.hpp>
#include <minimesh/core/util/numbers.hpp>

// gui
#include <minimesh/viz/mesh_viewer.hpp>
#include <minimesh/viz/opengl_headers.hpp>


using namespace minimesh;

// ======================================================
// Global variables
// ======================================================
namespace globalvars
{
Mesh_viewer viewer;
mohecore::Mesh_connectivity mesh;
mohecore::Mesh_connectivity mesh_backup; // Backup of original mesh before remeshing
mohecore::Mesh_modifier_uniform_remeshing modi_remesh(mesh);
mohecore::Mesh_modifier_edge_collapse modi_edge(mesh);
mohecore::Mesh_modifier_arap modi_arap(mesh);

//
int glut_main_window_id;
//
GLUI * glui;
//
int num_entities_to_simplify;
//
Eigen::Matrix3Xd deformed_vertex_positions; // Current displayed vertex positions (for mesh buffer)
//
bool show_collapse_overlay = false; // Toggle state for edge collapse visualization
int deform_mode = 0; // Toggle state for deform mode (enables anchor selection and ARAP)
int remesh_mode = 0; // Toggle state for remeshing mode (enables feature visualization)
//
int param_algorithm = 0; // 0 = Harmonic (Fixed), 1 = LSCM
//
float target_edge_length = 0.05f; // Target edge length for remeshing
//
// Remeshing parameters (with defaults matching the class)
int remesh_n_smoothing = 2;
float remesh_lambda = 0.4f;
float remesh_edge_flip_threshold = 0.8f;
float remesh_uncollapse_factor = 1.3f;
float remesh_feature_angle = 25.0f;
//
int const DEFORMATION_INTERVAL = 10; // Throttle interval in milliseconds
std::chrono::steady_clock::time_point last_deform_time = std::chrono::steady_clock::now();
}


// ======================================================
//              FREEGLUT CALL BACKS
// ======================================================
namespace freeglutcallback
{

// Forward declarations
void
update_anchor_visualization();

void
draw()
{
  globalvars::viewer.draw();
}


void
window_reshaped(int w, int h)
{
  bool should_redraw = false;
  should_redraw = should_redraw || globalvars::viewer.window_reshaped(w, h);

  if(should_redraw)
    glutPostRedisplay();
}


void
keyboard_pressed(unsigned char c, int x, int y)
{
  bool should_redraw = false;
  should_redraw = should_redraw || globalvars::viewer.keyboard_pressed(c, x, y);

  if(should_redraw)
    glutPostRedisplay();
}


void
keyboard_arrows_pressed(int c, int x, int y)
{
  bool should_redraw = false;
  should_redraw = should_redraw || globalvars::viewer.keyboard_pressed(c, x, y);

  if(should_redraw)
    glutPostRedisplay();
}


void
mouse_pushed(int button, int state, int x, int y)
{
  bool should_redraw = false;
  should_redraw = should_redraw || globalvars::viewer.mouse_pushed(button, state, x, y);

  //
  // NOTE: Sample of using Mesh_viewer for MESH DEFORMATION ASSINGMENT
  // Here is an example of how to use the selection feedback
  //
  {
    int clicked_on_vertex;
    bool did_user_click;
    globalvars::viewer.get_and_clear_vertex_selection(did_user_click, clicked_on_vertex);
    if(did_user_click)
    {
      printf("User just clicked on vertex %d \n", clicked_on_vertex);

      // Only toggle anchors when in deform mode AND in select mode
      if(globalvars::deform_mode && globalvars::viewer.get_mouse_function() == Mesh_viewer::MOUSE_SELECT)
      {
        // ARAP: Toggle anchor on/off for clicked vertex
        if(globalvars::modi_arap.is_anchor(clicked_on_vertex))
        {
          // Remove anchor
          if(globalvars::modi_arap.remove_anchor(clicked_on_vertex))
          {
            printf("Removed vertex %d as anchor\n", clicked_on_vertex);

            // Update visualization (automatically shown in deform mode)
            update_anchor_visualization();
            should_redraw = true;
          }
        }
        else
        {
          // Add anchor
          if(globalvars::modi_arap.add_anchor(clicked_on_vertex))
          {
            printf("Added vertex %d as anchor\n", clicked_on_vertex);

            // Update visualization (automatically shown in deform mode)
            update_anchor_visualization();
            should_redraw = true;
          }
        }
      }
    }
  }

  if(should_redraw)
    glutPostRedisplay();
}


void
mouse_moved(int x, int y)
{
  bool should_redraw = false;
  // Only process at a limited interval to avoid overloading
  auto now = std::chrono::steady_clock::now();
  auto duration_since_last = std::chrono::duration_cast<std::chrono::milliseconds>(now - globalvars::last_deform_time);
  if(duration_since_last.count() < globalvars::DEFORMATION_INTERVAL)
    return;
  globalvars::last_deform_time = now;
  should_redraw = should_redraw || globalvars::viewer.mouse_moved(x, y);
  {
    // Check for vertex pull and record new position and vertex
    bool has_pull_performed;
    Eigen::Vector3f pull_amount;
    int pulled_vert;
    globalvars::viewer.get_and_clear_vertex_displacement(has_pull_performed, pull_amount, pulled_vert);

    if(has_pull_performed && globalvars::viewer.get_mouse_function() == Mesh_viewer::MOUSE_MOVE_VERTEX)
    {
      // Update the position of the pulled vertex
      force_assert(pulled_vert != Mesh_viewer::invalid_index);
      globalvars::deformed_vertex_positions.col(pulled_vert) += pull_amount.cast<double>();

      // Check if we're in deform mode and if other vertices need to be pulled as well
      if(globalvars::deform_mode)
      {
        // ARAP: set the temporary anchor as an active anchor point

        Eigen::Vector3d pulled_position = globalvars::deformed_vertex_positions.col(pulled_vert);
        // update the entire mesh buffer based on arap deformation with anchors
        globalvars::modi_arap.deform_with_temp_anchor(
            pulled_vert, pulled_position, globalvars::deformed_vertex_positions);
      }
    }

    // Compute defragmentation maps to properly translate vertex positions
    // from old (sparse) indices to new (compact) mesh buffer indices
    mohecore::Mesh_connectivity::Defragmentation_maps defrag;
    globalvars::mesh.compute_defragmention_maps(defrag);

    // Create a compact position matrix for only active vertices
    Eigen::Matrix3Xf compact_positions(3, globalvars::mesh.n_active_vertices());
    for(int i = 0; i < globalvars::mesh.n_active_vertices(); ++i)
    {
      int old_vertex_idx = defrag.new2old_vertices[i];
      compact_positions.col(i) = globalvars::deformed_vertex_positions.col(old_vertex_idx).cast<float>();
    }

    // Update the mesh buffer with compact positions
    globalvars::viewer.get_mesh_buffer().set_vertex_positions(compact_positions);

    // Must rerender now.
    should_redraw = true;
  }


  if(should_redraw)
    glutPostRedisplay();
}

void
update_collapse_visualization()
{
  if(!globalvars::show_collapse_overlay)
  {
    // Clear the visualization by resetting all vertices to invisible
    int n_active_verts = globalvars::mesh.n_active_vertices();
    Eigen::Matrix4Xf vertex_colors(4, n_active_verts);
    vertex_colors.setConstant(0.7f); // invisible (R=0, G=0, B=0, A=0)
    globalvars::viewer.get_mesh_buffer().set_vertex_colors(vertex_colors);
    printf("Color overlay cleared\n");
    return;
  }

  // Show the top k candidates
  // printf("Updating visualization: showing top %d edge collapse candidates\n", globalvars::num_entities_to_simplify);

  // Get the top k candidates (returns half-edge indices)
  std::vector<int> top_k_he = globalvars::modi_edge.get_top_n_candidates(globalvars::num_entities_to_simplify);
  // printf("Found %zu candidates\n", top_k_he.size());

  // Compute defragmentation maps to handle vertex indexing
  mohecore::Mesh_connectivity::Defragmentation_maps defrag;
  globalvars::mesh.compute_defragmention_maps(defrag);

  // Prepare to color all vertices (default to white)
  int n_active_verts = globalvars::mesh.n_active_vertices();
  Eigen::Matrix4Xf vertex_colors(4, n_active_verts);
  vertex_colors.setConstant(0.7f); // Default: gray (R=0.7, G=0.7, B=0.7, A=0.7)

  // Color the vertices involved in top candidates
  // Top candidate: green
  // Rest: interpolate from yellow (rank 1) to red (last rank)
  int num_candidates = static_cast<int>(top_k_he.size());

  for(int rank = 0; rank < num_candidates; ++rank)
  {
    int he_idx = top_k_he[rank];
    auto he = globalvars::mesh.half_edge_at(he_idx);
    if(!he.is_active())
      continue;

    int origin_idx = he.origin().index();
    int dest_idx = he.dest().index();

    // Map to defragmented indices
    int origin_defrag = defrag.old2new_vertices[origin_idx];
    int dest_defrag = defrag.old2new_vertices[dest_idx];

    // Determine color based on rank
    float r, g, b, a;
    if(rank < 1)
    {
      // Top candidate: green
      r = 0.0f;
      g = 1.0f;
      b = 0.0f;
      a = 1.0f;

      // printf("Top candidate edge (%d, %d) \n", origin_idx, dest_idx);
    }
    else
    {
      // Rest of candidates: interpolate from yellow to red
      // Yellow (1, 1, 0) -> Red (1, 0, 0)
      float factor = static_cast<float>(rank - 1) / std::max(1.0f, static_cast<float>(num_candidates - 2));
      r = 1.0f;
      g = 1.0f; // Decreases from 1 to 0
      b = 0.0f;
      a = 1.0f - 0.8f * factor; // Decreases from 1 to 0.2
    }

    // Apply the same color to both vertices of the edge
    vertex_colors.col(origin_defrag) << r, g, b, a;
    vertex_colors.col(dest_defrag) << r, g, b, a;
  }

  // Apply vertex colors to mesh buffer
  globalvars::viewer.get_mesh_buffer().set_vertex_colors(vertex_colors);
}

void
subdivide_pressed(int)
{
  printf("Subdivide button was pressed \n");

  // Create a mesh modifier for the global mesh
  mohecore::Mesh_modifier_loop_subdivision modifier(globalvars::mesh);

  // Perform Loop subdivision (modifies mesh in-place)
  bool success = modifier.subdivide_loop();

  if(success)
  {
    printf("Loop subdivision completed successfully\n");

    // Rebuild the viewer with the modified mesh
    mohecore::Mesh_connectivity::Defragmentation_maps defrag;
    globalvars::mesh.compute_defragmention_maps(defrag);
    globalvars::viewer.get_mesh_buffer().rebuild(globalvars::mesh, defrag);

    // Reset deformed positions to match new mesh
    globalvars::deformed_vertex_positions.resize(3, globalvars::mesh.n_total_vertices());
    for(int i = 0; i < globalvars::mesh.n_total_vertices(); ++i)
    {
      globalvars::deformed_vertex_positions.col(i) = globalvars::mesh.vertex_at(i).xyz();
    }

    glutPostRedisplay();
  }
  else
  {
    printf("Loop subdivision failed\n");
  }
}

void
simplify_pressed(int)
{
  printf("Simplify button was pressed to remove %d entities \n", globalvars::num_entities_to_simplify);

  // Create a candidate structure to hold edge collapse information
  mohecore::Mesh_modifier_edge_collapse::MergeCandidate candidate{0.0, {0, 0}, Eigen::Vector3d::Zero(), 0};

  // Perform edge collapses up to the requested number
  int collapses = 0;
  while(collapses < globalvars::num_entities_to_simplify && globalvars::modi_edge.get_min_pair(candidate))
  {
    if(globalvars::modi_edge.collapse_edge(candidate))
    {
      printf("Collapsing edge (%d, %d) with error %.9f\n", candidate.pair.v1, candidate.pair.v2, candidate.error);
      collapses++;
    }
  }

  printf("Successfully collapsed %d edges\n", collapses);

  if(collapses > 0)
  {
    // Rebuild the viewer with the modified mesh
    mohecore::Mesh_connectivity::Defragmentation_maps defrag;
    globalvars::mesh.compute_defragmention_maps(defrag);
    globalvars::viewer.get_mesh_buffer().rebuild(globalvars::mesh, defrag);

    // Reset deformed positions to match new mesh
    globalvars::deformed_vertex_positions.resize(3, globalvars::mesh.n_total_vertices());
    for(int i = 0; i < globalvars::mesh.n_total_vertices(); ++i)
    {
      globalvars::deformed_vertex_positions.col(i) = globalvars::mesh.vertex_at(i).xyz();
    }

    // Update the visualization if overlay is enabled
    if(globalvars::show_collapse_overlay)
    {
      update_collapse_visualization();
    }

    glutPostRedisplay();
  }
  else
  {
    printf("No edges could be collapsed\n");
  }
}

void
show_top_candidates_pressed(int)
{
  // Toggle the color overlay on/off
  globalvars::show_collapse_overlay = !globalvars::show_collapse_overlay;

  if(globalvars::show_collapse_overlay)
  {
    printf("Color overlay enabled - showing top %d candidates\n", globalvars::num_entities_to_simplify);
  }
  else
  {
    printf("Color overlay disabled\n");
  }

  // Update the visualization based on new toggle state
  update_collapse_visualization();

  glutPostRedisplay();
}

void
show_spheres_pressed(int)
{
  //
  // Sample of using Mesh_viewer for MESH DEFORMATION ASSIGNMENT
  // Here I color the vertices (draw spheres on them)
  // Note that if you call rebuild, you have to redraw everything.
  //
  Eigen::VectorXi sphere_indices(3);
  sphere_indices << 0, 1, 2;
  Eigen::Matrix4Xf sphere_colors(4, 3);
  sphere_colors.col(0) << 1, 1, 0, 1;
  sphere_colors.col(1) << 0, 1, 1, 1;
  sphere_colors.col(2) << 0, 0, 1, 1;

  globalvars::viewer.get_mesh_buffer().set_colorful_spheres(sphere_indices, sphere_colors);

  glutPostRedisplay();
}

void
spinner_changed(int)
{
  // Update the visualization if overlay is enabled
  if(globalvars::show_collapse_overlay)
  {
    update_collapse_visualization();
    glutPostRedisplay();
  }
}

void
update_anchor_visualization()
{
  // Clear visualization if deform mode is off
  if(!globalvars::deform_mode)
  {
    Eigen::VectorXi empty_indices(0);
    Eigen::Matrix4Xf empty_colors(4, 0);
    globalvars::viewer.get_mesh_buffer().set_colorful_spheres(empty_indices, empty_colors);
    return;
  }

  // Get all anchor vertex indices from ARAP modifier
  std::vector<int> anchors = globalvars::modi_arap.get_static_anchors();

  // If no anchors, clear visualization
  if(anchors.empty())
  {
    Eigen::VectorXi empty_indices(0);
    Eigen::Matrix4Xf empty_colors(4, 0);
    globalvars::viewer.get_mesh_buffer().set_colorful_spheres(empty_indices, empty_colors);
    return;
  }

  // Convert vector to Eigen format
  Eigen::VectorXi sphere_indices(anchors.size());
  for(size_t i = 0; i < anchors.size(); ++i)
  {
    sphere_indices[i] = anchors[i];
  }

  // Create blue color for all anchors (RGBA: red=0, green=0, blue=1, alpha=1)
  Eigen::Matrix4Xf sphere_colors(4, anchors.size());
  for(size_t i = 0; i < anchors.size(); ++i)
  {
    sphere_colors.col(i) << 0.0f, 0.0f, 1.0f, 1.0f; // Blue
  }

  // Set the colored spheres
  globalvars::viewer.get_mesh_buffer().set_colorful_spheres(sphere_indices, sphere_colors);
}

Eigen::Matrix4Xf
angles_to_colors(const std::vector<double>& face_min_angles_radians,
                 mohecore::Mesh_connectivity& mesh,
                 const mohecore::Mesh_connectivity::Defragmentation_maps& defrag)
{
  int n_active_faces = mesh.n_active_faces();
  Eigen::Matrix4Xf colors(4, n_active_faces);

  // Tuning params
  constexpr double clamp_deg = 57.0;  
  constexpr float  power     = 3.0f;  

  for(int fid = 0; fid < mesh.n_total_faces(); ++fid)
  {
    if(!mesh.face_at(fid).is_active())
      continue;

    int active_face_idx = defrag.old2new_faces[fid];

    // Convert radians to degrees
    double angle_deg = face_min_angles_radians[fid] * 180.0 / M_PI;

    Eigen::Vector3f red(1.0f, 0.0f, 0.0f);
    Eigen::Vector3f gray(0.5f, 0.5f, 0.5f);

    Eigen::Vector3f rgb;

    if(angle_deg >= clamp_deg)
    {
      // Hard clamp to gray for angles in [58°, 60°] (and above, if any)
      rgb = gray;
    }
    else
    {
      // Map [0°, clamp_deg] -> s in [0, 1]
      float s = static_cast<float>(angle_deg / clamp_deg);
      s = std::min(1.0f, std::max(0.0f, s));


      float t = std::pow(s, 1/power);  

      // Interpolate red -> gray
      rgb = (1.0f - t) * red + t * gray;
    }

    // Alpha stays opaque
    colors.col(active_face_idx) = Eigen::Vector4f(rgb.x(), rgb.y(), rgb.z(), 1.0f);
  }

  return colors;
}

void
update_remesh_feature_visualization()
{
  // Clear visualization if remesh mode is off
  if(!globalvars::remesh_mode)
  {
    // Clear spheres
    Eigen::VectorXi empty_indices(0);
    Eigen::Matrix4Xf empty_colors(4, 0);
    globalvars::viewer.get_mesh_buffer().set_colorful_spheres(empty_indices, empty_colors);

    // Clear debug edges
    globalvars::viewer.get_mesh_buffer().set_debug_edge_colors(empty_indices, empty_colors);

    // Clear face colors (reset to default)
    int n_active_faces = globalvars::mesh.n_active_faces();
    Eigen::Matrix4Xf face_colors(4, n_active_faces);
    face_colors.setConstant(0.0f); // Transparent/default
    globalvars::viewer.get_mesh_buffer().set_face_colors(face_colors);

    return;
  }

  // Compute defragmentation maps to handle vertex indexing
  mohecore::Mesh_connectivity::Defragmentation_maps defrag;
  globalvars::mesh.compute_defragmention_maps(defrag);

  // Get feature edge and vertex type vectors from remesher via getters
  const std::vector<mohecore::Mesh_modifier_uniform_remeshing::VertexFeatureType>& vertex_feature_type =
      globalvars::modi_remesh.get_vertex_feature_types();

  // Collect feature vertices for sphere visualization
  std::vector<int> feature_vertex_indices;
  for(int v_idx = 0; v_idx < globalvars::mesh.n_total_vertices(); ++v_idx)
  {
    auto vert = globalvars::mesh.vertex_at(v_idx);
    if(!vert.is_active())
      continue;

    if(v_idx >= static_cast<int>(vertex_feature_type.size()))
      break;

    auto vtype = vertex_feature_type[v_idx];
    if(vtype == mohecore::Mesh_modifier_uniform_remeshing::EDGE ||
       vtype == mohecore::Mesh_modifier_uniform_remeshing::CORNER)
    {
      // Get the defragmented index for this vertex
      int active_idx = defrag.old2new_vertices[v_idx];
      feature_vertex_indices.push_back(active_idx);
    }
  }

  // Apply min angle visualization to faces
  const std::vector<double>& min_angles = globalvars::modi_remesh.get_face_min_angles();
  Eigen::Matrix4Xf face_colors = angles_to_colors(min_angles, globalvars::mesh, defrag);
  globalvars::viewer.get_mesh_buffer().set_face_colors(face_colors);

  // DEBUG: Visualize feature edges
  {
    const std::vector<bool>& is_feature_edge = globalvars::modi_remesh.get_feature_edges();
    std::vector<int> feature_edge_indices;

    // Build mapping from half-edges to edge_conn indices
    // This mirrors the logic in mesh_buffer.cpp rebuild()
    int edge_idx = 0;
    for(int i = 0; i < globalvars::mesh.n_active_half_edges(); ++i)
    {
      int old_he_idx = defrag.new2old_half_edges[i];
      auto he = globalvars::mesh.half_edge_at(old_he_idx);

      // Only process each edge once (same logic as mesh_buffer rebuild)
      if(he.index() > he.twin().index())
      {
        // Check if this half-edge OR its twin is marked as a feature edge
        if(old_he_idx < static_cast<int>(is_feature_edge.size()) &&
           (is_feature_edge[he.index()] || is_feature_edge[he.twin().index()]))
        {
          feature_edge_indices.push_back(edge_idx);
        }
        edge_idx++;
      }
    }

    // Set bright cyan color for feature edges
    if(!feature_edge_indices.empty())
    {
      Eigen::VectorXi edge_indices(feature_edge_indices.size());
      Eigen::Matrix4Xf edge_colors(4, feature_edge_indices.size());

      for(size_t i = 0; i < feature_edge_indices.size(); ++i)
      {
        edge_indices[i] = feature_edge_indices[i];
        // Bright cyan color (R=0, G=1, B=1, A=1)
        edge_colors.col(i) << 0.0f, 1.0f, 1.0f, 1.0f;
      }

      globalvars::viewer.get_mesh_buffer().set_debug_edge_colors(edge_indices, edge_colors);
    }
    else
    {
      // Clear debug edges
      Eigen::VectorXi empty_indices(0);
      Eigen::Matrix4Xf empty_colors(4, 0);
      globalvars::viewer.get_mesh_buffer().set_debug_edge_colors(empty_indices, empty_colors);
    }
  }

  // Draw spheres on feature vertices (light blue)
  if(!feature_vertex_indices.empty())
  {
    Eigen::VectorXi sphere_indices(feature_vertex_indices.size());
    Eigen::Matrix4Xf sphere_colors(4, feature_vertex_indices.size());

    for(size_t i = 0; i < feature_vertex_indices.size(); ++i)
    {
      sphere_indices[i] = feature_vertex_indices[i];
      // Light blue color (R=0.3, G=0.6, B=1.0, A=1)
      sphere_colors.col(i) << 0.3f, 0.6f, 1.0f, 1.0f;
    }

    // globalvars::viewer.get_mesh_buffer().set_colorful_spheres(sphere_indices, sphere_colors);
  }
  else
  {
    // Clear spheres if no feature vertices
    Eigen::VectorXi empty_indices(0);
    Eigen::Matrix4Xf empty_colors(4, 0);
    globalvars::viewer.get_mesh_buffer().set_colorful_spheres(empty_indices, empty_colors);
  }
}

void
deform_mode_changed(int)
{
  // When deform mode is toggled, update anchor visualization
  if(globalvars::deform_mode)
  {
    printf("Deform mode enabled - reinitializing ARAP with current mesh and resetting to rest pose\n");

    // Reinitialize modi_arap with the current mesh (which may have been subdivided/simplified)
    globalvars::modi_arap.initialize();

    // Reset deformed positions to match the current mesh (snap back to rest pose)
    globalvars::deformed_vertex_positions.resize(3, globalvars::mesh.n_total_vertices());
    for(int i = 0; i < globalvars::mesh.n_total_vertices(); ++i)
    {
      globalvars::deformed_vertex_positions.col(i) = globalvars::mesh.vertex_at(i).xyz();
    }

    // Update the viewer buffer with the reset positions
    globalvars::viewer.get_mesh_buffer().set_vertex_positions(globalvars::deformed_vertex_positions.cast<float>());

    update_anchor_visualization();
  }
  else
  {
    printf("Deform mode disabled - anchor visualization turned off\n");
    // Clear anchor visualization
    update_anchor_visualization();
  }

  glutPostRedisplay();
}

void
remesh_mode_changed(int)
{
  // When remesh mode is toggled, update visualization
  if(globalvars::remesh_mode)
  {
    printf("Remesh mode enabled - saving mesh backup and visualizing feature edges and vertices\n");

    // Save a backup of the current mesh before any remeshing operations
    globalvars::mesh_backup.copy(globalvars::mesh);
    printf("Mesh backup saved with %d vertices, %d faces\n",
           globalvars::mesh_backup.n_active_vertices(),
           globalvars::mesh_backup.n_active_faces());

    // Apply current parameter values to the remesher
    globalvars::modi_remesh.set_n_smoothing_iters(globalvars::remesh_n_smoothing);
    globalvars::modi_remesh.set_lambda_smoothing_damping(static_cast<double>(globalvars::remesh_lambda));
    globalvars::modi_remesh.set_edge_flip_threshold_dot(static_cast<double>(globalvars::remesh_edge_flip_threshold));
    globalvars::modi_remesh.set_uncollapse_threshold_factor(static_cast<double>(globalvars::remesh_uncollapse_factor));
    globalvars::modi_remesh.set_feature_angle_degrees(static_cast<double>(globalvars::remesh_feature_angle));

    // Re-initialize remesher with new feature angle (this marks features with the new threshold)
    globalvars::modi_remesh.initialize();
    globalvars::modi_remesh.compute_face_min_angles();

    // Update visualization (feature lists already built at initialization)
    update_remesh_feature_visualization();

    printf("Feature edges and vertices visualized\n");
  }
  else
  {
    printf("Remesh mode disabled - feature visualization turned off\n");
    // Clear feature visualization
    update_remesh_feature_visualization();
  }

  glutPostRedisplay();
}

void
reset_mesh_to_backup()
{
  if(!globalvars::remesh_mode)
  {
    printf("Cannot reset: not in remesh mode\n");
    return;
  }

  printf("Resetting mesh to backup...\n");

  // Restore mesh from backup
  globalvars::mesh.copy(globalvars::mesh_backup);

  // Rebuild the viewer with the restored mesh
  mohecore::Mesh_connectivity::Defragmentation_maps defrag;
  globalvars::mesh.compute_defragmention_maps(defrag);
  globalvars::viewer.get_mesh_buffer().rebuild(globalvars::mesh, defrag);

  // Reset deformed positions to match restored mesh
  globalvars::deformed_vertex_positions.resize(3, globalvars::mesh.n_total_vertices());
  for(int i = 0; i < globalvars::mesh.n_total_vertices(); ++i)
  {
    globalvars::deformed_vertex_positions.col(i) = globalvars::mesh.vertex_at(i).xyz();
  }

  // Re-initialize remesher with current parameters
  globalvars::modi_remesh.set_n_smoothing_iters(globalvars::remesh_n_smoothing);
  globalvars::modi_remesh.set_lambda_smoothing_damping(static_cast<double>(globalvars::remesh_lambda));
  globalvars::modi_remesh.set_edge_flip_threshold_dot(static_cast<double>(globalvars::remesh_edge_flip_threshold));
  globalvars::modi_remesh.set_uncollapse_threshold_factor(static_cast<double>(globalvars::remesh_uncollapse_factor));
  globalvars::modi_remesh.set_feature_angle_degrees(static_cast<double>(globalvars::remesh_feature_angle));
  globalvars::modi_remesh.initialize();
  globalvars::modi_remesh.compute_face_min_angles();

  // Update visualization
  update_remesh_feature_visualization();

  printf("Mesh restored to backup with %d vertices, %d faces\n",
         globalvars::mesh.n_active_vertices(),
         globalvars::mesh.n_active_faces());

  glutPostRedisplay();
}

void
remesh_reset_button_pressed(int)
{
  reset_mesh_to_backup();
}

void
remesh_parameter_changed(int)
{
  // When parameters change, reset mesh and reinitialize with new parameters
  if(!globalvars::remesh_mode)
  {
    // If not in remesh mode, just update the parameters without resetting
    return;
  }

  printf("Remesh parameters changed - resetting mesh and reinitializing...\n");
  reset_mesh_to_backup();
}

void
remesh_single_pass_pressed(int)
{
  printf("Remesh single pass button pressed with target edge length: %.4f\n", globalvars::target_edge_length);

  // Apply current parameter values to the remesher before running
  globalvars::modi_remesh.set_n_smoothing_iters(globalvars::remesh_n_smoothing);
  globalvars::modi_remesh.set_lambda_smoothing_damping(static_cast<double>(globalvars::remesh_lambda));
  globalvars::modi_remesh.set_edge_flip_threshold_dot(static_cast<double>(globalvars::remesh_edge_flip_threshold));
  globalvars::modi_remesh.set_uncollapse_threshold_factor(static_cast<double>(globalvars::remesh_uncollapse_factor));
  // Note: Feature angle not updated here as it requires re-initialization

  // Run a single remeshing pass (feature lists are managed internally by the remesher)
  globalvars::modi_remesh.run_single_pass(static_cast<double>(globalvars::target_edge_length), globalvars::remesh_n_smoothing);

  printf("Single remeshing pass completed\n");

  // Compute face minimal angles for quality visualization
  globalvars::modi_remesh.compute_face_min_angles();

  // Rebuild the viewer with the modified mesh
  mohecore::Mesh_connectivity::Defragmentation_maps defrag;
  globalvars::mesh.compute_defragmention_maps(defrag);
  globalvars::viewer.get_mesh_buffer().rebuild(globalvars::mesh, defrag);

  // Reset deformed positions to match new mesh
  globalvars::deformed_vertex_positions.resize(3, globalvars::mesh.n_total_vertices());
  for(int i = 0; i < globalvars::mesh.n_total_vertices(); ++i)
  {
    globalvars::deformed_vertex_positions.col(i) = globalvars::mesh.vertex_at(i).xyz();
  }

  // Update visualization (handles both feature vertices and min angle coloring based on remesh_mode)
  update_remesh_feature_visualization();

  glutPostRedisplay();
}

void
parameterization_pressed(int)
{
  printf("Parameterization button pressed with algorithm: %s\n",
         globalvars::param_algorithm == 0 ? "Harmonic (Fixed)" : "LSCM");

  bool success = false;

  if(globalvars::param_algorithm == 0)
  {
    // Harmonic parameterization (Fixed boundary)
    mohecore::Fixed_boundary_uv_param uv_param(globalvars::mesh);
    success = uv_param.compute_parameterization();

    if(success)
    {
      printf("Harmonic parameterization computed successfully\n");

      // Update deformed_vertex_positions with UV coordinates (z=0)
      for(int v = 0; v < globalvars::mesh.n_total_vertices(); ++v)
      {
        mohecore::Mesh_connectivity::Vertex_iterator vertex = globalvars::mesh.vertex_at(v);
        if(vertex.is_active())
        {
          Eigen::Vector2d uv = uv_param.get_uv_at_vertex(vertex.index());
          globalvars::deformed_vertex_positions.col(v) = Eigen::Vector3d(uv[0], uv[1], 0.0);
        }
      }
    }
  }
  else if(globalvars::param_algorithm == 1)
  {
    // LSCM parameterization
    mohecore::LSCM_uv_param uv_param(globalvars::mesh, mohecore::LSCM_uv_param::PinningStrategy::MAX_DISTANCE);
    success = uv_param.compute_parameterization();

    if(success)
    {
      printf("LSCM parameterization computed successfully\n");

      // Update deformed_vertex_positions with UV coordinates (z=0)
      for(int v = 0; v < globalvars::mesh.n_total_vertices(); ++v)
      {
        mohecore::Mesh_connectivity::Vertex_iterator vertex = globalvars::mesh.vertex_at(v);
        if(vertex.is_active())
        {
          Eigen::Vector2d uv = uv_param.get_uv_at_vertex(vertex.index());
          globalvars::deformed_vertex_positions.col(v) = Eigen::Vector3d(uv[0], uv[1], 0.0);
        }
      }
    }
  }

  if(success)
  {
    // Update the viewer buffer with the new flattened positions
    globalvars::viewer.get_mesh_buffer().set_vertex_positions(globalvars::deformed_vertex_positions.cast<float>());

    // Compute new bounding box from UV coordinates to rescale camera view
    Eigen::AlignedBox3f bbox;
    for(int v = 0; v < globalvars::mesh.n_total_vertices(); ++v)
    {
      mohecore::Mesh_connectivity::Vertex_iterator vertex = globalvars::mesh.vertex_at(v);
      if(vertex.is_active())
      {
        bbox.extend(globalvars::deformed_vertex_positions.col(v).cast<float>());
      }
    }

    // Re-initialize viewer with new bounding box to rescale camera
    globalvars::viewer.initialize(bbox);
    printf("Camera view rescaled to UV parameterization\n");

    // Redraw
    glutPostRedisplay();
  }
  else
  {
    printf("Parameterization failed\n");
  }
}
}

int
main(int argc, char * argv[])
{
  // Remember current folder
  foldertools::pushd();

  // Check for --isometric flag
  bool use_isometric_view = false;
  std::string mesh_filename;

  // Parse command line arguments
  for(int i = 1; i < argc; ++i)
  {
    std::string arg = argv[i];
    if(arg == "--isometric")
    {
      use_isometric_view = true;
      printf("Isometric view mode enabled\n");
    }
    else if(arg[0] != '-')
    {
      mesh_filename = arg;
    }
  }

  // If no command line argument is specified, load a hardcoded mesh.e
  // Useful when debugging with visual studio.
  // Change the hardcoded address to your needs.
  if(mesh_filename.empty())
  {
    // retrieve filepath of /home/sghys/projects/CPSC524-modeling/mesh/tetra_complex.obj

    mohecore::Mesh_io(globalvars::mesh).read_auto("/home/sghys/projects/CPSC524-modeling/hw4_mesh/cat.obj");
  }
  else // otherwise use the address specified in the command line
  {
    mohecore::Mesh_io(globalvars::mesh).read_auto(mesh_filename);
  }

  // Initialize GLUT window
  glutInit(&argc, argv);
  glutInitWindowSize(800, 600);
  glutInitDisplayMode(GLUT_STENCIL | GLUT_DEPTH | GLUT_RGBA | GLUT_DOUBLE);
  globalvars::glut_main_window_id = glutCreateWindow("Mesh Viewer");

  // Initialize GLUI window for buttons and ...
  globalvars::glui = GLUI_Master.create_glui("Controls");
  globalvars::glui->set_main_gfx_window(globalvars::glut_main_window_id);

  // Register callbacks
  glutDisplayFunc(freeglutcallback::draw);
  GLUI_Master.set_glutReshapeFunc(freeglutcallback::window_reshaped);
  GLUI_Master.set_glutKeyboardFunc(freeglutcallback::keyboard_pressed);
  GLUI_Master.set_glutSpecialFunc(freeglutcallback::keyboard_arrows_pressed);
  GLUI_Master.set_glutMouseFunc(freeglutcallback::mouse_pushed);
  glutMotionFunc(freeglutcallback::mouse_moved);
  GLUI_Master.set_glutIdleFunc(NULL);

  // Initialize the viewer (it needs the bounding box of the mesh)
  Eigen::AlignedBox3f bbox;
  for(int v = 0; v < globalvars::mesh.n_total_vertices(); ++v)
  {
    mohecore::Mesh_connectivity::Vertex_iterator vertex = globalvars::mesh.vertex_at(v);
    if(vertex.is_active())
    {
      bbox.extend(vertex.xyz().cast<float>());
    }
  }
  globalvars::viewer.initialize(bbox);

  // Set isometric view if requested
  if(use_isometric_view)
  {
    globalvars::viewer.set_isometric_view();
  }

  // Load the mesh in the viewer
  {
    mohecore::Mesh_connectivity::Defragmentation_maps defrag;
    globalvars::mesh.compute_defragmention_maps(defrag);
    globalvars::viewer.get_mesh_buffer().rebuild(globalvars::mesh, defrag);
  }

  // Setup background modifiers
  globalvars::modi_edge.initialize();
  globalvars::modi_remesh.initialize(); // Re-initialize after mesh is loaded

  // Compute initial quality data (visualization will be applied when remesh mode is enabled)
  globalvars::modi_remesh.compute_face_min_angles();

  //
  // Add radio buttons to see which mesh components to view
  // Please view GLUI's user manual to learn more.
  //

  GLUI_Panel * panel_view = globalvars::glui->add_panel("View mesh components");
  globalvars::glui->add_checkbox_to_panel(panel_view, "Show vertices", &globalvars::viewer.get_draw_vertices());
  globalvars::glui->add_checkbox_to_panel(panel_view, "Show edges", &globalvars::viewer.get_draw_edges());
  globalvars::glui->add_checkbox_to_panel(panel_view, "Show faces", &globalvars::viewer.get_draw_faces());
  globalvars::glui->add_checkbox_to_panel(panel_view, "Show axis", &globalvars::viewer.get_draw_axis());
  globalvars::glui->add_checkbox_to_panel(panel_view, "Show lighting", &globalvars::viewer.get_has_lighting());

  //
  // Add radio buttons to determine mouse left click functionality
  //
  GLUI_Panel * panel_mouse_func = globalvars::glui->add_panel("Mouse functionality");
  GLUI_RadioGroup * radio_group_mouse_func =
      globalvars::glui->add_radiogroup_to_panel(panel_mouse_func, &globalvars::viewer.get_mouse_function());
  for(int i = 0; i < Mesh_viewer::MOUSE_INVALID; ++i)
  {
    if(i == Mesh_viewer::MOUSE_VIEW)
      globalvars::glui->add_radiobutton_to_group(radio_group_mouse_func, "Pan and zoom");
    if(i == Mesh_viewer::MOUSE_SELECT)
      globalvars::glui->add_radiobutton_to_group(radio_group_mouse_func, "Select vertex");
    if(i == Mesh_viewer::MOUSE_MOVE_VERTEX)
      globalvars::glui->add_radiobutton_to_group(radio_group_mouse_func, "Move vertex");
  }

  //
  // Add subdivide button
  //
  GLUI_Button * button_subdivide =
      globalvars::glui->add_button("Subdivide Loop", -1, freeglutcallback::subdivide_pressed);
  button_subdivide->set_w(200);

  //
  // Add simplify button and a spinner to read how many entities to remove
  //
  globalvars::num_entities_to_simplify = 0;
  GLUI_Spinner * spinner_simplify = globalvars::glui->add_spinner("# of entities to simplify",
      GLUI_SPINNER_INT,
      &globalvars::num_entities_to_simplify,
      -1,
      freeglutcallback::spinner_changed);
  spinner_simplify->set_alignment(GLUI_ALIGN_CENTER);
  spinner_simplify->set_w(300);
  spinner_simplify->set_int_limits(0, 10000); // Prevent negative values
  spinner_simplify->set_speed(0.0001); // Increment by 1 when clicking arrows

  GLUI_Button * button_simplify = globalvars::glui->add_button("Simplify", -1, freeglutcallback::simplify_pressed);
  button_simplify->set_w(200);

  //
  // Add button to visualize top edge collapse candidates (toggle)
  //
  GLUI_Button * button_visualize =
      globalvars::glui->add_button("Visualize", -1, freeglutcallback::show_top_candidates_pressed);
  button_visualize->set_w(200);

  //
  // Add ARAP Deformation Panel
  //
  GLUI_Panel * panel_arap = globalvars::glui->add_panel("ARAP Deformation");

  // Add Deform checkbox that controls the entire deformation mode
  globalvars::glui->add_checkbox_to_panel(
      panel_arap, "Deform", &globalvars::deform_mode, -1, freeglutcallback::deform_mode_changed);

  //
  // Add UV Parameterization Panel
  //
  GLUI_Panel * panel_param = globalvars::glui->add_panel("UV Parameterization");

  // Add radio buttons for algorithm selection
  GLUI_RadioGroup * radio_group_param =
      globalvars::glui->add_radiogroup_to_panel(panel_param, &globalvars::param_algorithm);
  globalvars::glui->add_radiobutton_to_group(radio_group_param, "Harmonic");
  globalvars::glui->add_radiobutton_to_group(radio_group_param, "LSCM");

  // Add Parameterization button
  GLUI_Button * button_param =
      globalvars::glui->add_button_to_panel(panel_param, "Compute Parameterization", -1, freeglutcallback::parameterization_pressed);
  button_param->set_w(200);

  //
  // Add Remeshing Panel
  //
  GLUI_Panel * panel_remesh = globalvars::glui->add_panel("Remeshing");

  // Add Remesh Mode checkbox
  globalvars::glui->add_checkbox_to_panel(
      panel_remesh, "Remesh Mode", &globalvars::remesh_mode, -1, freeglutcallback::remesh_mode_changed);

  // Add spinner for target edge length
  GLUI_Spinner * spinner_edge_length = globalvars::glui->add_spinner_to_panel(panel_remesh,
      "Target Edge Length",
      GLUI_SPINNER_FLOAT,
      &globalvars::target_edge_length);
  spinner_edge_length->set_alignment(GLUI_ALIGN_CENTER);
  spinner_edge_length->set_w(300);
  spinner_edge_length->set_float_limits(0.001f, 100.0f);
  spinner_edge_length->set_speed(0.001f);

  // Add parameter spinners
  GLUI_Spinner * spinner_n_smoothing = globalvars::glui->add_spinner_to_panel(panel_remesh,
      "N Smoothing Iterations",
      GLUI_SPINNER_INT,
      &globalvars::remesh_n_smoothing,
      -1,
      freeglutcallback::remesh_parameter_changed);
  spinner_n_smoothing->set_alignment(GLUI_ALIGN_CENTER);
  spinner_n_smoothing->set_w(300);
  spinner_n_smoothing->set_int_limits(0, 20);

  GLUI_Spinner * spinner_lambda = globalvars::glui->add_spinner_to_panel(panel_remesh,
      "Lambda Smoothing",
      GLUI_SPINNER_FLOAT,
      &globalvars::remesh_lambda,
      -1,
      freeglutcallback::remesh_parameter_changed);
  spinner_lambda->set_alignment(GLUI_ALIGN_CENTER);
  spinner_lambda->set_w(300);
  spinner_lambda->set_float_limits(0.0f, 1.0f);
  spinner_lambda->set_speed(0.01f);

  GLUI_Spinner * spinner_edge_flip = globalvars::glui->add_spinner_to_panel(panel_remesh,
      "Edge Flip Threshold",
      GLUI_SPINNER_FLOAT,
      &globalvars::remesh_edge_flip_threshold,
      -1,
      freeglutcallback::remesh_parameter_changed);
  spinner_edge_flip->set_alignment(GLUI_ALIGN_CENTER);
  spinner_edge_flip->set_w(300);
  spinner_edge_flip->set_float_limits(0.0f, 1.0f);
  spinner_edge_flip->set_speed(0.01f);

  GLUI_Spinner * spinner_uncollapse = globalvars::glui->add_spinner_to_panel(panel_remesh,
      "Uncollapse Factor",
      GLUI_SPINNER_FLOAT,
      &globalvars::remesh_uncollapse_factor,
      -1,
      freeglutcallback::remesh_parameter_changed);
  spinner_uncollapse->set_alignment(GLUI_ALIGN_CENTER);
  spinner_uncollapse->set_w(300);
  spinner_uncollapse->set_float_limits(0.5f, 3.0f);
  spinner_uncollapse->set_speed(0.01f);

  GLUI_Spinner * spinner_feature_angle = globalvars::glui->add_spinner_to_panel(panel_remesh,
      "Feature Angle (degrees)",
      GLUI_SPINNER_FLOAT,
      &globalvars::remesh_feature_angle,
      -1,
      freeglutcallback::remesh_parameter_changed);
  spinner_feature_angle->set_alignment(GLUI_ALIGN_CENTER);
  spinner_feature_angle->set_w(300);
  spinner_feature_angle->set_float_limits(0.0f, 180.0f);
  spinner_feature_angle->set_speed(1.0f);

  // Add reset button
  GLUI_Button * button_remesh_reset =
      globalvars::glui->add_button_to_panel(panel_remesh, "Reset Mesh", -1, freeglutcallback::remesh_reset_button_pressed);
  button_remesh_reset->set_w(200);

  // Add button to run single remeshing pass
  GLUI_Button * button_remesh_pass =
      globalvars::glui->add_button_to_panel(panel_remesh, "Run Single Pass", -1, freeglutcallback::remesh_single_pass_pressed);
  button_remesh_pass->set_w(200);

  //
  // Add show spheres button to demo how to draw spheres on top of the vertices
  //
  globalvars::glui->add_button("Demo Showing Spheres", -1, freeglutcallback::show_spheres_pressed);

  //
  // Save the initial vertex positions
  //
  globalvars::deformed_vertex_positions.resize(3, globalvars::mesh.n_total_vertices());
  for(int i = 0; i < globalvars::mesh.n_total_vertices(); ++i)
  {
    globalvars::deformed_vertex_positions.col(i) = globalvars::mesh.vertex_at(i).xyz();
  }

  // Sync all glui variables
  globalvars::glui->sync_live();

  // Start main loop
  glutPostRedisplay(); // Draw everything again just for caution.
  glutMainLoop();

  // revert back to initial folder
  foldertools::popd();

  return 0;
}
