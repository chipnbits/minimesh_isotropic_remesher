// From standard library
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


    // If the mesh was defragmented, you should have done:
    // Mesh_connectivity::Defragmentation_maps defrag;
    // globalvars::mesh.compute_defragmention_maps(defrag);
    // deformed_vertex_positions.col(defrag.old2new_vertex[j]) += pull_amount.cast<double>();
    globalvars::viewer.get_mesh_buffer().set_vertex_positions(globalvars::deformed_vertex_positions.cast<float>());

    // update positions (only the viewer)


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
  std::vector<int> anchors = globalvars::modi_arap.get_anchors();

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
}


int
main(int argc, char * argv[])
{
  // Remember current folder
  foldertools::pushd();

  // If no command line argument is specified, load a hardcoded mesh.e
  // Useful when debugging with visual studio.
  // Change the hardcoded address to your needs.
  if(argc == 1)
  {
    // retrieve filepath of /home/sghys/projects/CPSC524-modeling/mesh/tetra_complex.obj

    mohecore::Mesh_io(globalvars::mesh).read_auto("/home/sghys/projects/CPSC524-modeling/mesh/camel_simple.obj");
  }
  else // otherwise use the address specified in the command line
  {
    mohecore::Mesh_io(globalvars::mesh).read_auto(argv[1]);
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

  // Load the mesh in the viewer
  {
    mohecore::Mesh_connectivity::Defragmentation_maps defrag;
    globalvars::mesh.compute_defragmention_maps(defrag);
    globalvars::viewer.get_mesh_buffer().rebuild(globalvars::mesh, defrag);
  }

  // Setup background modifier for edge collapse
  globalvars::modi_edge.initialize();

  // Setup ARAP modifier
  globalvars::modi_arap.initialize();

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
