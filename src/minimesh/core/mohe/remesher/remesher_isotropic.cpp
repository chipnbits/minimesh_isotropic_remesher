// Implementation of Uniform Explicit Isotropic Remeshing
// Based on Botsch and Kobbelt (2004)
//
// Project Part 1: Baseline remeshing algorithm
//

#include <Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <minimesh/core/mohe/mesh_analysis.hpp>
#include <minimesh/core/mohe/remesher/remesher_isotropic.hpp>
#include <minimesh/core/util/assert.hpp>
#include <set>
#include <vector>

namespace minimesh
{
namespace mohecore
{

// ============================================================
// Initialization
// ============================================================

void
Mesh_modifier_uniform_remeshing::initialize()
{
  // Mark feature edges and vertices on initialization
  printf("Marking feature edges and vertices...\n");
  mark_feature_edges_and_vertices();

  int total_feature_edges = 0;
  for(bool is_feature : _is_feature_edge)
  {
    if(is_feature)
      total_feature_edges++;
  }
  printf("Total feature edges marked: %d\n", total_feature_edges);
}

// ============================================================
// Main Entry Point
// ============================================================

void
Mesh_modifier_uniform_remeshing::remesh(double target_edge_length, int iterations)
{
  // Start time
  auto start_time = std::chrono::high_resolution_clock::now();
  for(int i = 0; i < iterations; ++i)
  {
    printf("Remeshing Iteration %d / %d\n", i + 1, iterations);
    run_single_pass(target_edge_length, N_SMOOTHING_ITERS);
  }
  // End time
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end_time - start_time;
  printf("Remeshing completed in %.3f seconds.\n", elapsed.count());

  // Compute face minimal angles for quality visualization
  compute_face_min_angles();
}

void
Mesh_modifier_uniform_remeshing::run_single_pass(double target_edge_length, int tangential_smoothing_iters)
{
  printf("Splitting long edges...\n");
  split_long_edges(target_edge_length, 4.0 / 3.0);
  printf("Collapsing short edges...\n");
  collapse_short_edges(target_edge_length, 4.0 / 5.0);
  printf("Flipping edges to optimize valence...\n");
  flip_edges_to_optimize_valence();
  printf("Performing tangential smoothing...\n");
  tangential_smoothing(tangential_smoothing_iters, SmoothingType::Barycenters);
}

void
Mesh_modifier_uniform_remeshing::remesh_to_target_edge_count(int target_edge_count,
    int iterations,
    const double rel_error_tol)
{
  double net_len = 0.0;
  int active_edges = 0;
  for(int e = 0; e < mesh().n_total_half_edges(); ++e)
  {
    auto he = mesh().half_edge_at(e);
    if(he.is_active() && e < he.twin().index())
    { // Count each edge once
      net_len += get_edge_length(e);
      active_edges++;
    }
  }
  double current_L = net_len / active_edges;
  double early_termination_seen = false;

  printf("Targeting %d edges. Starting with %d edges, Avg L = %f\n", target_edge_count, active_edges, current_L);

  for(int i = 0; i < iterations; ++i)
  {

    // Get edge count
    active_edges = 0;
    for(int he = 0; he < mesh().n_total_half_edges(); ++he)
    {
      if(mesh().half_edge_at(he).is_active() && he < mesh().half_edge_at(he).twin().index())
      {
        active_edges++;
      }
    }

    // Assume equilateral triangles and constant total mesh area A
    // A = N (sqrt(3)/4)*L^2 for A constant -> N_new *L_new^2 = N_cur*L_cur^2
    double scaling_factor = std::sqrt(static_cast<double>(active_edges) / target_edge_count);
    if(std::abs(1 - scaling_factor) < rel_error_tol)
    {
      if(!early_termination_seen)
      {
        early_termination_seen = true;
        printf("Close to target %i within tolerance, running one final iteration.\n", target_edge_count);
      }
      else
      {
        printf("Reached close to target %i within tolerance after %d iterations.\n", active_edges, i);
        break;
      }
    }
    else
    {
      early_termination_seen = false;
    }

    // Clamp to stop extreme changes
    scaling_factor = std::max(0.8, std::min(1.2, scaling_factor));
    double next_target_L = current_L * scaling_factor;

    printf("Iter %d: Edges %d -> Target %d. Updating L: %f -> %f\n",
        i + 1,
        active_edges,
        target_edge_count,
        current_L,
        next_target_L);

    run_single_pass(next_target_L, N_SMOOTHING_ITERS);

    for(int he = 0; he < mesh().n_total_half_edges(); ++he)
    {
      if(mesh().half_edge_at(he).is_active() && he < mesh().half_edge_at(he).twin().index())
      {
        active_edges++;
      }
    }


    current_L = next_target_L;
  }

  // Compute face minimal angles for quality visualization
  compute_face_min_angles();
}

void
Mesh_modifier_uniform_remeshing::remesh_to_target_vertex_count(int target_vertex_count,
    int iterations,
    const double rel_error_tol)
{
  double net_len = 0.0;
  int active_edges = 0;
  for(int e = 0; e < mesh().n_total_half_edges(); ++e)
  {
    auto he = mesh().half_edge_at(e);
    if(he.is_active() && e < he.twin().index())
    { // Count each edge once
      net_len += get_edge_length(e);
      active_edges++;
    }
  }
  double current_L = net_len / active_edges;

  int active_vertices = 0;
  for(int v = 0; v < mesh().n_total_vertices(); ++v)
  {
    if(mesh().vertex_at(v).is_active())
    {
      active_vertices++;
    }
  }

  double early_termination_seen = false;

  printf("Targeting %d vertices. Starting with %d vertices, Avg L = %f\n", target_vertex_count, active_vertices, current_L);

  for(int i = 0; i < iterations; ++i)
  {

    // Get vertex count
    active_vertices = 0;
    for(int v = 0; v < mesh().n_total_vertices(); ++v)
    {
      if(mesh().vertex_at(v).is_active())
      {
        active_vertices++;
      }
    }

    // Assume equilateral triangles with valence 6 and constant total mesh area A
    // A = (sqrt(3)/4)*L^2 * (N*6/2) for A constant -> N_new *L_new^2 = N_cur*L_cur^2
    double scaling_factor = std::sqrt(static_cast<double>(active_vertices) / target_vertex_count);
    if(std::abs(1 - scaling_factor) < rel_error_tol)
    {
      if(!early_termination_seen)
      {
        early_termination_seen = true;
        printf("Close to target %i within tolerance, running one final iteration.\n", target_vertex_count);
      }
      else
      {
        printf("Reached close to target %i within tolerance after %d iterations.\n", active_vertices, i);
        break;
      }
    }
    else
    {
      early_termination_seen = false;
    }

    // Clamp to stop extreme changes
    scaling_factor = std::max(0.8, std::min(1.2, scaling_factor));
    double next_target_L = current_L * scaling_factor;

    printf("Iter %d: Vertices %d -> Target %d. Updating L: %f -> %f\n",
        i + 1,
        active_vertices,
        target_vertex_count,
        current_L,
        next_target_L); 
    run_single_pass(next_target_L, N_SMOOTHING_ITERS);
    for(int v = 0; v < mesh().n_total_vertices(); ++v)
    {
      if(mesh().vertex_at(v).is_active())
      {
        active_vertices++;
      }
    }
    current_L = next_target_L;
  }

  // Compute face minimal angles for quality visualization
  compute_face_min_angles();
}


// ============================================================
// Core Remeshing Operations
// ============================================================

void
Mesh_modifier_uniform_remeshing::split_long_edges(double target_length, double alpha)
{
  double t = target_length * alpha;
  const int initial_he_count = mesh().n_total_half_edges();

  // Make a scratch list for half-edges to mark visted ones
  std::vector<bool> visited(initial_he_count, false);

  auto he = Mesh_connectivity::Half_edge_iterator();
  int twin_idx;
  double length;
  for(int i = 0; i < initial_he_count; ++i)
  {
    he = mesh().half_edge_at(i);
    twin_idx = he.twin().index();
    if(visited[i] || visited[twin_idx])
    {
      visited[i] = true;
      visited[twin_idx] = true;
      continue;
    }
    visited[i] = true;
    visited[twin_idx] = true;

    if(!he.is_active())
      continue;
    length = get_edge_length(he.index());
    if(length > t)
    {
      split_edge(he.index(), 0.5);
    }
  }
}

void
Mesh_modifier_uniform_remeshing::collapse_short_edges(double target_length, double beta)
{
  double t = target_length * beta;
  const int initial_he_count = mesh().n_total_half_edges();

  // Make a scratch list for half-edges to mark visted ones
  std::vector<bool> visited(initial_he_count, false);

  auto he = Mesh_connectivity::Half_edge_iterator();
  int twin_idx;
  double length;
  for(int i = 0; i < initial_he_count; ++i)
  {
    he = mesh().half_edge_at(i);
    twin_idx = he.twin().index();
    if(visited[i] || visited[twin_idx])
    {
      visited[i] = true;
      visited[twin_idx] = true;
      continue;
    }
    visited[i] = true;
    visited[twin_idx] = true;

    if(!he.is_active())
      continue;
    length = get_edge_length(he.index());
    if(length < t)
    {
      collapse_edge(he.index(), t);
    }
  }
}

void
Mesh_modifier_uniform_remeshing::flip_edges_to_optimize_valence()
{
  const int initial_he_count = mesh().n_total_half_edges();
  std::vector<bool> visited(initial_he_count, false);

  auto he = Mesh_connectivity::Half_edge_iterator();
  int twin_idx;
  for(int i = 0; i < initial_he_count; ++i)
  {
    he = mesh().half_edge_at(i);
    twin_idx = he.twin().index();
    if(visited[i] || visited[twin_idx])
    {
      visited[i] = true;
      visited[twin_idx] = true;
      continue;
    }
    visited[i] = true;
    visited[twin_idx] = true;

    if(!he.is_active())
      continue;

    if(should_flip_edge(he.index()) && is_legal_flip(he.index()))
    {
      if(flip_edge(he.index()))
      {
        // printf("Flipped edge %d to improve valence\n", he.index());
      }
    }
  }
}

void
Mesh_modifier_uniform_remeshing::tangential_smoothing(int smoothing_iters, SmoothingType type)
{
  const int num_vertices = mesh().n_total_vertices();
  // Placeholder for adaptive sizing field (not used in isotropic remeshing) default uniform
  std::vector<double> adaptive_sizing_field = std::vector<double>(num_vertices, 1.0);


  for(int iter = 0; iter < smoothing_iters; ++iter)
  {
    // Precompute mesh geometry information
    auto geometry_cache = compute_geometry_cache();
    // Storage for new vertex positions and if they need update
    auto new_pos = std::vector<Eigen::Vector3d>(num_vertices, Eigen::Vector3d(0, 0, 0));
    std::vector<bool> needs_update(num_vertices, false);

    for(int vi = 0; vi < num_vertices; ++vi)
    {
      auto vert = mesh().vertex_at(vi);
      VertexFeatureType vt = _vertex_feature_type[vi];
      if(!vert.is_active())
        continue;

      if(vt == CORNER)
      {
        new_pos[vi] = vert.xyz(); // Keep corner vertices fixed
        continue;
      }

      // if (vt == EDGE)
      // {
      //   new_pos[vi] = vert.xyz(); // Keep edge vertices fixed
      //   continue;
      //   // Skip movement along edge for now
      // }

      if (vt == EDGE)
      {
        // 1. Collect the two feature neighbors
        std::vector<Eigen::Vector3d> n_pos;
        std::vector<int> n_indices;

        auto ring = _m.vertex_ring_at(vi);
        do
        {
            int he_idx = ring.half_edge().index();
            if(_is_feature_edge[he_idx])
            {
                n_pos.push_back(ring.half_edge().origin().xyz());
                n_indices.push_back(ring.half_edge().origin().index());
            }
        } while(ring.advance());

        // We expect exactly 2 neighbors for a feature EDGE vertex
        if(n_pos.size() == 2)
        {
            Eigen::Vector3d p = vert.xyz();
            Eigen::Vector3d n0 = n_pos[0];
            Eigen::Vector3d n1 = n_pos[1];

            // 2. Calculate Weighted Lengths (Metric Space)
            // We divide Euclidean length by sizing field: Higher Sizing = Longer edges allowed.
            // We want: len0 / size0  ==  len1 / size1
            
            double size_here = adaptive_sizing_field[vi];
            double size0 = (adaptive_sizing_field[n_indices[0]] + size_here) * 0.5;
            double size1 = (adaptive_sizing_field[n_indices[1]] + size_here) * 0.5;

            double len0 = (n0 - p).norm();
            double len1 = (n1 - p).norm();

            double metric_len0 = len0 / size0;
            double metric_len1 = len1 / size1;

            // 3. Slide vertex along the incident edge to equalize metric lengths
            // If metric_len0 > metric_len1, n0 is "too far". Slide towards n0.
            
            Eigen::Vector3d move_dir(0,0,0);
            double imbalance = 0.0;

            if (metric_len0 > metric_len1)
            {
                move_dir = (n0 - p).normalized(); // Slide towards n0
                imbalance = metric_len0 - metric_len1;
            }
            else
            {
                move_dir = (n1 - p).normalized(); // Slide towards n1
                imbalance = metric_len1 - metric_len0;
            }

            // Apply movement
            // Scale the move so we don't jump too far (Relaxation step)
            // Converting metric imbalance back to Euclidean approximation for the step size
            double step_size = imbalance * size_here * 0.5 * LAMBDA_SMOOTHING_DAMPING;
            
            // Safety: Don't move past the neighbor (prevent tangling)
            double max_dist = (metric_len0 > metric_len1) ? len0 : len1;
            if(step_size > max_dist * 0.4) step_size = max_dist * 0.4;

            new_pos[vi] = p + move_dir * step_size;
            needs_update[vi] = true;
        }
        else
        {
            // Corner case (endpoints of wires) or non-manifold
            new_pos[vi] = vert.xyz();
        }
        
        continue; // Done with EDGE vertex
    }
      // ==========================
      // Simple non-FeatureCase: SMOOTH vertex
      // ==========================

      Eigen::Vector3d target_pos = vert.xyz();

      if(type == SmoothingType::Uniform)
      {
        // Move towards average of neighbors
        Eigen::Vector3d sum_neighbors(0, 0, 0);
        int neighbor_count = 0;
        auto ring = _m.vertex_ring_at(vi);
        do
        {
          auto nbr_p = ring.half_edge().origin().xyz();
          sum_neighbors += nbr_p;
          neighbor_count++;
        } while(ring.advance());

        assert(neighbor_count > 0);
        target_pos = sum_neighbors / static_cast<double>(neighbor_count);
      }
      else if(type == SmoothingType::Barycenters)
      {
        // Weighted average of face centroids
        Eigen::Vector3d sum_weighted_centroids(0, 0, 0);
        double weight_sum = 0.0;

        auto ring = _m.vertex_ring_at(vi);
        do
        {
          auto he = ring.half_edge();
          auto face = he.face();

          if(face.is_equal(mesh().hole()))
            continue;

          int v1 = he.origin().index();
          int v2 = he.next().origin().index();
          int v3 = he.next().next().origin().index();
          double sizing_at_barycenter =
              (adaptive_sizing_field[v1] + adaptive_sizing_field[v2] + adaptive_sizing_field[v3]) / 3.0;

          double weight = geometry_cache.face_areas[face.index()] * sizing_at_barycenter;
          sum_weighted_centroids += weight * geometry_cache.face_barycenters[face.index()];
          weight_sum += weight;
        } while(ring.advance());

        if(weight_sum > 1e-8)
        {
          target_pos = sum_weighted_centroids / weight_sum;
        }
        else
        {
          target_pos = vert.xyz(); // Fallback to current position
        }
      } // End of barycenters

      // Move vertex towards target position with damping
      Eigen::Vector3d old_pos = vert.xyz();
      Eigen::Vector3d move_vec = target_pos - old_pos;
      Eigen::Vector3d normal = geometry_cache.vertex_normals[vi];

      // Project move_vec onto tangent plane to stop volume changes
      Eigen::Vector3d tangential_move = move_vec - (move_vec.dot(normal)) * normal;
      new_pos[vi] = old_pos + LAMBDA_SMOOTHING_DAMPING * tangential_move;
      needs_update[vi] = true;
    } // End of vertex loop

    // Apply new positions
    for(int vi = 0; vi < num_vertices; ++vi)
    {
      if(needs_update[vi])
      {
        mesh().vertex_at(vi).data().xyz = new_pos[vi];
      }
    }

    printf("Total updated vertices in smoothing iteration %d: %i\n",
        iter + 1,
        static_cast<int>(std::count(needs_update.begin(), needs_update.end(), true)));
  } // End smoothing iters
}

// ============================================================
// Low-Level Mesh Modifiers
// ============================================================

int
Mesh_modifier_uniform_remeshing::split_edge(int he_index, double weight)
{
  // Validate input index
  if(he_index < 0 || he_index >= mesh().n_total_half_edges())
    return mesh().invalid_index;

  // Get the half-edge and its twin
  auto he_1 = mesh().half_edge_at(he_index);
  if(!he_1.is_active())
    return mesh().invalid_index;

  auto tw_1 = he_1.twin();
  if(!tw_1.is_active())
    return mesh().invalid_index;

  bool split_is_feature = _is_feature_edge[he_1.index()];

  // Get all involved vertices, half-edges, and faces
  auto he_2 = he_1.next();
  auto tw_2 = tw_1.next();
  auto he_3 = he_1.prev();
  auto tw_3 = tw_1.prev();
  auto F1 = he_1.face();
  auto F2 = tw_1.face();

  // Get the two vertices of the edge
  auto v1 = he_1.origin();
  auto v2 = tw_1.origin(); // Twin's origin is he's destination
  auto v3 = he_2.dest(); // Opposite vertex in face 1
  auto v4 = tw_2.dest(); // Opposite vertex in face 2

  bool he_face_is_hole = F1.is_equal(mesh().hole());
  bool tw_face_is_hole = F2.is_equal(mesh().hole());

  // Compute weighted position: weight * origin + (1-weight) * dest
  Eigen::Vector3d new_pos = weight * v1.xyz() + (1.0 - weight) * v2.xyz();

  // Create new vertex at weighted position
  auto v5 = mesh().add_vertex();
  v5.data().xyz = new_pos;
  tw_1.data().origin = v5.index(); // Update twin to point to new vertex

  // Create two new half-edges to split the original edge (dont recycle to avoid resplitting)
  auto he_4 = mesh().add_half_edge(false);
  auto tw_4 = mesh().add_half_edge(false);
  he_4.data().twin = tw_4.index();
  he_4.data().origin = v5.index();
  v5.data().half_edge = he_4.index();

  tw_4.data().twin = he_4.index();
  tw_4.data().origin = v2.index();
  v2.data().half_edge = tw_4.index();

  if(he_face_is_hole)
  {
    he_1.data().next = he_4.index();
    he_4.data().prev = he_1.index();
    he_4.data().next = he_2.index();
    he_2.data().prev = he_4.index();
    he_4.data().face = F1.index();
  }
  else
  {
    // Add new face, new twinned edge splitting F1
    auto F3 = mesh().add_face();
    auto he_5 = mesh().add_half_edge(false);
    auto he_6 = mesh().add_half_edge(false);
    he_5.data().twin = he_6.index();
    he_5.data().origin = v3.index();
    v3.data().half_edge = he_5.index();
    he_6.data().twin = he_5.index();
    he_6.data().origin = v5.index();
    v5.data().half_edge = he_6.index();

    // Set up new F1 cycle and face
    // Next cycle
    he_1.data().next = he_6.index();
    he_6.data().next = he_3.index();
    he_3.data().next = he_1.index(); // already set

    // Prev cycle
    he_3.data().prev = he_6.index();
    he_6.data().prev = he_1.index();
    he_1.data().prev = he_3.index(); // already set

    // Face assignments
    he_1.data().face = F1.index();
    he_6.data().face = F1.index();
    he_3.data().face = F1.index();
    F1.data().half_edge = he_1.index();

    // Set up new F3 cycle and face
    // Next cycle
    he_4.data().next = he_2.index();
    he_2.data().next = he_5.index();
    he_5.data().next = he_4.index();

    // Prev cycle
    he_5.data().prev = he_2.index();
    he_2.data().prev = he_4.index();
    he_4.data().prev = he_5.index();

    // Face assignments
    he_4.data().face = F3.index();
    he_2.data().face = F3.index();
    he_5.data().face = F3.index();
    F3.data().half_edge = he_4.index();
  }

  if(tw_face_is_hole)
  {
    tw_1.data().prev = tw_4.index();
    tw_4.data().next = tw_1.index();
    tw_4.data().prev = tw_3.index();
    tw_3.data().next = tw_4.index();
    tw_4.data().face = F2.index();
  }
  else
  {
    // Add new face, new twinned edge splitting F2
    auto F4 = mesh().add_face();
    auto tw_5 = mesh().add_half_edge(false);
    auto tw_6 = mesh().add_half_edge(false);
    tw_5.data().twin = tw_6.index();
    tw_5.data().origin = v5.index();
    v5.data().half_edge = tw_5.index();
    tw_6.data().twin = tw_5.index();
    tw_6.data().origin = v4.index();
    v4.data().half_edge = tw_6.index();

    // Set up new F2 cycle and face
    // Next cycle
    // tw_1.data().next = tw_2.index(); already set
    tw_2.data().next = tw_6.index();
    tw_6.data().next = tw_1.index();

    // Prev cycle
    tw_1.data().prev = tw_6.index();
    tw_6.data().prev = tw_2.index();
    // tw_2.data().prev = tw_1.index(); already set

    // Face assignments
    tw_1.data().face = F2.index();
    tw_2.data().face = F2.index();
    tw_6.data().face = F2.index();
    F2.data().half_edge = tw_1.index();

    // Set up new F4 cycle and face
    // Next cycle
    tw_4.data().next = tw_5.index();
    tw_5.data().next = tw_3.index();
    tw_3.data().next = tw_4.index();

    // Prev cycle
    tw_5.data().prev = tw_4.index();
    tw_4.data().prev = tw_3.index();
    tw_3.data().prev = tw_5.index();

    // Face assignments
    tw_4.data().face = F4.index();
    tw_3.data().face = F4.index();
    tw_5.data().face = F4.index();
    F4.data().half_edge = tw_4.index();
  }

  // Update all feature markings
  // Resize feature edge and vertex type containers if needed
  ensure_property_containers_size();
  _is_feature_edge[he_4.index()] = split_is_feature;
  _is_feature_edge[tw_4.index()] = split_is_feature;

  if(!he_face_is_hole)
  {
    // Update he5 and he6 if they were created
    _is_feature_edge[he_2.next().index()] = false;
    _is_feature_edge[he_1.next().index()] = false;
  }

  if(!tw_face_is_hole)
  {
    // Update tw5 and tw6 if they were created
    _is_feature_edge[tw_4.next().index()] = false;
    _is_feature_edge[tw_2.next().index()] = false;
  }
  if(split_is_feature)
  {
    // New vertex is feature if original
    _vertex_feature_type[v5.index()] = EDGE;

  }
  else
  {
    // New vertex is simple if edge is not feature
    _vertex_feature_type[v5.index()] = SIMPLE;
  }

  return v5.index();
}

bool
Mesh_modifier_uniform_remeshing::collapse_edge(int he_index, double threshold)
{
  int v1 = mesh().half_edge_at(he_index).origin().index();
  int v2 = mesh().half_edge_at(he_index).dest().index();
  auto v1_iter = mesh().vertex_at(v1);
  auto v2_iter = mesh().vertex_at(v2);

  bool v1_boundary = analysis::vertex_is_boundary(_m, v1);
  bool v2_boundary = analysis::vertex_is_boundary(_m, v2);

  // Determine conditional collapse position (do not move boundary vertices)
  Eigen::Vector3d new_pos;
  VertexFeatureType t1 = _vertex_feature_type[v1];
  VertexFeatureType t2 = _vertex_feature_type[v2];

  if (t1 == CORNER || t2 == CORNER)
  {
    // Do not collapse using corners
    return false;
  }

  if(v1_boundary && !v2_boundary)
  {
    if(t2 == EDGE) return false;
    new_pos = v1_iter.xyz(); // Keep v1 where it is
  }
  else if(!v1_boundary && v2_boundary)
  {
    if(t1 == EDGE) return false;
    new_pos = v2_iter.xyz(); // Snap v1 to v2's position
  }
  else
  {
    if(t1 == CORNER)
      new_pos = v1_iter.xyz();
    else if(t2 == CORNER)
      new_pos = v2_iter.xyz();
    else if(t1 == EDGE && t2 == SIMPLE)
      new_pos = v1_iter.xyz();
    else if(t1 == SIMPLE && t2 == EDGE)
      new_pos = v2_iter.xyz();
    else if(t1 == EDGE && t2 == EDGE)
      new_pos = 0.5 * (v1_iter.xyz() + v2_iter.xyz()); // Midpoint on feature edge
    else
      new_pos = 0.5 * (v2_iter.xyz() + v1_iter.xyz()); // Midpoint for simple-simple case
  }

  if(!is_legal_collapse(v1, v2, new_pos, UNCOLLAPSE_THRESHOLD_FACTOR * threshold))
    return false;


  auto he_v1_v2 = mesh().half_edge_at(he_index);
  auto he_v2_v1 = he_v1_v2.twin();
  auto u1_iter = he_v1_v2.next().dest();
  auto u2_iter = he_v2_v1.prev().origin();


  auto face_top = he_v1_v2.face();
  auto face_bottom = he_v2_v1.face();
  bool top_is_hole = face_top.is_equal(mesh().hole());
  bool bot_is_hole = face_bottom.is_equal(mesh().hole());

  // Abort if hole is complicating the collapse
  if(!top_is_hole && !bot_is_hole)
    {
        auto n1 = he_v1_v2.next().twin().face();
        auto n2 = he_v2_v1.prev().twin().face();

        if(n1.is_equal(mesh().hole()) || n2.is_equal(mesh().hole()))
        {
            return false; // ABORT: Neighbor is a hole, too risky.
        }
    }

    // Abort for other cases involving holes
    if(!top_is_hole && bot_is_hole)
    {
        auto n1 = he_v1_v2.next().twin().face();
        if(n1.is_equal(mesh().hole()))
        {
            return false; // ABORT: Neighbor is a hole, too risky.
        }
    }

    if(top_is_hole && !bot_is_hole)
    {
        auto n2 = he_v2_v1.prev().twin().face();
        if(n2.is_equal(mesh().hole()))
        {
            return false; // ABORT: Neighbor is a hole, too risky.
        }
    }


  // Remap all half-edges originating from v2 to originate from v1
  relabel_vertex(v2, v1);

  // Handle simple case of non-boundary faces
  if(!top_is_hole && !bot_is_hole)
  {
    auto face_f1 = he_v1_v2.next().twin().face();
    auto face_f2 = he_v2_v1.prev().twin().face();
    auto he_f1_associate = he_v1_v2.prev();
    auto he_f2_associate = he_v2_v1.next();
    // Get the half-edges that will need to be re-linked
    auto he_f1_bottom = he_v1_v2.next().twin().next();
    auto he_f2_top = he_v2_v1.prev().twin().prev();

    // outer half-edges of two adjacent v2 side faces
    // defined so that s1=s2 when collapsing edge of degenerate form
    auto he_s1 = he_v1_v2.next().twin().prev();
    auto he_s2 = he_v2_v1.prev().twin().next();

    // Top Face (f1): v1 -> v2 -> TopVert -> v1
    // Survivor: he_f1_associate (prev)
    // Victim:   he_v1_v2.next()
    bool top_victim_feature = _is_feature_edge[he_v1_v2.next().index()];
    if(top_victim_feature)
    {
        // Transfer feature status to the surviving edge and its twin
        _is_feature_edge[he_f1_associate.index()] = true;
        _is_feature_edge[he_f1_associate.twin().index()] = true;
    }

    // Bottom Face (f2): v2 -> v1 -> BotVert -> v2
    // Survivor: he_f2_associate (next)
    // Victim:   he_v2_v1.prev()
    bool bot_victim_feature = _is_feature_edge[he_v2_v1.prev().index()];
    if(bot_victim_feature)
    {
        // Transfer feature status to the surviving edge and its twin
        _is_feature_edge[he_f2_associate.index()] = true;
        _is_feature_edge[he_f2_associate.twin().index()] = true;
    }

    face_top.deactivate();
    face_bottom.deactivate();
    // Deactivate 6 half-edges
    he_v1_v2.next().twin().deactivate();
    he_v1_v2.next().deactivate();
    he_v1_v2.deactivate();
    he_v2_v1.prev().twin().deactivate();
    he_v2_v1.prev().deactivate();
    he_v2_v1.deactivate();
    // Deactivate vertex v2
    v2_iter.deactivate();

    // Update feature edge markings
    _is_feature_edge[he_v1_v2.next().twin().index()] = false;
    _is_feature_edge[he_v1_v2.next().index()] = false;
    _is_feature_edge[he_v1_v2.index()] = false;
    _is_feature_edge[he_v2_v1.prev().twin().index()] = false;
    _is_feature_edge[he_v2_v1.prev().index()] = false;
    _is_feature_edge[he_v2_v1.index()] = false;
    _vertex_feature_type[v2] = SIMPLE; // SIMPLE for safety

    // Check if this is a degenerate case where face_f1 == face_f2
    bool is_degenerate = he_s1.is_equal(he_s2);

    if(is_degenerate)
    {

      // f1 and f2 the same and s1, s2 the same
      he_f1_associate.data().prev = he_s1.index();
      he_s1.data().next = he_f1_associate.index();

      he_f1_associate.data().next = he_f2_associate.index();
      he_f2_associate.data().prev = he_f1_associate.index();

      he_f2_associate.data().next = he_s2.index();
      he_s2.data().prev = he_f2_associate.index();

      he_f1_associate.data().face = face_f1.index();
      he_f2_associate.data().face = face_f1.index();
      face_f1.data().half_edge = he_f1_associate.index();

      // // Run assertions to ensure connectivity is correct
      // // Check f1 = f2
      // assert(face_f1.is_equal(face_f2));
      // // Check s1 = s2
      // assert(he_s1.is_equal(he_s2));
      // // Check that cycle of 3 half-edges is correct and closed
      // const int he_start_index = he_f1_associate.index();
      // auto he_iter = he_f1_associate;
      // do{
      //   // Print out half edge twin, origin, prev, next, face for debugging
      //   std::cout << "Half-edge " << he_iter.index() << ": "
      //             << "twin=" << he_iter.twin().index() << ", "
      //             << "origin=" << he_iter.origin().index() << ", "
      //             << "prev=" << he_iter.prev().index() << ", "
      //             << "next=" << he_iter.next().index() << ", "
      //             << "face=" << he_iter.face().index() << std::endl;
      //   assert(he_iter.face().is_equal(face_f1));
      //   he_iter = he_iter.next();
      // }
      // while(he_iter.index() != he_start_index);
    }
    else
    {
      // Top edges re-linking
      he_f1_associate.data().next = he_f1_bottom.index();
      he_f1_bottom.data().prev = he_f1_associate.index();
      he_f1_associate.data().prev = he_f1_bottom.next().index();
      he_f1_bottom.next().data().next = he_f1_associate.index();
      // Relink face f1
      he_f1_associate.data().face = face_f1.index();
      face_f1.data().half_edge = he_f1_associate.index();

      // Bottom edges re-linking
      he_f2_associate.data().next = he_f2_top.prev().index();
      he_f2_top.prev().data().prev = he_f2_associate.index();
      he_f2_associate.data().prev = he_f2_top.index();
      he_f2_top.data().next = he_f2_associate.index();
      // Relink face f2
      he_f2_associate.data().face = face_f2.index();
      face_f2.data().half_edge = he_f2_associate.index();
    }

    // Update vertices to point to valid half-edges
    v1_iter.data().half_edge = he_f1_associate.twin().index();
    u1_iter.data().half_edge = he_f1_associate.index();
    u2_iter.data().half_edge = he_f2_associate.twin().index();
  }

  // ---------------------------------------------------------
  // CASE 2: Top is Triangle, Bottom is HOLE
  // ---------------------------------------------------------
  else if(!top_is_hole && bot_is_hole)
  {
    auto face_f1 = he_v1_v2.next().twin().face();
    auto he_f1_associate = he_v1_v2.prev();
    auto he_f1_bottom = he_v1_v2.next().twin().next();

    bool top_victim_feature = _is_feature_edge[he_v1_v2.next().index()];
    if(top_victim_feature)
    {
        _is_feature_edge[he_f1_associate.index()] = true;
        _is_feature_edge[he_f1_associate.twin().index()] = true;
    }

    face_top.deactivate();

    // Deactivate top triangle edges (F1)
    he_v1_v2.next().twin().deactivate();
    he_v1_v2.next().deactivate();
    he_v1_v2.deactivate();

    // Deactivate bottom edge (just the line)
    he_v2_v1.deactivate();
    v2_iter.deactivate();

    _is_feature_edge[he_v1_v2.next().twin().index()] = false;
    _is_feature_edge[he_v1_v2.next().index()] = false;
    _is_feature_edge[he_v1_v2.index()] = false;
    _is_feature_edge[he_v2_v1.index()] = false;
    _vertex_feature_type[v2] = SIMPLE; // SIMPLE for safety

    // --- Top Side: Standard Triangle Removal ---
    he_f1_associate.data().next = he_f1_bottom.index();
    he_f1_bottom.data().prev = he_f1_associate.index();
    he_f1_associate.data().prev = he_f1_bottom.next().index();
    he_f1_bottom.next().data().next = he_f1_associate.index();

    he_f1_associate.data().face = face_f1.index();
    face_f1.data().half_edge = he_f1_associate.index();

    // Update u1
    u1_iter.data().half_edge = he_f1_associate.index();

    // --- Bottom Side: Boundary Stitching ---
    auto prev_he = he_v2_v1.prev();
    auto next_he = he_v2_v1.next();

    // Close the loop
    prev_he.data().next = next_he.index();
    next_he.data().prev = prev_he.index();

    // v1 needs a valid pointer. he_f1_associate.twin() is a safe interior edge
    v1_iter.data().half_edge = he_f1_associate.twin().index();
  }

  // ---------------------------------------------------------
  // CASE 3: Top is HOLE, Bottom is Triangle
  // ---------------------------------------------------------
  else if(top_is_hole && !bot_is_hole)
  {
    // ONLY calculate Bottom Side variables
    auto face_f2 = he_v2_v1.prev().twin().face();
    auto he_f2_associate = he_v2_v1.next();
    auto he_f2_top = he_v2_v1.prev().twin().prev();

    bool bot_victim_feature = _is_feature_edge[he_v2_v1.prev().index()];
    if(bot_victim_feature)
    {
        _is_feature_edge[he_f2_associate.index()] = true;
        _is_feature_edge[he_f2_associate.twin().index()] = true;
    }

    face_bottom.deactivate();

    // Deactivate top edge (just the line)
    he_v1_v2.deactivate();

    // Deactivate bottom triangle edges
    he_v2_v1.prev().twin().deactivate();
    he_v2_v1.prev().deactivate();
    he_v2_v1.deactivate();
    v2_iter.deactivate();

    _is_feature_edge[he_v1_v2.index()] = false;
    _is_feature_edge[he_v2_v1.prev().twin().index()] = false;
    _is_feature_edge[he_v2_v1.prev().index()] = false;
    _is_feature_edge[he_v2_v1.index()] = false;
    _vertex_feature_type[v2] = SIMPLE; // SIMPLE for safety

    // --- Top Side: Boundary Stitching ---
    auto prev_he = he_v1_v2.prev();
    auto next_he = he_v1_v2.next();

    // Close the loop
    prev_he.data().next = next_he.index();
    next_he.data().prev = prev_he.index();

    // --- Bottom Side: Standard Triangle Removal ---
    he_f2_associate.data().next = he_f2_top.prev().index();
    he_f2_top.prev().data().prev = he_f2_associate.index();
    he_f2_associate.data().prev = he_f2_top.index();
    he_f2_top.data().next = he_f2_associate.index();

    he_f2_associate.data().face = face_f2.index();
    face_f2.data().half_edge = he_f2_associate.index();

    // Update u2
    u2_iter.data().half_edge = he_f2_associate.twin().index();

    // v1 needs a valid pointer. Top side is boundary, so pointing to 'next_he' is safe
    v1_iter.data().half_edge = next_he.index();
  }

  // ---------------------------------------------------------
  // CASE 4: Both are HOLES (Wire Edge / Strut)
  // ---------------------------------------------------------
  else
  {
    // Just stitch both boundary loops
    he_v1_v2.deactivate();
    he_v2_v1.deactivate();
    v2_iter.deactivate();

    _is_feature_edge[he_v1_v2.index()] = false;
    _is_feature_edge[he_v2_v1.index()] = false;
    _vertex_feature_type[v2] = SIMPLE; // SIMPLE for safety

    // Stitch Top
    auto t_prev = he_v1_v2.prev();
    auto t_next = he_v1_v2.next();
    t_prev.data().next = t_next.index();
    t_next.data().prev = t_prev.index();

    // Stitch Bottom
    auto b_prev = he_v2_v1.prev();
    auto b_next = he_v2_v1.next();
    b_prev.data().next = b_next.index();
    b_next.data().prev = b_prev.index();

    v1_iter.data().half_edge = t_next.index();
  }

  // Move vertex v1 to optimal position
  mesh().vertex_at(v1).data().xyz = new_pos;

  // Inherit feature type (most restrictive)
  if(t1 == CORNER || t2 == CORNER)
    _vertex_feature_type[v1] = CORNER;
  else if(t1 == EDGE || t2 == EDGE)
    _vertex_feature_type[v1] = EDGE;
  else
    _vertex_feature_type[v1] = SIMPLE;

  return true;
}

bool
Mesh_modifier_uniform_remeshing::flip_edge(int he_index)
{
  // HALF-EDGES
  Mesh_connectivity::Half_edge_iterator he0 = mesh().half_edge_at(he_index);
  Mesh_connectivity::Half_edge_iterator he1 = he0.twin();

  // meshes on the boundary are not flippable
  if(he0.face().is_equal(mesh().hole()) || he1.face().is_equal(mesh().hole()))
  {
    return false;
  }

  Mesh_connectivity::Half_edge_iterator he2 = he0.next();
  Mesh_connectivity::Half_edge_iterator he3 = he2.next();
  Mesh_connectivity::Half_edge_iterator he4 = he1.next();
  Mesh_connectivity::Half_edge_iterator he5 = he4.next();

  // VERTICES
  Mesh_connectivity::Vertex_iterator v0 = he1.origin();
  Mesh_connectivity::Vertex_iterator v1 = he0.origin();
  Mesh_connectivity::Vertex_iterator v2 = he3.origin();
  Mesh_connectivity::Vertex_iterator v3 = he5.origin();

  // FACES
  Mesh_connectivity::Face_iterator f0 = he0.face();
  Mesh_connectivity::Face_iterator f1 = he1.face();

  //
  // Now modify the connectivity
  //

  // HALF-EDGES
  he0.data().next = he3.index();
  he0.data().prev = he4.index();
  he0.data().origin = v3.index();
  //
  he1.data().next = he5.index();
  he1.data().prev = he2.index();
  he1.data().origin = v2.index();
  //
  he2.data().next = he1.index();
  he2.data().prev = he5.index();
  he2.data().face = f1.index();
  //
  he3.data().next = he4.index();
  he3.data().prev = he0.index();
  //
  he4.data().next = he0.index();
  he4.data().prev = he3.index();
  he4.data().face = f0.index();
  //
  he5.data().next = he2.index();
  he5.data().prev = he1.index();

  // VERTICES
  v0.data().half_edge = he2.index();
  v1.data().half_edge = he4.index();
  v2.data().half_edge = he1.index();
  v3.data().half_edge = he0.index();

  // FACES
  f0.data().half_edge = he0.index();
  f1.data().half_edge = he1.index();

  // operation successful
  return true;
}

// ============================================================
// Helpers & Topology Checks
// ============================================================


double
Mesh_modifier_uniform_remeshing::get_edge_length(int he_index) const
{
  Mesh_connectivity::Half_edge_iterator he = _m.half_edge_at(he_index);
  Mesh_connectivity::Vertex_iterator v0 = he.origin();
  Mesh_connectivity::Vertex_iterator v1 = he.dest();

  Eigen::Vector3d p0 = v0.xyz();
  Eigen::Vector3d p1 = v1.xyz();

  return (p1 - p0).norm();
}

int
Mesh_modifier_uniform_remeshing::get_vertex_valence(int v_index) const
{
  return analysis::vertex_valence(_m, v_index);
}

int
Mesh_modifier_uniform_remeshing::get_valence_deviation(int v_index) const
{
  int valence = get_vertex_valence(v_index);
  bool is_boundary = analysis::vertex_is_boundary(_m, v_index);

  int ideal_valence = is_boundary ? BOUNDARY_VALENCE : INTERIOR_VALENCE;
  return std::abs(valence - ideal_valence);
}

bool
Mesh_modifier_uniform_remeshing::should_flip_edge(int he_index)
{
  auto he = mesh().half_edge_at(he_index);
  auto twin = he.twin();

  if(_is_feature_edge[he.index()])
  {
    return false; // Dont flip feature edges
  }

  int v1 = he.origin().index();
  int v2 = twin.origin().index();
  int o1 = he.next().dest().index();
  int o2 = twin.next().dest().index();

  auto get_vertex_target = [&](int vid)
  { return analysis::vertex_is_boundary(_m, vid) ? BOUNDARY_VALENCE : INTERIOR_VALENCE; };

  int v1_val = get_vertex_valence(v1);
  int v2_val = get_vertex_valence(v2);
  int o1_val = get_vertex_valence(o1);
  int o2_val = get_vertex_valence(o2);

  int v1_tar = get_vertex_target(v1);
  int v2_tar = get_vertex_target(v2);
  int o1_tar = get_vertex_target(o1);
  int o2_tar = get_vertex_target(o2);

  auto abs_err = [](int val, int target) { return std::abs(val - target); };

  int current_deviation =
      abs_err(v1_val, v1_tar) + abs_err(v2_val, v2_tar) + abs_err(o1_val, o1_tar) + abs_err(o2_val, o2_tar);

  int new_deviation = abs_err(v1_val - 1, v1_tar) + abs_err(v2_val - 1, v2_tar) + abs_err(o1_val + 1, o1_tar) +
                      abs_err(o2_val + 1, o2_tar);

  // Only flip if it improves total deviation
  return new_deviation < current_deviation;
}

bool
Mesh_modifier_uniform_remeshing::is_legal_flip(int he_index)
{
  auto he = mesh().half_edge_at(he_index);
  auto twin = he.twin();
  int v1 = he.origin().index();
  int v2 = twin.origin().index();

  // Boundary edge is no flip
  if(he.face().is_equal(mesh().hole()) || twin.face().is_equal(mesh().hole()))
    return false;

  // o1, o2 are the opposite vertices
  int o1 = he.next().dest().index();
  int o2 = twin.next().dest().index();

  // Check for degenerate case where o1 == o2
  if(o1 == o2)
    return false;

  // Check if sides already connected
  if(get_halfedge_between_vertices(o1, o2) != -1)
    return false;

  // ==========================================
  // GEOMETRIC CHECKS (Convexity & Self-Intersection)
  // ==========================================

  // Compare change in normals before and after flip, threshold by dot product

  auto p_v1 = mesh().vertex_at(v1).xyz();
  auto p_v2 = mesh().vertex_at(v2).xyz();
  auto p_o1 = mesh().vertex_at(o1).xyz();
  auto p_o2 = mesh().vertex_at(o2).xyz();

  auto n_curr_1 = (p_v2 - p_v1).cross(p_o1 - p_v1).normalized();
  auto n_curr_2 = (p_v1 - p_v2).cross(p_o2 - p_v2).normalized();

  auto n_new_A = (p_v1 - p_o1).cross(p_o2 - p_o1).normalized();
  auto n_new_B = (p_v2 - p_o2).cross(p_o1 - p_o2).normalized();

  double min_dot = EDGE_FLIP_THRESHOLD_DOT;

  if(n_new_A.dot(n_curr_1) < min_dot)
    return false;
  if(n_new_B.dot(n_curr_2) < min_dot)
    return false;

  return true;
}


bool
Mesh_modifier_uniform_remeshing::is_legal_collapse(int v1, int v2, Eigen::Vector3d new_pos, double threshold)
{

  // Check for boundary conditions
  bool v1_boundary = analysis::vertex_is_boundary(_m, v1);
  bool v2_boundary = analysis::vertex_is_boundary(_m, v2);
  auto he = mesh().half_edge_at(get_halfedge_between_vertices(v1, v2));
  auto twin = he.twin();

  bool he_face_is_hole = he.face().is_equal(mesh().hole());
  bool tw_face_is_hole = twin.face().is_equal(mesh().hole());
  bool edge_is_boundary = (he_face_is_hole || tw_face_is_hole);

  // Boundary edge rule, do not collapse interior edge connecting two boundary vertices
  if(v1_boundary && v2_boundary && !edge_is_boundary)
    return false;

  // ==========================================
  // 1. FEATURE CHECKS
  // ==========================================
  bool is_feature_edge = _is_feature_edge[he.index()];
  VertexFeatureType v1_type = _vertex_feature_type[v1];
  VertexFeatureType v2_type = _vertex_feature_type[v2];

  if(is_feature_edge)
  {
    // No corners allowed for feature edge collapse
    if(v1_type == CORNER || v2_type == CORNER)
      return false;
  }
  else
  {
    // Stop collapse of parallel feature lines into a corner
    if(v1_type != SIMPLE && v2_type != SIMPLE)
      return false;
  }
  // At most one vertex can be a feature vertex or both are edge vertices and midpoint is okay now

  // Gather neighbors (excluding each other)
  std::set<int> n1 = get_all_neighbors_from_vertex(v1);
  std::set<int> n2 = get_all_neighbors_from_vertex(v2);
  n1.erase(v2);
  n2.erase(v1);

  // Opposite vertices across edge (the two triangle tips)
  std::set<int> allowed;
  if(!he_face_is_hole)
    allowed.insert(he.next().dest().index());
  if(!tw_face_is_hole)
    allowed.insert(twin.next().dest().index());


  // Common neighbors
  std::set<int> common;
  std::set_intersection(n1.begin(), n1.end(), n2.begin(), n2.end(), std::inserter(common, common.begin()));

  // Fail Condition A: the only common neighbors are u1 and u2 (link condition)
  std::set<int> common_minus_allowed;
  std::set_difference(common.begin(),
      common.end(),
      allowed.begin(),
      allowed.end(),
      std::inserter(common_minus_allowed, common_minus_allowed.begin()));
  if(!common_minus_allowed.empty())
    return false;

  // Condition B: there exists at least one neighbor NOT in common (tetrahedron check)
  std::set<int> symdiff;
  std::set_symmetric_difference(n1.begin(), n1.end(), n2.begin(), n2.end(), std::inserter(symdiff, symdiff.begin()));
  if(symdiff.empty())
    return false;

  // ==========================================
  // 2. GEOMETRIC CHECKS (Botsch et al. 2010 Requirements)
  // ==========================================

  // --- Requirement A: Edge Length Check, reject creating new long size edges ---

  // Combine all unique neighbors
  std::set<int> all_neighbors = n1;
  all_neighbors.insert(n2.begin(), n2.end());

  // Avoid sqrt by comparing squared lengths
  double threshold_sq = threshold * threshold;

  for(int nid : all_neighbors)
  {
    auto p_neighbor = mesh().vertex_at(nid).xyz();
    if((p_neighbor - new_pos).squaredNorm() > threshold_sq)
    {
      return false; // Resulting edge would be too long
    }
  }

  // --- Requirement B: Normal Flip / Intersection Check ---
  // Check face normals of the one-ring of v1 and v2 before and after the hypothetical collapse.
  
  // Lambda to check one vertex's neighbors
  auto simulate_normal_flip = [&](int v_idx, const Eigen::Vector3d& new_pos) -> bool {
      auto ring = _m.vertex_ring_at(v_idx);
      do {
          auto he = ring.half_edge();
          auto face = he.face();
          auto f1 = he.face();
          auto f2 = he.twin().face();
          
          // Skip the faces adjacent to the edge being collapsed (they will disappear)
          if (face.is_equal(f1) || face.is_equal(f2)) {
            continue;
          }

          // Get vertices of the triangle
          int v_curr = he.origin().index();
          int v_next = he.next().origin().index();
          int v_prev = he.prev().origin().index();

          // Compute Old Normal
          Eigen::Vector3d pA = mesh().vertex_at(v_curr).xyz();
          Eigen::Vector3d pB = mesh().vertex_at(v_next).xyz();
          Eigen::Vector3d pC = mesh().vertex_at(v_prev).xyz();
          Eigen::Vector3d n_old = (pB - pA).cross(pC - pA).normalized();

          // Compute New Normal: Replace v1 or v2 with p_new
          if (v_curr == v1 || v_curr == v2) pA = new_pos;
          if (v_next == v1 || v_next == v2) pB = new_pos;
          if (v_prev == v1 || v_prev == v2) pC = new_pos;

          Eigen::Vector3d n_new = (pB - pA).cross(pC - pA);
          
          // Check for degenerate triangle
          if(n_new.norm() < 1e-12) return false; 
          n_new.normalize();

          // Botsch et al. suggest checking that the normal doesn't flip more than 90 degrees (dot < 0)
          // A safer threshold is usually 0.2 to keep triangles somewhat regular.
          if (n_old.dot(n_new) < 0.2) { 
              return false; 
          }

      } while (ring.advance());
      return true;
  };

  if (!simulate_normal_flip(v1, new_pos)) return false;
  if (!simulate_normal_flip(v2, new_pos)) return false;

  return true;
}

void
Mesh_modifier_uniform_remeshing::relabel_vertex(int old_id, int new_id)
{
  auto he = mesh().vertex_at(old_id).half_edge();
  const int he_start_index = he.index();
  do
  {
    assert(he.origin().index() == old_id);
    he.data().origin = new_id;
    he = he.twin().next();
  } while(he.index() != he_start_index);
}

std::set<int>
Mesh_modifier_uniform_remeshing::get_all_neighbors_from_vertex(int vertex_id) const
{
  std::set<int> neighbors;
  auto ring_iter = _m.vertex_ring_at(vertex_id);
  do
  {
    // he always points to vertex, so origin is the neighbor
    neighbors.insert(ring_iter.half_edge().origin().index());
  } while(ring_iter.advance());
  return neighbors;
}

int
Mesh_modifier_uniform_remeshing::get_halfedge_between_vertices(const int v0, const int v1)
{
  // Traverse the one-ring of v0 to find the half-edge pointing to v1
  Mesh_connectivity::Vertex_ring_iterator ring = _m.vertex_ring_at(v0);
  do
  {
    Mesh_connectivity::Half_edge_iterator he = ring.half_edge();
    if(he.origin().index() == v1)
    {
      // Found half-edge from v1 to v0, return its twin (v0 to v1)
      return he.twin().index();
    }
  } while(ring.advance());

  return -1; // No such half-edge found
}

std::vector<int>
Mesh_modifier_uniform_remeshing::get_one_ring_neighbors(int v_index) const
{
  std::vector<int> neighbors;
  Mesh_connectivity::Vertex_ring_iterator ring = _m.vertex_ring_at(v_index);

  do
  {
    Mesh_connectivity::Half_edge_iterator he = ring.half_edge();
    neighbors.push_back(he.origin().index());
  } while(ring.advance());

  return neighbors;
}

Eigen::Vector3d
Mesh_modifier_uniform_remeshing::calculate_normal(const Eigen::Vector3d & p0,
    const Eigen::Vector3d & p1,
    const Eigen::Vector3d & p2) const
{
  Eigen::Vector3d u = p1 - p0;
  Eigen::Vector3d v = p2 - p0;
  Eigen::Vector3d n = u.cross(v);
  return n.normalized();
}

Mesh_modifier_uniform_remeshing::GeometryCache
Mesh_modifier_uniform_remeshing::compute_geometry_cache()
{
  GeometryCache cache;
  const int n_vertices = mesh().n_total_vertices();
  const int n_faces = mesh().n_total_faces();
  const int n_half_edges = mesh().n_total_half_edges();
  // Preset all values to zero with memory alloc
  cache.vertex_normals.resize(n_vertices, Eigen::Vector3d(0.0, 0.0, 0.0));
  cache.face_areas.resize(n_faces, 0.0);
  cache.face_barycenters.resize(n_faces, Eigen::Vector3d(0.0, 0.0, 0.0));
  cache.half_edge_angles.resize(n_half_edges, 0.0);

  for(int f = 0; f < n_faces; ++f)
  {
    auto face = mesh().face_at(f);
    if(!face.is_active() || face.is_equal(mesh().hole()))
      continue;

    // Get the three half edges and vertices of the triangle
    auto he1 = face.half_edge();
    auto he2 = he1.next();
    auto he3 = he2.next();
    auto v1 = he1.origin();
    auto v2 = he2.origin();
    auto v3 = he3.origin();
    // Cycle v1 -> v2 -> v3 -> v1 half edges: he1, he2, he3

    // Compute face normal
    Eigen::Vector3d p1 = v1.xyz();
    Eigen::Vector3d p2 = v2.xyz();
    Eigen::Vector3d p3 = v3.xyz();

    Eigen::Vector3d u = p2 - p1;
    Eigen::Vector3d v = p3 - p1;
    Eigen::Vector3d cross_prod = u.cross(v);

    // Compute face area (0.5 * |cross product|)
    double double_area = cross_prod.norm();
    cache.face_areas[f] = 0.5 * double_area;

    Eigen::Vector3d face_normal;
    if(double_area > 1e-12)
    {
      face_normal = cross_prod / double_area; // Normalized face normal
    }
    else
    {
      face_normal = Eigen::Vector3d(0, 0, 1); // Degenerate fallback
    }

    // Compute face barycenter
    cache.face_barycenters[f] = (p1 + p2 + p3) / 3.0;

    // Angle computations
    Eigen::Vector3d e12 = (p2 - p1).normalized(); // v1 -> v2
    Eigen::Vector3d e13 = (p3 - p1).normalized(); // v1 -> v3
    Eigen::Vector3d e23 = (p3 - p2).normalized(); // v2 -> v3
    Eigen::Vector3d e21 = -e12; // v2 -> v1
    Eigen::Vector3d e31 = -e13; // v3 -> v1
    Eigen::Vector3d e32 = -e23; // v3 -> v2

    double dot1 = std::max(-1.0, std::min(1.0, e12.dot(e13)));
    double angle1 = std::acos(dot1);
    double dot2 = std::max(-1.0, std::min(1.0, e23.dot(e21)));
    double angle2 = std::acos(dot2);
    double dot3 = std::max(-1.0, std::min(1.0, e31.dot(e32)));
    double angle3 = std::acos(dot3);

    // Store half-edge angles
    cache.half_edge_angles[he1.index()] = angle1;
    cache.half_edge_angles[he2.index()] = angle2;
    cache.half_edge_angles[he3.index()] = angle3;

    // Add face normal contribution to each vertex normal, weight by angle
    cache.vertex_normals[v1.index()] += face_normal * angle1;
    cache.vertex_normals[v2.index()] += face_normal * angle2;
    cache.vertex_normals[v3.index()] += face_normal * angle3;
  }

  // Normalize vertex normals
  for(int v = 0; v < n_vertices; ++v)
  {
    if(mesh().vertex_at(v).is_active())
    {
      cache.vertex_normals[v].normalize();
    }
  }

  return cache;
}

void
Mesh_modifier_uniform_remeshing::compute_face_min_angles()
{
  const int n_faces = mesh().n_total_faces();
  _face_min_angles.resize(n_faces, 0.0);

  for(int f = 0; f < n_faces; ++f)
  {
    auto face = mesh().face_at(f);
    if(!face.is_active() || face.is_equal(mesh().hole()))
    {
      _face_min_angles[f] = 0.0;
      continue;
    }

    // Get the three half edges and vertices of the triangle
    auto he1 = face.half_edge();
    auto he2 = he1.next();
    auto he3 = he2.next();
    auto v1 = he1.origin();
    auto v2 = he2.origin();
    auto v3 = he3.origin();

    // Get vertex positions
    Eigen::Vector3d p1 = v1.xyz();
    Eigen::Vector3d p2 = v2.xyz();
    Eigen::Vector3d p3 = v3.xyz();

    // Compute angles at each vertex
    // Angle at v1 (between edges v1->v2 and v1->v3)
    Eigen::Vector3d e12 = (p2 - p1).normalized();
    Eigen::Vector3d e13 = (p3 - p1).normalized();
    double dot1 = std::max(-1.0, std::min(1.0, e12.dot(e13)));
    double angle1 = std::acos(dot1);

    // Angle at v2 (between edges v2->v3 and v2->v1)
    Eigen::Vector3d e23 = (p3 - p2).normalized();
    Eigen::Vector3d e21 = -e12;
    double dot2 = std::max(-1.0, std::min(1.0, e23.dot(e21)));
    double angle2 = std::acos(dot2);

    // Angle at v3 (between edges v3->v1 and v3->v2)
    Eigen::Vector3d e31 = -e13;
    Eigen::Vector3d e32 = -e23;
    double dot3 = std::max(-1.0, std::min(1.0, e31.dot(e32)));
    double angle3 = std::acos(dot3);

    // Store the minimum angle for this face
    _face_min_angles[f] = std::min({angle1, angle2, angle3});
  }
}

void
Mesh_modifier_uniform_remeshing::mark_feature_edges_and_vertices()
{
  // Precompute all the face normals
  std::vector<Eigen::Vector3d> face_normals(mesh().n_total_faces(), Eigen::Vector3d(0.0, 0.0, 0.0));
  for(int f = 0; f < mesh().n_total_faces(); ++f)
  {
    auto face = mesh().face_at(f);
    if(!face.is_active() || face.is_equal(mesh().hole()))
      continue;

    // Get the three half edges and vertices of the triangle
    auto he1 = face.half_edge();
    auto he2 = he1.next();
    auto he3 = he2.next();
    auto v1 = he1.origin();
    auto v2 = he2.origin();
    auto v3 = he3.origin();
    // Cycle v1 -> v2 -> v3 -> v1 half edges: he1, he2, he3

    // Compute face normal
    Eigen::Vector3d p1 = v1.xyz();
    Eigen::Vector3d p2 = v2.xyz();
    Eigen::Vector3d p3 = v3.xyz();

    Eigen::Vector3d u = p2 - p1;
    Eigen::Vector3d v = p3 - p1;
    Eigen::Vector3d cross_prod = u.cross(v);

    double area_double = cross_prod.norm();
    if(area_double > 1e-12)
    {
      face_normals[f] = cross_prod / area_double; // Normalized face normal
    }
    else
    {
      face_normals[f] = Eigen::Vector3d(0, 0, 1); // Degenerate fallback
    }
  }


  const int n_half_edges = mesh().n_total_half_edges();
  const int n_vertices = mesh().n_total_vertices();
  const double cos_feature_angle_rad = std::cos(FEATURE_ANGLE_DEGREES * (M_PI / 180.0));

  printf("Marking feature edges with angle threshold: %f degrees (cosine: %f)\n",
      FEATURE_ANGLE_DEGREES,
      cos_feature_angle_rad);
  printf("Total half-edges: %d, vertices: %d\n", n_half_edges, n_vertices);

  _is_feature_edge.assign(n_half_edges, false);
  _vertex_feature_type.assign(n_vertices, SIMPLE);

  // Mark feature edges
  for(int he_idx = 0; he_idx < n_half_edges; ++he_idx)
  {
    auto he = mesh().half_edge_at(he_idx);
    if(!he.is_active())
      continue; // Skip inactive half-edges

    auto tw_index = he.twin().index();
    if(he_idx > tw_index)
      continue; // Process each edge only once

    auto f1 = he.face().index();
    auto f2 = he.twin().face().index();

    bool is_boundary = (f1 == mesh().hole().index() || f2 == mesh().hole().index());

    if(is_boundary)
    {
      _is_feature_edge[he_idx] = true; // Boundary edge is feature edge
      _is_feature_edge[tw_index] = true;
    }
    else
    {
      // Both faces are valid, check angle between normals
      double cos_beta = face_normals[f1].dot(face_normals[f2]);
      if(cos_beta < cos_feature_angle_rad)
      {
        // Mark as feature edge if angle exceeds threshold
        _is_feature_edge[he_idx] = true;
        _is_feature_edge[tw_index] = true;
      }
    }
  }

  // Classify vertices based on incident feature edges
  for(int v_idx = 0; v_idx < n_vertices; ++v_idx)
  {
    if(!mesh().vertex_at(v_idx).is_active())
      continue; // Skip inactive vertices
    _vertex_feature_type[v_idx] = classify_vertex_feature_type(v_idx);
  }
  return;
}

Mesh_modifier_uniform_remeshing::VertexFeatureType
Mesh_modifier_uniform_remeshing::classify_vertex_feature_type(int v_idx)
{
  // Count incoming feature edges
  auto ring = _m.vertex_ring_at(v_idx);
  int feature_edge_count = 0;
  do
  {
    auto he = ring.half_edge();
    if(_is_feature_edge[he.index()])
      feature_edge_count++;
  } while(ring.advance());

  if(feature_edge_count == 0)
  {
    return SIMPLE;
  }
  else if(feature_edge_count == 1)
  {
    return CORNER;
  }
  else if(feature_edge_count == 2)
  {
    return EDGE;
  }
  else
  {
    return CORNER;
  }
}

} // end of mohecore
} // end of minimesh
