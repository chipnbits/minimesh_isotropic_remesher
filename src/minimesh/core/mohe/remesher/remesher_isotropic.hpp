// Implementation of Uniform Explicit Isotropic Remeshing
// Based on Botsch and Kobbelt (2004)
//
// Project Part 1: Baseline remeshing algorithm
//

#pragma once

#include <Eigen/Core>
#include <minimesh/core/mohe/mesh_connectivity.hpp>
#include <set>
#include <vector>

namespace minimesh
{
namespace mohecore
{

class Mesh_modifier_uniform_remeshing
{
public:
  Mesh_modifier_uniform_remeshing(Mesh_connectivity & mesh_in)
  : _m(mesh_in)
  {
    // run initialization
    initialize();
  }

  enum class SmoothingType
  {
    Uniform, // Simple barycenter of neighbors (Isotropic)
    Barycenters // Area-weighted or CVT (Centroidal Voronoi) approximation
  };

  // Figure 3.14 https://nccastaff.bournemouth.ac.uk/jmacey/MastersProject/MSc15/08Tanja/
  enum VertexFeatureType
  {
    SIMPLE = 0,
    EDGE = 1, // Two incident feature edges
    CORNER = 2 // 1 or more than 2 incident feature edges
  };

  struct GeometryCache
  {
    // Indexed by vertex index
    std::vector<Eigen::Vector3d> vertex_normals;
    // Indexed by face index
    std::vector<Eigen::Vector3d> face_barycenters;
    // Indexed by face index
    std::vector<double> face_areas;
    // Indexed by half-edge index originating from corner
    std::vector<double> half_edge_angles;
  };

  // Get the underlying mesh
  Mesh_connectivity & mesh() { return _m; }
  const Mesh_connectivity & mesh() const { return _m; }

  // Initialize the remesher - marks feature edges and vertices
  void initialize();

  // Getters for feature lists
  const std::vector<bool>& get_feature_edges() const { return _is_feature_edge; }
  const std::vector<VertexFeatureType>& get_vertex_feature_types() const { return _vertex_feature_type; }

  //
  // Main Entry Point
  // Executes the iterative remeshing algorithm described by
  //
  // Params:
  //   target_edge_length: The uniform global target length (L)
  //   iterations: Number of times to run the split/collapse/flip/smooth loop
  //
  void remesh(double target_edge_length, int iterations = 5);
  void run_single_pass(double target_edge_length, int tangential_smoothing_iters);

public:
  // pointer to the mesh that we are working on.
  Mesh_connectivity & _m;

  const int INTERIOR_VALENCE = 6;
  const int BOUNDARY_VALENCE = 4;
  const double EDGE_FLIP_THRESHOLD_DOT = 0.7; // Cosine of angle threshold for normal deviation check
  const double LAMBDA_SMOOTHING_DAMPING = 0.7; // Damping factor for vertex smoothing
  const int N_SMOOTHING_ITERS = 5; // Number of smoothing iterations

  // Feature edges
  const double FEATURE_ANGLE_DEGREES = 45.0; // Degrees
  const double FEATURE_ANGLE_COSINE = cos(FEATURE_ANGLE_DEGREES * M_PI / 180.0);

  std::vector<bool> _is_feature_edge;
  std::vector<VertexFeatureType> _vertex_feature_type;
  void ensure_property_containers_size() {
    const std::size_t n_vertices =
      static_cast<std::size_t>(mesh().n_total_vertices());
    const std::size_t n_half_edges =
      static_cast<std::size_t>(mesh().n_total_half_edges());
    if (n_vertices > _vertex_feature_type.size()) {
      _vertex_feature_type.resize(mesh().n_total_vertices(), SIMPLE);
    }
    if (n_half_edges > _is_feature_edge.size()) {
      _is_feature_edge.resize(mesh().n_total_half_edges(), false);
    }
  }

  // ============================================================
  // Core Remeshing Operations
  // ============================================================

  //
  // 1. Edge Splitting
  // Iterates over all edges. If edge_length > target_length*alpha the edge is split at its midpoint. Lit value is 4/3.
  void split_long_edges(double target_length, double alpha = 4.0 / 3.0);

  //
  // 2. Edge Collapsing
  // Iterates over all edges. If edge_length < beta * target_length, the edge is collapsed. Lit value is 4/5.
  void collapse_short_edges(double target_length, double beta = 4.0 / 5.0);

  //
  // 3. Edge Flipping
  // Flips edges to improve vertex valence towards the ideal value 4-boundary, 6-internal.
  void flip_edges_to_optimize_valence();

  //
  // 4. Vertex Smoothing
  // Applies tangential smoothing to redistribute vertices for better mesh quality.
  // Usually involves moving vertex toward the centroid of its neighbors,
  // then projecting back to the tangent plane (or original surface if available).
  //
  void tangential_smoothing(int smoothing_iters = 1, SmoothingType type = SmoothingType::Barycenters);

  // ============================================================
  // Low-Level Mesh Modifiers
  // ============================================================

  // Splits a specific edge (given by half-edge index) at a weighted position along the edge.
  // weight = 0.0 gives the origin vertex, weight = 1.0 gives the destination vertex, weight = 0.5 gives the midpoint.
  // Returns the index of the new vertex, or -1 on failure.
  int split_edge(int he_index, double weight = 0.5);

  // Collapses the edge defined by he_index (collapses origin into destination).
  // Returns true if successful.
  bool collapse_edge(int he_index, double threshold);

  // Flips the edge shared by the two triangles adjacent to he_index.
  // Returns true if successful.
  bool flip_edge(int he_index);
  // ============================================================
  // Helpers & Topology Checks
  // ============================================================

  // Returns the Euclidean length of the edge associated with the half-edge
  double get_edge_length(int he_index) const;

  // Calculates the valence (degree) of a vertex
  int get_vertex_valence(int v_index) const;

  // Calculates the deviation of valence from ideal (6 for internal, 4 for boundary)
  // Used by flip_edges_to_optimize_valence
  int get_valence_deviation(int v_index) const;

  // Determines if flipping the edge would improve valence deviation
  bool should_flip_edge(int he_index);

  // Checks edge collapse legality according to Botsch & Kobbelt criteri and length threshold
  bool is_legal_collapse(int v1, int v2, double threshold);

  // Check edge flip legality according to Botsch & Kobbelt criteria
  bool is_legal_flip(int he_index);


  // Helper to get all neighbors of a vertex as a set (for collapse legality checks)
  std::set<int> get_all_neighbors_from_vertex(int v_index) const;

  // Helper to relabel all half-edges originating from old_id to new_id
  void relabel_vertex(int old_id, int new_id);

  // Helper to find specific half-edge between vertices (reused from your sample)
  int get_halfedge_between_vertices(const int v0, const int v1);

  // Helper to get all neighbors of a vertex (for smoothing and checks)
  std::vector<int> get_one_ring_neighbors(int v_index) const;

  // Compute a surface normal given three points in ccw cycle order
  Eigen::Vector3d calculate_normal(const Eigen::Vector3d & p0,
      const Eigen::Vector3d & p1,
      const Eigen::Vector3d & p2) const;

  // Precompute geometry information for the current mesh state
  GeometryCache compute_geometry_cache();

  // Classify vertices based on incident feature edges, assumes feature edges are marked
  VertexFeatureType classify_vertex_feature_type(int v_index);

  // initialize feature edges and vertices
  void mark_feature_edges_and_vertices();
};
} // end of mohecore
} // end of minimesh