#pragma once

//
// mesh_modifier_edge_collapse.hpp
//
// Edge collapse mesh modifier.
//
// Author: Shayan Hoshyari
//


#include <string>
#include <unordered_set>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <minimesh/core/mohe/mesh_connectivity.hpp>

namespace minimesh
{
namespace mohecore
{


class Mesh_modifier_arap
{
public:
  struct Neighbor
  {
    int index; // Neighbor vertex index
    double weight; // Cotangent weight w_ij
  };

  struct Anchors
  {
    // --- public API ---
    void init(int n)
    {
      _n = n;
      _is_static.assign(n, 0);
      _static_list.clear();
    }

    // mark/unmark static anchors (bumps version only on real change)
    bool addStatic(int v)
    {
      if(!_valid(v))
        return false;
      if(_is_static[v])
        return false; // already static

      // If here, we are adding a new static anchor
      _is_static[v] = true;
      _static_list.push_back(v);
      _version++;
      return true;
    }

    bool removeStatic(int v)
    {
      if(!_valid(v))
        return false;
      if(!_is_static[v])
        return false; // not static

      // If here, we are removing a static anchor
      _is_static[v] = false;
      _static_list.erase(std::remove(_static_list.begin(), _static_list.end(), v), _static_list.end());
      _version++;
      return true;
    }

    void clearStatic()
    {
      for(int v : _static_list)
      {
        _is_static[v] = false;
      }
      _static_list.clear();
      _version++;
    }

    bool isStaticAnchor(int v) const
    {
      if(!_valid(v))
        return false;
      return _is_static[v];
    }

    const std::vector<int> & getStaticList() const { return _static_list; }

    // set the current temporary (moving) anchor id, return bool if changed
    bool setTemp(int v)
    {
      if(_temp_now == v)
      {
        return false;
      }
      else
      {
        _temp_now = v;
        // Check if total set of anchors changed
        if(!isStaticAnchor(v))
        {
          _version++;
        }
        return true;
      }
    }
    bool hasTemp() const { return _temp_now != -1; }
    int getTemp() const { return _temp_now; }
    bool removeTemp()
    {
      if(_temp_now == -1)
      {
        return false;
      }
      else
      {
        _version++;
        _temp_now = -1;
        _temp_pos = Eigen::Vector3d::Zero();
        return true;
      }
    }
    Eigen::Vector3d getTempPosition() const { return _temp_pos; }
    void setTempPosition(const Eigen::Vector3d & pos) { _temp_pos = pos; }


    int version() const { return _version; }

  private:
    bool _valid(int v) const { return v >= 0 && v < _n; }
    int _n = 0;
    std::vector<bool> _is_static; // Check if vertex is static anchor
    std::vector<int> _static_list; // List of static anchor indices
    int _temp_now = -1; // Current temporary (moving) anchor index
    Eigen::Vector3d _temp_pos; // Current temporary (moving) anchor position

    uint64_t _version = 0; // b std::vector<int> _pos_in_free_base; // vertex -> index in _free_base (or -1)
  };


  // Trivial constructor
  Mesh_modifier_arap(Mesh_connectivity & mesh_in)
  : _m(mesh_in)
  {
    _anchors.init(mesh_in.n_total_vertices());
  }

  // Get the underlying mesh
  Mesh_connectivity & mesh() { return _m; }
  const Mesh_connectivity & mesh() const { return _m; }

  //
  // Given two vertices, this function return the index of the half-edge going from v0 to v1.
  // Returns mesh::invalid_index if no half-edge exists between the two vertices.
  //
  int get_halfedge_between_vertices(const int v0, const int v1);

  //
  // Flip an edge in a mesh
  // Input: The mesh, and the index of a half-edge on the edge we wish to flip
  // Return true if operation successful, and false if operation not possible
  //
  // Assumption: mesh is all triangles
  //
  // NOTE: To see how this method works, take a look at edge-flip.svg
  //
  bool flip_edge(const int he_index);

  //
  // ARAP Deformation Methods (stubs for now)
  //

  // Initialize ARAP structures (call this after mesh is loaded)
  void initialize();

  // Build adjacency list with cotangent weights for all vertices
  void build_adjacency_list();

  // Build the Laplacian matrix for the mesh
  void build_laplacian_matrix();

  // Build the free/free and free/constraint blocks of the Laplacian matrix from the constraint indices
  void build_blocks_from_constraints();

  void build_rhs_bf_direct();


  // Add a vertex as an anchor point
  // Returns true if successful, false if vertex is already an anchor
  bool add_anchor(const int vertex_index);

  // Remove a vertex from anchors
  // Returns true if successful, false if vertex wasn't an anchor
  bool remove_anchor(const int vertex_index);

  // Check if a vertex is an anchor
  bool is_anchor(const int vertex_index) const;

  // Clear all anchors
  void clear_anchors();

  // Get all anchor vertex indices (returns vector for convenience)
  std::vector<int> get_static_anchors();

  // Perform ARAP deformation by moving an anchor (in-place version for GUI performance)
  // vertex_index: the anchor to move
  // new_position: the new position for the anchor
  // deformed_positions: output matrix (3 x n_vertices) with deformed positions
  // Returns true if deformation was successful
  // NOTE: This is a stub - actual ARAP implementation to be added later
  //       For now, just copies original positions and updates the anchor
  bool deform_with_temp_anchor(const int vertex_index,
      const Eigen::Vector3d & new_position,
      Eigen::Matrix3Xd & deformed_positions);

  // Compute ARAP deformation and return new positions (convenience version for testing)
  // vertex_index: the anchor to move
  // new_position: the new position for the anchor
  // Returns: Matrix (3 x n_vertices) with deformed positions
  // Throws: std::runtime_error if deformation fails
  Eigen::Matrix3Xd compute_deformation(const int vertex_index, const Eigen::Vector3d & new_position);

  // Apply ARAP deformation directly to mesh vertex positions (for CLI & file output)
  // vertex_index: the anchor to move
  // new_position: the new position for the anchor
  // Returns true if deformation was successful
  // NOTE: This modifies the actual mesh vertex positions in-place
  bool apply_deformation_to_mesh(const int vertex_index, const Eigen::Vector3d & new_position);


private:
  // pointer to the mesh that we are working on.
  Mesh_connectivity & _m;

  Anchors _anchors; // Anchor management
  int _last_anchor_version = -1; // Last known version of anchors

  // Adjacency list with cotangent weights
  // adj[i] holds (neighbor_index, cotangent_weight) for each neighbor of vertex i
  std::vector<std::vector<Neighbor>> _adj;

  Eigen::Matrix3Xd _vertices_rest; // Rest positions of vertices
  Eigen::Matrix3Xd _vertices_deformed; // Deformed positions of vertices

  Eigen::SparseMatrix<double> _L; // Laplacian matrix across all params
  Eigen::SparseMatrix<double> _Lff; // Free-free Laplacian matrix
  Eigen::SparseMatrix<double> _Lfc; // Free-anchor Laplacian matrix

  std::vector<Eigen::Matrix3d> _R; // Rotation matrices for each vertex
  Eigen::MatrixXd _Bf; // Reduced right-hand side matrix for free vertices
  Eigen::MatrixXd _C; // Constraint matrix
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> _solver;

  // Freed index vertices, constrained index vertices, and previous constrained index vertices
  std::vector<int> _free_idx, _cons_idx, _last_cons_idx;


  // Core ARAP solver - all public deformation methods delegate to this
  // vertex_index: the temporary anchor to move
  // new_position: the new position for the temporary anchor
  // output: output matrix (3 x n_vertices) with deformed positions
  // Returns true if deformation was successful
  bool _solve_arap(const int temp_anchor, const Eigen::Vector3d & pulled, Eigen::Matrix3Xd & output);

  Eigen::MatrixXd _gather_constraint_matrix();

  // Perform an iteration to update the internal deformation state using Lp = B and rotation estimation
  bool _solve_deformation_with_anchors();

  // Perform an iteration of R_i estimation
  bool _solve_rotation_estimation();
};

} // end of mohecore
} // end of minimesh
