#include <Eigen/Geometry>
#include <algorithm>
#include <minimesh/core/mohe/mesh_modifier_arap.hpp>
#include <minimesh/core/util/assert.hpp>

namespace minimesh
{
namespace mohecore
{


//
// Given two vertices, this function return the index of the half-edge going from v0 to v1.
// Returns -1 if no half-edge exists between the two vertices.
//
int
Mesh_modifier_arap::get_halfedge_between_vertices(const int v0, const int v1)
{
  // Get a ring iterator for v0
  Mesh_connectivity::Vertex_ring_iterator ring_iter = mesh().vertex_ring_at(v0);

  int answer = mesh().invalid_index;

  // Loop over all half-edges that end at v0.
  do
  {
    // Make sure that the half-edge does end at v0
    assert(ring_iter.half_edge().dest().index() == v0);

    // If the half-edge also starts and v1, then it's twin
    // goes from v0 to v1. This would be the half-edge that
    // we were looking for
    if(ring_iter.half_edge().origin().index() == v1)
    {
      answer = ring_iter.half_edge().twin().index();
    }
  } while(ring_iter.advance());

  if(answer != mesh().invalid_index)
  {
    assert(mesh().half_edge_at(answer).origin().index() == v0);
    assert(mesh().half_edge_at(answer).dest().index() == v1);
  }

  return answer;
}

// ARAP deformation methods

void
Mesh_modifier_arap::initialize()
{
  // Defragment the mesh first to ensure compact, sequential vertex indices
  Mesh_connectivity::Defragmentation_maps maps;
  mesh().defragment_in_place(maps);
  int n_vertices = mesh().n_total_vertices();

  // Initialize anchors
  _anchors.init(n_vertices);
  _anchors.clearStatic();

  // Initialize rest positions matrix (3 x n_vertices)
  _vertices_rest.resize(3, n_vertices);
  // Store rest positions for all vertices (all are active after defragmentation)
  for(int i = 0; i < n_vertices; ++i)
  {
    Mesh_connectivity::Vertex_iterator vi = mesh().vertex_at(i);
    _vertices_rest.col(i) = vi.xyz();
  }

  // Initialize deformed positions matrix (3 x n_vertices) to rest positions
  _vertices_deformed.resize(3, n_vertices);
  _vertices_deformed = _vertices_rest;

  // Initialize rotation matrices
  _R.resize(n_vertices);
  for(auto & R_i : _R)
  {
    R_i = Eigen::Matrix3d::Identity();
  }

  _Bf = Eigen::Matrix3Xd::Zero(3, n_vertices);

  // Build cotangent weight adjacency list
  build_adjacency_list();
  build_laplacian_matrix();
  // Memorize the sparsity structure of L for the solver
  _solver.analyzePattern(_L);
}

///
// Weights and neighbors builders
///

static inline double
_safe_cot(const Eigen::Vector3d & a, const Eigen::Vector3d & b, double eps = 1e-12)
{
  const double denom = a.cross(b).norm();
  if(denom <= eps)
    return 0.0; // robust: treat degenerate as 0
  return a.dot(b) / denom;
}

static inline double
_opposite_cot_in_face(Mesh_connectivity & M, Mesh_connectivity::Half_edge_iterator he)
{
  if(he.face().index() == M.hole_index)
    return 0.0;

  // triangle i -> j -> k; opposite angle is at vertex k
  auto he_next = he.next();
  int k = he_next.dest().index();

  Eigen::Vector3d vk = M.vertex_at(k).xyz();
  Eigen::Vector3d vi = he.origin().xyz();
  Eigen::Vector3d vj = he.dest().xyz();

  // vectors from k to i and k to j
  return _safe_cot(vi - vk, vj - vk);
}

std::vector<double>
_compute_opposite_cot_per_halfedge(Mesh_connectivity & M)
{
  std::vector<double> opp_cot(M.n_total_half_edges(), 0.0);
  for(int h = 0; h < M.n_total_half_edges(); ++h)
  {
    auto he = M.half_edge_at(h);
    if(!he.is_active())
      throw std::runtime_error("inactive half-edge (defragment first)");
    opp_cot[h] = _opposite_cot_in_face(M, he); // 0 on boundary
  }
  return opp_cot;
}

void
Mesh_modifier_arap::build_adjacency_list()
{

  // Helped by https://rodolphe-vaillant.fr/entry/69/c-code-for-cotangent-weights-over-a-triangular-mesh
  // After defragmentation in initialize(), n_total == n_active and all vertices are active
  int n_verts = mesh().n_active_vertices();
  // Reset adjacency list
  _adj.clear();
  _adj.resize(n_verts);

  // Compute cotangent weights for each half-edge
  auto opp_cot = _compute_opposite_cot_per_halfedge(mesh());

  // Loop over all vertices and store neighbor list with cotangent weights
  for(int i = 0; i < n_verts; ++i)
  {
    Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(i);
    if(!v.is_active())
    {
      throw std::runtime_error("Mesh not properly defragmented on initialize()");
    }
    auto v_ring = mesh().vertex_ring_at(i);
    do
    {
      Mesh_connectivity::Half_edge_iterator he_in = v_ring.half_edge();
      // Get the neighbor vertex (origin of the half-edge pointing to i)
      int j = he_in.origin().index();
      // w_ij = (cot(alpha) + cot(beta)) / 2
      double w_ij = (opp_cot[he_in.index()] + opp_cot[he_in.twin().index()]) / 2.0;
      // Store the neighbor and weight in adjacency list of vertex i
      _adj[i].push_back(Neighbor{j, w_ij});

    } while(v_ring.advance()); // Continue around the vertex until return to start
  }
}

// Assumes that the w_ij adjacency list is already built
void
Mesh_modifier_arap::build_laplacian_matrix()
{
  const int n = mesh().n_total_vertices();
  _L.resize(n, n);

  // Triplet is a (row, column, value) entry for sparse matrix construction
  std::vector<Eigen::Triplet<double>> tripletList;
  tripletList.reserve(n * 7); // Prealloc memory (average ~6 neighbors + diagonal)

  // Build the entries of L using (i,j, -w_ij) for off-diagonal and (i,i, sum_j w_ij) for diagonal
  for(int i = 0; i < n; ++i)
  {
    double sum_wij = 0.0;
    for(const Neighbor & neighbor : _adj[i])
    {
      int j = neighbor.index;
      double w_ij = neighbor.weight;
      sum_wij += w_ij;
      tripletList.emplace_back(i, j, -w_ij); // Off-diagonal entries
    }
    tripletList.emplace_back(i, i, sum_wij); // Diagonal entry
  }

  // Construct the sparse matrix from triplets
  _L.setFromTriplets(tripletList.begin(), tripletList.end());
}

void
Mesh_modifier_arap::build_blocks_from_constraints()
{
  const int n = _L.rows();

  // maps: global -> local
  std::vector<int> fmap(n, -1), cmap(n, -1);
  for(int a = 0; a < (int)_free_idx.size(); ++a)
    fmap[_free_idx[a]] = a;
  for(int a = 0; a < (int)_cons_idx.size(); ++a)
    cmap[_cons_idx[a]] = a;

  // Build Lff and Lfc by filtering triplets from full L
  std::vector<Eigen::Triplet<double>> Tf, Tfc;
  Tf.reserve(_L.nonZeros());
  Tfc.reserve(_L.nonZeros() / 4);

  for(int k = 0; k < _L.outerSize(); ++k)
  {
    for(Eigen::SparseMatrix<double>::InnerIterator it(_L, k); it; ++it)
    {
      int i = it.row(); // row
      int j = it.col(); // col
      double v = it.value();

      if(fmap[i] != -1)
      { // only rows of free vertices survive
        if(fmap[j] != -1)
        {
          Tf.emplace_back(fmap[i], fmap[j], v); // Lff
        }
        else if(cmap[j] != -1)
        {
          Tfc.emplace_back(fmap[i], cmap[j], v); // Lfc
        }
      }
    }
  }

  _Lff.resize(_free_idx.size(), _free_idx.size());
  _Lfc.resize(_free_idx.size(), _cons_idx.size());
  _Lff.setFromTriplets(Tf.begin(), Tf.end());
  _Lfc.setFromTriplets(Tfc.begin(), Tfc.end());
}

// Fill C (|cons| x 3) with the current target positions of the anchors including temporary moved one
Eigen::MatrixXd
Mesh_modifier_arap::_gather_constraint_matrix()
{
  // Allocate a matrix C of size (m x 3) for m constrained vertices
  const int m = (int)_cons_idx.size();
  Eigen::MatrixXd C(m, 3);
  for(int a = 0; a < m; ++a)
  {
    int idx = _cons_idx[a];
    if(idx == _anchors.getTemp())
    {
      C.row(a) = _anchors.getTempPosition().transpose();
    }
    else
    {
      // static anchors: keep at their current (or rest) position
      C.row(a) = _vertices_rest.col(idx).transpose(); // or current position
    }
  }

  return C;
}

void
Mesh_modifier_arap::build_rhs_bf_direct()
{
  // Gather constraint matrix
  _C = _gather_constraint_matrix();

  const int nf = (int)_free_idx.size();
  _Bf.resize(nf, 3);

  // Compute RHS for each free vertex only (avoid computing for constrained vertices)
  for(int a = 0; a < nf; ++a)
  {
    int i = _free_idx[a];
    Eigen::Vector3d b_i = Eigen::Vector3d::Zero();

    // Sum over all neighbors of this free vertex
    for(const Neighbor & nb : _adj[i])
    {
      int j = nb.index;
      double w_ij = nb.weight;
      Eigen::Vector3d p_ij = _vertices_rest.col(i) - _vertices_rest.col(j);
      b_i += w_ij * 0.5 * (_R[i] + _R[j]) * p_ij;
    }

    _Bf.row(a) = b_i.transpose();
  }

  // Subtract constraint contribution
  _Bf -= _Lfc * _C;
}

///
// Anchor management methods
///

bool
Mesh_modifier_arap::add_anchor(const int vertex_index)
{
  return _anchors.addStatic(vertex_index);
}


bool
Mesh_modifier_arap::remove_anchor(const int vertex_index)
{
  return _anchors.removeStatic(vertex_index);
}


bool
Mesh_modifier_arap::is_anchor(const int vertex_index) const
{
  return _anchors.isStaticAnchor(vertex_index);
}


void
Mesh_modifier_arap::clear_anchors()
{
  _anchors.clearStatic();
}


std::vector<int>
Mesh_modifier_arap::get_static_anchors()
{
  return _anchors.getStaticList();
}

/* Core ARAP Deformation Solver (Private)
Uses the set anchors in original static position with one temporary anchor moved to a new position.
Assumes output matrix is already properly sized and initialized with positions.
All public deformation methods delegate to this.
*/
bool
Mesh_modifier_arap::_solve_arap(const int temp_anchor, const Eigen::Vector3d & pulled, Eigen::Matrix3Xd & output)
{
  // Check if vertex is active
  Mesh_connectivity::Vertex_iterator v_pulled = mesh().vertex_at(temp_anchor);
  if(!v_pulled.is_active())
  {
    return false;
  }

  // Assign the temp anchor and check for rebuilding the constraint structures
  _anchors.setTemp(temp_anchor);
  _anchors.setTempPosition(pulled);
  int new_anchor_version = _anchors.version();
  if(_last_anchor_version != new_anchor_version)
  {
    // Rebuild constraint related structures
    _last_anchor_version = new_anchor_version;
    _cons_idx.clear();
    _free_idx.clear();
    for(int i = 0; i < mesh().n_total_vertices(); ++i)
    {
      if(_anchors.isStaticAnchor(i) or (i == temp_anchor))
      {
        _cons_idx.push_back(i);
      }
      else
      {
        _free_idx.push_back(i);
      }
    }
    build_blocks_from_constraints();
    // Factorize Lff for solving
    _solver.compute(_Lff);
    if(_solver.info() != Eigen::Success)
    {
      return false; // Factorization failed
    }
  }

  printf("Deformed mesh with anchor vertex %d \n", temp_anchor);
  printf(" - New anchor position: (%f, %f, %f) \n", pulled.x(), pulled.y(), pulled.z());

  _solve_deformation_with_anchors();
  // Copy deformed positions to output
  output = _vertices_deformed;

  return true;
}



// Solve the Lp = B constrained system with current Ri and given anchors
bool
Mesh_modifier_arap::_solve_deformation_with_anchors()
{
  // Verify L is precomputed
  if(_L.rows() == 0 || _L.cols() == 0)
  {
    return false;
  }

  // Verify vertices_rest is initialized
  if(_vertices_rest.cols() == 0)
  {
    return false;
  }
  // Build reduced RHS Bf and store in _Bf
  build_rhs_bf_direct();
  // Solve for free vertex positions: Lff * Pf = Bf
  Eigen::MatrixXd Pf = _solver.solve(_Bf);
  if(_solver.info() != Eigen::Success)
  {
    return false; // Solve failed
  }

  // Assemble full output positions
  for(int a = 0; a < (int)_free_idx.size(); ++a)
  {
    int i = _free_idx[a];
    _vertices_deformed.col(i) = Pf.row(a).transpose();
  }
  for(int a = 0; a < (int)_cons_idx.size(); ++a)
  {
    int i = _cons_idx[a];
    _vertices_deformed.col(i) = _C.row(a).transpose();
  }
  return true;
}


/* ARAP Deformation - In-place version (for GUI performance)
Uses the set anchors in original static position with one temporary anchor moved to a new position.
Modifies the passed set of deformed positions in-place (no allocation).
*/
bool
Mesh_modifier_arap::deform_with_temp_anchor(const int vertex_index,
    const Eigen::Vector3d & new_position,
    Eigen::Matrix3Xd & deformed_positions)
{
  // Direct pass-through to solver (GUI provides pre-allocated matrix)
  return _solve_arap(vertex_index, new_position, deformed_positions);
}


/* ARAP Deformation - Returning version (for testing & convenience)
Allocates a new matrix, initializes with current mesh positions, computes deformation, and returns it.
*/
Eigen::Matrix3Xd
Mesh_modifier_arap::compute_deformation(const int vertex_index, const Eigen::Vector3d & new_position)
{
  // Allocate new matrix
  Eigen::Matrix3Xd result(3, mesh().n_total_vertices());

  // Initialize with current mesh positions
  for(int i = 0; i < mesh().n_total_vertices(); ++i)
  {
    Mesh_connectivity::Vertex_iterator vi = mesh().vertex_at(i);
    if(vi.is_active())
    {
      result.col(i) = vi.xyz();
    }
  }

  // Solve
  bool success = _solve_arap(vertex_index, new_position, result);

  if(!success)
  {
    throw std::runtime_error("ARAP deformation failed for vertex " + std::to_string(vertex_index));
  }

  return result;
}

/* ARAP Deformation - Mesh-modifying version (for CLI & file output)
Computes deformation in a temporary matrix and applies directly to mesh vertex positions.
*/
bool
Mesh_modifier_arap::apply_deformation_to_mesh(const int vertex_index, const Eigen::Vector3d & new_position)
{

  Eigen::Matrix3Xd deformed_positions = compute_deformation(vertex_index, new_position);

  // Apply deformed positions back to mesh
  for(int i = 0; i < mesh().n_total_vertices(); ++i)
  {
    Mesh_connectivity::Vertex_iterator vi = mesh().vertex_at(i);
    if(vi.is_active())
    {
      vi.data().xyz = deformed_positions.col(i);
    }
  }

  return true;
}


} // end of mohecore
} // end of minimesh
