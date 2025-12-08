#include <minimesh/core/mohe/fixed_uv_param.hpp>
#include <minimesh/core/util/assert.hpp>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <set>
#include <cmath>

namespace minimesh
{
namespace mohecore
{

//
// Given two vertices, this function return the index of the half-edge going from v0 to v1.
// Returns -1 if no half-edge exists between the two vertices.
//
int
Fixed_boundary_uv_param::get_halfedge_between_vertices(const int v0, const int v1)
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


bool
Fixed_boundary_uv_param::flip_edge(const int he_index)
{
  //
  // Take a reference to all involved entities
  //

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
} // All done


//
// Get UV coordinate for a specific vertex
//
Eigen::Vector2d
Fixed_boundary_uv_param::get_uv_at_vertex(int vertex_index) const
{
  assert(vertex_index >= 0 && vertex_index < (int)_uv_coords.size());
  return _uv_coords[vertex_index];
}


//
// Initialize boundary and interior vertex lists
//
void
Fixed_boundary_uv_param::initialize_vertex_lists()
{
  _boundary_loop_indices.clear();
  _interior_vertex_indices.clear();

  const int nV = mesh().n_total_vertices();
  std::vector<bool> is_boundary(nV, false);

  // Find a boundary half-edge by iterating through all half-edges
  Mesh_connectivity::Half_edge_iterator start_he;
  bool found = false;
  for(int i = 0; i < mesh().n_total_half_edges(); ++i)
  {
    auto he = mesh().half_edge_at(i);
    if(he.is_active() && he.face().is_equal(mesh().hole()))
    {
      start_he = he;
      found = true;
      break;
    }
  }

  if(found){
    Mesh_connectivity::Half_edge_iterator cur_he = start_he;
    do
    {
      auto v = cur_he.origin();
      int idx = v.index();

      if (v.is_active())         // <-- only track active boundary vertices
      {
        _boundary_loop_indices.push_back(idx);
        is_boundary[idx] = true;    // mark as boundary
      }

      cur_he = cur_he.next();
    } while (!cur_he.is_equal(start_he));
  }

  // Collect active interior vertices (active and not marked as boundary)
  for (int i = 0; i < nV; ++i)
  {
    auto v = mesh().vertex_at(i);
    if (v.is_active() && !is_boundary[i])
    {
      _interior_vertex_indices.push_back(i);
    }
  }
}


//
// Map boundary vertices to UV space (e.g., a circle)
//
void
Fixed_boundary_uv_param::map_boundary_to_uv()
{
  //Map based on arc length: preserve boundary shape better
  const auto& v = _boundary_loop_indices;
  const int n_boundary = static_cast<int>(v.size());
  assert(n_boundary > 0);
  // use min valued index as starting point
  int minVal = *std::min_element(v.begin(), v.end());
  
  // ring iterate back to boundary from starting vertex
  auto ring_iter = mesh().vertex_ring_at(minVal);
  assert(ring_iter.reset_boundary()); // should be on boundary
  // on edge pointing to minVal on boundary, so next is the exiting half-edge
  auto start_he = ring_iter.half_edge().next();

  // Compute vN-1 edge lengths, gather vertex indices excluding start vertex
  std::vector<double> edge_lengths;
  std::vector<int> vertex_indices;
  double total_length = 0.0;
  auto cur_he = start_he;
  do
  {
    double edge_length = (cur_he.origin().xyz() - cur_he.dest().xyz()).norm();
    edge_lengths.push_back(edge_length);
    vertex_indices.push_back(cur_he.dest().index());
    total_length += edge_length;
    cur_he = cur_he.next();
  } while (!cur_he.is_equal(start_he));
  assert(edge_lengths.size() == static_cast<std::size_t>(n_boundary));

  // assign UVs based on cumulative length on unit circle
  double cumulative_length = 0.0;
  cur_he = start_he;
  for (std::size_t i = 0; i < edge_lengths.size(); ++i)
  {
    cumulative_length += edge_lengths[i];
    double angle = 2.0 * M_PI * cumulative_length / total_length;
    _uv_coords[vertex_indices[i]] = Eigen::Vector2d(std::cos(angle), std::sin(angle));
  }
}

//
// Compute cotangent weight for each half-edge
// Returns a vector indexed by half-edge index
// Twin half-edges will have the same weight value
//
bool
Fixed_boundary_uv_param::compute_edge_weights()
{
  const int n_halfedges = mesh().n_total_half_edges();
  _Dij.resize(n_halfedges, 0.0);
  _lambda_ij.resize(n_halfedges, 0.0);

  // HARVEST COTANGENT WEIGHTS
  for(int h = 0; h < n_halfedges; ++h)
  {
    auto he = mesh().half_edge_at(h);
    int twin_idx = he.twin().index();

    if(!he.is_active())
      continue;  // Skip inactive half-edges

    // Case of boundary half-edge, cotangent is 0.0
    if (he.face().index() == mesh().hole_index)
    {
      continue;
    }

    // Compute cotangent of angle opposite to half-edge h
    // The opposite vertex is at the next().dest() position
    int k = he.next().dest().index();
    Eigen::Vector3d xyz_i = he.origin().xyz();
    Eigen::Vector3d xyz_j = he.dest().xyz();
    Eigen::Vector3d xyz_k = mesh().vertex_at(k).xyz();

    // Get the cotangent of angle at k between edges (k->i) and (k->j)
    Eigen::Vector3d v_ki = xyz_i - xyz_k;
    Eigen::Vector3d v_kj = xyz_j - xyz_k;
    const double denom = v_ki.cross(v_kj).norm();

    // Robust approach, only accumulate for non-degenerate angles
    if(denom > 1e-12)
    {
      double cot_value = v_ki.dot(v_kj) / denom;
      _Dij[h] += cot_value;
      _Dij[twin_idx] += cot_value;  // Same weight for twin
    }
  }

  // NORMALIZE WEIGHTS
  // sum of outgoing weights at each vertex is 1.0
  const int n_vertices = mesh().n_total_vertices();
  for (int v=0; v < n_vertices; ++v)
  {
    // Sum of outgoing weights
    double weight_sum = 0.0;

    // Iterate over outgoing half-edges to get sum of weights
    auto ring_iter = mesh().vertex_ring_at(v);
    do
    {
      auto he = ring_iter.half_edge().twin(); // outgoing half-edge (iterator points to incoming)
      if(he.is_active())
      {
        weight_sum += _Dij[he.index()];
      }
    } while(ring_iter.advance());

    // Normalize outgoing weights using the sum
    if(weight_sum > 1e-12)
    {
      ring_iter = mesh().vertex_ring_at(v);
      do
      {
        auto he = ring_iter.half_edge().twin();
        if(he.is_active())
        {
          int he_idx = he.index();
          _lambda_ij[he_idx] = _Dij[he_idx] / weight_sum;
        }
      } while(ring_iter.advance());
    }
  }

  // PRint summary statistics
  // double min_weight = std::numeric_limits<double>::max();;
  // double max_weight = std::numeric_limits<double>::lowest();
  // double avg_weight = 0.0;
  // int count = 0;
  // for(int h = 0; h < n_halfedges; ++h)
  // {
  //   if(mesh().half_edge_at(h).is_active())
  //   {
  //     double w = _lambda_ij[h];
  //     if(w < min_weight) min_weight = w;
  //     if(w > max_weight) max_weight = w;
  //     avg_weight += w;
  //     count++;
  //   }
  // }
  // if(count > 0)
  //   avg_weight /= static_cast<double>(count);
  // printf("Cotangent Weights: min = %.6f, max = %.6f, avg = %.6f\n", min_weight, max_weight, avg_weight);

  return true;
}

//
// Solve for interior vertex UV positions
//
bool
Fixed_boundary_uv_param::solve_interior_uvs()
{
  const int n_interior = static_cast<int>(_interior_vertex_indices.size());
  if(n_interior == 0) return true;

  std::vector<int> mesh_to_solver_idx(mesh().n_total_vertices(), -1);
  for(int i = 0; i < n_interior; ++i)
  {
    mesh_to_solver_idx[_interior_vertex_indices[i]] = i;
  }

  // System: (I - Lambda_interior) * U_interior = Lambda_boundary * U_boundary
  Eigen::SparseMatrix<double> A(n_interior, n_interior);
  std::vector<Eigen::Triplet<double>> coefficients;
  
  Eigen::VectorXd b_u = Eigen::VectorXd::Zero(n_interior);
  Eigen::VectorXd b_v = Eigen::VectorXd::Zero(n_interior);

  // 3. Fill the Matrix and RHS
  for(int i = 0; i < n_interior; ++i)
  {
    int v_idx = _interior_vertex_indices[i];

    // Diagonal entry: A(i,i) = 1
    coefficients.push_back(Eigen::Triplet<double>(i, i, 1.0));

    // Iterate over neighbors of the interior vertex v_idx
    auto ring_iter = mesh().vertex_ring_at(v_idx);
    do
    {
      auto he_outgoing = ring_iter.half_edge().twin();
      int neighbor_idx = he_outgoing.dest().index();
      
      // Get normalized weight lambda_ij (stored on the outgoing half-edge index)
      double weight = _lambda_ij[he_outgoing.index()];

      // Check if neighbor is Interior or Boundary
      int neighbor_solver_idx = mesh_to_solver_idx[neighbor_idx];

      if(neighbor_solver_idx != -1) 
      {
        // CASE: Neighbor is Interior
        coefficients.push_back(Eigen::Triplet<double>(i, neighbor_solver_idx, -weight));
      }
      else
      {
        // CASE: Neighbor is Boundary (Fixed)
        // Move to RHS: +lambda_ij * u_fixed 
        Eigen::Vector2d fixed_uv = _uv_coords[neighbor_idx];
        b_u[i] += weight * fixed_uv.x();
        b_v[i] += weight * fixed_uv.y();
      }

    } while(ring_iter.advance());
  }

  A.setFromTriplets(coefficients.begin(), coefficients.end());
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);

  if(solver.info() != Eigen::Success)
  {
    return false;
  }

  Eigen::VectorXd x_u = solver.solve(b_u);
  Eigen::VectorXd x_v = solver.solve(b_v);

  if(solver.info() != Eigen::Success)
  {
    return false;
  }

  // Store results back to _uv_coords
  for(int i = 0; i < n_interior; ++i)
  {
    int v_idx = _interior_vertex_indices[i];
    _uv_coords[v_idx] = Eigen::Vector2d(x_u[i], x_v[i]);
  }

  return true;
}


//
// Compute the UV parameterization
//
bool
Fixed_boundary_uv_param::compute_parameterization()
{
  _uv_coords.resize(mesh().n_total_vertices());
  initialize_vertex_lists();

  // Check that we have a boundary
  if(_boundary_loop_indices.empty())
  {
    return false;
  }

  // Step 2: Map boundary vertices to UV space
  map_boundary_to_uv();
  // Step 2.5: Compute edge weights to populate _Dij and _lambda_ij
  if (!compute_edge_weights())
  {
    return false;
  }

  // Step 3: Solve for interior vertex positions using _lambda_ij weights
  bool success = solve_interior_uvs();

  if(success)
  {
    _is_computed = true;
  }

  return success;
}




} // end of mohecore
} // end of minimesh
