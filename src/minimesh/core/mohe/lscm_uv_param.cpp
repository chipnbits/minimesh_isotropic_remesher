#include <minimesh/core/mohe/lscm_uv_param.hpp>
#include <minimesh/core/util/assert.hpp>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/Dense>
#include <vector>
#include <set>
#include <cmath>

namespace minimesh
{
namespace mohecore
{

bool
LSCM_uv_param::choose_pinned_vertices()
{
  // find one boundary half-edge
  Mesh_connectivity::Half_edge_iterator start_he;
  bool found = false;
  for(int h = 0; h < mesh().n_total_half_edges(); ++h)
  {
    auto he = mesh().half_edge_at(h);
    if(he.is_active() && he.face().is_equal(mesh().hole()))
    {
      start_he = he;
      found = true;
      break;
    }
  }
  if(!found) return false; // closed mesh – out of scope

  int best_i = -1, best_j = -1;

  if(_pinning_strategy == PinningStrategy::LONGEST_EDGE)
  {
    // Find the longest edge on the boundary
    double best_d2 = -1.0;
    auto cur_he = start_he;
    do
    {
      if(cur_he.is_active())
      {
        auto v0 = cur_he.origin();
        auto v1 = cur_he.dest();
        if(v0.is_active() && v1.is_active())
        {
          Eigen::Vector3d p0 = v0.xyz();
          Eigen::Vector3d p1 = v1.xyz();
          double d2 = (p1 - p0).squaredNorm();
          if(d2 > best_d2)
          {
            best_d2 = d2;
            best_i = v0.index();
            best_j = v1.index();
          }
        }
      }
      cur_he = cur_he.next();
    } while(!cur_he.is_equal(start_he));
  }
  else // PinningStrategy::MAX_DISTANCE
  {
    // collect boundary vertices
    std::vector<int> boundary;
    const int nV = mesh().n_total_vertices();
    std::vector<bool> is_boundary(nV, false);

    // walk boundary loop
    auto cur_he = start_he;
    do
    {
      auto v = cur_he.origin();
      if(v.is_active())
      {
        int idx = v.index();
        if(!is_boundary[idx])
        {
          is_boundary[idx] = true;
          boundary.push_back(idx);
        }
      }
      cur_he = cur_he.next();
    } while(!cur_he.is_equal(start_he));

    if(boundary.size() < 2) return false;

    // brute-force max Euclidean distance on boundary vertices
    double best_d2 = -1.0;

    for(size_t a = 0; a < boundary.size(); ++a)
    {
      Eigen::Vector3d pa = mesh().vertex_at(boundary[a]).xyz();
      for(size_t b = a + 1; b < boundary.size(); ++b)
      {
        Eigen::Vector3d pb = mesh().vertex_at(boundary[b]).xyz();
        double d2 = (pa - pb).squaredNorm();
        if(d2 > best_d2)
        {
          best_d2 = d2;
          best_i = boundary[a];
          best_j = boundary[b];
        }
      }
    }
  }

  if(best_i < 0 || best_j < 0) return false;

  _pin0 = best_i;
  _pin1 = best_j;

  _uv_coords[_pin0] = Eigen::Vector2d(0.0, 1.0);
  _uv_coords[_pin1] = Eigen::Vector2d(0.0, 0.0);

  return true;
}

bool
LSCM_uv_param::build_lscm_system(Eigen::SparseMatrix<double>& A,
                                 Eigen::VectorXd& b)
{
  const int nV = _mesh.n_total_vertices();
  const int nF = _mesh.n_total_faces();

  // map mesh vertex -> compact index for free vertices (converts mesh indices to linear system indices)
  std::vector<int> mesh_to_free(nV, -1);
  int free_count = 0;
  for(int v = 0; v < nV; ++v)
  {
    auto vert = _mesh.vertex_at(v);
    if(!vert.is_active()) continue;
    if(v == _pin0 || v == _pin1) continue;
    mesh_to_free[v] = free_count++;
  }

  if(free_count == 0) return false;

  // System size: 2 rows per face (X and Y constraints), 2 cols per free vertex
  const int n_rows = 2 * nF; 
  const int n_cols = 2 * free_count;

  A.resize(n_rows, n_cols);
  b = Eigen::VectorXd::Zero(n_rows);

  std::vector<Eigen::Triplet<double>> triples;
  triples.reserve(n_rows * 6); // Approx 6 entries per row

  int current_face = 0;
  for(int f = 0; f < nF; ++f)
  {
    auto face = _mesh.face_at(f);
    if(!face.is_active() || face.is_equal(_mesh.hole())) continue;

    // --- 1. Get Triangle Vertices, edges of face ---
    auto he0 = face.half_edge();
    auto he1 = he0.next();
    auto he2 = he1.next();

    int va = he0.origin().index();
    int vb = he1.origin().index();
    int vc = he2.origin().index(); 

    Eigen::Vector3d p1 = he0.origin().xyz();
    Eigen::Vector3d p2 = he1.origin().xyz();
    Eigen::Vector3d p3 = he2.origin().xyz();

    // --- 2. Geometric Constants---
    // We need to match P3 - P1 = s * Rot * (P2 - P1)
    Eigen::Vector3d e12 = p2 - p1; // Vector P1->P2
    Eigen::Vector3d e13 = p3 - p1; // Vector P1->P3

    double len12_sq = e12.squaredNorm();
    double cross_norm = e12.cross(e13).norm(); // 2 * Area
    double dot_prod   = e12.dot(e13);

    if(len12_sq < 1e-12) continue; // Degenerate edge

    // Calculate the coefficients directly without trig functions
    // C_sim = s * cos(alpha) 
    // -> sin(a3)/sin(a2) x cos(a1) = |cxa|/|c||a| x |a.b|/|a||b| x |b||c|/|bxc|= |a.b|/|a|^2
    // S_sim = s * sin(alpha)
    // -> |a x b| = |a||b|sin(alpha) but |axb| is double  triangle area so |axb| = |axc| = |bxc|
    // -> sin(a3)/sin(a2) x sin(a1) = |axb|/|a||b| x  |cxa|/|c||a| x |b||c|/|bxc| = |axb|/|a|^2
    double C_sim = dot_prod / len12_sq;
    double S_sim = cross_norm / len12_sq;

    // --- 3. Build Linear Rows ---
    // We minimize || (u3,v3) - (u1,v1) - R * ((u2,v2) - (u1,v1)) ||^2
    // R is scaled rotation matrix R=
    //   C_sim  * (u2-u1) S_sim * (v2-v1)
    //   -S_sim * (u2-u1) C_sim * (v2-v1)
    // So we have two equations per face:

    // Row 1 (X-residual): u3 - u1 - Rot_x = 0
    //   => u3 - u1 - [ C(u2-u1) + S(v2-v1) ] = 0
    //   => u3 + (C - 1)u1 - C*u2 - S*v2 + S*v1 = 0
    // Coefficients for X row:
    //   u1: (C_sim - 1)
    //   v1: S_sim
    //   u2: -C_sim
    //   v2: -S_sim
    //   u3: 1
    //   v3: 0

    // Row 2 (Y-residual): v3 - v1 - Rot_y = 0
    //   => v3 - v1 - [ -S(u2-u1) + C(v2-v1) ] = 0
    //   => v3 + (C - 1)v1 - C*v2 + S*u2 - S*u1 = 0
    // Coefficients for Y row:
    //   u1: -S_sim
    //   v1: (C_sim - 1)
    //   u2: S_sim
    //   v2: -C_sim
    //   u3: 0
    //   v3: 1
    
    // Row indices in matrix for two equations per face
    int row_x = 2 * current_face;
    int row_y = 2 * current_face + 1;

    // Helper lambda to clean up triplet insertion
    auto add_triplet = [&](int v_idx, double val_u, double val_v, int row) {
      if(v_idx == _pin0 || v_idx == _pin1) {
        // Pinned: Move to RHS (b vector)
        double u_pinned = _uv_coords[v_idx].x();
        double v_pinned = _uv_coords[v_idx].y();
        b[row] -= (val_u * u_pinned + val_v * v_pinned);
      } else {
        // Free: Add to Matrix A with augmented U vector [u..., v...]
        int col_u = mesh_to_free[v_idx];
        int col_v = col_u + free_count;
        if(std::abs(val_u) > 1e-9) triples.emplace_back(row, col_u, val_u);
        if(std::abs(val_v) > 1e-9) triples.emplace_back(row, col_v, val_v);
      }
    };

    // --- Fill Row X ---
    int idx1 = va;
    int idx2 = vb;
    int idx3 = vc;

    add_triplet(idx1, (C_sim - 1.0),  S_sim,     row_x);
    add_triplet(idx2, -C_sim,        -S_sim,     row_x);
    add_triplet(idx3, 1.0,            0.0,       row_x);

    // --- Fill Row Y ---
    add_triplet(idx1, -S_sim,        (C_sim - 1.0), row_y);
    add_triplet(idx2, S_sim,         -C_sim,        row_y);
    add_triplet(idx3, 0.0,            1.0,          row_y);

    current_face++;
  }


  int used_rows = 2 * current_face;
  if(used_rows == 0) return false;

  A.conservativeResize(used_rows, n_cols);
  b.conservativeResize(used_rows);
  A.setFromTriplets(triples.begin(), triples.end());

  return true;
}

bool
LSCM_uv_param::solve_system(const Eigen::SparseMatrix<double>& A,
                            const Eigen::VectorXd& b)
{
  Eigen::SparseMatrix<double> AtA = A.transpose() * A;
  Eigen::VectorXd Atb = A.transpose() * b;

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(AtA);
  if(solver.info() != Eigen::Success)
    return false;

  Eigen::VectorXd x = solver.solve(Atb);
  if(solver.info() != Eigen::Success)
    return false;

  // write back to _uv_coords
  const int nV = mesh().n_total_vertices();

  // rebuild mapping to free indices
  std::vector<int> mesh_to_free(nV, -1);
  int free_count = 0;
  for(int v = 0; v < nV; ++v)
  {
    auto vert = mesh().vertex_at(v);
    if(!vert.is_active()) continue;
    if(v == _pin0 || v == _pin1) continue;
    mesh_to_free[v] = free_count++;
  }

  for(int v = 0; v < nV; ++v)
  {
    auto vert = mesh().vertex_at(v);
    if(!vert.is_active()) continue;

    if(v == _pin0 || v == _pin1)
      continue; // already set

    int free_id = mesh_to_free[v];
    if(free_id < 0) continue;

    double u = x[free_id];
    double vcoord = x[free_id + free_count];

    _uv_coords[v] = Eigen::Vector2d(u, vcoord);
  }

  return true;
}

bool
LSCM_uv_param::compute_parameterization()
{
  _uv_coords.assign(mesh().n_total_vertices(), Eigen::Vector2d::Zero());

  if(!choose_pinned_vertices())
    return false;

  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd b;

  if(!build_lscm_system(A, b))
    return false;

  if(!solve_system(A, b))
    return false;

  _is_computed = true;
  return true;
}


} // end of mohecore
} // end of minimesh
