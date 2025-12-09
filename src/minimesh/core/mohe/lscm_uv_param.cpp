#include <minimesh/core/mohe/lscm_uv_param.hpp>
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

bool
LSCM_uv_param::choose_pinned_vertices()
{
  // collect boundary vertices
  std::vector<int> boundary;
  const int nV = mesh().n_total_vertices();

  std::vector<bool> is_boundary(nV, false);

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
  int best_i = -1, best_j = -1;

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

  _pin0 = best_i;
  _pin1 = best_j;

  // Fix them in texture space, e.g. along x-axis
  _uv_coords[_pin0] = Eigen::Vector2d(0.0, 0.0);
  _uv_coords[_pin1] = Eigen::Vector2d(1.0, 0.0);

  return true;
}

static void
triangle_3d_to_2d(const Eigen::Vector3d& p0,
                  const Eigen::Vector3d& p1,
                  const Eigen::Vector3d& p2,
                  double& x0, double& y0,
                  double& x1, double& y1,
                  double& x2, double& y2)
{
  Eigen::Vector3d e1 = p1 - p0;
  Eigen::Vector3d e2 = p2 - p0;

  double len_e1 = e1.norm();
  if(len_e1 < 1e-12)
  {
    x0=y0=x1=y1=x2=y2=0.0;
    return;
  }

  Eigen::Vector3d ex = e1 / len_e1;
  Eigen::Vector3d n  = e1.cross(e2);
  double n_norm = n.norm();
  if(n_norm < 1e-12)
  {
    // Degenerate triangle – just project onto e1
    Eigen::Vector3d ey = Eigen::Vector3d::UnitY();
    if(fabs(ex.dot(ey)) > 0.9) ey = Eigen::Vector3d::UnitZ();
    ey = (ey - ey.dot(ex)*ex).normalized();

    x0 = y0 = 0.0;
    x1 = len_e1; y1 = 0.0;
    Eigen::Vector3d v2 = p2 - p0;
    x2 = v2.dot(ex);
    y2 = v2.dot(ey);
    return;
  }
  n.normalize();
  Eigen::Vector3d ey = n.cross(ex);

  x0 = 0.0; y0 = 0.0;
  x1 = len_e1; y1 = 0.0;

  Eigen::Vector3d v2 = p2 - p0;
  x2 = v2.dot(ex);
  y2 = v2.dot(ey);
}

bool
LSCM_uv_param::build_lscm_system(Eigen::SparseMatrix<double>& A,
                                 Eigen::VectorXd& b)
{
  const int nV = _mesh.n_total_vertices();
  const int nF = _mesh.n_total_faces();

  // map mesh vertex -> compact index for free vertices
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

  //https://www.cs.jhu.edu/~misha/Fall09/Levy02.pdf
  // A has 4 M_f blocks that are Fx(free) each
  const int max_rows = 2 * nF; // 2 rows per face (for x and y components)
  const int n_cols   = 2 * free_count;

  A.resize(max_rows, n_cols);
  b = Eigen::VectorXd::Zero(max_rows);

  std::vector<Eigen::Triplet<double>> triples;
  triples.reserve(max_rows * 10);

  int current_face = 0;
  for(int f = 0; f < nF; ++f)
  {
    auto face = _mesh.face_at(f);
    if(!face.is_active()) continue;
    if(face.is_equal(_mesh.hole())) continue;

    auto he0 = face.half_edge();
    auto he1 = he0.next();
    auto he2 = he1.next();

    int va = he0.origin().index();
    int vb = he1.origin().index();
    int vc = he2.origin().index(); // usce c as the pivot

    Eigen::Vector3d pa = he0.origin().xyz();
    Eigen::Vector3d pb = he1.origin().xyz();
    Eigen::Vector3d pc = he2.origin().xyz();

    // triangle area
    Eigen::Vector3d e0 = pb - pa;
    Eigen::Vector3d e1 = pc - pa;
    double dT_f = e0.cross(e1).norm();
    if(dT_f < 1e-16) continue;

    double w = std::sqrt(dT_f); // row weight, matches energy sum A_T*|...|^2

    // angle at vertex c
    Eigen::Vector3d v1 = pa - pc;
    Eigen::Vector3d v2 = pb - pc;
    double len1 = v1.norm();
    double len2 = v2.norm();
    if(len1 < 1e-16 || len2 < 1e-16) continue;

    double cos_theta = v1.dot(v2) / (len1 * len2);
    cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
    double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));

    int row_real = 2 * current_face; 
    int row_imag = row_real + 1;

    auto add_coeff = [&](int mesh_v,
                         double coeff_u,
                         double coeff_v,
                         int row)
    {
      // x = [u_free(0..free_count-1), v_free(0..free_count-1)]
      if(mesh_v == _pin0 || mesh_v == _pin1)
      {
        // move to RHS: sum_free(...) = - sum_pinned(coeff * value)
        Eigen::Vector2d uv = _uv_coords[mesh_v];
        double u_p = uv.x();
        double v_p = uv.y();

        b[row] -= w * (coeff_u * u_p + coeff_v * v_p);
      }
      else
      {
        int free_id = mesh_to_free[mesh_v];
        if(free_id < 0) return; 

        int col_u = free_id;
        int col_v = free_id + free_count;

        if(std::abs(coeff_u) > 0.0)
          triples.emplace_back(row, col_u, w * coeff_u);
        if(std::abs(coeff_v) > 0.0)
          triples.emplace_back(row, col_v, w * coeff_v);
      }
    };

    // === REAL ROW (from real part of complex equation) ===
    // u_a - u_c - (cosθ (u_b - u_c) - sinθ (v_b - v_c)) = 0
    // -> coefficients:
    //   u_a:  +1
    //   u_b:  -cosθ
    //   u_c:  cosθ - 1
    //   v_b:  +sinθ
    //   v_c:  -sinθ

    add_coeff(va, /*coeff_u=*/ 1.0,      /*coeff_v=*/ 0.0,      row_real);
    add_coeff(vb, /*coeff_u=*/-cos_theta,/*coeff_v=*/ sin_theta,row_real);
    add_coeff(vc, /*coeff_u=*/ cos_theta - 1.0,
                 /*coeff_v=*/-sin_theta, row_real);

    // === IMAG ROW (from imaginary part) ===
    // v_a - v_c - (sinθ (u_b - u_c) + cosθ (v_b - v_c)) = 0
    // -> coefficients:
    //   u_b:  -sinθ
    //   u_c:  +sinθ
    //   v_a:  +1
    //   v_b:  -cosθ
    //   v_c:  cosθ - 1

    add_coeff(va, /*coeff_u=*/ 0.0,      /*coeff_v=*/ 1.0,      row_imag);
    add_coeff(vb, /*coeff_u=*/-sin_theta,/*coeff_v=*/-cos_theta,row_imag);
    add_coeff(vc, /*coeff_u=*/ sin_theta,
                 /*coeff_v=*/ cos_theta - 1.0, row_imag);

    ++current_face;
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
