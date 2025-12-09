#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <minimesh/core/mohe/mesh_connectivity.hpp>
#include <string>

namespace minimesh
{
namespace mohecore
{

class LSCM_uv_param
{
public:
  LSCM_uv_param(Mesh_connectivity & m)
  : _mesh(m)
  , _is_computed(false)
  {
  }

  bool compute_parameterization();

  Eigen::Vector2d get_uv_at_vertex(int v) const
  {
    assert(v >= 0 && v < (int)_uv_coords.size());
    return _uv_coords[v];
  }

private:
  Mesh_connectivity & mesh() { return _mesh; }
  const Mesh_connectivity & mesh() const { return _mesh; }

  bool choose_pinned_vertices(); // pick the 2 anchors
  bool build_lscm_system(Eigen::SparseMatrix<double> & A, Eigen::VectorXd & b);
  bool solve_system(const Eigen::SparseMatrix<double> & A, const Eigen::VectorXd & b);

private:
  Mesh_connectivity & _mesh;
  std::vector<Eigen::Vector2d> _uv_coords;
  bool _is_computed;

  int _pin0 = -1;
  int _pin1 = -1;
};
}
}