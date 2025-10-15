#include <Eigen/Core>
#include <Eigen/Dense>
#include <minimesh/core/mohe/mesh_io.hpp>
#include <minimesh/core/mohe/mesh_modifier_edge_collapse.hpp>

#include <iostream>

#include <doctest.h>

using namespace minimesh::mohecore;

namespace
{
// Helper function to compare doubles with tolerance
const double EPSILON = 1e-10;

// Helper to create a plane from normal and distance
Eigen::Vector4d
make_plane(double nx, double ny, double nz, double d)
{
  return Eigen::Vector4d(nx, ny, nz, d);
}

// Helper to create a full 4x4 matrix from outer product of plane vector
Eigen::Matrix4d
plane_outer_product(const Eigen::Vector4d & p)
{
  return p * p.transpose();
}

// Helper to load mesh from file
bool
load_test_mesh(Mesh_connectivity & mesh, const std::string & filename)
{
  Mesh_io io(mesh);
  try
  {
    io.read_obj_general(filename);
    return mesh.check_sanity_slowly(false);
  }
  catch(...)
  {
    std::cerr << "Failed to load mesh: " << filename << std::endl;
    return false;
  }
}

} // anonymous namespace

TEST_CASE("SymQuadric: Default constructor initializes to zero")
{
  SymQuadric Q;

  for(int i = 0; i < 10; ++i)
  {
    CHECK(Q.q[i] == 0.0);
  }
}

TEST_CASE("SymQuadric: setZero properly zeros all entries")
{
  SymQuadric Q;

  // Set some non-zero values
  for(int i = 0; i < 10; ++i)
  {
    Q.q[i] = static_cast<double>(i + 1);
  }

  // Reset to zero
  Q.setZero();

  for(int i = 0; i < 10; ++i)
  {
    CHECK(Q.q[i] == 0.0);
  }
}

TEST_CASE("SymQuadric: operator+= adds quadrics correctly")
{
  SymQuadric Q1, Q2;

  // Set Q1 to some values
  for(int i = 0; i < 10; ++i)
  {
    Q1.q[i] = static_cast<double>(i);
  }

  // Set Q2 to some values
  for(int i = 0; i < 10; ++i)
  {
    Q2.q[i] = static_cast<double>(i * 2);
  }

  // Add Q2 to Q1
  Q1 += Q2;

  // Check results
  for(int i = 0; i < 10; ++i)
  {
    CHECK(Q1.q[i] == doctest::Approx(i + i * 2));
  }
}

TEST_CASE("SymQuadric: operator+ creates new quadric with sum")
{
  SymQuadric Q1, Q2;

  for(int i = 0; i < 10; ++i)
  {
    Q1.q[i] = 1.0;
    Q2.q[i] = 2.0;
  }

  SymQuadric Q3 = Q1 + Q2;

  // Check Q3 has sum
  for(int i = 0; i < 10; ++i)
  {
    CHECK(Q3.q[i] == doctest::Approx(3.0));
  }

  // Check Q1 and Q2 unchanged
  for(int i = 0; i < 10; ++i)
  {
    CHECK(Q1.q[i] == doctest::Approx(1.0));
    CHECK(Q2.q[i] == doctest::Approx(2.0));
  }
}

TEST_CASE("SymQuadric: addPlane with simple plane")
{
  SymQuadric Q;

  // Add a plane with normal [1, 0, 0] and d = 0
  Eigen::Vector4d plane = make_plane(1.0, 0.0, 0.0, 0.0);
  Q.addPlane(plane);

  // Expected: outer product of [1, 0, 0, 0]
  // Only a00 should be 1.0, rest should be 0
  CHECK(Q.q[0] == doctest::Approx(1.0)); // a00
  CHECK(Q.q[1] == doctest::Approx(0.0)); // a01
  CHECK(Q.q[2] == doctest::Approx(0.0)); // a02
  CHECK(Q.q[3] == doctest::Approx(0.0)); // a03
  CHECK(Q.q[4] == doctest::Approx(0.0)); // a11
  CHECK(Q.q[5] == doctest::Approx(0.0)); // a12
  CHECK(Q.q[6] == doctest::Approx(0.0)); // a13
  CHECK(Q.q[7] == doctest::Approx(0.0)); // a22
  CHECK(Q.q[8] == doctest::Approx(0.0)); // a23
  CHECK(Q.q[9] == doctest::Approx(0.0)); // a33
}

TEST_CASE("SymQuadric: addPlane accumulates multiple planes")
{
  SymQuadric Q;

  // Add first plane: [1, 0, 0, 0]
  Q.addPlane(make_plane(1.0, 0.0, 0.0, 0.0));

  // Add second plane: [0, 1, 0, 0]
  Q.addPlane(make_plane(0.0, 1.0, 0.0, 0.0));

  // Expected: sum of two outer products
  CHECK(Q.q[0] == doctest::Approx(1.0)); // a00 = 1*1
  CHECK(Q.q[1] == doctest::Approx(0.0)); // a01 = 1*0 + 0*1
  CHECK(Q.q[4] == doctest::Approx(1.0)); // a11 = 0*0 + 1*1
}

TEST_CASE("SymQuadric: addPlane matches full matrix outer product")
{
  SymQuadric Q;

  // Test with a general plane
  Eigen::Vector4d plane(0.6, 0.8, 0.0, -2.0);
  Q.addPlane(plane);

  // Compute expected full matrix
  Eigen::Matrix4d expected = plane_outer_product(plane);
  Eigen::Matrix4d actual = Q.toMatrix();

  // Compare all entries
  for(int i = 0; i < 4; ++i)
  {
    for(int j = 0; j < 4; ++j)
    {
      CHECK(actual(i, j) == doctest::Approx(expected(i, j)));
    }
  }
}

TEST_CASE("SymQuadric: toMatrix reconstructs full symmetric matrix")
{
  SymQuadric Q;

  // Set specific values in compact form
  Q.q[0] = 1.0; // a00
  Q.q[1] = 2.0; // a01
  Q.q[2] = 3.0; // a02
  Q.q[3] = 4.0; // a03
  Q.q[4] = 5.0; // a11
  Q.q[5] = 6.0; // a12
  Q.q[6] = 7.0; // a13
  Q.q[7] = 8.0; // a22
  Q.q[8] = 9.0; // a23
  Q.q[9] = 10.0; // a33

  Eigen::Matrix4d M = Q.toMatrix();

  // Check diagonal
  CHECK(M(0, 0) == doctest::Approx(1.0));
  CHECK(M(1, 1) == doctest::Approx(5.0));
  CHECK(M(2, 2) == doctest::Approx(8.0));
  CHECK(M(3, 3) == doctest::Approx(10.0));

  // Check upper triangle
  CHECK(M(0, 1) == doctest::Approx(2.0));
  CHECK(M(0, 2) == doctest::Approx(3.0));
  CHECK(M(0, 3) == doctest::Approx(4.0));
  CHECK(M(1, 2) == doctest::Approx(6.0));
  CHECK(M(1, 3) == doctest::Approx(7.0));
  CHECK(M(2, 3) == doctest::Approx(9.0));

  // Check symmetry
  CHECK(M(1, 0) == doctest::Approx(2.0));
  CHECK(M(2, 0) == doctest::Approx(3.0));
  CHECK(M(3, 0) == doctest::Approx(4.0));
  CHECK(M(2, 1) == doctest::Approx(6.0));
  CHECK(M(3, 1) == doctest::Approx(7.0));
  CHECK(M(3, 2) == doctest::Approx(9.0));
}

TEST_CASE("SymQuadric: evalMul_xt_Q_x for zero quadric")
{
  SymQuadric Q; // Zero quadric
  Eigen::Vector3d x(1.0, 2.0, 3.0);

  double result = Q.evalMul_xt_Q_x(x);
  CHECK(result == doctest::Approx(0.0));
}

TEST_CASE("SymQuadric: evalMul_xt_Q_x matches full matrix computation")
{
  SymQuadric Q;

  // Build quadric from a plane
  Eigen::Vector4d plane(1.0, 0.0, 0.0, -5.0);
  Q.addPlane(plane);

  // Test point
  Eigen::Vector3d x(2.0, 3.0, 4.0);

  // Compute using compact form
  double compact_result = Q.evalMul_xt_Q_x(x);

  // Compute using full matrix
  Eigen::Matrix4d M = Q.toMatrix();
  Eigen::Vector4d x_homog(x[0], x[1], x[2], 1.0);
  double matrix_result = x_homog.transpose() * M * x_homog;

  CHECK(compact_result == doctest::Approx(matrix_result));
}

TEST_CASE("SymQuadric: evalMul_xt_Q_x with multiple planes")
{
  SymQuadric Q;

  // Add multiple planes
  Q.addPlane(make_plane(1.0, 0.0, 0.0, -1.0));
  Q.addPlane(make_plane(0.0, 1.0, 0.0, -2.0));
  Q.addPlane(make_plane(0.0, 0.0, 1.0, -3.0));

  Eigen::Vector3d x(1.0, 2.0, 3.0);

  // Compute using compact form
  double compact_result = Q.evalMul_xt_Q_x(x);

  // Compute using full matrix
  Eigen::Matrix4d M = Q.toMatrix();
  Eigen::Vector4d x_homog(x[0], x[1], x[2], 1.0);
  double matrix_result = x_homog.transpose() * M * x_homog;

  CHECK(compact_result == doctest::Approx(matrix_result));
}

TEST_CASE("SymQuadric: evalMul_xt_Q_x measures distance to plane")
{
  SymQuadric Q;

  // Plane: z = 0 (normal [0, 0, 1], d = 0)
  Eigen::Vector4d plane(0.0, 0.0, 1.0, 0.0);
  Q.addPlane(plane);

  // Point on plane should have zero error
  Eigen::Vector3d on_plane(1.0, 2.0, 0.0);
  CHECK(Q.evalMul_xt_Q_x(on_plane) == doctest::Approx(0.0));

  // Point at distance d from plane should have error d^2
  Eigen::Vector3d off_plane(1.0, 2.0, 3.0);
  double error = Q.evalMul_xt_Q_x(off_plane);
  CHECK(error == doctest::Approx(9.0)); // 3^2 = 9
}

TEST_CASE("SymQuadric: solveMinimizer with positive definite matrix")
{
  SymQuadric Q;

  // Create a well-conditioned positive definite quadric
  // Add three orthogonal planes meeting at origin
  Q.addPlane(make_plane(1.0, 0.0, 0.0, 0.0));
  Q.addPlane(make_plane(0.0, 1.0, 0.0, 0.0));
  Q.addPlane(make_plane(0.0, 0.0, 1.0, 0.0));

  Eigen::Vector3d x_opt;
  bool success = Q.solveMinimizer(x_opt);

  CHECK(success);
  CHECK(x_opt[0] == doctest::Approx(0.0));
  CHECK(x_opt[1] == doctest::Approx(0.0));
  CHECK(x_opt[2] == doctest::Approx(0.0));
}

TEST_CASE("SymQuadric: solveMinimizer finds correct minimizer")
{
  SymQuadric Q;

  // Create planes that should meet at a specific point
  // Three planes: x=1, y=2, z=3
  Q.addPlane(make_plane(1.0, 0.0, 0.0, -1.0));
  Q.addPlane(make_plane(0.0, 1.0, 0.0, -2.0));
  Q.addPlane(make_plane(0.0, 0.0, 1.0, -3.0));

  Eigen::Vector3d x_opt;
  bool success = Q.solveMinimizer(x_opt);

  CHECK(success);
  CHECK(x_opt[0] == doctest::Approx(1.0).epsilon(0.001));
  CHECK(x_opt[1] == doctest::Approx(2.0).epsilon(0.001));
  CHECK(x_opt[2] == doctest::Approx(3.0).epsilon(0.001));

  // Verify that this is indeed the minimizer
  double error_at_min = Q.evalMul_xt_Q_x(x_opt);

  // Check nearby points have higher error
  Eigen::Vector3d nearby = x_opt + Eigen::Vector3d(0.1, 0.1, 0.1);
  double error_nearby = Q.evalMul_xt_Q_x(nearby);
  CHECK(error_nearby > error_at_min);
}

TEST_CASE("SymQuadric: solveMinimizer matches full matrix solution")
{
  SymQuadric Q;

  // Build a general quadric
  Q.addPlane(make_plane(0.6, 0.8, 0.0, -2.0));
  Q.addPlane(make_plane(0.0, 0.6, 0.8, -1.0));
  Q.addPlane(make_plane(0.8, 0.0, 0.6, -3.0));

  // Solve using compact form
  Eigen::Vector3d x_compact;
  bool success_compact = Q.solveMinimizer(x_compact);

  // Solve using full matrix
  Eigen::Matrix4d M = Q.toMatrix();
  Eigen::Matrix3d A = M.block<3, 3>(0, 0);
  Eigen::Vector3d b = M.block<3, 1>(0, 3);

  Eigen::FullPivLU<Eigen::Matrix3d> lu(A);
  bool success_matrix = (lu.rcond() >= 1e-12);
  Eigen::Vector3d x_matrix = lu.solve(-b);

  CHECK(success_compact == success_matrix);
  if(success_compact && success_matrix)
  {
    CHECK(x_compact[0] == doctest::Approx(x_matrix[0]));
    CHECK(x_compact[1] == doctest::Approx(x_matrix[1]));
    CHECK(x_compact[2] == doctest::Approx(x_matrix[2]));
  }
}

TEST_CASE("SymQuadric: solveMinimizer fails for singular matrix")
{
  SymQuadric Q;

  // Add only one plane - creates rank-1 matrix (singular)
  Q.addPlane(make_plane(1.0, 0.0, 0.0, -1.0));

  Eigen::Vector3d x_opt;
  bool success = Q.solveMinimizer(x_opt);

  CHECK_FALSE(success);
}

TEST_CASE("SymQuadric: solveMinimizer fails for singular matrix")
{
  SymQuadric Q;

  // Add only two planes, will be singular in 3D
  Q.addPlane(make_plane(1.0, 0.0, 0.0, -1.0));
  Q.addPlane(make_plane(0.0, 1.0, 0.0, -1.0));

  Eigen::Vector3d x_opt;
  bool success = Q.solveMinimizer(x_opt);

  // Check rcond directly
  Eigen::Matrix4d M = Q.toMatrix();
  Eigen::Matrix3d A = M.block<3, 3>(0, 0);
  Eigen::FullPivLU<Eigen::Matrix3d> lu(A);

  INFO("rcond: " << lu.rcond());
  // Expect rcond to be very small
  INFO("isInvertible: " << lu.isInvertible());
  INFO("rank: " << lu.rank());

  // Should fail due to near-singularity
  CHECK_FALSE(success);
}

TEST_CASE("SymQuadric: solveMinimizer fails for singular matrix")
{
  SymQuadric Q;

  // Add only two planes, will be singular in 3D
  Q.addPlane(make_plane(1.0, 0.0, 0.0, 0.0));
  Q.addPlane(make_plane(0.0, 1.0, 0.0, 0.0));
  Q.addPlane(make_plane(0.0, 1.0, 1e-8, 0.0)); // Nearly singular
  Eigen::Vector3d x_opt;
  bool success = Q.solveMinimizer(x_opt);

  // Check rcond directly
  Eigen::Matrix4d M = Q.toMatrix();
  Eigen::Matrix3d A = M.block<3, 3>(0, 0);
  Eigen::FullPivLU<Eigen::Matrix3d> lu(A);

  INFO("rcond: " << lu.rcond());
  // Expect rcond to be very small
  INFO("isInvertible: " << lu.isInvertible());
  INFO("rank: " << lu.rank());
  INFO("Solved for x_opt: " << x_opt.transpose()); // For debugging

  // Should fail due to near-singularity
  CHECK_FALSE(success);
}


TEST_CASE("SymQuadric: Integration test with plane distance formula")
{
  SymQuadric Q;

  // Create a horizontal plane at z = 5: normal [0, 0, 1], d = -5
  // Plane equation: 0*x + 0*y + 1*z - 5 = 0  =>  z = 5
  Eigen::Vector4d plane(0.0, 0.0, 1.0, -5.0);
  Q.addPlane(plane);

  // Test various points
  Eigen::Vector3d p1(0.0, 0.0, 5.0); // On plane
  Eigen::Vector3d p2(0.0, 0.0, 8.0); // Distance 3 above
  Eigen::Vector3d p3(1.0, 2.0, 2.0); // Distance 3 below

  double e1 = Q.evalMul_xt_Q_x(p1);
  double e2 = Q.evalMul_xt_Q_x(p2);
  double e3 = Q.evalMul_xt_Q_x(p3);

  CHECK(e1 == doctest::Approx(0.0));
  CHECK(e2 == doctest::Approx(9.0)); // (8-5)^2 = 9
  CHECK(e3 == doctest::Approx(9.0)); // (2-5)^2 = 9
}

TEST_CASE("SymQuadric: Stress test with random planes")
{
  SymQuadric Q;

  // Add several random unit-normal planes
  Q.addPlane(make_plane(0.577, 0.577, 0.577, -1.0));
  Q.addPlane(make_plane(0.707, 0.707, 0.0, -2.0));
  Q.addPlane(make_plane(0.0, 0.707, 0.707, -3.0));
  Q.addPlane(make_plane(0.816, 0.408, 0.408, -4.0));

  // Verify compact and full matrix give same evaluation
  Eigen::Vector3d test_point(1.5, 2.5, 3.5);

  double compact_result = Q.evalMul_xt_Q_x(test_point);

  Eigen::Matrix4d M = Q.toMatrix();
  Eigen::Vector4d x_homog(test_point[0], test_point[1], test_point[2], 1.0);
  double matrix_result = x_homog.transpose() * M * x_homog;

  CHECK(compact_result == doctest::Approx(matrix_result));

  // Try to solve for minimizer
  Eigen::Vector3d x_opt;
  bool success = Q.solveMinimizer(x_opt);

  if(success)
  {
    // Verify x_opt is a local minimum
    double error_at_opt = Q.evalMul_xt_Q_x(x_opt);

    // Test several nearby points
    for(int i = -2; i <= 2; ++i)
    {
      for(int j = -2; j <= 2; ++j)
      {
        for(int k = -2; k <= 2; ++k)
        {
          if(i == 0 && j == 0 && k == 0)
            continue;

          Eigen::Vector3d perturbed = x_opt + 0.01 * Eigen::Vector3d(i, j, k);
          double error_perturbed = Q.evalMul_xt_Q_x(perturbed);
          CHECK(error_perturbed >= error_at_opt);
        }
      }
    }
  }
}

TEST_CASE("SymQuadric: Commutative property of addition")
{
  SymQuadric Q1, Q2;

  Q1.addPlane(make_plane(1.0, 0.0, 0.0, -1.0));
  Q2.addPlane(make_plane(0.0, 1.0, 0.0, -2.0));

  SymQuadric Qa = Q1 + Q2;
  SymQuadric Qb = Q2 + Q1;

  for(int i = 0; i < 10; ++i)
  {
    CHECK(Qa.q[i] == doctest::Approx(Qb.q[i]));
  }
}

TEST_CASE("SymQuadric: Associative property of addition")
{
  SymQuadric Q1, Q2, Q3;

  Q1.addPlane(make_plane(1.0, 0.0, 0.0, -1.0));
  Q2.addPlane(make_plane(0.0, 1.0, 0.0, -2.0));
  Q3.addPlane(make_plane(0.0, 0.0, 1.0, -3.0));

  SymQuadric Qa = (Q1 + Q2) + Q3;
  SymQuadric Qb = Q1 + (Q2 + Q3);

  for(int i = 0; i < 10; ++i)
  {
    CHECK(Qa.q[i] == doctest::Approx(Qb.q[i]));
  }
}

TEST_CASE("Mesh_modifier_edge_collapse: initialize on tetrahedron")
{
  // Load the tetrahedron mesh
  Mesh_connectivity mesh;
  REQUIRE(load_test_mesh(mesh, "mesh/tetra.obj"));

  // Verify basic mesh properties
  INFO("Loaded tetrahedron: V=",
      mesh.n_active_vertices(),
      " F=",
      mesh.n_active_faces(),
      " HE=",
      mesh.n_active_half_edges());

  CHECK(mesh.n_active_vertices() == 4);
  CHECK(mesh.n_active_faces() == 4);
  CHECK(mesh.check_sanity_slowly(false));

  // Create modifier and initialize (computes quadrics and builds valid pairs)
  Mesh_modifier_edge_collapse modifier(mesh);

  // This should not throw an exception
  REQUIRE_NOTHROW(modifier.initialize());

  // After initialization, mesh should still be valid
  CHECK(mesh.check_sanity_slowly(false));
  CHECK(mesh.n_active_vertices() == 4);
  CHECK(mesh.n_active_faces() == 4);

  // Print some information about the mesh for verification
  INFO("Tetrahedron vertices:");
  for(int i = 0; i < mesh.n_total_vertices(); ++i)
  {
    auto v = mesh.vertex_at(i);
    if(v.is_active())
    {
      auto xyz = v.xyz();
      INFO("  v", i, ": (", xyz[0], ", ", xyz[1], ", ", xyz[2], ")");
    }
  }

  // Print out all the results of the quadrics in matrix form for each vertex
  INFO("Vertex quadrics:");
  for(int i = 0; i < mesh.n_total_vertices(); ++i)
  {
    auto v = mesh.vertex_at(i);
    if(v.is_active())
    {
      auto Q = modifier.vertex_quadric(i);
      Eigen::Matrix4d M = Q.toMatrix();
      INFO("  Q for v", i, ":\n", M);
      INFO(" Vertex location: ", v.xyz().transpose());
      INFO("Vertex quadric evaluation at vertex position: ", Q.evalMul_xt_Q_x(v.xyz()));
    }
  }
}

TEST_CASE("Mesh_modifier_edge_collapse: heap versioning with reinserted pairs")
{
  // Test that re-inserting pairs creates duplicates in heap but validation filters them
  // The number of valid candidates should remain the same, but version counts should increment

  Mesh_connectivity mesh;
  REQUIRE(load_test_mesh(mesh, "mesh/tetra.obj"));

  Mesh_modifier_edge_collapse modifier(mesh);
  REQUIRE_NOTHROW(modifier.initialize());

  // Get all pairs containing vertex 0
  auto v0_pairs = modifier.get_all_pairs_from_vertex(0);

  INFO("Vertex 0 has " << v0_pairs.size() << " pairs");
  REQUIRE(v0_pairs.size() > 0);

  // Count initial valid candidates
  std::vector<Mesh_modifier_edge_collapse::MergeCandidate> initial_candidates;
  Mesh_modifier_edge_collapse::MergeCandidate candidate;
  while(modifier.get_min_pair(candidate))
  {
    initial_candidates.push_back(candidate);
  }
  int initial_count = initial_candidates.size();

  INFO("Initial candidate count: " << initial_count);
  REQUIRE(initial_count > 0);

  // Re-initialize to reset heap
  modifier.initialize();

  // Re-insert all pairs from vertex 0 (this will increment their versions and add duplicates)
  for(const auto & pair : v0_pairs)
  {
    modifier.add_or_update_pair(pair.v1, pair.v2);
  }

  // Now when we get min pairs, we should get the same count but with incremented versions
  // The heap contains duplicates but get_min_pair should filter out stale entries
  std::vector<Mesh_modifier_edge_collapse::MergeCandidate> new_candidates;
  while(modifier.get_min_pair(candidate))
  {
    new_candidates.push_back(candidate);
  }

  INFO("New candidate count: " << new_candidates.size());

  // Should have same number of valid candidates despite heap duplicates
  CHECK(new_candidates.size() == initial_count);

  // Check that pairs from vertex 0 have incremented versions (version 2)
  // Other pairs should still have version 1
  int v0_pairs_found = 0;
  int other_pairs_found = 0;

  for(const auto & cand : new_candidates)
  {
    if(cand.pair.v1 == 0 || cand.pair.v2 == 0)
    {
      INFO("Vertex 0 pair (" << cand.pair.v1 << ", " << cand.pair.v2 << ") has version " << cand.version);
      CHECK(cand.version == 2);
      v0_pairs_found++;
    }
    else
    {
      CHECK(cand.version == 1);
      other_pairs_found++;
    }
  }

  CHECK(v0_pairs_found == v0_pairs.size());
  INFO("Found " << v0_pairs_found << " vertex 0 pairs with incremented versions");
  INFO("Found " << other_pairs_found << " other pairs with original versions");
}

TEST_CASE("Mesh_modifier_edge_collapse: heap versioning with manual invalidation")
{
  // Test that manually incrementing version counters invalidates pairs
  // without reinserting them, effectively removing them from valid results

  Mesh_connectivity mesh;
  REQUIRE(load_test_mesh(mesh, "mesh/tetra.obj"));

  Mesh_modifier_edge_collapse modifier(mesh);
  REQUIRE_NOTHROW(modifier.initialize());

  // Get all pairs containing vertex 0
  auto v0_pairs = modifier.get_all_pairs_from_vertex(0);

  INFO("Vertex 0 has " << v0_pairs.size() << " pairs");
  REQUIRE(v0_pairs.size() > 0);

  // Count initial valid candidates
  std::vector<Mesh_modifier_edge_collapse::MergeCandidate> initial_candidates;
  Mesh_modifier_edge_collapse::MergeCandidate candidate;
  while(modifier.get_min_pair(candidate))
  {
    initial_candidates.push_back(candidate);
  }
  int initial_count = initial_candidates.size();

  INFO("Initial candidate count: " << initial_count);
  REQUIRE(initial_count > 0);

  // Re-initialize to get a fresh heap
  modifier.initialize();

  // Now manually invalidate all pairs from vertex 0 by incrementing their version
  // WITHOUT adding them to the heap. This makes all heap entries for these pairs stale.
  for(const auto & pair : v0_pairs)
  {
    modifier.invalidate_pair(pair);
  }

  // Count how many valid candidates remain after invalidation
  std::vector<Mesh_modifier_edge_collapse::MergeCandidate> remaining_candidates;
  while(modifier.get_min_pair(candidate))
  {
    remaining_candidates.push_back(candidate);
  }
  int remaining_count = remaining_candidates.size();

  INFO("Remaining candidate count after invalidation: " << remaining_count);

  // The count should have decreased by exactly the number of vertex 0 pairs
  int expected_remaining = initial_count - v0_pairs.size();
  CHECK(remaining_count == expected_remaining);

  // Verify that none of the remaining candidates involve vertex 0
  for(const auto & cand : remaining_candidates)
  {
    CHECK(cand.pair.v1 != 0);
    CHECK(cand.pair.v2 != 0);
    INFO("Remaining pair: (" << cand.pair.v1 << ", " << cand.pair.v2 << ")");
  }

  INFO("Successfully invalidated " << v0_pairs.size() << " pairs without reinsertion");
  INFO("Verified that invalidated pairs do not appear in results");
}

TEST_CASE("Mesh_modifier_edge_collapse: get_top_n_candidates doesn't modify heap")
{
  // Test that get_top_n_candidates truly behaves like a peek - heap state is unchanged

  Mesh_connectivity mesh;
  REQUIRE(load_test_mesh(mesh, "mesh/camel_simple.obj"));

  Mesh_modifier_edge_collapse modifier(mesh);
  REQUIRE_NOTHROW(modifier.initialize());

  // First pass: drain the heap to get baseline
  std::vector<int> baseline_half_edges;
  Mesh_modifier_edge_collapse::MergeCandidate candidate;
  while(modifier.get_min_pair(candidate))
  {
    int he = modifier.get_halfedge_between_vertices(candidate.pair.v1, candidate.pair.v2);
    baseline_half_edges.push_back(he);
  }

  INFO("Baseline: heap contains " << baseline_half_edges.size() << " candidates");
  REQUIRE(baseline_half_edges.size() > 0);

  // Re-initialize to reset heap to same state
  modifier.initialize();

  // Call get_top_n_candidates (peek operation)
  std::vector<int> top_3 = modifier.get_top_n_candidates(3);
  INFO("get_top_n_candidates returned " << top_3.size() << " candidates");

  // Second pass: drain the heap again - should match baseline exactly
  std::vector<int> after_peek_half_edges;
  while(modifier.get_min_pair(candidate))
  {
    int he = modifier.get_halfedge_between_vertices(candidate.pair.v1, candidate.pair.v2);
    after_peek_half_edges.push_back(he);
  }

  INFO("After peek: heap contains " << after_peek_half_edges.size() << " candidates");

  // Verify the heap contents are identical (same count and same order)
  CHECK(after_peek_half_edges.size() == baseline_half_edges.size());
  for(size_t i = 0; i < std::min(baseline_half_edges.size(), after_peek_half_edges.size()); ++i)
  {
    CHECK(after_peek_half_edges[i] == baseline_half_edges[i]);
  }

  // Verify top_3 matches the first 3 from baseline
  int expected_count = std::min(3, static_cast<int>(baseline_half_edges.size()));
  CHECK(top_3.size() == expected_count);
  for(int i = 0; i < expected_count; ++i)
  {
    CHECK(top_3[i] == baseline_half_edges[i]);
  }
}

TEST_CASE("Mesh_modifier_edge_collapse: get_top_n_candidates cleans up stale entries")
{
  // Test that get_top_n_candidates removes stale entries while peeking
  // Even with stale entries present, heap behavior after peek should match clean heap

  Mesh_connectivity mesh;
  REQUIRE(load_test_mesh(mesh, "mesh/camel_simple.obj"));

  Mesh_modifier_edge_collapse modifier(mesh);
  REQUIRE_NOTHROW(modifier.initialize());

  // First pass: drain clean heap to get baseline (what we expect after cleanup)
  std::vector<int> baseline_half_edges;
  Mesh_modifier_edge_collapse::MergeCandidate candidate;
  while(modifier.get_min_pair(candidate))
  {
    int he = modifier.get_halfedge_between_vertices(candidate.pair.v1, candidate.pair.v2);
    baseline_half_edges.push_back(he);
  }
  INFO("Baseline (clean heap): " << baseline_half_edges.size() << " candidates");

  // Re-initialize and add stale entries
  modifier.initialize();

  // Get all pairs from vertex 0 and re-insert them to create stale entries
  auto v0_pairs = modifier.get_all_pairs_from_vertex(0);
  INFO("Adding stale entries for " << v0_pairs.size() << " vertex 0 pairs");
  for(const auto& pair : v0_pairs)
  {
    modifier.add_or_update_pair(pair.v1, pair.v2);  // Creates version 2, version 1 becomes stale
  }

  // Now heap has duplicate entries: version 1 (stale) and version 2 (valid) for v0 pairs
  // Call get_top_n_candidates - should clean up stale entries
  std::vector<int> top_candidates = modifier.get_top_n_candidates(10);
  INFO("get_top_n_candidates returned " << top_candidates.size() << " candidates");

  // Drain heap after peek - should match baseline (stale entries cleaned up)
  std::vector<int> after_cleanup_half_edges;
  while(modifier.get_min_pair(candidate))
  {
    int he = modifier.get_halfedge_between_vertices(candidate.pair.v1, candidate.pair.v2);
    after_cleanup_half_edges.push_back(he);
  }
  INFO("After cleanup: " << after_cleanup_half_edges.size() << " candidates");

  // Should have same count as baseline (all stale entries removed)
  CHECK(after_cleanup_half_edges.size() == baseline_half_edges.size());

  // Should match baseline exactly (same order, same half-edges)
  for(size_t i = 0; i < std::min(baseline_half_edges.size(), after_cleanup_half_edges.size()); ++i)
  {
    CHECK(after_cleanup_half_edges[i] == baseline_half_edges[i]);
  }
}

TEST_CASE("Mesh_modifier_edge_collapse: get_top_n_candidates with n larger than heap")
{
  // Test that requesting more candidates than available works correctly

  Mesh_connectivity mesh;
  REQUIRE(load_test_mesh(mesh, "mesh/tetra.obj"));

  Mesh_modifier_edge_collapse modifier(mesh);
  REQUIRE_NOTHROW(modifier.initialize());

  // Count total candidates first
  Mesh_modifier_edge_collapse::MergeCandidate candidate;
  int total_count = 0;
  while(modifier.get_min_pair(candidate))
  {
    total_count++;
  }
  INFO("Mesh has " << total_count << " total candidates");

  // Re-initialize and request more than available
  modifier.initialize();
  std::vector<int> all_candidates = modifier.get_top_n_candidates(100);

  INFO("Requested 100, got " << all_candidates.size());

  // Should return all available candidates
  CHECK(all_candidates.size() == total_count);

  // Verify heap is still intact
  int count_after = 0;
  while(modifier.get_min_pair(candidate))
  {
    count_after++;
  }

  CHECK(count_after == total_count);
  INFO("Heap still contains " << count_after << " valid candidates after peek");
}

TEST_CASE("Mesh_modifier_edge_collapse: tetrahedron all edges illegal collapse")
{
  // Test that all edges in a tetrahedron fail the is_legal_collapse test
  // A tetrahedron is the minimal 3D simplicial complex - no edge can be collapsed
  // without destroying the 3D structure

  Mesh_connectivity mesh;
  REQUIRE(load_test_mesh(mesh, "mesh/tetra.obj"));

  Mesh_modifier_edge_collapse modifier(mesh);
  REQUIRE_NOTHROW(modifier.initialize());

  // Verify it's a tetrahedron
  CHECK(mesh.n_active_vertices() == 4);
  CHECK(mesh.n_active_faces() == 4);

  INFO("Testing tetrahedron edge collapse legality");

  // Collect all unique edges by iterating through all half-edges
  std::set<std::pair<int, int>> unique_edges;
  for(int he_idx = 0; he_idx < mesh.n_total_half_edges(); ++he_idx)
  {
    auto he = mesh.half_edge_at(he_idx);
    if(he.is_active())
    {
      int v1 = he.origin().index();
      int v2 = he.dest().index();
      // Store edge in canonical form (smaller index first)
      if(v1 > v2)
        std::swap(v1, v2);
      unique_edges.insert({v1, v2});
    }
  }

  INFO("Found " << unique_edges.size() << " unique edges in tetrahedron");

  // A tetrahedron has 6 edges (complete graph on 4 vertices)
  CHECK(unique_edges.size() == 6);

  // Check that ALL edges are illegal to collapse
  int illegal_count = 0;
  for(const auto & edge : unique_edges)
  {
    bool is_legal = modifier.is_legal_collapse(edge.first, edge.second);
    INFO("Edge (" << edge.first << ", " << edge.second << ") is_legal: " << is_legal);

    if(!is_legal)
    {
      illegal_count++;
    }

    // All edges should be illegal
    CHECK_FALSE(is_legal);
  }

  CHECK(illegal_count == 6);
  INFO("All " << illegal_count << " edges in tetrahedron are correctly marked as illegal to collapse");
}

TEST_CASE("Mesh_modifier_edge_collapse: pyramid_square base diagonal not collapsable")
{
  // Test the square pyramid mesh
  // The base has 4 vertices forming a square that's triangulated with a diagonal
  // According to the requirement, there should be exactly one edge on the base
  // that is not collapsable (the diagonal edge used for triangulation)

  Mesh_connectivity mesh;
  REQUIRE(load_test_mesh(mesh, "mesh/pyramid_square.obj"));

  Mesh_modifier_edge_collapse modifier(mesh);
  REQUIRE_NOTHROW(modifier.initialize());

  // Verify it's a square pyramid
  INFO("Pyramid: V=" << mesh.n_active_vertices()
       << " F=" << mesh.n_active_faces()
       << " HE=" << mesh.n_active_half_edges());

  CHECK(mesh.n_active_vertices() == 5);
  CHECK(mesh.n_active_faces() == 6); // 4 sides + 2 base triangles

  // Identify the apex and base vertices by z-coordinate
  int apex_vertex = -1;
  std::vector<int> base_vertices;

  for(int i = 0; i < mesh.n_total_vertices(); ++i)
  {
    auto v = mesh.vertex_at(i);
    if(v.is_active())
    {
      auto xyz = v.xyz();
      INFO("v" << i << ": (" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << ")");

      // Apex should be at z=1, base vertices at z=0
      if(std::abs(xyz[2] - 1.0) < 1e-6)
      {
        apex_vertex = i;
      }
      else if(std::abs(xyz[2] - 0.0) < 1e-6)
      {
        base_vertices.push_back(i);
      }
    }
  }

  REQUIRE(apex_vertex >= 0);
  REQUIRE(base_vertices.size() == 4);
  INFO("Apex vertex: " << apex_vertex);
  INFO("Base vertices: " << base_vertices[0] << ", " << base_vertices[1]
       << ", " << base_vertices[2] << ", " << base_vertices[3]);

  // Collect all edges on the base (edges between base vertices only)
  std::set<std::pair<int, int>> base_edges;
  for(int he_idx = 0; he_idx < mesh.n_total_half_edges(); ++he_idx)
  {
    auto he = mesh.half_edge_at(he_idx);
    if(he.is_active())
    {
      int v1 = he.origin().index();
      int v2 = he.dest().index();

      // Check if both vertices are in the base
      bool v1_is_base = std::find(base_vertices.begin(), base_vertices.end(), v1) != base_vertices.end();
      bool v2_is_base = std::find(base_vertices.begin(), base_vertices.end(), v2) != base_vertices.end();

      if(v1_is_base && v2_is_base)
      {
        // Store edge in canonical form
        if(v1 > v2)
          std::swap(v1, v2);
        base_edges.insert({v1, v2});
      }
    }
  }

  INFO("Found " << base_edges.size() << " edges on the base");

  // The base is triangulated, so it should have 5 edges:
  // 4 perimeter edges + 1 diagonal
  CHECK(base_edges.size() == 5);

  // Check collapse legality for all base edges
  int legal_count = 0;
  int illegal_count = 0;

  for(const auto & edge : base_edges)
  {
    bool is_legal = modifier.is_legal_collapse(edge.first, edge.second);
    INFO("Base edge (" << edge.first << ", " << edge.second << ") is_legal: " << is_legal);

    if(is_legal)
    {
      legal_count++;
    }
    else
    {
      illegal_count++;
    }
  }

  INFO("Base edges: " << legal_count << " legal, " << illegal_count << " illegal");

  // According to the requirement, exactly one edge should be illegal (not collapsable)
  CHECK(illegal_count == 1);
  CHECK(legal_count == 4);
}

TEST_CASE("Mesh_modifier_edge_collapse: collapse edge on pyramid_square")
{
  // Test that collapse_edge works correctly on the pyramid_square mesh
  // Load the mesh, get the minimum error pair, and attempt to collapse it

  Mesh_connectivity mesh;
  REQUIRE(load_test_mesh(mesh, "mesh/pyramid_square.obj"));

  Mesh_modifier_edge_collapse modifier(mesh);
  REQUIRE_NOTHROW(modifier.initialize());

  // Record initial mesh state
  int initial_vertices = mesh.n_active_vertices();
  int initial_faces = mesh.n_active_faces();
  int initial_half_edges = mesh.n_active_half_edges();

  INFO("Initial mesh state: V=" << initial_vertices
       << " F=" << initial_faces
       << " HE=" << initial_half_edges);

  CHECK(initial_vertices == 5);
  CHECK(initial_faces == 6);

  // Get the minimum error pair
  Mesh_modifier_edge_collapse::MergeCandidate candidate;
  bool found = modifier.get_min_pair(candidate);

  REQUIRE(found);
  INFO("Min pair: (" << candidate.pair.v1 << ", " << candidate.pair.v2 << ")");
  INFO("Error: " << candidate.error);
  INFO("Optimal position: (" << candidate.x_opt[0] << ", "
       << candidate.x_opt[1] << ", " << candidate.x_opt[2] << ")");

  // Attempt to collapse the edge
  bool collapse_success = modifier.collapse_edge(candidate);

  INFO("Collapse " << (collapse_success ? "succeeded" : "failed"));

  if(collapse_success)
  {
    // Verify mesh is still valid after collapse
    CHECK(mesh.check_sanity_slowly(false));

    // Check that vertex count decreased by 1
    int final_vertices = mesh.n_active_vertices();
    INFO("Final vertices: " << final_vertices);
    CHECK(final_vertices == initial_vertices - 1);

    // Check that face count decreased (should lose at least 2 faces for interior edge)
    int final_faces = mesh.n_active_faces();
    INFO("Final faces: " << final_faces);
    CHECK(final_faces < initial_faces);

    // Print final mesh state
    INFO("Final mesh state: V=" << final_vertices
         << " F=" << final_faces
         << " HE=" << mesh.n_active_half_edges());
  }
  else
  {
    // If collapse failed, mesh should remain unchanged
    INFO("Collapse was illegal, mesh should remain unchanged");
    CHECK(mesh.n_active_vertices() == initial_vertices);
    CHECK(mesh.n_active_faces() == initial_faces);
  }
}
