#include <doctest.h>

#include <Eigen/Core>
#include <minimesh/core/mohe/mesh_connectivity.hpp>
#include <minimesh/core/mohe/mesh_io.hpp>
#include <minimesh/core/mohe/mesh_modifier.hpp>

#include <iostream>
#include <string>

using namespace minimesh;

namespace {
static bool load_test_mesh(mohecore::Mesh_connectivity& mesh, const std::string& filename) {
    mohecore::Mesh_io io(mesh);
    try {
        io.read_obj_general(filename);
        return mesh.check_sanity_slowly(false);
    } catch (...) {
        std::cerr << "Failed to load mesh: " << filename << std::endl;
        return false;
    }
}

// Helper function to validate basic subdivision properties
[[maybe_unused]] static bool validate_subdivision_properties(mohecore::Mesh_connectivity& original,
                                           mohecore::Mesh_connectivity& subdivided) {
    // Basic sanity checks
    if (!subdivided.check_sanity_slowly(false)) return false;

    // Subdivision should increase vertex and face count
    if (subdivided.n_active_vertices() <= original.n_active_vertices()) return false;
    if (subdivided.n_active_faces() <= original.n_active_faces()) return false;

    return true;
}
} // namespace

TEST_CASE("Subdivision basic setup") {
    mohecore::Mesh_connectivity mesh;
    REQUIRE(load_test_mesh(mesh, "mesh/cube.obj"));

    // Store original counts for comparison
    const int original_vertices = mesh.n_active_vertices();
    const int original_faces = mesh.n_active_faces();
    const int original_edges = mesh.n_active_half_edges() / 2;

    CHECK(original_vertices == 8);
    CHECK(original_faces == 12);
    CHECK(original_edges == 18);

    CAPTURE(original_vertices);

    // TODO: Add actual subdivision operation here
    // Example placeholder:
    // mohecore::Mesh_modifier modifier(mesh);
    // modifier.subdivide_loop(); // or whatever subdivision method you implement

    INFO("Original mesh loaded successfully for subdivision testing");
}

TEST_CASE("Triangle edge division test") {
    mohecore::Mesh_connectivity mesh;
    mohecore::Mesh_modifier modifier(mesh);

    // Create a simple triangle mesh programmatically
    std::vector<double> xyz = {
        0.0, 0.0, 0.0,  // vertex 0
        1.0, 0.0, 0.0,  // vertex 1
        0.5, 1.0, 0.0   // vertex 2
    };
    std::vector<int> triangle_verts = {0, 1, 2};

    mesh.build_from_triangles(xyz, triangle_verts);

    // Verify initial mesh is valid
    REQUIRE(mesh.check_sanity_slowly(false));
    REQUIRE(mesh.n_active_vertices() == 3);
    REQUIRE(mesh.n_active_faces() == 1);
    REQUIRE(mesh.n_active_half_edges() == 6); // 3 edges * 2 half-edges each

    // Get a half-edge on the triangle to divide
    auto face = mesh.face_at(0);
    auto he_to_divide = face.half_edge();
    int he_index = he_to_divide.index();

    // Store original vertex positions for validation
    auto origin_pos = he_to_divide.origin().xyz();
    auto dest_pos = he_to_divide.dest().xyz();

    // Divide the edge at midpoint (default weight = 0.5)
    bool success = modifier.divide_edge(he_index);
    REQUIRE(success);

    // Verify mesh is still valid after division
    CHECK(mesh.check_sanity_slowly(false));

    // Check that we have one more vertex
    CHECK(mesh.n_active_vertices() == 4);

    // Check that we have two more half-edges (edge was split into two)
    CHECK(mesh.n_active_half_edges() == 8);

    // Find the new vertex and verify its position
    bool found_new_vertex = false;
    Eigen::Vector3d expected_midpoint = (origin_pos + dest_pos) * 0.5;

    for (int v = 0; v < mesh.n_total_vertices(); ++v) {
        auto vertex = mesh.vertex_at(v);
        if (vertex.is_active()) {
            Eigen::Vector3d pos = vertex.xyz();
            if ((pos - expected_midpoint).norm() < 1e-10) {
                found_new_vertex = true;
                break;
            }
        }
    }

    CHECK(found_new_vertex);
}

TEST_CASE("Weighted edge division test") {
    mohecore::Mesh_connectivity mesh;
    mohecore::Mesh_modifier modifier(mesh);

    // Create a simple triangle mesh
    std::vector<double> xyz = {
        0.0, 0.0, 0.0,  // vertex 0
        2.0, 0.0, 0.0,  // vertex 1
        1.0, 2.0, 0.0   // vertex 2
    };
    std::vector<int> triangle_verts = {0, 1, 2};

    mesh.build_from_triangles(xyz, triangle_verts);
    REQUIRE(mesh.check_sanity_slowly(false));

    // Get a half-edge to divide
    auto face = mesh.face_at(0);
    auto he_to_divide = face.half_edge();
    int he_index = he_to_divide.index();

    auto origin_pos = he_to_divide.origin().xyz();
    auto dest_pos = he_to_divide.dest().xyz();

    // Divide edge with weight 0.25 (closer to origin)
    bool success = modifier.divide_edge(he_index, 0.25);
    REQUIRE(success);

    // Verify mesh is still valid
    CHECK(mesh.check_sanity_slowly(false));
    CHECK(mesh.n_active_vertices() == 4);

    // Find the new vertex and verify its weighted position
    bool found_weighted_vertex = false;
    Eigen::Vector3d expected_pos = 0.25 * origin_pos + 0.75 * dest_pos;

    for (int v = 0; v < mesh.n_total_vertices(); ++v) {
        auto vertex = mesh.vertex_at(v);
        if (vertex.is_active()) {
            Eigen::Vector3d pos = vertex.xyz();
            if ((pos - expected_pos).norm() < 1e-10) {
                found_weighted_vertex = true;
                break;
            }
        }
    }

    CHECK(found_weighted_vertex);
}

TEST_CASE("Edge division edge cases") {
    mohecore::Mesh_connectivity mesh;
    mohecore::Mesh_modifier modifier(mesh);

    // Create a simple triangle mesh
    std::vector<double> xyz = {
        0.0, 0.0, 0.0,  // vertex 0
        1.0, 0.0, 0.0,  // vertex 1
        0.5, 1.0, 0.0   // vertex 2
    };
    std::vector<int> triangle_verts = {0, 1, 2};

    mesh.build_from_triangles(xyz, triangle_verts);
    REQUIRE(mesh.check_sanity_slowly(false));

    auto face = mesh.face_at(0);
    auto he_to_divide = face.half_edge();
    int he_index = he_to_divide.index();

    // Test edge cases for weight parameter
    SUBCASE("Weight = 0.0 (at destination)") {
        bool success = modifier.divide_edge(he_index, 0.0);
        REQUIRE(success);
        CHECK(mesh.check_sanity_slowly(false));
        CHECK(mesh.n_active_vertices() == 4);
    }

    SUBCASE("Weight = 1.0 (at origin)") {
        bool success = modifier.divide_edge(he_index, 1.0);
        REQUIRE(success);
        CHECK(mesh.check_sanity_slowly(false));
        CHECK(mesh.n_active_vertices() == 4);
    }

    SUBCASE("Invalid half-edge index") {
        bool success = modifier.divide_edge(-1);
        CHECK_FALSE(success);
        // Mesh should remain unchanged
        CHECK(mesh.n_active_vertices() == 3);
    }
}

TEST_CASE("Loop subdivision mesh properties") {
    // Test Loop subdivision properties on different meshes
    const char* test_meshes[] = {"mesh/cube.obj", "mesh/tetra.obj", "mesh/camel_simple.obj"};

    for (const char* mesh_file : test_meshes) {
        SUBCASE(mesh_file) {
            mohecore::Mesh_connectivity mesh;
            REQUIRE(load_test_mesh(mesh, mesh_file));

            // Store original mesh properties
            const int V_orig = mesh.n_active_vertices();
            const int E_orig = mesh.n_active_half_edges() / 2;  // Half-edges to edges
            const int F_orig = mesh.n_active_faces();

            CAPTURE(V_orig);
            CAPTURE(E_orig);
            CAPTURE(F_orig);

            // Apply Loop subdivision
            mohecore::Mesh_modifier modifier(mesh);
            bool success = modifier.subdivide_loop();
            REQUIRE(success);

            // Verify mesh is still valid after subdivision
            CHECK(mesh.check_sanity_slowly(false));

            // Get new mesh properties
            const int V_new = mesh.n_active_vertices();
            const int E_new = mesh.n_active_half_edges() / 2;
            const int F_new = mesh.n_active_faces();

            CAPTURE(V_new);
            CAPTURE(E_new);
            CAPTURE(F_new);

            // Verify Loop subdivision properties:
            // - faces is n × 4 (each triangle becomes 4 triangles)
            CHECK(F_new == F_orig * 4);

            // - vertices is V + (1/2)|E| (original vertices + edge midpoints)
            CHECK(V_new == V_orig + E_orig);

            // - edges: Each original edge becomes 2 edges, plus 3 new interior edges per face
            // E_new = 2*E_orig + 3*F_orig
            CHECK(E_new == 2 * E_orig + 3 * F_orig);

            // Additional sanity checks
            CHECK(V_new > V_orig);
            CHECK(E_new > E_orig);
            CHECK(F_new > F_orig);
        }
    }
}

TEST_CASE("Loop subdivision properties - simple triangle") {
    mohecore::Mesh_connectivity mesh;
    mohecore::Mesh_modifier modifier(mesh);

    // Create a simple triangle mesh
    std::vector<double> xyz = {
        0.0, 0.0, 0.0,  // vertex 0
        1.0, 0.0, 0.0,  // vertex 1
        0.5, 1.0, 0.0   // vertex 2
    };
    std::vector<int> triangle_verts = {0, 1, 2};

    mesh.build_from_triangles(xyz, triangle_verts);
    REQUIRE(mesh.check_sanity_slowly(false));

    // Store original properties for single triangle
    const int V_orig = 3;  // vertices
    const int E_orig = 3;  // edges
    const int F_orig = 1;  // faces

    CHECK(mesh.n_active_vertices() == V_orig);
    CHECK(mesh.n_active_half_edges() / 2 == E_orig);
    CHECK(mesh.n_active_faces() == F_orig);

    // Apply Loop subdivision
    bool success = modifier.subdivide_loop();
    REQUIRE(success);

    // Verify mesh is still valid
    CHECK(mesh.check_sanity_slowly(false));

    // Get new properties
    const int V_new = mesh.n_active_vertices();
    const int E_new = mesh.n_active_half_edges() / 2;
    const int F_new = mesh.n_active_faces();

    // Check Loop subdivision formulas for single triangle:
    // - faces: 1 × 4 = 4
    CHECK(F_new == 4);

    // - vertices: 3 + 3 = 6 (original vertices + edge midpoints)
    CHECK(V_new == 6);

    // - edges: Each original edge becomes 2 edges, plus 3 new interior edges per face
    // E_new = 2*3 + 3*1 = 9
    CHECK(E_new == 9);
}
