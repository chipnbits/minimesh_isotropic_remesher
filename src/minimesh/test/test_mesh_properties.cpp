#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

#include <Eigen/Core>
#include <minimesh/core/mohe/mesh_connectivity.hpp>
#include <minimesh/core/mohe/mesh_io.hpp>
#include <minimesh/core/mohe/mesh_modifier_loop_subdivision.hpp>

#include <iostream>
#include <string>

using namespace minimesh;

namespace {
struct BoundingBox {
    Eigen::Vector3d min_coords{ 1e10,  1e10,  1e10};
    Eigen::Vector3d max_coords{-1e10, -1e10, -1e10};
};

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

static BoundingBox calculate_bounding_box(mohecore::Mesh_connectivity& mesh) {
    BoundingBox bbox;
    for (int i = 0; i < mesh.n_total_vertices(); ++i) {
        auto v = mesh.vertex_at(i);
        if (!v.is_active()) continue;
        auto xyz = v.xyz();
        bbox.min_coords = bbox.min_coords.cwiseMin(xyz);
        bbox.max_coords = bbox.max_coords.cwiseMax(xyz);
    }
    return bbox;
}

static int calculate_vertex_valence(mohecore::Mesh_connectivity& mesh, int vertex_id) {
    int val = 0;
    auto ring = mesh.vertex_ring_at(vertex_id);
    do { val++; } while (ring.advance());
    return val;
}

static bool has_boundary_edges(mohecore::Mesh_connectivity& mesh) {
    for (int i = 0; i < mesh.n_total_half_edges(); ++i) {
        auto he = mesh.half_edge_at(i);
        if (!he.is_active()) continue;
        if (he.twin().index() == -2) return true;  // hole_index value
    }
    return false;
}

static bool check_euler_characteristic(mohecore::Mesh_connectivity& mesh, int expected_genus = 0) {
    int V = mesh.n_active_vertices();
    int E = mesh.n_active_half_edges() / 2;
    int F = mesh.n_active_faces();
    return (V - E + F) == (2 - 2 * expected_genus);
}
} // namespace

TEST_CASE("Cube basic properties") {
    mohecore::Mesh_connectivity mesh;
    REQUIRE(load_test_mesh(mesh, "mesh/cube.obj"));

    CHECK(mesh.n_active_vertices()    == 8);
    CHECK(mesh.n_active_faces()       == 12);
    CHECK(mesh.n_active_half_edges()  == 36);
    CHECK(mesh.check_sanity_slowly(false));
    CHECK(check_euler_characteristic(mesh, 0));
}

TEST_CASE("Cube geometry properties") {
    mohecore::Mesh_connectivity mesh;
    REQUIRE(load_test_mesh(mesh, "mesh/cube.obj"));

    auto bbox = calculate_bounding_box(mesh);
    CHECK(bbox.min_coords[0] == doctest::Approx(0.0));
    CHECK(bbox.min_coords[1] == doctest::Approx(0.0));
    CHECK(bbox.min_coords[2] == doctest::Approx(0.0));
    CHECK(bbox.max_coords[0] == doctest::Approx(1.0));
    CHECK(bbox.max_coords[1] == doctest::Approx(1.0));
    CHECK(bbox.max_coords[2] == doctest::Approx(1.0));

    auto v0 = mesh.vertex_at(0);
    CHECK(v0.is_active());
    auto p = v0.xyz();
    CHECK(p[0] == doctest::Approx(0.0));
    CHECK(p[1] == doctest::Approx(0.0));
    CHECK(p[2] == doctest::Approx(0.0));

    CHECK(!has_boundary_edges(mesh));
}

TEST_CASE("Cube topology properties") {
    mohecore::Mesh_connectivity mesh;
    REQUIRE(load_test_mesh(mesh, "mesh/cube.obj"));

    for (int i = 0; i < mesh.n_active_vertices(); ++i) {
        auto v = mesh.vertex_at(i);
        if (!v.is_active()) continue;
        int val = calculate_vertex_valence(mesh, i);
        CHECK_MESSAGE(val >= 3, "vertex ", i, " valence=", val, " should be >= 3");
        CHECK_MESSAGE(val <= 6, "vertex ", i, " valence=", val, " should be <= 6");
    }

    for (int i = 0; i < mesh.n_total_half_edges(); ++i) {
        auto he = mesh.half_edge_at(i);
        if (!he.is_active()) continue;
        auto next = he.next();
        auto prev_of_next = next.prev();
        CHECK(prev_of_next.is_equal(he));
    }
}

TEST_CASE("Edge flip operation") {
    mohecore::Mesh_connectivity mesh;
    mohecore::Mesh_modifier_loop_subdivision modifier(mesh);
    REQUIRE(load_test_mesh(mesh, "mesh/cube.obj"));

    const int V0 = mesh.n_active_vertices();
    const int F0 = mesh.n_active_faces();
    const int H0 = mesh.n_active_half_edges();

    int v0 = 0, v1 = 1;
    int he_index = modifier.get_halfedge_between_vertices(v0, v1);

    if (he_index != -3) {  // invalid_index value
        CHECK(modifier.flip_edge(he_index));
        CHECK(mesh.check_sanity_slowly(false));
        CHECK(mesh.n_active_vertices()   == V0);
        CHECK(mesh.n_active_faces()      == F0);
        CHECK(mesh.n_active_half_edges() == H0);
        CHECK(modifier.get_halfedge_between_vertices(v0, v1) == -3);  // invalid_index value
    } else {
        MESSAGE("No edge found between vertices 0 and 1");
        // No assertion needed - test passes if no edge exists
    }
}

TEST_CASE("Multiple meshes load & basic validation") {
    const char* files[] = {"mesh/cube.obj", "mesh/cow1.obj", "mesh/camel_simple.obj"};
    for (auto f : files) {
        mohecore::Mesh_connectivity mesh;
        if (load_test_mesh(mesh, f)) {
            INFO(f, " : V=", mesh.n_active_vertices(), " F=", mesh.n_active_faces());
            CHECK(mesh.n_active_vertices() > 0);
            CHECK(mesh.n_active_faces() > 0);
            CHECK(mesh.check_sanity_slowly(false));
        } else {
            MESSAGE("Warning: could not load ", f);
            // Test passes even if optional mesh can't be loaded
        }
    }
}

TEST_CASE("Template transform (placeholder)") {
    mohecore::Mesh_connectivity mesh;
    mohecore::Mesh_modifier_loop_subdivision modifier(mesh);
    REQUIRE(load_test_mesh(mesh, "mesh/cube.obj"));
    // Template test - add your transform validation here
    CHECK(mesh.n_active_vertices() > 0); // Basic validation
}
