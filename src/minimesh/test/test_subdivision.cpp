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

TEST_CASE("Triangle subdivision test") {
    mohecore::Mesh_connectivity mesh;
    mohecore::Mesh_modifier modifier(mesh);

    // Create a simple triangle mesh programmatically for testing
    // TODO: Implement triangle creation or load a simple triangle mesh

    WARN("Triangle subdivision test not implemented yet");
    // Placeholder test that passes
    CHECK(true);
}

TEST_CASE("Quad subdivision test") {
    mohecore::Mesh_connectivity mesh;
    mohecore::Mesh_modifier modifier(mesh);

    // TODO: Test subdivision on quad-based meshes

    WARN("Quad subdivision test not implemented yet");
    // Placeholder test that passes
    CHECK(true);
}

TEST_CASE("Subdivision edge cases") {
    mohecore::Mesh_connectivity mesh;
    mohecore::Mesh_modifier modifier(mesh);

    // TODO: Test subdivision on edge cases like:
    // - Single triangle
    // - Non-manifold geometry
    // - Meshes with boundaries

    WARN("Subdivision edge cases test not implemented yet");
    // Placeholder test that passes
    CHECK(true);
}

TEST_CASE("Subdivision preserves topology") {
    mohecore::Mesh_connectivity mesh;
    REQUIRE(load_test_mesh(mesh, "mesh/cube.obj"));

    // TODO: Verify that subdivision preserves topological properties
    // - Genus should remain the same
    // - No holes should be introduced
    // - Mesh should remain manifold

    WARN("Subdivision topology preservation test not implemented yet");
    // Placeholder test that passes
    CHECK(true);
}