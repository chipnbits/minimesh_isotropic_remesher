#include <doctest.h>

#include <Eigen/Core>
#include <minimesh/core/mohe/mesh_connectivity.hpp>
#include <minimesh/core/mohe/mesh_io.hpp>
#include <minimesh/core/mohe/mesh_analysis.hpp>

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
} // namespace

TEST_CASE("Mesh utils - is_triangular_mesh") {
    mohecore::Mesh_connectivity mesh;

    // Cube should be triangular (made of triangles)
    REQUIRE(load_test_mesh(mesh, "mesh/cube.obj"));
    CHECK(mohecore::analysis::is_triangular_mesh(mesh));

    // Camel should also be triangular (made of triangles)
    REQUIRE(load_test_mesh(mesh, "mesh/camel.obj"));
    CHECK(mohecore::analysis::is_triangular_mesh(mesh));

    INFO("Triangular mesh check completed");
}

TEST_CASE("Mesh utils - count_connected_components") {
    mohecore::Mesh_connectivity mesh;

    // Single connected cube should have 1 component
    REQUIRE(load_test_mesh(mesh, "mesh/cube.obj"));
    CHECK(mohecore::analysis::count_connected_components(mesh) == 1);

    // Octopus mesh has 4 connected components (body + 3 separate pieces)
    REQUIRE(load_test_mesh(mesh, "mesh/octopus.obj"));
    int octopus_components = mohecore::analysis::count_connected_components(mesh);
    INFO("Expected 4, got ", octopus_components);
    CHECK(octopus_components == 4);
}