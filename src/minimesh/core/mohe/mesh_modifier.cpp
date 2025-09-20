#include <minimesh/core/mohe/mesh_modifier.hpp>
#include <minimesh/core/util/assert.hpp>
#include <minimesh/core/mohe/mesh_analysis.hpp>
#include <queue>

namespace minimesh
{
namespace mohecore
{


//
// Given two vertices, this function return the index of the half-edge going from v0 to v1.
// Returns -1 if no half-edge exists between the two vertices.
//
int Mesh_modifier::get_halfedge_between_vertices(const int v0, const int v1)
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

bool Mesh_modifier::flip_edge(const int he_index)
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

bool Mesh_modifier::divide_edge(const int he_index, const double weight)
{
	// Validate input index
	if (he_index < 0 || he_index >= mesh().n_total_half_edges()) return false;

	// Get the half-edge and its twin
	auto he = mesh().half_edge_at(he_index);
	if (!he.is_active()) return false;

	auto twin = he.twin();
	if (!twin.is_active()) return false;

	// Get the two vertices of the edge
	auto origin_vertex = he.origin();
	auto dest_vertex = twin.origin();  // Twin's origin is he's destination

	// Compute weighted position: weight * origin + (1-weight) * dest
	Eigen::Vector3d new_pos = weight * origin_vertex.xyz() + (1.0 - weight) * dest_vertex.xyz();

	// Create new vertex at weighted position
	auto new_vertex = mesh().add_vertex();
	new_vertex.data().xyz = new_pos;

	// Create two new half-edges to split the original edge
	auto new_he = mesh().add_half_edge();
	auto new_twin = mesh().add_half_edge();

	// Store original connectivity info before modification
	int he_next = he.data().next;
	int he_face = he.data().face;

	int twin_next = twin.data().next;
	int twin_face = twin.data().face;

	// Update original half-edge: now goes from origin to new vertex
	he.data().next = new_he.index();
	// prev stays the same
	// face stays the same
	// origin stays the same
	he.data().twin = new_twin.index();

	// Update new half-edge: goes from new vertex to destination
	new_he.data().next = he_next;
	new_he.data().prev = he.index();
	new_he.data().twin = twin.index();
	new_he.data().face = he_face;
	new_he.data().origin = new_vertex.index();

	// Update original twin: now goes from dest to new vertex
	twin.data().next = new_twin.index();
	// prev stays the same
	// face stays the same
	// origin stays the same (dest_vertex)
	twin.data().twin = new_he.index();

	// Update new twin: goes from new vertex to origin
	new_twin.data().next = twin_next;
	new_twin.data().prev = twin.index();
	new_twin.data().twin = he.index();
	new_twin.data().face = twin_face;
	new_twin.data().origin = new_vertex.index();

	// Update next/prev pointers of adjacent half-edges
	auto new_he_next_idx = new_he.data().next;
	if (new_he_next_idx != mesh().invalid_index) {
		new_he.next().data().prev = new_he.index();
	}
	auto new_twin_next_idx = new_twin.data().next;
	if (new_twin_next_idx != mesh().invalid_index) {
		mesh().half_edge_at(new_twin_next_idx).data().prev = new_twin.index();
	}

	// Update vertex half-edge pointers
	new_vertex.data().half_edge = new_he.index();

	// Origin and dest vertices retained the same half-edge pointers
	// Face retained the same half-edge pointers

	return true;
}

bool Mesh_modifier::subdivide_loop()
{
	
	// Placeholder
	return true;
}


} // end of mohecore
} // end of minimesh
