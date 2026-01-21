#include <minimesh/core/mohe/mesh_analysis.hpp>
#include <queue>
#include <vector>

namespace minimesh
{
namespace mohecore
{
namespace analysis
{

bool is_triangular_mesh(Mesh_connectivity& mesh)
{
	for(int face_id = 0; face_id < mesh.n_total_faces(); ++face_id)
	{
		auto face = mesh.face_at(face_id);
		if(face.is_active() && face.n_vertices() != 3)
		{
			return false;
		}
	}
	return true;
}

int count_connected_components(Mesh_connectivity& mesh)
{
	if(mesh.n_active_vertices() == 0) return 0;

	// Track visited vertices
	std::vector<bool> visited_vertices(mesh.n_total_vertices(), false);
	int components = 0;

	// Traverse all vertices
	for(int v = 0; v < mesh.n_total_vertices(); ++v)
	{
		// Bypass if not active or already visited
		auto vertex = mesh.vertex_at(v);
		if(!vertex.is_active() || visited_vertices[v]) continue;

		// Start new component and perform BFS
		components++;
		std::queue<int> queue;
		queue.push(v);
		visited_vertices[v] = true;

		while(!queue.empty())
		{
			int current_v = queue.front();
			queue.pop();

			// Visit all neighbors through vertex ring
			auto ring = mesh.vertex_ring_at(current_v);
			do {
				auto neighbor_vertex = ring.half_edge().origin();
				int neighbor_id = neighbor_vertex.index();

				if(neighbor_vertex.is_active() && !visited_vertices[neighbor_id])
				{
					visited_vertices[neighbor_id] = true;
					queue.push(neighbor_id);
				}
			} while(ring.advance());
		}
	}

	return components;
}

	int vertex_valence(Mesh_connectivity& mesh, int vertex_id)
	{
		int valence = 0;
		auto ring = mesh.vertex_ring_at(vertex_id);
		do { valence++; } while (ring.advance());
		return valence;
	}

	bool vertex_is_boundary(Mesh_connectivity& mesh, int vertex_id)
	{
		auto ring = mesh.vertex_ring_at(vertex_id);
		return ring.reset_boundary();
	}
} // end of analysis
} // end of mohecore
} // end of minimesh