```c++
Mesh_connectivity::Defragmentation_maps maps;

// 1) Stop using any live iterators/indices into the mesh.
//    (They will be invalid after defrag.)

// 2) Defragment in place (build maps, rebuild mesh, swap in)
mesh.defragment_in_place(maps);

// 3) Remap ALL external arrays/sets that you own using maps.
//    Example helpers shown below.

// Per-vertex arrays:
remap_dense_vector(vertex_quadrics, maps.old2new_vertices);
remap_dense_vector(vertex_flags,   maps.old2new_vertices);

// Per-half-edge arrays:
remap_dense_vector(he_flags,       maps.old2new_half_edges);

// Per-face arrays:
remap_dense_vector(face_labels,    maps.old2new_faces);

// 4) Any structures that store indices inside their payloads
//    (like a queue of edge ids) must be rebuilt or remapped.
//    Easiest & safest: rebuild from the current mesh.
rebuild_edge_priority_queue_from_mesh(mesh);

// 5) Resume work. All indices you use from now on are compact.
```

For the case of coloring in the mesh, the compact form must be used for the coloring matrix to be valid:

```c++
Mesh_connectivity::Defragmentation_maps maps;
mesh.compute_defragmention_maps(maps); // leaves mesh as-is

// Vertices 
// 3 x n_active_vertices, columns are RGB in [0,1]
Eigen::MatrixXf Vcols(3, mesh.n_active_vertices());
Vcols.setZero();

for (int vid = 0; vid < mesh.n_total_vertices(); ++vid) {
    int nv = maps.old2new_vertices[vid];
    if (nv == Mesh_connectivity::invalid_index) continue;  // skip inactive
    Eigen::Vector3f c = color_for_vertex_old_id(vid);      // your logic
    Vcols.col(nv) = c;
}
mesh_buffer.set_vertex_colors(Vcols);

//Faces 
Eigen::MatrixXf Fcols(3, mesh.n_active_faces());
Fcols.setZero();

for (int fid = 0; fid < mesh.n_total_faces(); ++fid) {
    int nf = maps.old2new_faces[fid];
    if (nf == Mesh_connectivity::invalid_index) continue;
    Eigen::Vector3f c = color_for_face_old_id(fid);        // your logic
    Fcols.col(nf) = c;
}
mesh_buffer.set_face_colors(Fcols);
mesh_buffer.update_opengl_buffers();
```
