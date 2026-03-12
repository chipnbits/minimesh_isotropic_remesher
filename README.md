# 3D Mesh Uniform Isotropic Remeshing

**Author:** Simon Ghyselincks

## Overview Summary

This project implements a uniform isotropic remeshing algorithm for 3D meshes based upon a synthesis of techniques described by Botsch and Kobbelt (2004) and Tanja Munz (2015). The remeshing is a parameterization-free approach that relies solely upon local mesh operations to iteratively adjust a mesh towards a targeted edge length with desirable isotropic properties, while preserving important features.

### Operation of Remesher
The remeshing algorithm can be viewed and operated either through the GUI application or the CLI application included in the MiniMesh framework. The GUI contains new radio buttons to toggle the remeshing mode to `on`. When activated, the important feature edges are identified with blue highlights, and a target edge length can be set in the numerical field. The remeshing is done iteratively, pressing the `Run Single Pass` button with the desired target edge length.

In the CLI application, the remesher is set up to iteratively remesh towards either a target edge or vertex count. The relevant function calls for the remesher are:

  ```cpp
  // Remesh to target edge count (approximate)
  void remesh_to_target_edge_count(int target_edge_count, int iterations, const double rel_error_tol = .02);

  // Remesh to target vertex count (approximate)
  void remesh_to_target_vertex_count(int target_vertex_count, int iterations, const double rel_error_tol = .02);
  ```

  These calls will iteratively adjust the mesh to reach the desired edge or vertex count, within a specified relative error tolerance or maximum number of iterations.


### References

- Mario Botsch and Leif Kobbelt. **A Remeshing Approach to Multiresolution Modeling.**  
  *Proceedings of the Eurographics/ACM SIGGRAPH Symposium on Geometry Processing*, 2004.

- Munz Tanja. **Curvature Adaptive Remeshing.**  
  Master’s Thesis, Bournemouth University, 2015.  
  https://nccastaff.bournemouth.ac.uk/jmacey/MastersProject/MSc15/08Tanja/report.pdf


