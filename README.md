# 2D Mesh Parameterization 

**Author:** Simon Ghyselincks
**Student Number:** 12145124

## Overview

This project implements mesh parameterization. Two algorithms are provided:
1. Fixed Boundary Parameterization using cotangent weights.
2. Least Squares Conformal Maps (LSCM).

The easiest way to get started is through the GUI, wich has radio buttons to select the desired algorithm and initiate the parameterization process. The CLI version runs the same parameterization with save to file.
## Running the Applications

The CLI and GUI can be run using the following aliased commands:

- **CLI**: `meshcli <filepath> <algorithm>`
  - Example: `meshcli ./mesh/cat.obj lscm`
  - The `<algorithm>` argument can be either `fixed` (for Fixed Boundary Parameterization) or `lscm` (for Least Squares Conformal Maps)

- **GUI**: `meshgui <filepath>`
  - Example: `meshgui ./hw4_mesh/camel_head.obj`

## CLI Usage
The CLI will read the mesh file and then apply either the fixed boundary or the lscm version. Fixed boundary is using the cotan weight scheme, while lscm is using the least squares conformal maps approach. The UV coordiantes are written back into a 2D mesh that can be inspected using the GUI. CLI renders are found in the `exports/` folder.

## GUI Usage
Run the GUI with the path to mesh file. Select the desired parameterization algorithm using the radio buttons then click "Parameterize Mesh" to execute. The resulting UV parameterization will be displayed in the viewport.