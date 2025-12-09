# 2D Mesh Parameterization 

**Author:** Simon Ghyselincks
**Student Number:** 12145124

## Overview

This project implements mesh parameterization

## Running the Applications

The CLI and GUI can be run using the following aliased commands:

- **CLI**: `meshcli <filepath> <algorithm>`
  - Example: `meshcli ./mesh/cat.obj lscm`
  - The `<algorithm>` argument can be either `fixed` (for Fixed Boundary Parameterization) or `lscm` (for Least Squares Conformal Maps)

- **GUI**: `meshgui <filepath>`
  - Example: `meshgui ./exports/cat_cli.obj`

## CLI Usage
The CLI will read the mesh file and then apply either the fixed boundary or the lscm version. Fixed boundary is using the cotan weight scheme, while lscm is using the least squares conformal maps approach. The UV coordiantes are written back into a 2D mesh that can be inspected using the GUI. CLI renders are found in the `exports/` folder.

## GUI Usage
Simply run the GUI with the path to the processed exported mesh file. The GUI will visualize the 2D parameterized mesh. The boundary vertices are counted and highlighted in red to ensure that they are being mapped correctly.