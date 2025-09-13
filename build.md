Latest version of CMake extension and some build configurations are probvided through the `CMakePresets.json` and the `launch.json` files.

There is a slower debug build and a faster optimized build. Since this is a development project, the debug build is the default. An automatic build and run with the cube.obj can be done using the CMake extension `Ctrl+Shift+B` from anywhere in the project. The default rebuild method is incremental and not a full clean rebuild.

New presets accessed through the bottom menu bar in VSCode can be added through the `CMakePresets.json` file.

### Command Line Interface
```bash
cd build-dbg
./bin/minimeshcli ../mesh/cube.obj
```

### Graphical User Interface
```bash
cd build-dbg
./bin/minimeshgui ../mesh/cube.obj
```

## Build Outputs

- **Executables**: `build-dbg/bin/` or `build-opt/bin/`
  - `minimeshcli` - Command line mesh processing
  - `minimeshgui` - Interactive GUI application

- **Libraries**:
  - `libminimeshcore.a` - Core mesh data structures
  - `libminimeshviz.a` - Visualization components

## File Structure

- `build-dbg/`, `build-opt/` - Build directories
- `src/minimesh/` - Source code
- `third-party/` - Dependencies
- `mesh/` - Example mesh files
- `.vscode/` - VSCode configuration