# Minimesh Build System

## Quick Start

### Option 1: Using Robust Build Scripts (Recommended)
```bash
# Configure and build debug version
./configure.sh debug
./build.sh debug

# Or configure and build release version
./configure.sh release
./build.sh release

# Clean rebuild
./build.sh debug clean
```

### Option 2: Manual CMake (For Reference)
```bash
# Debug build
cd build-dbg
cmake .. \
    -DCMAKE_BUILD_TYPE=debug \
    -DEIGEN3_DIR="$(pwd)/../third-party/eigen" \
    -DFREEGLUT_DIR="$(pwd)/../third-party/freeglut/bin-opensuse-gcc731/debug" \
    -DGLUI_DIR="$(pwd)/../third-party/glui/bin-opensuse-gcc731/debug"
make -j

# Release build
cd build-opt
cmake .. \
    -DCMAKE_BUILD_TYPE=relwithdebinfo \
    -DEIGEN3_DIR="$(pwd)/../third-party/eigen" \
    -DFREEGLUT_DIR="$(pwd)/../third-party/freeglut/bin-opensuse-gcc731/release" \
    -DGLUI_DIR="$(pwd)/../third-party/glui/bin-opensuse-gcc731/release"
make -j
```

## Running the Applications

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

## Dependencies

### System Requirements (Linux/WSL)
- CMake 3.11+
- GCC/G++ compiler
- OpenGL development libraries
- FreeGLUT development libraries

### Third-party Libraries (Included)
- Eigen 3.x (Linear algebra)
- FreeGLUT (OpenGL utilities)
- GLUI (GUI toolkit)

## VSCode Integration

The repository includes pre-configured VSCode settings:
- **IntelliSense**: C++ code navigation and autocomplete
- **Debugging**: GDB integration for both CLI and GUI
- **Build Tasks**: CMake integration

## Troubleshooting

### GUI Display Issues (WSL)
If the GUI shows menus but no graphics viewport:
```bash
# Try software rendering
LIBGL_ALWAYS_SOFTWARE=1 ./bin/minimeshgui ../mesh/cube.obj

# Install mesa utilities
sudo apt install mesa-utils
```

### Build Errors
1. Ensure third-party dependencies are extracted: `unzip third-party.zip`
2. Clean and reconfigure: `rm -rf build-dbg/* && ./configure.sh debug`
3. Check system OpenGL packages: `dpkg -l | grep -E "(libgl|freeglut)"`

## Development Workflow

1. **Setup**: Run `./configure.sh debug` once
2. **Development**: Edit source files in `src/`
3. **Build**: Run `./build.sh debug` to compile changes
4. **Debug**: Use VSCode debugger or `gdb bin/minimeshcli`
5. **Test**: Run executables from `build-dbg/bin/`

## File Structure

- `configure.sh`, `build.sh` - Robust build scripts
- `build-dbg/`, `build-opt/` - Build directories
- `src/minimesh/` - Source code
- `third-party/` - Dependencies
- `mesh/` - Example mesh files
- `.vscode/` - VSCode configuration