#!/bin/bash

# Robust CMake configuration script for minimesh
# Usage: ./configure.sh [debug|release]

set -e  # Exit on error

# Get the directory where this script is located (project root)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Default to debug build
BUILD_TYPE="${1:-debug}"

if [ "$BUILD_TYPE" = "debug" ]; then
    BUILD_DIR="build-dbg"
    CMAKE_BUILD_TYPE="debug"
elif [ "$BUILD_TYPE" = "release" ]; then
    BUILD_DIR="build-opt"
    CMAKE_BUILD_TYPE="relwithdebinfo"
else
    echo "Usage: $0 [debug|release]"
    exit 1
fi

echo "=== Configuring minimesh for $BUILD_TYPE build ==="
echo "Project root: $SCRIPT_DIR"
echo "Build directory: $BUILD_DIR"

# Create build directory if it doesn't exist
mkdir -p "$BUILD_DIR"

# Change to build directory
cd "$BUILD_DIR"

# Configure with CMake using absolute paths to avoid path resolution issues
cmake .. \
    -DCMAKE_BUILD_TYPE="$CMAKE_BUILD_TYPE" \
    -DEIGEN3_DIR="$SCRIPT_DIR/third-party/eigen" \
    -DFREEGLUT_DIR="$SCRIPT_DIR/third-party/freeglut/bin-opensuse-gcc731/$BUILD_TYPE" \
    -DGLUI_DIR="$SCRIPT_DIR/third-party/glui/bin-opensuse-gcc731/$BUILD_TYPE"

echo ""
echo "=== Configuration complete! ==="
echo "To build: cd $BUILD_DIR && make -j"
echo "To run CLI: cd $BUILD_DIR && ./bin/minimeshcli ../mesh/cube.obj"
echo "To run GUI: cd $BUILD_DIR && ./bin/minimeshgui ../mesh/cube.obj"