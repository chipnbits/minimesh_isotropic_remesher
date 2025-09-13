#!/bin/bash

# Robust build script for minimesh
# Usage: ./build.sh [debug|release] [clean]

set -e  # Exit on error

# Get the directory where this script is located (project root)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Parse arguments
BUILD_TYPE="${1:-debug}"
CLEAN_BUILD="${2}"

if [ "$BUILD_TYPE" = "debug" ]; then
    BUILD_DIR="build-dbg"
elif [ "$BUILD_TYPE" = "release" ]; then
    BUILD_DIR="build-opt"
else
    echo "Usage: $0 [debug|release] [clean]"
    exit 1
fi

echo "=== Building minimesh ($BUILD_TYPE) ==="

# Check if build directory exists
if [ ! -d "$BUILD_DIR" ]; then
    echo "Build directory $BUILD_DIR does not exist. Running configure first..."
    ./configure.sh "$BUILD_TYPE"
fi

cd "$BUILD_DIR"

# Clean build if requested
if [ "$CLEAN_BUILD" = "clean" ]; then
    echo "Cleaning previous build..."
    make clean || true
fi

# Build the project
echo "Building with $(nproc) parallel jobs..."
make -j$(nproc)

echo ""
echo "=== Build complete! ==="
echo "Executables are in: $BUILD_DIR/bin/"
ls -la bin/

echo ""
echo "To run CLI: cd $BUILD_DIR && ./bin/minimeshcli ../mesh/cube.obj"
echo "To run GUI: cd $BUILD_DIR && ./bin/minimeshgui ../mesh/cube.obj"