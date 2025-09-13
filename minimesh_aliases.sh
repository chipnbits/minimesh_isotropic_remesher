#!/bin/bash

# Minimesh Project Aliases
# Source this file to add convenient aliases for minimesh CLI and GUI
# Usage: source ./minimesh_aliases.sh

# Get the absolute path to the project directory
MINIMESH_PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Define paths to executables (prioritize debug build, fallback to release)
if [ -f "$MINIMESH_PROJECT_DIR/build-dbg/bin/minimeshcli" ]; then
    MINIMESH_CLI="$MINIMESH_PROJECT_DIR/build-dbg/bin/minimeshcli"
    MINIMESH_GUI="$MINIMESH_PROJECT_DIR/build-dbg/bin/minimeshgui"
    echo "Using debug build executables"
elif [ -f "$MINIMESH_PROJECT_DIR/build-opt/bin/minimeshcli" ]; then
    MINIMESH_CLI="$MINIMESH_PROJECT_DIR/build-opt/bin/minimeshcli"
    MINIMESH_GUI="$MINIMESH_PROJECT_DIR/build-opt/bin/minimeshgui"
    echo "Using release build executables"
else
    echo "Error: No minimesh executables found. Please build the project first."
    return 1
fi

# CLI alias - runs from project directory so relative paths work
alias minimeshcli='cd "$MINIMESH_PROJECT_DIR" && "$MINIMESH_CLI"'

# GUI alias with optional filename parameter
minimeshgui() {
    cd "$MINIMESH_PROJECT_DIR"
    if [ $# -eq 0 ]; then
        # No arguments - just run GUI
        "$MINIMESH_GUI"
    else
        # With filename argument
        "$MINIMESH_GUI" "$@"
    fi
}

# Export the function so it's available in subshells
export -f minimeshgui


echo "Minimesh aliases loaded successfully!"
echo ""
echo "Available commands:"
echo "  minimeshcli           - Run the CLI tool"
echo "  minimeshgui [file]    - Run the GUI tool (optionally with a mesh file)"
echo ""
echo "Example usage:"
echo "  minimeshcli"
echo "  minimeshgui mesh/bunny.obj"
echo "  minimeshgui"