#!/bin/bash
# Build the project

set -e  # Exit on error

cd "$(dirname "$0")/.."  # Go to project root

# Detect number of CPU cores
NCPU=$(nproc)

echo "Building with $NCPU cores..."
cd build
cmake --build . -j$NCPU

echo "Build complete!"

