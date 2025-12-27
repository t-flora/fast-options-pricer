#!/bin/bash
# Configure the project with CMake
# Usage: ./configure.sh [BUILD_TYPE]
# BUILD_TYPE: Debug, Release, RelWithDebInfo, or MinSizeRel (default: Release)

set -e  # Exit on error

cd "$(dirname "$0")/.."  # Go to project root

# Get build type from argument or default to Release
BUILD_TYPE="${1:-Release}"

echo "Cleaning build directory..."
rm -rf build/*

echo "Configuring with CMake (BUILD_TYPE=$BUILD_TYPE)..."
cd build
cmake .. -DCMAKE_BUILD_TYPE="$BUILD_TYPE"

echo "Configuration complete!"