#!/bin/bash
# Configure, build, and run benchmarks

set -e  # Exit on error

SCRIPT_DIR="$(dirname "$0")"

echo "=== Step 1: Configure ==="
"$SCRIPT_DIR/configure.sh"

echo ""
echo "=== Step 2: Build ==="
"$SCRIPT_DIR/build.sh"

echo ""
echo "=== Step 3: Run Benchmarks ==="
"$SCRIPT_DIR/run_benchmarks.sh" "$@"

