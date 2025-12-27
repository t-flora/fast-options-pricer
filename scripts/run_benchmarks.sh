#!/bin/bash
# Run benchmarks

set -e  # Exit on error

cd "$(dirname "$0")/.."  # Go to project root

echo "Running benchmarks..."
./build/benchmarks/RunBenchmarks.out "$@"