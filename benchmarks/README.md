# Benchmarks

This folder contains performance benchmarks for the fast-options-pricer project using [Google Benchmark](https://github.com/google/benchmark).

## Overview

The benchmarks measure the performance of various option pricing implementations to:
- Establish baseline performance metrics
- Compare different optimization strategies
- Track performance regressions
- Guide optimization efforts

## Current Benchmarks

### `bench_baseline.cpp`
Benchmarks the legacy Monte Carlo option pricing algorithm. This serves as the baseline for measuring improvements in future optimizations.

**What it measures:**
- Monte Carlo simulation with 1,000 paths per iteration
- 100 time steps per path
- European call option pricing using geometric Brownian motion
- Uses standard normal random number generation (C++ `<random>` library)

**Test parameters:**
- Spot price (S₀): 100.0
- Strike price (K): 100.0
- Time to maturity (T): 30.0 years
- Risk-free rate (r): 0.08
- Volatility (σ): 0.3
- Option type: Call (type = 1)

## Building and Running

### Quick Start

From the project root:

```bash
# Full build (configure, build, and run benchmarks)
./scripts/full_build.sh

# Or run individually:
./scripts/configure.sh        # Configure with CMake
./scripts/build.sh            # Build the project
./scripts/run_benchmarks.sh   # Run benchmarks
```

### Build with Debug Symbols

If you want to profile or debug the benchmarks:

```bash
./scripts/configure.sh RelWithDebInfo
./scripts/build.sh
./scripts/run_benchmarks.sh
```

### Running Specific Benchmarks

Pass Google Benchmark flags to the run script:

```bash
# Filter specific benchmarks
./scripts/run_benchmarks.sh --benchmark_filter=Legacy

# Run for a specific time
./scripts/run_benchmarks.sh --benchmark_min_time=5.0

# Show additional statistics
./scripts/run_benchmarks.sh --benchmark_repetitions=10

# JSON output for analysis
./scripts/run_benchmarks.sh --benchmark_format=json --benchmark_out=results.json
```

## Adding New Benchmarks

1. Create a new `.cpp` file in this folder (e.g., `bench_optimized.cpp`)

2. Include the benchmark header and implement your benchmark:

```cpp
#include <benchmark/benchmark.h>

static void BM_YourBenchmark(benchmark::State& state) {
    // Setup (outside the loop)
    
    for (auto _ : state) {
        // Code to benchmark
        benchmark::DoNotOptimize(result);
    }
}

BENCHMARK(BM_YourBenchmark);
BENCHMARK_MAIN();
```

3. Update `CMakeLists.txt` to add the new benchmark:

```cmake
add_executable(YourBenchmark.out bench_optimized.cpp)
target_link_libraries(YourBenchmark.out PRIVATE benchmark::benchmark)
```

4. Build and run:

```bash
./scripts/build.sh
./build/benchmarks/YourBenchmark.out
```

## Understanding the Output

Google Benchmark provides output like this:

```
-------------------------------------------------------------------
Benchmark                         Time             CPU   Iterations
-------------------------------------------------------------------
BM_LegacyMonteCarlo           12345 ns        12340 ns        56789
```

- **Time**: Wall clock time per iteration
- **CPU**: CPU time per iteration
- **Iterations**: Number of times the benchmark was run to get stable results

## Tips for Good Benchmarks

1. **Isolate the code**: Only benchmark the critical path you want to measure
2. **Avoid optimization**: Use `benchmark::DoNotOptimize()` to prevent compiler optimizations from skipping your code
3. **Setup outside the loop**: Move initialization code outside `for (auto _ : state)`
4. **Realistic workloads**: Use realistic parameters that match production use cases
5. **Multiple iterations**: Run with `--benchmark_repetitions` to account for variance
6. **Consistent environment**: Close other applications and disable CPU frequency scaling for more stable results

## Dependencies

- **Google Benchmark**: Automatically fetched via CMake's FetchContent
- **C++17 or later**: Required for the codebase

## References

- [Google Benchmark User Guide](https://github.com/google/benchmark/blob/main/docs/user_guide.md)
- [Reducing Benchmark Variance](https://github.com/google/benchmark/blob/main/docs/reducing_variance.md)

