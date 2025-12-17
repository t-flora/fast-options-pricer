add_executable(RunBenchmarks.out bench_baseline.cpp)

# Link against Google Benchmark
target_link_libraries(RunBenchmarks.out PRIVATE benchmark::benchmark)

# Allow us to include the legacy headers
target_include_directories(RunBenchmarks.out PRIVATE ${CMAKE_SOURCE_DIR}/src/legacy)