#include <benchmark/benchmark.h>
#include <cmath>
#include <vector>
#include <random> // Replacing Boost for immediate dependency-free running

// --- HARDCOPY OF YOUR LEGACY STRUCTS FOR ISOLATION ---
struct OptionData {
    double K; double T; double r; double sig;
    int type; // 1 == call, -1 == put
    double myPayOffFunction(double S) {
        if (type == 1) return std::max(S - K, 0.0);
        else return std::max(K - S, 0.0);
    }
};

// Simple standard normal generator to replace BoostNormal temporarily
double get_standard_normal() {
    static std::mt19937 generator(42);
    static std::normal_distribution<double> dist(0.0, 1.0);
    return dist(generator);
}

// --- THE LOGIC FROM TestMC.cpp ADAPTED FOR BENCHMARK ---
static void BM_LegacyMonteCarlo(benchmark::State& state) {
    // Setup Data (Moved out of the timed loop)
    OptionData myOption;
    myOption.T = 30.0;
    myOption.K = 100.0;
    myOption.sig = 0.3;
    myOption.r = 0.08;
    myOption.type = 1;
    double S_0 = 100.0;
    
    long N = 100; // Time steps
    
    // The Benchmark Loop
    for (auto _ : state) {
        // We simulate ONE path per iteration for micro-benchmarking, 
        // OR we can simulate a batch. Let's simulate a small batch (e.g. 1000)
        // so the timer isn't dominated by overhead.
        
        double price = 0.0;
        long NSim = 1000; 
        double k = myOption.T / double(N);
        double sqrk = std::sqrt(k);
        double VOld, VNew;

        for (long i = 1; i <= NSim; ++i) {
            VOld = S_0;
            for (long index = 1; index < N; ++index) {
                double dW = get_standard_normal();
                
                // Your Legacy Drift/Diffusion Logic inline
                // drift = (data->r)*X
                // diffusion = data->sig * pow(X, betaCEV) (assuming beta=1.0)
                
                double drift = myOption.r * VOld;
                double diffusion = myOption.sig * VOld; // Optimized pow(X, 1) to X
                
                VNew = VOld + (k * drift) + (sqrk * diffusion * dW);
                VOld = VNew;
            }
            price += myOption.myPayOffFunction(VNew);
        }
        
        // Prevent compiler optimization from skipping the work
        benchmark::DoNotOptimize(price);
    }
}

// Register the function as a benchmark
BENCHMARK(BM_LegacyMonteCarlo);

BENCHMARK_MAIN();