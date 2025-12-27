#include "option_data.hpp" 
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

template <class T> void print(const std::vector<T>& myList)
{  // A generic print function for vectors
	
	std::cout << std::endl << "Size of vector is " << myList.size() << "\n[";

	// We must use a const iterator here, otherwise we get a compiler error.
	typename std::vector<T>::const_iterator i;
	for (i = myList.begin(); i != myList.end(); ++i)
	{
			std::cout << *i << ",";

	}

	std::cout << "]\n";
}

// generic function for standard deviation of call output prices
template <typename U>
U std_dev(const std::vector<U>& prices, const double r, const double T) {

	// compute sample standard deviation using the computational formula:
	// SD = sqrt((1/(n-1)) * (sum(x_i^2) - (sum(x_i))^2/n))
	U sum_of_squares = 0;
	U sum = 0;
	for (size_t i = 0; i < prices.size(); i++) {
		sum_of_squares += prices[i]*prices[i];
		sum += prices[i];
	}
	
	U n = prices.size();
	U variance = (sum_of_squares - (sum * sum / n)) / (n - 1);
	
	// Standard deviation of undiscounted payoffs
	U sd = sqrt(variance);
	
	// Apply discounting: SD(c*X) = |c| * SD(X) where c = exp(-r*T)
	sd *= exp(-r*T);

    return sd;
}

// generic function for standard error of call output prices
template <typename U>
U std_err(const std::vector<U>& prices, const double r, const double T) {
	return std_dev(prices, r, T) / sqrt(prices.size());
}

namespace SDEDefinition
{ // Defines drift + diffusion + data
	OptionData* data;				// The data for the option MC

	double drift(double t, double X)
	{ // Drift term
		return (data->r)*X; // r - D
	}
	
	double diffusion(double t, double X)
	{ // Diffusion term
	
		double betaCEV = 1.0;
		return data->sig * pow(X, betaCEV);
		
	}

	double diffusionDerivative(double t, double X)
	{ // Diffusion term, needed for the Milstein method
	
		double betaCEV = 1.0;
		return 0.5 * (data->sig) * (betaCEV) * pow(X, 2.0 * betaCEV - 1.0);
	}
} // End of namespace

// Function to run a single MC pricing test
double runMCPricing(OptionData& option, double S_0, long N, long NSim, std::mt19937_64& rng) {
	std::vector<double> x(N + 1);
	double dt = option.T / N;
	for(int i = 0; i <= N; ++i) {
		x[i] = i * dt;
	}

	double k = option.T / double(N);
	double sqrk = sqrt(k);
	
	std::normal_distribution<double> myNormal(0.0, 1.0);
	
	using namespace SDEDefinition;
	SDEDefinition::data = &option;

	double price = 0.0;
	int coun = 0; // Number of times S hits origin

	for (long i = 1; i <= NSim; ++i) {
		if ((i/50000) * 50000 == i) {
			std::cout << "  Simulation " << i << "/" << NSim << std::endl;
		}

		double VOld = S_0;
		double VNew;
		
		for (unsigned long index = 1; index < x.size(); ++index) {
			double dW = myNormal(rng);
			VNew = VOld + (k * drift(x[index-1], VOld))
						+ (sqrk * diffusion(x[index-1], VOld) * dW);
			VOld = VNew;

			if (VNew <= 0.0) coun++;
		}
		
		double tmp = option.myPayOffFunction(VNew);
		price += tmp / double(NSim);
	}
	
	// Discount the average price
	price *= exp(-option.r * option.T);
	
	if (coun > 0) {
		std::cout << "  Warning: Origin hit " << coun << " times" << std::endl;
	}
	
	return price;
}

int main()
{
	std::cout << "Monte Carlo Option Pricer - Testing Accuracy\n";
	std::cout << std::string(60, '=') << "\n\n";
	
	// Test parameters
	long N = 1000;      // Time steps
	long NSim = 500000; // Simulations (increased for better accuracy)
	
	std::cout << "Parameters:" << std::endl;
	std::cout << "  Time steps (N): " << N << std::endl;
	std::cout << "  Simulations (NSim): " << NSim << std::endl;
	std::cout << std::endl;
	
	std::mt19937_64 rng(42); // Fixed seed for reproducibility
	
	// Test cases: T, K, sig, r, S_0, expected_price, type_name
	struct TestCase {
		double T, K, sig, r, S_0, expected;
		int type; // 1 = Call, -1 = Put
		std::string name;
	};
	
	std::vector<TestCase> testCases = {
		{0.25, 65.0, 0.30, 0.08, 60.0, 2.13337, 1, "Batch 1 Call"},
		{0.25, 65.0, 0.30, 0.08, 60.0, 5.84628, -1, "Batch 1 Put"},
		{1.0, 100.0, 0.20, 0.00, 100.0, 7.96557, 1, "Batch 2 Call"},
		{1.0, 100.0, 0.20, 0.00, 100.0, 7.96557, -1, "Batch 2 Put"}
	};
	
	// Run tests
	for (const auto& test : testCases) {
		std::cout << "Testing " << test.name << std::endl;
		std::cout << "  Parameters: S0=" << test.S_0 << ", K=" << test.K 
				  << ", T=" << test.T << ", r=" << test.r << ", σ=" << test.sig << std::endl;
		std::cout << "  Expected: " << test.expected << std::endl;
		
		OptionData option;
		option.K = test.K;
		option.T = test.T;
		option.r = test.r;
		option.sig = test.sig;
		option.type = test.type;
		
		double computed = runMCPricing(option, test.S_0, N, NSim, rng);
		double error = std::abs(computed - test.expected);
		double rel_error = error / test.expected * 100.0;
		
		std::cout << "  Computed: " << computed << std::endl;
		std::cout << "  Absolute Error: " << error << std::endl;
		std::cout << "  Relative Error: " << rel_error << "%" << std::endl;
		
		// Check 2 decimal accuracy (0.005 tolerance)
		bool passed = (error < 0.005);
		std::cout << "  Status: " << (passed ? "PASS ✓" : "FAIL ✗") 
				  << " (2 decimal accuracy)" << std::endl;
		std::cout << std::endl;
	}
	
	std::cout << std::string(60, '=') << std::endl;

	return 0;
}