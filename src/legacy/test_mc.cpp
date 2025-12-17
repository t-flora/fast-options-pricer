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

	// compute difference between square of prices and mean square of prices
	U sum_of_squares = 0;
	U mean_square = 0; // compute mean square of prices - first sum, then divide by size
	for (size_t i = 0; i < prices.size(); i++) {
		sum_of_squares += prices[i]*prices[i];
		mean_square += prices[i];
	}
	
	mean_square = mean_square*mean_square / prices.size(); // divide by size

	// compute standard deviation
	U sd = sqrt(((sum_of_squares-mean_square)/(prices.size()-1)) * exp(-r*T));

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

int main()
{
	std::cout <<  "1 factor MC with explicit Euler\n";
	// OptionData myOption;
	OptionData myOption, option2;
	// T K sig r S b, call, put
	// { {0.25, 65, 0.30, 0.08, 60}, 2.13337, 5.84628 },	// batch 1
	// { {1.0, 100, 0.2, 0.0, 100}, 7.96557, 7.96557 },		// batch 2
	// batch 1
	// myOption.K = 65.0;
	// myOption.T = 0.25;
	// myOption.r = 0.08;
	// myOption.sig = 0.3;
	// myOption.type = -1;	// Put -1, Call +1
	// double S_0 = 60;

	// batch 2
	myOption.K = 100.0;
	myOption.T = 1.0;
	myOption.r = 0.0;
	myOption.sig = 0.2;
	myOption.type = 1;	// Put -1, Call +1
	double S_0 = 100;

	long N = 1000;
	std::cout << "Number of subintervals in time: " << N << std::endl;
	// std::cin >> N;

    std::vector<double> x(N + 1);
    double dt = myOption.T / N;
    for(int i = 0; i <= N; ++i) {
        x[i] = i * dt;
    }

	// Create the basic SDE (Context class)
	double VOld = S_0;
	double VNew;


	// V2 mediator stuff
	long NSim = 100000;
	std::cout << "Number of simulations: " << NSim << std::endl;
	// std::cin >> NSim;

	double k = myOption.T / double (N);
	double sqrk = sqrt(k);

	// Normal random number
	double dW;
	double price = 0.0;	// Option price

	// NormalGenerator is a base class
    std::mt19937_64 rng(42);
	std::normal_distribution<double> myNormal(0.0, 1.0);

	using namespace SDEDefinition;
	SDEDefinition::data = &myOption;

	std::vector<double> prices(NSim);
	int coun = 0; // Number of times S hits origin

	// A.
	for (long i = 1; i <= NSim; ++i)
	{ // Calculate a path at each iteration
			
		if ((i/10000) * 10000 == i)
		{// Give status after each 1000th iteration
				std::cout << i << std::endl;
		}

		VOld = S_0;
		for (unsigned long index = 1; index < x.size(); ++index)
		{

			// Create a random number
			dW = myNormal(rng);
				
			// The FDM (in this case explicit Euler)
			VNew = VOld  + (k * drift(x[index-1], VOld))
						+ (sqrk * diffusion(x[index-1], VOld) * dW);

			VOld = VNew;

			// Spurious values
			if (VNew <= 0.0) coun++;
		}
			
		double tmp = myOption.myPayOffFunction(VNew);
		prices[i-1] = tmp;
		price += (tmp)/double(NSim);
	}
	
	// D. Finally, discounting the average price
	price *= exp(-myOption.r * myOption.T);

	// Cleanup not needed; std::normal_distribution is not allocated with new

	std::cout << "Values for batch 2 call" << std::endl;
	std::cout << "Price, after discounting: " << price << ", " << std::endl;
	std::cout << "Number of times origin is hit: " << coun << std::endl;
	std::cout << "Standard deviation: " << std_dev(prices, myOption.r, myOption.T) << endl;
	std::cout << "Standard error: " << std_err(prices, myOption.r, myOption.T) << endl;

	// file output handling - append result to file (this was changed between runs of different batches)
	std::ofstream out_file("out_b2_call.txt", std::ios::app);  // Open in append mode
	
	// check if file is empty
	out_file.seekp(0, std::ios::end);
	if (out_file.tellp() == 0) { // if file is empty, add headers
		out_file << "N,NSIM,SD,SE" << std::endl;
	}
	
	// append the data
	out_file << N << ","
			<< NSim << ","
			<< std_dev(prices, myOption.r, myOption.T) << "," 
			<< std_err(prices, myOption.r, myOption.T) << std::endl;

	return 0;
}