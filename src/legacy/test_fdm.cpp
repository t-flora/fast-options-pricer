// A simple mesher on a 1d domain. We divide
// an interval into J+1 mesh points, J-1 of which
// are internal mesh points.

#include <vector>
#include <iostream>
#include <string>
#include "option_data.hpp"
#include <cmath>
#include <algorithm>
using namespace std;


class Mesher
{
public:
		double a, b;	// In space

		Mesher()
		{
			a =0.0; b = 1.0;
		}

		Mesher (double A, double B)
		{ // Describe the domain of integration

			a = A;
			b = B;
		}

		std::vector<double> xarr(int J)
		{
			// NB Full array (includes end points)

			double h = (b - a) / double (J);
			
			int size = J+1;
			int start = 1;

			std::vector<double> result(size, start);
			result[0] = a;

			for (unsigned int j = 1; j < result.size(); j++)
			{
				result[j] = result[j-1] + h;
			}

			return result;
		}

		std::vector<double> Xarr(int J)
		{ // Return as an STL vector

			// NB Full array (includes end points)

			double h = (b - a) / double (J);
			
			int size = J+1;
			int start = 1;

			std::vector<double> result(size, start);
			result[0] = a;

			for (unsigned int j = 1; j < result.size(); j++)
			{
				result[j] = result[j-1] + h;
			}

			return result;
		}

};

namespace ParabolicIBVP
{

	// Coefficients of PDE equation
	double (*sigma)(double x, double t); // Diffusion term
	double (*mu)(double x, double t);	 // Convection term
	double (*b)(double x, double t);	 // Free term
	double (*f)(double x, double t);	 // The forcing term term


	// (Dirichlet) boundary conditions
	double (*BCL)(double t);	 // The left-hand boundary condition
	double (*BCR)(double t);	 // The right-hand boundary condition

	// Initial condition
	double (*IC)(double x);		// The condition at time t = 0

}

using namespace ParabolicIBVP;

class FDM
{
public:
		// Solver solver;

	//	Vector<double, long> A, B, C;	// LHS coefficients at level n+1
		std::vector<double> a, bb, c;	// LHS coefficients at level n
		std::vector<double> RHS;		// Inhomogeneous term
		
		std::vector<double> vecOld; // Sol at n
		std::vector<double> vecNew; // Sol at n+1

		FDM()
		{

		}

		void initIC(const std::vector<double>& xarr)
		{ // Initialise the solutin at time zero. This occurs only 
		  // at the interior mesh points of xarr (and there are J-1 
		  // of them).
		  

			vecOld = std::vector<double>(xarr.size());

			// Initialise at the boundaries
			vecOld[0] = BCL(0.0);
			vecOld[vecOld.size()-1] = BCR(0.0);

			// Now initialise values in interior of interval using
			// the initial function 'IC' from the PDE
			for (unsigned int j = 1; j < xarr.size()-1; j++)
			{
				vecOld[j] = IC(xarr[j]);
			}

			//print(vecOld);

			vecNew = vecOld; // V2 optimise
		
		}

		const std::vector<double>& current() const
		{
			return vecNew;
		}

		void calculateCoefficients(const std::vector<double>& xarr, double tprev, double tnow)
		{ // Calculate the coefficients for the solver

			// Explicit method
		//	A = Vector<double, long> (xarr.Size(), xarr.MinIndex(), 0.0);
		//	C = A;
		//	B = Vector<double, long> (xarr.Size(), xarr.MinIndex(), 1.0);

			a = std::vector<double> (xarr.size()-2);
			bb = std::vector<double> (xarr.size()-2);
			c = std::vector<double> (xarr.size()-2);
			RHS = std::vector<double> (xarr.size()-2);

			double tmp1, tmp2;
			double k = tnow - tprev;
			double h = xarr[1] - xarr[0];

			for (unsigned int j = 1; j < xarr.size()-1; j++)
			{

				tmp1 = k * ((sigma)(xarr[j], tprev)/(h*h));
				tmp2 = k * (((mu)(xarr[j], tprev)* 0.5)/h);
	
				a[j-1] = tmp1 - tmp2;
				bb[j-1] = 1.0 - (2.0 * tmp1) + (k * (b)(xarr[j], tprev));
				c[j-1] = tmp1 + tmp2;
				RHS[j-1] = k * f(xarr[j], tprev);
			}

		}

		void solve (double tnow)
		{
			// Explicit method

			vecNew[0] = BCL(tnow);
			vecNew[vecNew.size()-1] = BCR(tnow);

			for (unsigned int i = 1; i < vecNew.size()-1; i++)
			{
				vecNew[i] = (a[i-1] * vecOld[i-1])
									+ (bb[i-1] * vecOld[i])
									+ (c[i-1] * vecOld[i+1]) - RHS[i-1];
			}
			vecOld = vecNew; // Not the most efficient, V2 can optimise it
	
		}

};


class FDMDirector
{

private:
	FDMDirector () {}

	double T;
	double Xmax;

	double k;
	long J, N;
	double tprev, tnow;
	FDM fdm;

public:
	std::vector<double> xarr; // Mesh array in space S
	std::vector<double> tarr; // Mesh array in time 

public:
	FDMDirector (double XM, double TM, long J, long NT)
	{

		T = TM;
		J = J;
		N = NT;
		Xmax = XM;
		fdm = FDM();

		// Create meshes in S and t
		Mesher mx(0.0, Xmax);
		xarr = mx.xarr(J);

		Mesher mt(0.0, T);
		tarr = mt.xarr(NT);

		Start();
	}

	
	const std::vector<double>& current() const
	{
		return fdm.current();
	}

	void Start() // Calculate next level
	{
		// Steps 1, 2: Get stuff from Mesher
		tprev = tnow = 0.0;
		k = T/N;
	
		// Step 3: Update new mesh array in FDM scheme
		fdm.initIC(xarr);

	}

	void doit()
	{
		// Step 4, 5: Get new coefficient arrays + solve
		
		for (unsigned int n = 1; n < tarr.size(); ++n)
		{
				tnow = tarr[n]; // n+1
				fdm.calculateCoefficients(xarr, tprev, tnow);
				fdm.solve(tnow);
				tprev = tnow; // n becomes n+1
		}
		
	}
};

// Global option data for PDE coefficients
OptionData* globalOptionData = nullptr;

// PDE coefficient functions for Black-Scholes
double sigmaPDE(double x, double t) {
	return 0.5 * globalOptionData->sig * globalOptionData->sig * x * x;
}

double muPDE(double x, double t) {
	return globalOptionData->r * x;
}

double bPDE(double x, double t) {
	return -globalOptionData->r;
}

double fPDE(double x, double t) {
	return 0.0;
}

// Boundary conditions
double BCLPDE(double t) {
	// For a call: V(0,t) = 0
	// For a put: V(0,t) = K*exp(-r*(T-t))
	if (globalOptionData->type == 1) {
		return 0.0; // Call
	} else {
		return globalOptionData->K * exp(-globalOptionData->r * (globalOptionData->T - t)); // Put
	}
}

double BCRPDE(double t) {
	// For a call: V(S_max,t) = S_max - K*exp(-r*(T-t))
	// For a put: V(S_max,t) = 0
	if (globalOptionData->type == 1) {
		// We need S_max, which we'll store temporarily
		// This will be set dynamically
		return 0.0; // Placeholder, will be updated in runFDMPricing
	} else {
		return 0.0; // Put
	}
}

// Helper function to set right boundary for calls
double S_max_global = 0.0;
double BCRPDE_call(double t) {
	if (globalOptionData->type == 1) {
		return S_max_global - globalOptionData->K * exp(-globalOptionData->r * (globalOptionData->T - t));
	}
	return 0.0;
}

// Initial condition (payoff at maturity)
double ICPDE(double x) {
	return globalOptionData->myPayOffFunction(x);
}

// Function to run a single FDM pricing test
double runFDMPricing(OptionData& option, double S_0, long J, long N) {
	// Set up global option data
	globalOptionData = &option;
	
	// Set spatial domain: [0, S_max]
	// S_max should be large enough to capture the relevant domain
	// A common choice is 3-5 times the strike price
	S_max_global = std::max(5.0 * option.K, 3.0 * S_0);
	
	double h = S_max_global / J;  // Space step
	double k = option.T / N;      // Time step
	
	// Check stability condition for explicit FDM
	// For Black-Scholes: k <= h^2 / (sigma^2 * S_max^2)
	double sigma_max = 0.5 * option.sig * option.sig * S_max_global * S_max_global;
	double stability_factor = k * sigma_max / (h * h);
	
	std::cout << "  Grid info: h=" << h << ", k=" << k << ", stability=" << stability_factor << std::endl;
	
	// If unstable, adjust N to satisfy stability
	if (stability_factor > 0.5) {
		N = static_cast<long>(std::ceil(option.T * sigma_max / (0.4 * h * h)));
		k = option.T / N;
		stability_factor = k * sigma_max / (h * h);
		std::cout << "  Adjusted N=" << N << " for stability (factor=" << stability_factor << ")" << std::endl;
	}
	
	// Set up PDE coefficients
	ParabolicIBVP::sigma = sigmaPDE;
	ParabolicIBVP::mu = muPDE;
	ParabolicIBVP::b = bPDE;
	ParabolicIBVP::f = fPDE;
	ParabolicIBVP::BCL = BCLPDE;
	ParabolicIBVP::BCR = BCRPDE_call;
	ParabolicIBVP::IC = ICPDE;
	
	// Create FDM solver
	FDMDirector fdmDirector(S_max_global, option.T, J, N);
	
	// Solve the PDE
	fdmDirector.doit();
	
	// Get the solution
	const std::vector<double>& solution = fdmDirector.current();
	
	// Interpolate to find the value at S_0
	// Find the grid points surrounding S_0
	const std::vector<double>& xarr = fdmDirector.xarr;
	
	// Linear interpolation
	double price = 0.0;
	for (size_t i = 0; i < xarr.size() - 1; ++i) {
		if (xarr[i] <= S_0 && S_0 <= xarr[i + 1]) {
			// Linear interpolation
			double weight = (S_0 - xarr[i]) / (xarr[i + 1] - xarr[i]);
			price = solution[i] * (1.0 - weight) + solution[i + 1] * weight;
			break;
		}
	}
	
	return price;
}

int main()
{
	std::cout << "Finite Difference Method Option Pricer - Testing Accuracy\n";
	std::cout << std::string(60, '=') << "\n\n";
	
	// Test parameters
	long J = 500;   // Space steps
	long N = 1000;  // Time steps
	
	std::cout << "Parameters:" << std::endl;
	std::cout << "  Space steps (J): " << J << std::endl;
	std::cout << "  Time steps (N): " << N << std::endl;
	std::cout << std::endl;
	
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
		
		double computed = runFDMPricing(option, test.S_0, J, N);
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