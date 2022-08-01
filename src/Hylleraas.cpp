#include <iostream>
#include <vector>
#include <chrono>

#include "integrator_test.h"
#include "energies.h"

using namespace std ;

double input_d() {
	double d ;
	cin >> d ;
	return d;
}

unsigned int input_ui() {
	unsigned int ui ;
	cin >> ui ;
	return ui;
}

bool input_b() {
	bool b ;
	cin >> b ;
	return b ;
}

int main(){
	// Available minimizers (to be converted to enum later)
	constexpr size_t do_wolfe = 1 ;
	constexpr size_t do_poly2 = 2 ;
	constexpr size_t do_barzilai_borwein = 3 ;

	// Get required parameters from user
	const double alpha0 = input_d() ;
	size_t minimizer = input_ui() ;
	const size_t Z = input_ui() ;
	const size_t n = input_ui() ;
	const size_t m = input_ui() ;
	const size_t k = input_ui() ;

	const double max_step_size = 0.1 ;
	
	// Create working arrays
	const size_t dim = (n+1)*(m+1)*(k/2+1) ;
	const size_t dim2 = dim*dim ;
	vector<double> H(dim2), S(dim2) ;

	test_integrator();

	vector<double> coefficients(dim) ;

	size_t num_iter = 20 ;
	double alpha = alpha0 ;
	double step_size = max_step_size ;
	double min_step_size = 1.0e-5 ;
	bool converged = false ;
	const double eps_gradient = 1.0e-5 ;

	// Set parameters here to make the compiler happy
	double denergy_dalpha_old = 0.0 ;
	double alpha_old = alpha ;
	double energy = 0.0 ;
	double denergy_dalpha = 0.0 ;
	double d2energy_dalpha2 = 0.0 ;

	// Parameters of Wolfe condition
	const double c1 = 0.0001 ;
	const double c2 = 0.9 ;

	if (minimizer==do_wolfe) {
		printf("Do Wolfe update\n") ;
	} else if (minimizer==do_poly2) {
		printf("Do Taylor update\n") ;
	} else if (minimizer==do_barzilai_borwein) {
		printf("Do Barzilai-Borwein update\n") ;
	} else {
		printf("Unknown minimizer! Switch to Wolfe update!") ;
		minimizer = do_wolfe ;
	}

	printf("\niter time step_size alpha norm grad energy\n") ;

	for (size_t iter = 1 ; iter <= num_iter ; iter++) {

                double tstart = double(clock()) ;

		// Solve the Schrödinger equation for the given value of alpha
		calc_energy(alpha, n, m, k, Z, coefficients, energy, denergy_dalpha, d2energy_dalpha2) ;

		// Determine Gamma
		if (minimizer==do_poly2) {
			// Determine minimum of 2nd order Taylor polynomial
			step_size = 1.0/d2energy_dalpha2 ;
		} else if (iter > 1 && minimizer==do_barzilai_borwein) {
			// Here, it reduces to a 1st order forward estimate of the second derivative
			step_size = abs((alpha-alpha_old)/(denergy_dalpha-denergy_dalpha_old)) ;
			denergy_dalpha_old = denergy_dalpha ;
			alpha_old = alpha ;
		} else {
			vector<double> coefficients2(dim) ;
			double denergy_dalpha2, energy2, d2energy_dalpha22 ;
			step_size = -1.0 ;
			do {
				if (step_size < 0.0) {
					step_size = max_step_size ;
				} else {
					step_size /= 2.0 ;
				}
				if (step_size == 0.0) {
					converged = true ;
					break ;
				}
				calc_energy(alpha-step_size*denergy_dalpha, n, m, k, Z, coefficients2, energy2, denergy_dalpha2, d2energy_dalpha22) ;
			}
			while (energy2 >= energy-c1*step_size*denergy_dalpha*denergy_dalpha && denergy_dalpha2 >= c2*denergy_dalpha && step_size >= min_step_size) ;
		}

		// Enforce a certain step size
		if (minimizer!=do_poly2) { step_size = max(step_size, min_step_size) ; }

		// Update alpha
		alpha -= step_size*denergy_dalpha ;

		double tend = double(clock()) ;

		// Print statistics
                printf("%lu %f %f %f %f %f\n", iter, (tend-tstart)/CLOCKS_PER_SEC, step_size, alpha, denergy_dalpha, energy) ;

		if (abs(denergy_dalpha) < eps_gradient) converged = true ;

		if (converged) break ;

	}

	return 0;
}
