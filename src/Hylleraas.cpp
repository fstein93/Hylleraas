#include <iostream>
#include <vector>
#include <numeric>
#include <chrono>

#include "integrator.h"
#include "integrator_test.h"
#include "blas_util.h"

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

void calc_H(vector<double>& H, vector<double>& dH_dalpha, vector<double>& d2H_dalpha2, const size_t Z, const integrator& Integrator, const size_t n, const size_t m, const size_t k){
	// Kinetic energy terms
	size_t idx = 0 ;
	for (size_t n1 = 0 ; n1 <= n ; n1++) {
		for (size_t m1 = 0 ; m1 <= m ; m1++){
			for (size_t k1 = 0 ; k1 <= k ; k1+=2) {
				for (size_t n2 = 0 ; n2 <= n ; n2++) {
					for (size_t m2 = 0 ; m2 <= m ; m2++) {
						for (size_t k2 = 0 ; k2 <= k ; k2+=2) {
							if (n1 < n2 || ((n1 == n2) && (m1 < m2 || (m1 == m2 && k1 <= k2)))) {
								const double element = Integrator.integral_kinetic(n1, m1, k1, n2, m2, k2) ;
								H[idx] = element ;
                                                                const double deriv = Integrator.fac_dalpha_kinetic(n1+n2, m1+m2, k1+k2)*element ;
                                                                dH_dalpha[idx] += deriv ;
                                                                d2H_dalpha2[idx] = Integrator.fac_d2alpha_kinetic(n1+n2, m1+m2, k1+k2)*deriv ;
							}
							idx++ ;
						}
					}
				}
			}
		}
	}
	// Nuclear attraction terms
	idx = 0 ;
        const double Z_double = (double) Z ;
	for (size_t n1 = 0 ; n1 <= n ; n1++) {
		for (size_t m1 = 0 ; m1 <= m ; m1++){
			for (size_t k1 = 0 ; k1 <= k ; k1+=2) {
				for (size_t n2 = 0 ; n2 <= n ; n2++) {
					for (size_t m2 = 0 ; m2 <= m ; m2++) {
						for (size_t k2 = 0 ; k2 <= k ; k2+=2) {
							if (n1 < n2 || ((n1 == n2) && (m1 < m2 || (m1 == m2 && k1 <= k2)))) {
								const double element = -Z_double*Integrator.integral_nuclear(n1+n2, m1+m2, k1+k2) ;
								H[idx] += element ;
                                                                const double deriv = Integrator.fac_dalpha_nuclear(n1+n2, m1+m2, k1+k2)*element ;
                                                                dH_dalpha[idx] += deriv ;
                                                                d2H_dalpha2[idx] += Integrator.fac_d2alpha_nuclear(n1+n2, m1+m2, k1+k2)*deriv ;
							}
							idx++ ;
						}
					}
				}
			}
		}
	}
	// Electron-repulsion term
	idx = 0 ;
	for (size_t n1 = 0 ; n1 <= n ; n1++) {
		for (size_t m1 = 0 ; m1 <= m ; m1++){
			for (size_t k1 = 0 ; k1 <= k ; k1+=2) {
				for (size_t n2 = 0 ; n2 <= n ; n2++) {
					for (size_t m2 = 0 ; m2 <= m ; m2++) {
						for (size_t k2 = 0 ; k2 <= k ; k2+=2) {
							if (n1 < n2 || ((n1 == n2) && (m1 < m2 || (m1 == m2 && k1 <= k2)))) {
								const double element = Integrator.integral_repulsion(n1+n2, m1+m2, k1+k2) ;
								H[idx] += element ;
								const double deriv = Integrator.fac_dalpha_repulsion(n1+n2, m1+m2, k1+k2)*element ;
								dH_dalpha[idx] += deriv ;
								d2H_dalpha2[idx] += Integrator.fac_d2alpha_repulsion(n1+n2, m1+m2, k1+k2)*deriv ;
							}
							idx++ ;
						}
					}
				}
			}
		}
	}
}

void calc_S(vector<double>& S, vector<double>& dS_dalpha, vector<double>& d2S_dalpha2, const integrator & Integrator, const size_t n, const size_t m, const size_t k){
	size_t idx = 0 ;
	for (size_t n1 = 0 ; n1 <= n ; n1++) {
		for (size_t m1 = 0 ; m1 <= m ; m1++){
			for (size_t k1 = 0 ; k1 <= k ; k1+=2) {
				for (size_t n2 = 0 ; n2 <= n ; n2++) {
					for (size_t m2 = 0 ; m2 <= m ; m2++) {
						for (size_t k2 = 0 ; k2 <= k ; k2+=2) {
							if (n1 < n2 || ((n1 == n2) && (m1 < m2 || (m1 == m2 && k1 <= k2)))) {
								const double element = Integrator.integral_overlap(n1+n2, m1+m2, k1+k2) ;
								S[idx] = element ;
								dS_dalpha[idx] = Integrator.fac_dalpha_overlap(n1+n2, m1+m2, k1+k2)*element ;
								d2S_dalpha2[idx] = Integrator.fac_d2alpha_overlap(n1+n2, m1+m2, k1+k2)*dS_dalpha[idx] ;
							}
							idx++ ;
						}
					}
				}
			}
		}
	}
}

void print_vector(const vector<double>& a) {
	for (const auto& x: a) {
		printf("%f ", x) ;
	}
	printf("\n") ;
}

void calc_energy(const double alpha, const size_t n, const size_t m, const size_t k, const size_t Z, vector<double>& coefficients, double& energy, double& denergy_dalpha, double& d2energy_dalpha2) {
	constexpr bool do_debug = false ;
	constexpr bool do_debug_summary = true ;

	const size_t dim = (n+1)*(m+1)*(k/2+1) ;
	const size_t dim2 = dim*dim ;
        
	const integrator Integrator(alpha, n, m, k);
                
	// Calculate arrays
	vector<double> H(dim2), dH_dalpha(dim2), d2H_dalpha2(dim2), h_coeff(dim), S(dim2), dS_dalpha(dim2), d2S_dalpha2(dim2) ;
	calc_H(H, dH_dalpha, d2H_dalpha2, Z, Integrator, n, m, k) ;
        calc_S(S, dS_dalpha, d2S_dalpha2, Integrator, n, m, k) ;

        if (do_debug) {
                printf("\nH ") ;
                print_vector(H) ;
                printf("S ") ;
                print_vector(S) ;
                printf("dH ") ;
                print_vector(dH_dalpha) ;
                printf("dS ") ;
                print_vector(dS_dalpha) ;
                printf("d2H ") ;
                print_vector(d2H_dalpha2) ;
                printf("d2S ") ;
                print_vector(d2S_dalpha2) ;
        }

	// Copy the matrices H and S because we will need them later
	vector<double> H_copy(H) ;
	vector<double> S_copy(S) ;

        // Determine the coefficients and the energy
        calc_first_eig(H, S, coefficients, energy) ;

	if (do_debug) {
		printf("coeff ") ;
		print_vector(coefficients) ;
	}

        // Calculate C^T*H*C and 1.0/C^T*S*C because H and S will change in calc_first_eig
        // We do not need the case A=H without debugging
        double coeff_h_coeff = 0.0 ;
        if (do_debug) {
                matrix_vector_prod(0.0, h_coeff, 1.0, H_copy, coefficients) ;
                coeff_h_coeff = inner_product(coefficients.begin(), coefficients.end(), h_coeff.begin(), 0.0) ;
        }
        matrix_vector_prod(0.0, h_coeff, 1.0, S_copy, coefficients) ;
        const double inv_norm = 1.0/inner_product(coefficients.begin(), coefficients.end(), h_coeff.begin(), 0.0) ;

	matrix_vector_prod(0.0, h_coeff, 1.0, dH_dalpha, coefficients) ;
	const double coeff_dh_coeff = inner_product(coefficients.begin(), coefficients.end(), h_coeff.begin(), 0.0) ;
	matrix_vector_prod(0.0, h_coeff, 1.0, dS_dalpha, coefficients) ;
	const double coeff_ds_coeff = inner_product(coefficients.begin(), coefficients.end(), h_coeff.begin(), 0.0) ;

        matrix_vector_prod(0.0, h_coeff, 1.0, d2H_dalpha2, coefficients) ;
        const double coeff_d2h_coeff = inner_product(coefficients.begin(), coefficients.end(), h_coeff.begin(), 0.0) ;
        matrix_vector_prod(0.0, h_coeff, 1.0, d2S_dalpha2, coefficients) ;
        const double coeff_d2s_coeff = inner_product(coefficients.begin(), coefficients.end(), h_coeff.begin(), 0.0) ;

	denergy_dalpha = (coeff_dh_coeff-energy*coeff_ds_coeff)*inv_norm ;
	d2energy_dalpha2 = (coeff_d2h_coeff-energy*coeff_d2s_coeff-2.0*denergy_dalpha*coeff_ds_coeff)*inv_norm ;

	if (do_debug_summary) { printf("Debug %f %f %f %f %f\n", inv_norm, energy, coeff_h_coeff*inv_norm, denergy_dalpha, d2energy_dalpha2) ; }
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
