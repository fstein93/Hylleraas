#include <iostream>
#include <vector>
#include <numeric>

#include "integrator.h"
#include "blas_util.h"
#include "energies.h"

using namespace std ;

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
