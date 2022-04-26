#include <iostream>
#include <vector>
#include <stdio.h>
#include <numeric>
#include <cmath>
#include <chrono>

#include "basic_integrals.h"
#include "integrals.h"

using namespace std ;

double input_d(){
	double d ;
	cin >> d ;
	return d;
}

unsigned int input_ui(){
	unsigned int ui ;
	cin >> ui ;
	return ui;
}

void calc_H(vector<double>& H, vector<double>& dH_dalpha, const size_t Z, const integrator Integrator, const double alpha, const size_t n, const size_t m, const size_t k){
	// Kinetic energy terms
	size_t idx = 0 ;
	for (size_t n1 = 0 ; n1 <= n ; n1++) {
		for (size_t m1 = 0 ; m1 <= m ; m1++){
			for (size_t k1 = 0 ; k1 <= k ; k1++) {
				for (size_t n2 = 0 ; n2 <= n ; n2++) {
					for (size_t m2 = 0 ; m2 <= m ; m2++) {
						for (size_t k2 = 0 ; k2 <= k ; k2++) {
							const double element = integral_kinetic(Integrator, alpha, n1, m1, k1, n2, m2, k2) ;
							H[idx] += element ;
							dH_dalpha[idx] = fma(fac_dalpha_kinetic(Integrator, n1+n2, m1+m2, k1+k2), element, dH_dalpha[idx]) ;
							idx++ ;
						}
					}
				}
			}
		}
	}
	// Nuclear attraction terms
	idx = 0 ;
	for (size_t n1 = 0 ; n1 <= n ; n1++) {
		for (size_t m1 = 0 ; m1 <= m ; m1++){
			for (size_t k1 = 0 ; k1 <= k ; k1++) {
				for (size_t n2 = 0 ; n2 <= n ; n2++) {
					for (size_t m2 = 0 ; m2 <= m ; m2++) {
						for (size_t k2 = 0 ; k2 <= k ; k2++) {
							const double element = -integral_nuclear(Integrator, Z, n1+n2, m1+m2, k1+k2) ;
							H[idx] += element ;
							dH_dalpha[idx] += fma(fac_dalpha_nuclear(Integrator, n1+n2, m1+m2, k1+k2), element, dH_dalpha[idx]) ;
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
			for (size_t k1 = 0 ; k1 <= k ; k1++) {
				for (size_t n2 = 0 ; n2 <= n ; n2++) {
					for (size_t m2 = 0 ; m2 <= m ; m2++) {
						for (size_t k2 = 0 ; k2 <= k ; k2++) {
							const double element = integral_repulsion(Integrator, n1+n2, m1+m2, k1+k2) ;
							H[idx] += element ;
							dH_dalpha[idx] += fma(fac_dalpha_repulsion(Integrator, n1+n2, m1+m2, k1+k2), element, dH_dalpha[idx]) ;
							idx++ ;
						}
					}
				}
			}
		}
	}
}

void calc_S(vector<double>& S, vector<double>& dS_dalpha, const integrator Integrator, const size_t n, const size_t m, const size_t k){
	size_t idx = 0 ;
	for (size_t n1 = 0 ; n1 <= n ; n1++) {
		for (size_t m1 = 0 ; m1 <= m ; m1++){
			for (size_t k1 = 0 ; k1 <= k ; k1++) {
				for (size_t n2 = 0 ; n2 <= n ; n2++) {
					for (size_t m2 = 0 ; m2 <= m ; m2++) {
						for (size_t k2 = 0 ; k2 <= k ; k2++) {
							const double element = integral_overlap(Integrator, n1+n2, m1+m2, k1+k2) ;
							S[idx] += element ;
							dS_dalpha[idx] = fma(fac_dalpha_overlap(Integrator, n1+n2, m1+m2, k1+k2), element, dS_dalpha[idx]) ;
							idx++ ;
						}
					}
				}
			}
		}
	}
}

void matrix_vector_prod(const double alpha, vector<double> & y, const double beta, vector<double> const & A, vector<double> const & x) {
	if (alpha != 1.0) {
		if (alpha == 0.0) {
			for (size_t i = 0 ; i < y.size() ; i++) {
				y[i] = 0.0 ;
			}
		} else {
			for (size_t i = 0 ; i < y.size() ; i++) {
				y[i] *= alpha ;
			}
		}
	}
	if (beta != 0.0) {
		if (beta == 1.0) {
			size_t idx = 0 ;
			for (size_t i = 0 ; i < y.size() ; i++) {
				double & yi = y[i] ;
				for (size_t j = 0 ; j < x.size() ; j++) {
					yi = fma(A[idx], x[j], yi) ;
					idx++ ;
				}
			}
		} else {
			size_t idx = 0 ;
			for (size_t i = 0 ; i < y.size() ; i++) {
				double yi = 0.0 ;
				for (size_t j = 0 ; j < x.size() ; j++) {
					yi = fma(A[idx], x[j], yi) ;
					idx++ ;
				}
				y[i] = fma(beta, yi, y[i]) ;
				y[i] += beta*yi ;
			}
		}
	}
}

void vector_add(const double alpha, vector<double> & y, const double beta, vector<double> const & x) {
	if (alpha != 1.0) {
		if (alpha == 0.0) {
			for (size_t i = 0 ; i < y.size() ; i++) {
				y[i] = 0.0 ;
			}
		} else {
			for (size_t i = 0 ; i < y.size() ; i++) {
				y[i] *= alpha ;
			}
		}
	}
	if (beta != 0.0) {
		if (beta == 1.0) {
			for (size_t i = 0 ; i < y.size() ; i++) {
				y[i] += x[i] ;
			}
		} else {
			for (size_t i = 0 ; i < y.size() ; i++) {
				y[i] = fma(beta, x[i], y[i]) ;
			}
		}
	}
}

int main(){
	// Get required parameters from user
	const double alpha = input_d() ;
	const size_t Z = input_ui() ;
	const size_t n = input_ui() ;
	const size_t m = input_ui() ;
	const size_t k = input_ui() ;
	
	// Create working arrays
	const size_t dim = (n+1)*(m+1)*(k+1) ;
	const size_t dim2 = dim*dim ;
	vector<double> H(dim2), S(dim2), dH_dalpha(dim2), dS_dalpha(dim2) ;

        const integrator Integrator(alpha, n, m, k);
	
	if (true) {
		printf("exp_integral %zu %f\n", n+m+2*k+5, 2.0*alpha);
		for (size_t i = 0 ; i < n+m+2*k+6 ; i++) {
			printf("%zu %f\n", i, Integrator.exp_integral(i));
		}
		printf("integral_basic1 %zu %zu %zu %f\n", n, m, k, 2.0*alpha);
		for (size_t n1 = 0 ; n1 <= n+2 ; n1++) {
			for (size_t m1 = 0 ; m1 <= m+1 ; m1++) {
				for (size_t k1 = 0 ; k1 <= k+1 ; k1++) {
					printf("%zu %zu %zu %f\n", n1, m1, k1, Integrator.integral_plain(n1, m1, k1));
				}
			}
		}
		printf("integral_basic2 %zu %zu %zu %f\n", n, m, k, 2.0*alpha);
		for (size_t n1 = 0 ; n1 <= n ; n1++) {
			for (size_t m1 = 0 ; m1 <= m ; m1++) {
				for (size_t k1 = 0 ; k1 <= k ; k1++) {
					printf("%zu %zu %zu %f\n", n1, m1, k1, Integrator.integral_st(n1, m1, k1));
				}
			}
		}
	}
	
	// Calculate arrays
	if (true) {
		double tstart, tend ;
		tstart = clock() ;
		calc_H(H, dH_dalpha, Z, Integrator, alpha, n, m, k) ;
		tend = clock() ;
		printf("Time calc_H : %f\n", (tstart-tend)/CLOCKS_PER_SEC);
		tstart = clock() ;
		calc_S(S, dS_dalpha, Integrator, n, m, k) ;
		tend = clock() ;
		printf("Time calc_S : %f\n", (tstart-tend)/CLOCKS_PER_SEC);
		
		vector<double> coefficients(dim) ;
		
		coefficients[0] = 1.0 ;
		
		tstart = clock() ;
		
		vector<double> h_coeff(dim) ;
		matrix_vector_prod(0.0, h_coeff, 1.0, H, coefficients);
		
		vector<double> s_coeff(dim) ;
		matrix_vector_prod(0.0, s_coeff, 1.0, S, coefficients);
		
		const double coeff_h_coeff = inner_product(coefficients.begin(), coefficients.end(), h_coeff.begin(), 0.0) ;
		const double coeff_s_coeff = inner_product(coefficients.begin(), coefficients.end(), s_coeff.begin(), 0.0) ;
		const double inv_norm = 1.0/coeff_s_coeff ;
		
		const double energy = coeff_h_coeff/coeff_s_coeff ;
		
		vector<double> nabla_energy(dim) ;
		vector_add(0.0, nabla_energy, 1.0, h_coeff);
		vector_add(1.0, nabla_energy, -energy, s_coeff);
		vector_add(2.0*inv_norm, nabla_energy, 0.0, h_coeff);
		
		vector<double> dh_dalpha_coeff(dim) ;
		matrix_vector_prod(0.0, dh_dalpha_coeff, 1.0, dH_dalpha, coefficients);
		
		vector<double> ds_dalpha_coeff(dim) ;
		matrix_vector_prod(0.0, ds_dalpha_coeff, 1.0, dS_dalpha, coefficients);
		
		const double coeff_dh_dalpha_coeff = inner_product(coefficients.begin(), coefficients.end(), dh_dalpha_coeff.begin(), 0.0) ;
		const double coeff_ds_dalpha_coeff = inner_product(coefficients.begin(), coefficients.end(), ds_dalpha_coeff.begin(), 0.0) ;
		
		double denergy_dalpha = inv_norm*(coeff_dh_dalpha_coeff-energy*coeff_ds_dalpha_coeff) ;
		
		printf("Energy: %f\n", energy);
		printf("Norm gradient: %f\n", sqrt(inner_product(nabla_energy.begin(), nabla_energy.end(), nabla_energy.begin(), 0.0)+denergy_dalpha*denergy_dalpha)) ;
		
		tend = clock() ;
		printf("Time energy : %f\n", (tstart-tend)/CLOCKS_PER_SEC);
	}
	
	// Set values of arrays
	return 1;
}
