#include <iostream>
#include <vector>
#include <stdio.h>
#include <numeric>
#include <cmath>
#include <chrono>
#include <cblas.h>

#include "integrator.h"
//#include "integrator_test.h"

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

void calc_H(vector<double>& H, vector<double>& dH_dalpha, const size_t Z, const integrator & Integrator, const size_t n, const size_t m, const size_t k){
	// Kinetic energy terms
	size_t idx = 0 ;
	for (size_t n1 = 0 ; n1 <= n ; n1++) {
		for (size_t m1 = 0 ; m1 <= m ; m1++){
			for (size_t k1 = 0 ; k1 <= k ; k1++) {
				for (size_t n2 = n1 ; n2 <= n ; n2++) {
					for (size_t m2 = (n1>n2 ? 0 : m1) ; m2 <= m ; m2++) {
						for (size_t k2 = ((n1 > n2 || m1 > m2) ? 0 : k1) ; k2 <= k ; k2++) {
							const double element = Integrator.integral_kinetic(n1, m1, k1, n2, m2, k2) ;
							H[idx] += element ;
							dH_dalpha[idx] = fma(Integrator.fac_dalpha_kinetic(n1+n2, m1+m2, k1+k2), element, dH_dalpha[idx]) ;
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
			for (size_t k1 = 0 ; k1 <= k ; k1++) {
				for (size_t n2 = n1 ; n2 <= n ; n2++) {
					for (size_t m2 = (n1>n2 ? 0 : m1) ; m2 <= m ; m2++) {
						for (size_t k2 = ((n1 > n2 || m1 > m2) ? 0 : k1) ; k2 <= k ; k2++) {
							const double element = -Z_double*Integrator.integral_nuclear(n1+n2, m1+m2, k1+k2) ;
							H[idx] += element ;
							dH_dalpha[idx] = fma(Integrator.fac_dalpha_nuclear(n1+n2, m1+m2, k1+k2), element, dH_dalpha[idx]) ;
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
				for (size_t n2 = n1 ; n2 <= n ; n2++) {
					for (size_t m2 = (n1>n2 ? 0 : m1) ; m2 <= m ; m2++) {
						for (size_t k2 = ((n1 > n2 || m1 > m2) ? 0 : k1) ; k2 <= k ; k2++) {
							const double element = Integrator.integral_repulsion(n1+n2, m1+m2, k1+k2) ;
							H[idx] += element ;
							dH_dalpha[idx] = fma(Integrator.fac_dalpha_repulsion(n1+n2, m1+m2, k1+k2), element, dH_dalpha[idx]) ;
							idx++ ;
						}
					}
				}
			}
		}
	}
}

void calc_S(vector<double>& S, vector<double>& dS_dalpha, const integrator & Integrator, const size_t n, const size_t m, const size_t k){
	size_t idx = 0 ;
	for (size_t n1 = 0 ; n1 <= n ; n1++) {
		for (size_t m1 = 0 ; m1 <= m ; m1++){
			for (size_t k1 = 0 ; k1 <= k ; k1++) {
				for (size_t n2 = n1 ; n2 <= n ; n2++) {
					for (size_t m2 = (n1>n2 ? 0 : m1) ; m2 <= m ; m2++) {
						for (size_t k2 = 0 ; k2 <= k ; k2++) {
							const double element = Integrator.integral_overlap(n1+n2, m1+m2, k1+k2) ;
							S[idx] += element ;
							dS_dalpha[idx] = fma(Integrator.fac_dalpha_overlap(n1+n2, m1+m2, k1+k2), element, dS_dalpha[idx]) ;
							idx++ ;
						}
					}
				}
			}
		}
	}
}

void matrix_vector_prod(const double alpha, vector<double> & y, const double beta, vector<double> const & A, vector<double> const & x) {
  cblas_dsymv(CblasRowMajor, CblasUpper, (blasint) (y.size()), beta, &A[0], (blasint) (x.size()), &x[0], 1, alpha, &y[0], 1);
}

void vector_add(const double alpha, vector<double> & y, const double beta, vector<double> & x) {
  cblas_dscal((blasint) (y.size()), alpha, &y[0], 1);
  cblas_daxpy((blasint) (y.size()), beta, &y[0], 1, &x[0], 1);
}

void test_integrator() {
	constexpr double alpha1 = 0.5 ;
	const integrator Integrator1(alpha1, 3, 3, 3);

        constexpr double factorials[] = {1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0, 39916800.0, 479001600.0} ;

	double error_exp_integral = 0.0 ;
        for (size_t i = 0 ; i < 13 ; i++) {
		error_exp_integral += abs(Integrator1.exp_integral(i)/factorials[i]-1.0) ;
	}

        double error_kinetic = 0.0 ;
        error_kinetic += abs(Integrator1.integral_kinetic(0, 0, 0, 0, 0, 0)*pow(alpha1, 4)-1.0) ;
	error_kinetic += abs(Integrator1.integral_kinetic(0, 0, 0, 0, 1, 0)*pow(alpha1, 5)*(16.0/25.0)-1.0) ;
        error_kinetic += abs(Integrator1.integral_kinetic(0, 1, 0, 0, 1, 0)*pow(alpha1, 6)/4.0-1.0) ;

        double error_nuclear = 0.0 ;
        error_nuclear += abs(0.5*Integrator1.integral_nuclear(0, 0, 0)*pow(alpha1, 5)-1.0) ;
        error_nuclear += abs(Integrator1.integral_nuclear(0, 1, 0)*pow(alpha1, 6)*(4.0/15.0)-1.0) ;
        error_nuclear += abs(Integrator1.integral_nuclear(0, 2, 0)*pow(alpha1, 7)/9.0-1.0) ;

        double error_repulsion = 0.0 ;
        error_repulsion += abs(1.6*Integrator1.integral_repulsion(0, 0, 0)*pow(alpha1, 5)-1.0) ;
        error_repulsion += abs(Integrator1.integral_repulsion(0, 1, 0)*pow(alpha1, 6)-1.0) ;
        error_repulsion += abs(Integrator1.integral_repulsion(0, 2, 0)*pow(alpha1, 7)*(16.0/35.0)-1.0) ;

        double error_overlap = 0.0 ;
        error_overlap += abs(Integrator1.integral_overlap(0, 0, 0)*pow(alpha1, 6)-1.0) ;
        error_overlap += abs(Integrator1.integral_overlap(0, 1, 0)*pow(alpha1, 7)*(16.0/35.0)-1.0) ;
        error_overlap += abs(Integrator1.integral_overlap(0, 2, 0)*pow(alpha1, 8)/6.0-1.0) ;

	constexpr double alpha2 = 1.0 ;
        const integrator Integrator2(alpha2, 3, 3, 3);

        constexpr double factorials2[] = {0.5, 0.25, 0.25, 0.375, 0.75, 1.875, 5.625, 19.6875, 78.75, 354.375, 1771.875, 9745.3125, 58471.875} ;

        for (size_t i = 0 ; i < 13 ; i++) {
                error_exp_integral += abs(Integrator2.exp_integral(i)/factorials2[i]-1.0) ;
        }

        error_kinetic += abs(Integrator2.integral_kinetic(0, 0, 0, 0, 0, 0)*pow(alpha2, 4)-1.0) ;
        error_kinetic += abs(Integrator2.integral_kinetic(0, 0, 0, 0, 1, 0)*pow(alpha2, 5)*(16.0/25.0)-1.0) ;
        error_kinetic += abs(Integrator2.integral_kinetic(0, 1, 0, 0, 1, 0)*pow(alpha2, 6)/4.0-1.0) ;
        error_nuclear += abs(0.5*Integrator2.integral_nuclear(0, 0, 0)*pow(alpha2, 5)-1.0) ;
        error_nuclear += abs(Integrator2.integral_nuclear(0, 1, 0)*pow(alpha2, 6)*(4.0/15.0)-1.0) ;
        error_nuclear += abs(Integrator2.integral_nuclear(0, 2, 0)*pow(alpha2, 7)/9.0-1.0) ;
        error_repulsion += abs(1.6*Integrator2.integral_repulsion(0, 0, 0)*pow(alpha2, 5)-1.0) ;
        error_repulsion += abs(Integrator2.integral_repulsion(0, 1, 0)*pow(alpha2, 6)-1.0) ;
        error_repulsion += abs(Integrator2.integral_repulsion(0, 2, 0)*pow(alpha2, 7)*(16.0/35.0)-1.0) ;
        error_overlap += abs(Integrator2.integral_overlap(0, 0, 0)*pow(alpha2, 6)-1.0) ;
	error_overlap += abs(Integrator2.integral_overlap(0, 1, 0)*pow(alpha2, 7)*(16.0/35.0)-1.0) ;
	error_overlap += abs(Integrator2.integral_overlap(0, 2, 0)*pow(alpha2, 8)/6.0-1.0) ;

	cout << "Error exp_integral " << error_exp_integral << endl ;

	cout << "Error kinetic " << error_kinetic << endl ;

        cout << "Error nuclear " << error_nuclear << endl ;

        cout << "Error repulsion " << error_repulsion << endl ;

        cout << "Error overlap " << error_overlap << endl ;
}

int main(){
	// Get required parameters from user
	const double alpha0 = input_d() ;
	const size_t Z = input_ui() ;
	const size_t n = input_ui() ;
	const size_t m = input_ui() ;
	const size_t k = input_ui() ;
	
	// Create working arrays
	const size_t dim = (n+1)*(m+1)*(k+1) ;
	const size_t dim2 = dim*dim ;
	vector<double> H(dim2), S(dim2), dH_dalpha(dim2), dS_dalpha(dim2) ;

	test_integrator();

	vector<double> coefficients(dim) ;
	
	coefficients[0] = 1.0 ;

	size_t num_iter = 10 ;
	double alpha = alpha0 ;

	for (size_t iter = 1 ; iter <= num_iter ; iter++) {

                double tstart = double(clock()) ;

                printf("\nIteration %lu\n", iter) ;

	        const integrator Integrator(alpha, n, m, k);

        	// Calculate arrays
	        double tstart2, tend2 ;
        	tstart2 = double(clock()) ;
	        calc_H(H, dH_dalpha, Z, Integrator, n, m, k) ;
        	tend2 = double(clock()) ;
	        printf("Time calc_H : %f\n", (tend2-tstart2)/CLOCKS_PER_SEC);
        	tstart2 = double(clock()) ;
	        calc_S(S, dS_dalpha, Integrator, n, m, k) ;
        	tend2 = double(clock()) ;
	        printf("Time calc_S : %f\n", (tend2-tstart2)/CLOCKS_PER_SEC);

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
	
		double tend = double(clock()) ;
		printf("Time energy : %f\n", (tend-tstart)/CLOCKS_PER_SEC);

	}
	
	return 0;
}
