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

extern "C" {
extern void dsygvx_(int*, char*, char*, char*, int*, double*, int*, double*, int*, double*, double*, int*, int*, double*, int*, double*, double*, int*, double*, int*, int*, int*, int*) ;

extern double dlamch_(char*) ;

extern int ilaenv_(int*, char*, char*, int*, int*, int*, int*) ;
}

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

void calc_H(vector<double>& H, vector<double>& dH_dalpha, const size_t Z, const integrator & Integrator, const size_t n, const size_t m, const size_t k){
	// Kinetic energy terms
	size_t idx = 0 ;
	for (size_t n1 = 0 ; n1 <= n ; n1++) {
		for (size_t m1 = 0 ; m1 <= m ; m1++){
			for (size_t k1 = 0 ; k1 <= k ; k1++) {
				for (size_t n2 = 0 ; n2 <= n ; n2++) {
					for (size_t m2 = 0 ; m2 <= m ; m2++) {
						for (size_t k2 = 0 ; k2 <= k ; k2++) {
							if (n1 < n2 || ((n1 == n2) && (m1 < m2 || (m1 == m2 && k1 <= k2)))) {
								const double element = Integrator.integral_kinetic(n1, m1, k1, n2, m2, k2) ;
								H[idx] = element ;
								dH_dalpha[idx] = fma(Integrator.fac_dalpha_kinetic(n1+n2, m1+m2, k1+k2), element, dH_dalpha[idx]) ;
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
			for (size_t k1 = 0 ; k1 <= k ; k1++) {
				for (size_t n2 = 0 ; n2 <= n ; n2++) {
					for (size_t m2 = 0 ; m2 <= m ; m2++) {
						for (size_t k2 = 0 ; k2 <= k ; k2++) {
							if (n1 < n2 || ((n1 == n2) && (m1 < m2 || (m1 == m2 && k1 <= k2)))) {
								const double element = -Z_double*Integrator.integral_nuclear(n1+n2, m1+m2, k1+k2) ;
								H[idx] += element ;
								dH_dalpha[idx] = fma(Integrator.fac_dalpha_nuclear(n1+n2, m1+m2, k1+k2), element, dH_dalpha[idx]) ;
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
			for (size_t k1 = 0 ; k1 <= k ; k1++) {
				for (size_t n2 = 0 ; n2 <= n ; n2++) {
					for (size_t m2 = 0 ; m2 <= m ; m2++) {
						for (size_t k2 = 0 ; k2 <= k ; k2++) {
							if (n1 < n2 || ((n1 == n2) && (m1 < m2 || (m1 == m2 && k1 <= k2)))) {
								const double element = Integrator.integral_repulsion(n1+n2, m1+m2, k1+k2) ;
								H[idx] += element ;
								dH_dalpha[idx] = fma(Integrator.fac_dalpha_repulsion(n1+n2, m1+m2, k1+k2), element, dH_dalpha[idx]) ;
							}
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
				for (size_t n2 = 0 ; n2 <= n ; n2++) {
					for (size_t m2 = 0 ; m2 <= m ; m2++) {
						for (size_t k2 = 0 ; k2 <= k ; k2++) {
							if (n1 < n2 || ((n1 == n2) && (m1 < m2 || (m1 == m2 && k1 <= k2)))) {
								const double element = Integrator.integral_overlap(n1+n2, m1+m2, k1+k2) ;
								S[idx] = element ;
								dS_dalpha[idx] = fma(Integrator.fac_dalpha_overlap(n1+n2, m1+m2, k1+k2), element, dS_dalpha[idx]) ;
							}
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

void vector_add(const double alpha, vector<double> & y, const double beta, const vector<double> & x) {
  cblas_dscal((blasint) (y.size()), alpha, &y[0], 1);
  cblas_daxpy((blasint) (y.size()), beta, &x[0], 1, &y[0], 1);
}

void print_vector(const vector<double>& a) {
	for (const auto& x: a) {
		printf("%f ", x) ;
	}
	printf("\n") ;
}

void test_integrator() {
	constexpr size_t size_array = 2 ;
	constexpr double alpha[2] = {0.5, 1.0} ;
	constexpr size_t size_factorials = 13 ;
	constexpr double factorials[size_factorials] = {1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0, 39916800.0, 479001600.0} ;

	double error_exp_integral = 0.0 ;
	double error_kinetic = 0.0 ;
        double error_nuclear = 0.0 ;
        double error_repulsion = 0.0 ;
        double error_overlap = 0.0 ;

	for (size_t i = 0 ; i < size_array ; i++) {
		const double& alpha1 = alpha[i] ;
		const integrator Integrator(alpha1, 3, 3, 3);

        	for (size_t j = 0 ; j < size_factorials ; j++) {
			error_exp_integral += abs(Integrator.exp_integral(j)*pow(2.0*alpha1, j+1)/factorials[j]-1.0) ;
		}

	        error_kinetic += abs(Integrator.integral_kinetic(0, 0, 0, 0, 0, 0)*pow(alpha1, 4)-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(0, 0, 0, 0, 1, 0)*pow(alpha1, 5)*(16.0/25.0)-1.0) ;
	        error_kinetic += abs(Integrator.integral_kinetic(0, 1, 0, 0, 1, 0)*pow(alpha1, 6)/4.0-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(0, 0, 0, 0, 0, 1)*pow(alpha1, 5)-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(0, 0, 1, 0, 0, 1)*0.4*pow(alpha1, 6)-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(0, 0, 0, 0, 1, 1)*(16.0/35.0)*pow(alpha1, 6)-1.0) ;
                error_kinetic += abs(Integrator.integral_kinetic(0, 0, 0, 1, 0, 0)*0.5*pow(alpha1, 5)-1.0) ;
                error_kinetic += abs(Integrator.integral_kinetic(0, 0, 0, 1, 0, 1)*(1.0/3.0)*pow(alpha1, 6)-1.0) ;
                error_kinetic += abs(Integrator.integral_kinetic(0, 0, 0, 1, 1, 0)*(32.0/115.0)*pow(alpha1, 6)-1.0) ;
                error_kinetic += abs(Integrator.integral_kinetic(0, 0, 0, 1, 1, 1)*(32.0/245.0)*pow(alpha1, 7)-1.0) ;
                error_kinetic += abs(Integrator.integral_kinetic(0, 0, 1, 0, 1, 0)*(16.0/35.0)*pow(alpha1, 6)-1.0) ;
                error_kinetic += abs(Integrator.integral_kinetic(0, 0, 1, 0, 1, 1)*(16.0/103.0)*pow(alpha1, 7)-1.0) ;
                error_kinetic += abs(Integrator.integral_kinetic(0, 0, 1, 1, 0, 0)*0.5*pow(alpha1, 6)-1.0) ;
                error_kinetic += abs(Integrator.integral_kinetic(0, 0, 1, 1, 0, 1)*(2.0/15.0)*pow(alpha1, 7)-1.0) ;
                error_kinetic += abs(Integrator.integral_kinetic(0, 0, 1, 1, 1, 0)*(32.0/175.0)*pow(alpha1, 7)-1.0) ;
                error_kinetic += abs(Integrator.integral_kinetic(0, 0, 1, 1, 1, 1)*(32.0/707.0)*pow(alpha1, 8)-1.0) ;

	        error_nuclear += abs(Integrator.integral_nuclear(0, 0, 0)*0.5*pow(alpha1, 5)-1.0) ;
        	error_nuclear += abs(Integrator.integral_nuclear(0, 1, 0)*pow(alpha1, 6)*(4.0/15.0)-1.0) ;
	        error_nuclear += abs(Integrator.integral_nuclear(0, 2, 0)*pow(alpha1, 7)/9.0-1.0) ;

	        error_repulsion += abs(Integrator.integral_repulsion(0, 0, 0)*1.6*pow(alpha1, 5)-1.0) ;
        	error_repulsion += abs(Integrator.integral_repulsion(0, 1, 0)*pow(alpha1, 6)-1.0) ;
	        error_repulsion += abs(Integrator.integral_repulsion(0, 2, 0)*pow(alpha1, 7)*(16.0/35.0)-1.0) ;

	        error_overlap += abs(Integrator.integral_overlap(0, 0, 0)*pow(alpha1, 6)-1.0) ;
        	error_overlap += abs(Integrator.integral_overlap(0, 1, 0)*pow(alpha1, 7)*(16.0/35.0)-1.0) ;
	        error_overlap += abs(Integrator.integral_overlap(0, 2, 0)*pow(alpha1, 8)/6.0-1.0) ;

	}

	cout << "Error exp_integral " << error_exp_integral << endl ;

	cout << "Error kinetic " << error_kinetic << endl ;

        cout << "Error nuclear " << error_nuclear << endl ;

        cout << "Error repulsion " << error_repulsion << endl ;

        cout << "Error overlap " << error_overlap << endl ;
}

void calc_first_eig(vector<double>& H, vector<double>& S, vector<double>& coefficients, double& energy) {
	// Parameters for dlamch
	char cmach = 'S' ;

	// Parameters for ilaenv
	int ispec = 1 ;
	char name[] = "dsytrd" ;
	char opts[] = "VIU" ;
	int n1 = (int) coefficients.size() ;
        int n2 = 1 ;
        int n3 = -1 ;
        int n4 = -1 ;

	// Finally, set up dsygvx
	int itype = 1 ;
	char jobz = 'V' ;
	char range = 'I' ;
	char uplo = 'U' ;
	double vl = 0.0 ;
	double vu = 0.0 ;
	int il = 1 ;
	int iu = 1 ;
	double abstol = 2.0*dlamch_(&cmach) ;
	int found_evals, info ;
	int dim = (int) coefficients.size() ;
	int lwork = max(ilaenv_(&ispec, &name[0], &opts[0], &n1, &n2, &n3, &n4)+3, 8)*dim ;
	vector<int> iwork(2*((size_t) dim)), ifail((size_t) dim) ;
	vector<double> work((size_t) lwork), evals((size_t) dim) ;
	dsygvx_(&itype, &jobz, &range, &uplo, &dim, &H[0], &dim, &S[0], &dim, &vl, &vu, &il, &iu, &abstol, &found_evals, &evals[0], &coefficients[0], &dim, &work[0], &lwork, &iwork[0], &ifail[0], &info) ;
	energy = evals[0] ;
}

void calc_energy(const double alpha, const size_t n, const size_t m, const size_t k, const size_t Z, vector<double>& coefficients, double& energy, double& denergy_dalpha) {
	const size_t dim = (n+1)*(m+1)*(k+1) ;
	const size_t dim2 = dim*dim ;
        
	const integrator Integrator(alpha, n, m, k);
                
	// Calculate arrays
	vector<double> H(dim2), dH_dalpha(dim2), h_coeff(dim), S(dim2), dS_dalpha(dim2) ;
	calc_H(H, dH_dalpha, Z, Integrator, n, m, k) ;
        calc_S(S, dS_dalpha, Integrator, n, m, k) ;

	// Determine the coefficients and the energy
	calc_first_eig(H, S, coefficients, energy) ;

	matrix_vector_prod(0.0, h_coeff, 1.0, dH_dalpha, coefficients) ;
	matrix_vector_prod(1.0, h_coeff, -energy, dS_dalpha, coefficients) ;
	denergy_dalpha = inner_product(coefficients.begin(), coefficients.end(), h_coeff.begin(), 0.0) ;
}

int main(){
	// Get required parameters from user
	const double alpha0 = input_d() ;
	const bool do_wolfe = input_b() ;
	const size_t Z = input_ui() ;
	const size_t n = input_ui() ;
	const size_t m = input_ui() ;
	const size_t k = input_ui() ;

	const double gamma0 = 0.5 ;
	
	// Create working arrays
	const size_t dim = (n+1)*(m+1)*(k+1) ;
	const size_t dim2 = dim*dim ;
	vector<double> H(dim2), S(dim2), dH_dalpha(dim2), dS_dalpha(dim2) ;

	test_integrator();

	vector<double> coefficients(dim) ;

	size_t num_iter = 10 ;
	double alpha = alpha0 ;
	double gamma = gamma0 ;
	bool converged = false ;
	const double eps_gradient = 1.0e-5 ;

	double denergy_dalpha_old = 0.0 ;
	double alpha_old = alpha ;
	double energy = 0.0 ;
	double denergy_dalpha = 0.0 ;

	// Parameters of Wolfe condition
	const double c1 = 0.0001 ;
	const double c2 = 0.9 ;
	const double rho0 = 0.1 ;

	if (do_wolfe) {
		printf("Do Wolfe update: true") ;
	} else {
		printf("Do Wolfe update: false") ;
	}

	printf("iter time gamma alpha norm grad energy\n") ;

	for (size_t iter = 1 ; iter <= num_iter ; iter++) {

                double tstart = double(clock()) ;

		// Solve the Schrödinger equation for the given value of alpha
		calc_energy(alpha, n, m, k, Z, coefficients, energy, denergy_dalpha) ;

		// Determine Gamma
		if (iter > 1 && !do_wolfe) {
			denergy_dalpha_old = denergy_dalpha-denergy_dalpha_old ;
			const double inv_norm = 1.0/abs(denergy_dalpha_old) ;
			alpha_old = alpha-alpha_old ;
			gamma = inv_norm*alpha_old ;
		} else {
			vector<double> coefficients2(dim) ;
			double denergy_dalpha2, energy2 ;
			gamma = -1.0 ;
			do {
				if (gamma < 0.0) {
					gamma = rho0 ;
				} else {
					gamma /= 2.0 ;
				}
				if (gamma == 0.0) {
					converged = true ;
					break ;
				}
				calc_energy(alpha-gamma*denergy_dalpha, n, m, k, Z, coefficients2, energy2, denergy_dalpha2) ;
			}
			while (energy2 >= energy-c1*gamma*denergy_dalpha*denergy_dalpha && denergy_dalpha2 >= c2*denergy_dalpha) ;
		}

		denergy_dalpha_old = denergy_dalpha ;
		alpha_old = alpha ;
		alpha -= gamma*denergy_dalpha ;

		double tend = double(clock()) ;

                printf("%lu %f %f %f %f %f\n", iter, (tend-tstart)/CLOCKS_PER_SEC, gamma, alpha, denergy_dalpha, energy) ;

		if (abs(denergy_dalpha) < eps_gradient) converged = true ;

		if (converged) break ;

	}

	return 0;
}
