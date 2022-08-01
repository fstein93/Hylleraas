#include <iostream>
#include <vector>
#include <cblas.h>
#include <cmath>

#include "blas_util.h"

using namespace std ;

extern "C" {
extern void dsygvx_(int*, char*, char*, char*, int*, double*, int*, double*, int*, double*, double*, int*, int*, double*, int*, double*, double*, int*, double*, int*, int*, int*, int*) ;

extern double dlamch_(char*) ;

extern int ilaenv_(int*, char*, char*, int*, int*, int*, int*) ;
}

void matrix_vector_prod(const double alpha, vector<double> & y, const double beta, vector<double> const & A, vector<double> const & x) {
  cblas_dsymv(CblasRowMajor, CblasUpper, static_cast<blasint>(y.size()), beta, &A[0], static_cast<blasint>(x.size()), &x[0], 1, alpha, &y[0], 1);
}

void vector_add(const double alpha, vector<double> & y, const double beta, const vector<double> & x) {
  cblas_dscal(static_cast<blasint>(y.size()), alpha, &y[0], 1);
  cblas_daxpy(static_cast<blasint>(y.size()), beta, &x[0], 1, &y[0], 1);
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
	// We have the lower half because in Fortran, we have column-major
	char uplo = 'L' ;
	double vl = 0.0 ;
	double vu = 0.0 ;
	int il = 1 ;
	int iu = 1 ;
	double abstol = 2.0*dlamch_(&cmach) ;
	int found_evals, info ;
	size_t dim = coefficients.size() ;
	int dim_int = static_cast<int>(dim) ;
	int lwork = max(ilaenv_(&ispec, &name[0], &opts[0], &n1, &n2, &n3, &n4)+3, 8)*dim_int ;
	vector<int> iwork(2*dim), ifail(dim) ;
	vector<double> work(static_cast<size_t>(lwork)), evals(dim) ;
	dsygvx_(&itype, &jobz, &range, &uplo, &dim_int, &H[0], &dim_int, &S[0], &dim_int, &vl, &vu, &il, &iu, &abstol, &found_evals, &evals[0], &coefficients[0], &dim_int, &work[0], &lwork, &iwork[0], &ifail[0], &info) ;
	if (info == 0) {
		energy = evals[0] ;
	} else {
		printf("DSYGVX failed\n") ;
		energy = H[0]/S[0] ;
		coefficients[0] = pow(S[0], 0.5) ;
		for (size_t i = 1 ; i < dim ; i++) {
			coefficients[i] = 0.0 ;
		}
	}
}

