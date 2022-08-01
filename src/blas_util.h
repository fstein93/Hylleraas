#ifndef BLAS_UTIL_H
#define BLAS_UTIL_H

void matrix_vector_prod(const double, std::vector<double>&, const double, std::vector<double> const &, std::vector<double> const &) ;

void vector_add(const double, std::vector<double> &, const double, const std::vector<double>&) ;

void calc_first_eig(std::vector<double>& H, std::vector<double>& S, std::vector<double>& coefficients, double& energy) ;

#endif
