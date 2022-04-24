#ifndef BASIC_INTEGRALS_H
#define BASIC_INTEGRALS_H

constexpr inline double exp_integral(const double alpha, const size_t n) {
	return (n>0 ? (((double) n)*(1.0/alpha)*exp_integral(alpha, n-1)) : 1.0/alpha);
}

constexpr inline double integral_basic1(const double alpha, const size_t n, const size_t m, const size_t k) {
	return 2.0/(((double) (2*k+1))*((double) (2*k+m+2)))*exp_integral(alpha, 2*k+m+n+2);
}

constexpr inline double fac_dalpha_basic1(const double alpha, const size_t n, const size_t m, const size_t k) {
	return -((double) (2*k+m+n+3))/alpha ;
}

constexpr inline double integral_basic2(const double alpha, const size_t n, const size_t m, const size_t k) {
	return 4.0*((double) 4*k+m+5)/(((double) (2*k+1))*((double) (2*k+3))*((double) (2*k+m+2))*((double) (2*k+m+4)))*exp_integral(alpha, 2*k+m+n+4);
}

constexpr inline double fac_dalpha_basic2(const double alpha, const size_t n, const size_t m, const size_t k) {
	return -((double) (2*k+m+n+5))/alpha ;
}

constexpr inline double integral_basic3(const double alpha, const size_t n, const size_t m, const size_t k) {
	return 4.0/(((double) (2*k+1))*((double) (2*k+3))*((double) (2*k+m+4)))*exp_integral(alpha, 2*k+m+n+4);
}

constexpr inline double fac_dalpha_basic3(const double alpha, const size_t n, const size_t m, const size_t k) {
	return -((double) (2*k+m+n+5))/alpha ;
}

constexpr inline double integral_basic4(const double alpha, const size_t n, const size_t m, const size_t k) {
	return -4.0/(((double) (2*k+1))*((double) (2*k+m+2))*((double) (2*k+m+4)))*exp_integral(alpha, 2*k+m+n+4);
}

constexpr inline double fac_dalpha_basic4(const double alpha, const size_t n, const size_t m, const size_t k) {
	return -((double) (2*k+m+n+5))/alpha ;
}

#endif
