#ifndef INTEGRALS_H
#define INTEGRRALS_H

constexpr inline double integral_overlap(const double alpha, const size_t n, const size_t m, const size_t k) {
	return integral_basic2(2.0*alpha, n, m+1, k);
}

constexpr inline double fac_dalpha_overlap(const double alpha, const size_t n, const size_t m, const size_t k) {
	return -((double) (2*k+m+n+5))/alpha ;
}

constexpr inline double integral_nuclear(const double alpha, const size_t Z, const size_t n, const size_t m, const size_t k) {
	return 4.0*((double) Z)*integral_basic1(2.0*alpha, n+1, m+1, k);
}

constexpr inline double fac_dalpha_nuclear(const double alpha, const size_t n, const size_t m, const size_t k) {
	return -((double) (2*k+m+n+4))/alpha ;
}

constexpr inline double integral_kinetic(const double alpha, const size_t n1, const size_t m1, const size_t k1, const size_t n2, const size_t m2, const size_t k2) {
	return alpha*alpha*integral_basic2(2.0*alpha, n1+n2, m1+m2+1, k1+k2)
	+(n1+n2>0 ? -alpha*((double) n1+n2)*integral_basic2(2.0*alpha, n1+n2-1, m1+m2+1, k1+k2) : 0.0)
	+((n1>0 && n2>0) ? ((double) n1)*((double) n2)*integral_basic2(2.0*alpha, n1+n2-2, m1+m2+1, k1+k2) : 0.0)
	+((k1>0 && k2>0) ? 4.0*((double) k1)*((double) k2)*integral_basic2(2.0*alpha, n1+n2, m1+m2+1, k1+k2-1) : 0.0)
	+((m1>0 && m2>0) ? ((double) m1)*((double) m2)*integral_basic2(2.0*alpha, n1+n2, m1+m2-1, k1+k2) : 0.0)
	+((m1>0) ? -((double) m1)*(alpha-((n2>0) ? ((double) n2) : 0.0))*integral_basic3(2.0*alpha, n1+n2+1, m1+m2-1, k1+k2) : 0.0)
	+((m2>0) ? -((double) m2)*(alpha-((n1>0) ? ((double) n1) : 0.0))*integral_basic3(2.0*alpha, n1+n2+1, m1+m2-1, k1+k2) : 0.0)
	+((m1>0 && k2>0) ? -2.0*((double) m1)*((double) k2)*integral_basic4(2.0*alpha, n1+n2, m1+m2-1, k1+k2) : 0.0)
	+((m2>0 && k1>0) ? -2.0*((double) m2)*((double) k1)*integral_basic4(2.0*alpha, n1+n2, m1+m2-1, k1+k2) : 0.0);
}

constexpr inline double fac_dalpha_kinetic(const double alpha, const size_t n, const size_t m, const size_t k) {
	return -((double) (2*k+m+n+3))/alpha ;
}

constexpr inline double integral_repulsion(const double alpha, const size_t n, const size_t m, const size_t k) {
	return integral_basic2(2.0*alpha, n, m, k);
}

constexpr inline double fac_dalpha_repulsion(const double alpha, const size_t n, const size_t m, const size_t k) {
	return -((double) (2*k+m+n+4))/alpha ;
}

#endif
