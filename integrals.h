#ifndef INTEGRALS_H
#define INTEGRALS_H

inline double integral_overlap(const integrator Integrator, const size_t n, const size_t m, const size_t k) {
	return Integrator.integral_st(n, m+1, k);
}

inline double fac_dalpha_overlap(const double alpha, const size_t n, const size_t m, const size_t k) {
	return -((double) (2*k+m+n+5))/alpha ;
}

inline double integral_nuclear(const integrator Integrator, const size_t Z, const size_t n, const size_t m, const size_t k) {
	return 4.0*((double) Z)*Integrator.integral_plain(n+1, m+1, k);
}

inline double fac_dalpha_nuclear(const double alpha, const size_t n, const size_t m, const size_t k) {
	return -((double) (2*k+m+n+4))/alpha ;
}

inline double integral_kinetic(const integrator Integrator, const double alpha, const size_t n1, const size_t m1, const size_t k1, const size_t n2, const size_t m2, const size_t k2) {
	return alpha*alpha*Integrator.integral_st(n1+n2, m1+m2+1, k1+k2)
	+(n1+n2>0 ? -alpha*((double) n1+n2)*Integrator.integral_st(n1+n2-1, m1+m2+1, k1+k2) : 0.0)
	+((n1>0 && n2>0) ? ((double) n1)*((double) n2)*Integrator.integral_st(n1+n2-2, m1+m2+1, k1+k2) : 0.0)
	+((k1>0 && k2>0) ? 4.0*((double) k1)*((double) k2)*Integrator.integral_st(n1+n2, m1+m2+1, k1+k2-1) : 0.0)
	+((m1>0 && m2>0) ? ((double) m1)*((double) m2)*Integrator.integral_st(n1+n2, m1+m2-1, k1+k2) : 0.0)
	+((m1>0) ? -((double) m1)*(alpha-((n2>0) ? ((double) n2) : 0.0))*Integrator.integral_ut(n1+n2+1, m1+m2-1, k1+k2) : 0.0)
	+((m2>0) ? -((double) m2)*(alpha-((n1>0) ? ((double) n1) : 0.0))*Integrator.integral_ut(n1+n2+1, m1+m2-1, k1+k2) : 0.0)
	+((m1>0 && k2>0) ? -2.0*((double) m1)*((double) k2)*Integrator.integral_su(n1+n2, m1+m2-1, k1+k2) : 0.0)
	+((m2>0 && k1>0) ? -2.0*((double) m2)*((double) k1)*Integrator.integral_su(n1+n2, m1+m2-1, k1+k2) : 0.0);
}

inline double fac_dalpha_kinetic(const double alpha, const size_t n, const size_t m, const size_t k) {
	return -((double) (2*k+m+n+3))/alpha ;
}

inline double integral_repulsion(const integrator Integrator, const size_t n, const size_t m, const size_t k) {
	return Integrator.integral_st(n, m, k);
}

inline double fac_dalpha_repulsion(const double alpha, const size_t n, const size_t m, const size_t k) {
	return -((double) (2*k+m+n+4))/alpha ;
}

#endif
