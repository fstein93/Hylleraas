#ifndef INTEGRATOR_H
#define INTEGRATOR_H

class integrator {
	public:
		integrator(const double alpha_, const size_t n, const size_t m, const size_t k) {
                        alpha = alpha_;
			if (basic_integrals != NULL) delete basic_integrals ;
			size_basic_integrals = n+m+2*k+7 ;
			basic_integrals = new double[n+m+2*k+6] ;
			basic_integrals[0] = 0.5/alpha;
			for (size_t i = 1 ; i < n+m+2*k+6 ; i++ ) {
				basic_integrals[i] = ((double) i)*(0.5/alpha)*basic_integrals[i-1];
			}
		}
		~integrator() {if (basic_integrals != NULL) delete basic_integrals ;}
		double exp_integral(const size_t n) const {
			if (n >= size_basic_integrals) return 0.0 ;
			return basic_integrals[n] ;
		}
		double integral_plain(const size_t n, const size_t m, const size_t k) const {
			return 2.0/(((double) (2*k+1))*((double) (2*k+m+2)))*basic_integrals[2*k+m+n+2];
		}
		double fac_dintegral_plain(const size_t n, const size_t m, const size_t k) const {
			return -((double) (2*k+m+n+3))/alpha ;
		}
		double integral_st(const size_t n, const size_t m, const size_t k) const {
			return 4.0*((double) (4*k+m+5))/(((double) (2*k+1))*((double) (2*k+3))*((double) (2*k+m+2))*((double) (2*k+m+4)))*exp_integral(2*k+m+n+4);
		}
		double fac_dintegral_st(const size_t n, const size_t m, const size_t k) const {
			return -((double) (2*k+m+n+5))/alpha ;
		}
		double integral_ut(const size_t n, const size_t m, const size_t k) const {
			return 4.0/(((double) (2*k+1))*((double) (2*k+3))*((double) (2*k+m+4)))*exp_integral(2*k+m+n+4);
		}
		double fac_dintegral_ut(const size_t n, const size_t m, const size_t k) const {
			return -((double) (2*k+m+n+5))/alpha ;
		}
		double integral_su(const size_t n, const size_t m, const size_t k) const {
			return -4.0/(((double) (2*k+1))*((double) (2*k+m+2))*((double) (2*k+m+4)))*exp_integral(2*k+m+n+4);
		}
		double fac_dintegral_su(const size_t n, const size_t m, const size_t k) const {
			return -((double) (2*k+m+n+5))/alpha ;
		}

		double integral_overlap(const size_t n, const size_t m, const size_t k) const {
                        return 4.0*((double) (4*k+m+6))/(((double) (2*k+1))*((double) (2*k+3))*((double) (2*k+m+3))*((double) (2*k+m+5)))*exp_integral(2*k+m+n+5);
		}

		double fac_dalpha_overlap(const size_t n, const size_t m, const size_t k) const {
                        return -((double) (2*k+m+n+6))/alpha ;
		}

		double integral_nuclear(const size_t n, const size_t m, const size_t k) const {
                        return 8.0/(((double) (2*k+1))*((double) (2*k+m+3)))*exp_integral(2*k+m+n+4);
		}

		double fac_dalpha_nuclear(const size_t n, const size_t m, const size_t k) const {
                        return -((double) (2*k+m+n+5))/alpha ;
		}

		double integral_kinetic(const size_t n1, const size_t m1, const size_t k1, const size_t n2, const size_t m2, const size_t k2) const {
		        return alpha*alpha*integral_st(n1+n2, m1+m2+1, k1+k2)
		        +(n1+n2>0 ? -alpha*((double) (n1+n2))*integral_st(n1+n2-1, m1+m2+1, k1+k2) : 0.0)
		        +((n1>0 && n2>0) ? ((double) n1)*((double) n2)*integral_st(n1+n2-2, m1+m2+1, k1+k2) : 0.0)
		        +((k1>0 && k2>0) ? 4.0*((double) k1)*((double) k2)*integral_st(n1+n2, m1+m2+1, k1+k2-1) : 0.0)
		        +((m1>0 && m2>0) ? ((double) m1)*((double) m2)*integral_st(n1+n2, m1+m2-1, k1+k2) : 0.0)
		        +((m1>0) ? -((double) m1)*(alpha-((n2>0) ? ((double) n2) : 0.0))*integral_ut(n1+n2+1, m1+m2-1, k1+k2) : 0.0)
		        +((m2>0) ? -((double) m2)*(alpha-((n1>0) ? ((double) n1) : 0.0))*integral_ut(n1+n2+1, m1+m2-1, k1+k2) : 0.0)
		        +((m1>0 && k2>0) ? -2.0*((double) m1)*((double) k2)*integral_su(n1+n2, m1+m2-1, k1+k2) : 0.0)
		        +((m2>0 && k1>0) ? -2.0*((double) m2)*((double) k1)*integral_su(n1+n2, m1+m2-1, k1+k2) : 0.0);
		}

		double fac_dalpha_kinetic(const size_t n, const size_t m, const size_t k) const {
                        return -((double) (2*k+m+n+4))/alpha ;
		}

		double integral_repulsion(const size_t n, const size_t m, const size_t k) const {
                        return 4.0*((double) (4*k+m+5))/(((double) (2*k+1))*((double) (2*k+3))*((double) (2*k+m+2))*((double) (2*k+m+4)))*exp_integral(2*k+m+n+4);
		}

		double fac_dalpha_repulsion(const size_t n, const size_t m, const size_t k) const {
                        return -((double) (2*k+m+n+5))/alpha ;
		}
	private:
		double* basic_integrals = NULL ;
		size_t size_basic_integrals = 0 ;
		double alpha ;
};

#endif
