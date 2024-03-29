#ifndef INTEGRATOR_H
#define INTEGRATOR_H

class integrator {
	public:
		integrator(const double alpha_, const size_t n, const size_t m, const size_t k) {
                        alpha = alpha_;
			if (basic_integrals != NULL) delete basic_integrals ;
			size_basic_integrals = 2*(n+m+k)+7 ;
			basic_integrals = new double[size_basic_integrals] ;
			basic_integrals[0] = 0.5/alpha;
			for (size_t i = 1 ; i < size_basic_integrals ; i++ ) {
				basic_integrals[i] = ((double) i)*(0.5/alpha)*basic_integrals[i-1];
			}
		}
		~integrator() {if (basic_integrals != NULL) delete basic_integrals ;}
		double exp_integral(const size_t n) const {
			if (n >= size_basic_integrals) return 0.0 ;
			return basic_integrals[n] ;
		}
		double integral_st(const size_t n, const size_t m, const size_t k) const {
			return (k%2==0 ? 4.0*((double) (2*k+m+5))/(((double) (k+1))*((double) (k+3))*((double) (k+m+2))*((double) (k+m+4)))*exp_integral(k+m+n+4) : 0.0) ;
		}
		double fac_dintegral_st(const size_t n, const size_t m, const size_t k) const {
			return -static_cast<double>(k+m+n+5)/alpha ;
		}
                double fac_d2integral_st(const size_t n, const size_t m, const size_t k) const {
                        return -static_cast<double>(k+m+n+6)/alpha ;
                }
		double integral_ut(const size_t n, const size_t m, const size_t k) const {
			return (k%2==0 ? 4.0/(static_cast<double>(k+1)*static_cast<double>(k+3)*static_cast<double>(k+m+4))*exp_integral(k+m+n+4) : 0.0) ;
		}
		double fac_dintegral_ut(const size_t n, const size_t m, const size_t k) const {
			return -static_cast<double>(k+m+n+5)/alpha ;
		}
                double fac_d2integral_ut(const size_t n, const size_t m, const size_t k) const {
                        return -static_cast<double>(k+m+n+6)/alpha ;
                }
		double integral_su(const size_t n, const size_t m, const size_t k) const {
			return (k%2==0 ? -4.0/(static_cast<double>(k+1)*static_cast<double>(k+m+2)*static_cast<double>(k+m+4))*exp_integral(k+m+n+4) : 0.0) ;
		}
		double fac_dintegral_su(const size_t n, const size_t m, const size_t k) const {
			return -static_cast<double>(k+m+n+5)/alpha ;
		}
                double fac_d2integral_su(const size_t n, const size_t m, const size_t k) const {
                        return -static_cast<double>(k+m+n+6)/alpha ;
                }

		double integral_overlap(const size_t n, const size_t m, const size_t k) const {
                        return (k%2 == 0 ? 4.0*static_cast<double>(2*k+m+6)/(static_cast<double>(k+1)*static_cast<double>(k+3)*static_cast<double>(k+m+3)*static_cast<double>(k+m+5))*exp_integral(k+m+n+5) : 0.0);
		}

		double fac_dalpha_overlap(const size_t n, const size_t m, const size_t k) const {
                        return -static_cast<double>(k+m+n+6)/alpha ;
		}

                double fac_d2alpha_overlap(const size_t n, const size_t m, const size_t k) const {
                        return -static_cast<double>(k+m+n+7)/alpha ;
                }

		double integral_nuclear(const size_t n, const size_t m, const size_t k) const {
                        return (k%2 == 0 ? 8.0/(static_cast<double>(k+1)*static_cast<double>(k+m+3))*exp_integral(k+m+n+4) : 0.0);
		}

		double fac_dalpha_nuclear(const size_t n, const size_t m, const size_t k) const {
                        return -static_cast<double>(k+m+n+5)/alpha ;
		}

                double fac_d2alpha_nuclear(const size_t n, const size_t m, const size_t k) const {
                        return -static_cast<double>(k+m+n+6)/alpha ;
                }

		double integral_kinetic(const size_t n1, const size_t m1, const size_t k1, const size_t n2, const size_t m2, const size_t k2) const {
		        return ((k1+k2)%2 == 0 ? alpha*alpha*integral_st(n1+n2, m1+m2+1, k1+k2) // d2/ds2 (1)
		        +(n1+n2>0 ? -alpha*static_cast<double>(n1+n2)*integral_st(n1+n2-1, m1+m2+1, k1+k2) : 0.0) // d2/ds2 (2)
		        +((n1>0 && n2>0) ? static_cast<double>(n1)*static_cast<double>(n2)*integral_st(n1+n2-2, m1+m2+1, k1+k2) : 0.0) // d2/ds2 (3)
		        +((k1>0 && k2>0) ? static_cast<double>(k1)*static_cast<double>(k2)*integral_st(n1+n2, m1+m2+1, k1+k2-2) : 0.0) // d2/dt2
		        +((m1>0 && m2>0) ? static_cast<double>(m1)*static_cast<double>(m2)*integral_st(n1+n2, m1+m2-1, k1+k2) : 0.0) // d2/du2
		        +((m1>0) ? static_cast<double>(m1)*(-alpha*integral_ut(n1+n2+1, m1+m2-1, k1+k2)+((n2>0) ? static_cast<double>(n2)*integral_ut(n1+n2, m1+m2-1, k1+k2) : 0.0)) : 0.0) // d/ds*d/du
		        +((m2>0) ? static_cast<double>(m2)*(-alpha*integral_ut(n1+n2+1, m1+m2-1, k1+k2)+((n1>0) ? static_cast<double>(n1)*integral_ut(n1+n2, m1+m2-1, k1+k2) : 0.0)) : 0.0)
		        +((m1>0 && k2>0) ? -static_cast<double>(m1)*static_cast<double>(k2)*integral_su(n1+n2, m1+m2-1, k1+k2) : 0.0) // d/dt*d/du
		        +((m2>0 && k1>0) ? -static_cast<double>(m2)*static_cast<double>(k1)*integral_su(n1+n2, m1+m2-1, k1+k2) : 0.0) : 
			(k2>0 ? static_cast<double>(k2)*(alpha*integral_st(n1+n2, m1+m2+1, k1+k2-1)+(n1>0 ? -static_cast<double>(n1)*integral_st(n1+n2-1, m1+m2+1, k1+k2-1) : 0.0)) : 0.0)
                        +(k1>0 ? static_cast<double>(k1)*(alpha*integral_st(n1+n2, m1+m2+1, k1+k2-1)+(n2>0 ? -static_cast<double>(n2)*integral_st(n1+n2-1, m1+m2+1, k1+k2-1) : 0.0)) : 0.0) );
		}

		double fac_dalpha_kinetic(const size_t n, const size_t m, const size_t k) const {
                        return -static_cast<double>(k+m+n+4)/alpha ;
		}

                double fac_d2alpha_kinetic(const size_t n, const size_t m, const size_t k) const {
                        return -static_cast<double>(k+m+n+5)/alpha ;
                }

		double integral_repulsion(const size_t n, const size_t m, const size_t k) const {
                        return (k%2==0 ? 4.0*static_cast<double>(2*k+m+5)/(static_cast<double>(k+1)*static_cast<double>(k+3)*static_cast<double>(k+m+2)*static_cast<double>(k+m+4))*exp_integral(k+m+n+4) : 0.0);
		}

		double fac_dalpha_repulsion(const size_t n, const size_t m, const size_t k) const {
                        return -static_cast<double>(k+m+n+5)/alpha ;
		}

                double fac_d2alpha_repulsion(const size_t n, const size_t m, const size_t k) const {
                        return -static_cast<double>(k+m+n+6)/alpha ;
                }
	private:
		double* basic_integrals = NULL ;
		size_t size_basic_integrals = 0 ;
		double alpha = 0.0 ;
};

#endif
