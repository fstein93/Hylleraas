#ifndef BASIC_INTEGRALS_H
#define BASIC_INTEGRALS_H

#include <vector>

class integrator {
	public:
		integrator(const double alpha_, const size_t n, const size_t m, const size_t k) {
                        alpha = alpha_;
			basic_integrals.resize(n+m+2*k+6);
			basic_integrals[0] = 0.5/alpha;
			for (size_t i = 1 ; i < basic_integrals.size() ; i++ ) {
				basic_integrals[i] = ((double) i)*(0.5/alpha)*basic_integrals[i-1];
			}
		}
		double exp_integral(const size_t n) const {
			return basic_integrals[n] ;
		}
		double integral_plain(const size_t n, const size_t m, const size_t k) const {
			return 2.0/(((double) (2*k+1))*((double) (2*k+m+2)))*exp_integral(2*k+m+n+2);
		}
		double fac_dintegral_plain(const size_t n, const size_t m, const size_t k) const {
			return -((double) (2*k+m+n+3))/alpha ;
		}
		double integral_st(const size_t n, const size_t m, const size_t k) const {
			return 4.0*((double) 4*k+m+5)/(((double) (2*k+1))*((double) (2*k+3))*((double) (2*k+m+2))*((double) (2*k+m+4)))*exp_integral(2*k+m+n+4);
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
	private:
		std::vector<double> basic_integrals ;
		double alpha ;
};

#endif
