#include <iostream>
#include <cmath>

#include "integrator.h"
#include "integrator_test.h"

using namespace std ;

void test_integrator() {
	constexpr size_t size_array = 5 ;
	constexpr double alpha[5] = {0.5, 1.0, 2.0, 4.0, 10.0} ;
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
		error_kinetic += abs(Integrator.integral_kinetic(1, 1, 1, 1, 1, 1)/387.0*pow(alpha1, 10)-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(1, 1, 1, 2, 2, 2)*(64.0/2513511.0)*pow(alpha1, 13)-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(1, 1, 1, 3, 3, 3)/10135125.0*pow(alpha1, 16)-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(2, 2, 2, 1, 1, 1)*(64.0/2513511.0)*pow(alpha1, 13)-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(3, 3, 3, 1, 1, 1)/10135125.0*pow(alpha1, 16)-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(2, 2, 2, 3, 3, 3)*(32.0/179302233825.0)*pow(alpha1, 19)-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(3, 3, 3, 2, 2, 2)*(32.0/179302233825.0)*pow(alpha1, 19)-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(3, 3, 3, 3, 3, 3)/5296906345500.0*pow(alpha1, 22)-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(1, 2, 3, 3, 2, 1)/11240775.0*pow(alpha1, 16)-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(2, 1, 3, 0, 2, 1)*(32.0/1107513.0)*pow(alpha1, 13)-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(1, 1, 0, 3, 0, 1)*(4.0/2205.0)*pow(alpha1, 10)-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(3, 1, 3, 0, 2, 1)*(64.0/14142843.0)*pow(alpha1, 14)-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(1, 3, 0, 2, 1, 0)/3495.0*pow(alpha1, 11)-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(2, 2, 3, 0, 3, 1)*(16.0/23501907.0)*pow(alpha1, 15)-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(0, 3, 2, 3, 0, 3)*(64.0/86837751.0)*pow(alpha1, 15)-1.0) ;
		error_kinetic += abs(Integrator.integral_kinetic(2, 1, 2, 3, 2, 3)*(16.0/1261214955.0)*pow(alpha1, 17)-1.0) ;

	        error_nuclear += abs(Integrator.integral_nuclear(0, 0, 0)*0.5*pow(alpha1, 5)-1.0) ;
        	error_nuclear += abs(Integrator.integral_nuclear(0, 1, 0)*pow(alpha1, 6)*(4.0/15.0)-1.0) ;
	        error_nuclear += abs(Integrator.integral_nuclear(0, 2, 0)*pow(alpha1, 7)/9.0-1.0) ;
		error_nuclear += abs(Integrator.integral_nuclear(0, 0, 5)) ;
		error_nuclear += abs(Integrator.integral_nuclear(2, 3, 6)*(16.0/30405375.0)*pow(alpha1, 16)-1.0) ;
		error_nuclear += abs(Integrator.integral_nuclear(3, 6, 3)) ;
		error_nuclear += abs(Integrator.integral_nuclear(6, 5, 4)*(8.0/123743795175.0)*pow(alpha1, 20)-1.0) ;
		error_nuclear += abs(Integrator.integral_nuclear(5, 3, 5)) ;
		error_nuclear += abs(Integrator.integral_nuclear(6, 3, 2)*(32.0/212837625.0)*pow(alpha1, 16)-1.0) ;
		error_nuclear += abs(Integrator.integral_nuclear(3, 3, 5)) ;
		error_nuclear += abs(Integrator.integral_nuclear(3, 0, 1)) ;
		error_nuclear += abs(Integrator.integral_nuclear(4, 4, 3)) ;
		error_nuclear += abs(Integrator.integral_nuclear(3, 3, 6)*(2.0/30405375.0)*pow(alpha1, 17)-1.0) ;
		error_nuclear += abs(Integrator.integral_nuclear(4, 1, 6)*(8.0/18243225.0)*pow(alpha1, 16)-1.0) ;

	        error_repulsion += abs(Integrator.integral_repulsion(0, 0, 0)*1.6*pow(alpha1, 5)-1.0) ;
        	error_repulsion += abs(Integrator.integral_repulsion(0, 1, 0)*pow(alpha1, 6)-1.0) ;
	        error_repulsion += abs(Integrator.integral_repulsion(0, 2, 0)*pow(alpha1, 7)*(16.0/35.0)-1.0) ;
		error_repulsion += abs(Integrator.integral_repulsion(1, 0, 3)) ;
		error_repulsion += abs(Integrator.integral_repulsion(1, 4, 5)) ;
		error_repulsion += abs(Integrator.integral_repulsion(6, 5, 0)*(4.0/16891875.0)*pow(alpha1, 16)-1.0) ;
		error_repulsion += abs(Integrator.integral_repulsion(5, 4, 0)*(64.0/6081075.0)*pow(alpha1, 14)-1.0) ;
		error_repulsion += abs(Integrator.integral_repulsion(1, 3, 4)/1080.0*pow(alpha1, 13)-1.0) ;
		error_repulsion += abs(Integrator.integral_repulsion(0, 3, 3)) ;
		error_repulsion += abs(Integrator.integral_repulsion(6, 2, 1)) ;
		error_repulsion += abs(Integrator.integral_repulsion(4, 3, 4)/368550.0*pow(alpha1, 16)-1.0) ;
		error_repulsion += abs(Integrator.integral_repulsion(3, 3, 0)/540.0*pow(alpha1, 11)-1.0) ;
		error_repulsion += abs(Integrator.integral_repulsion(1, 5, 6)/1143450.0*pow(alpha1, 17)-1.0) ;

	        error_overlap += abs(Integrator.integral_overlap(0, 0, 0)*pow(alpha1, 6)-1.0) ;
        	error_overlap += abs(Integrator.integral_overlap(0, 1, 0)*pow(alpha1, 7)*(16.0/35.0)-1.0) ;
	        error_overlap += abs(Integrator.integral_overlap(0, 2, 0)*pow(alpha1, 8)/6.0-1.0) ;
		error_overlap += abs(Integrator.integral_overlap(2, 2, 4)/7020.0*pow(alpha1, 14)-1.0) ;
		error_overlap += abs(Integrator.integral_overlap(5, 5, 2)*(16.0/723647925.0)*pow(alpha1, 18)-1.0) ;
		error_overlap += abs(Integrator.integral_overlap(5, 5, 1)) ;
		error_overlap += abs(Integrator.integral_overlap(3, 0, 5)) ;
		error_overlap += abs(Integrator.integral_overlap(6, 2, 5)) ;
		error_overlap += abs(Integrator.integral_overlap(1, 2, 0)/24.0*pow(alpha1, 9)-1.0) ;
		error_overlap += abs(Integrator.integral_overlap(2, 6, 3)) ;
		error_overlap += abs(Integrator.integral_overlap(0, 4, 3)) ;
		error_overlap += abs(Integrator.integral_overlap(5, 0, 1)) ;
		error_overlap += abs(Integrator.integral_overlap(1, 3, 2)*(64.0/27027.0)*pow(alpha1, 12)-1.0) ;

	}

	cout << "Error exp_integral " << error_exp_integral << endl ;

	cout << "Error kinetic " << error_kinetic << endl ;

        cout << "Error nuclear " << error_nuclear << endl ;

        cout << "Error repulsion " << error_repulsion << endl ;

        cout << "Error overlap " << error_overlap << endl ;
}
