#include "RedshiftBiasEFT.h"

double IntegrandBispectrumAP (const int & id, const double & k1, const double & k2, const double & k3, const double & mu, const double & phi, const ParamsP11 & p, const double & aperp, const double & apar) {
	
	double f1 = p.f ; 
	double APfactor = 1./(2.*pow(apar,2)*pow(aperp,4)) ;

	double q1, q2, q3, nu1, nu2, nu3 ;

	APtransform (aperp, apar, k1, k2, k3, mu, phi, q1, q2, q3, nu1, nu2, nu3) ;

	switch (id) {
		case 0: // 1*
			return APfactor* (f1*pow(q1,-1)*pow(q2,-1)*pow(q3,-1)*(q1*P11(q2,p)*P11(q3,p)*pow(nu1,2)*(7*pow(q2,2)*(nu2*nu3 + pow(1 - pow(nu2,2),0.5)*pow(1 - pow(nu3,2),0.5)) + 2*q2*q3*(7 - 4*pow(nu3,2) + pow(nu2,2)*(-4 + 8*pow(nu3,2)) + 8*nu2*nu3*pow(1 - pow(nu2,2),0.5)*pow(1 - pow(nu3,2),0.5)) + 7*(f1*nu2*nu3*pow(q1,2) + pow(q3,2)*(nu2*nu3 + pow(1 - pow(nu2,2),0.5)*pow(1 - pow(nu3,2),0.5)))) + P11(q1,p)*(q3*P11(q2,p)*pow(nu3,2)*(7*nu1*nu2*(pow(q2,2) + f1*pow(q3,2)) + 7*pow(q2,2)*pow(1 - pow(nu1,2),0.5)*pow(1 - pow(nu2,2),0.5) + 7*pow(q1,2)*(nu1*nu2 + pow(1 - pow(nu1,2),0.5)*pow(1 - pow(nu2,2),0.5)) + 2*q1*q2*(7 - 4*pow(nu2,2) + pow(nu1,2)*(-4 + 8*pow(nu2,2)) + 8*nu1*nu2*pow(1 - pow(nu1,2),0.5)*pow(1 - pow(nu2,2),0.5))) + q2*P11(q3,p)*pow(nu2,2)*(7*pow(q1,2)*(nu1*nu3 + pow(1 - pow(nu1,2),0.5)*pow(1 - pow(nu3,2),0.5)) + 2*q1*q3*(7 - 4*pow(nu3,2) + pow(nu1,2)*(-4 + 8*pow(nu3,2)) + 8*nu1*nu3*pow(1 - pow(nu1,2),0.5)*pow(1 - pow(nu3,2),0.5)) + 7*(f1*nu1*nu3*pow(q2,2) + pow(q3,2)*(nu1*nu3 + pow(1 - pow(nu1,2),0.5)*pow(1 - pow(nu3,2),0.5)))))))/7. ;
			break ;

		case 1: // b1*
			return APfactor* pow(q1,-1)*pow(q2,-1)*pow(q3,-1)*(q1*P11(q2,p)*P11(q3,p)*((pow(q2,2) + pow(q3,2))*(nu2*nu3 + pow(1 - pow(nu2,2),0.5)*pow(1 - pow(nu3,2),0.5)) + f1*nu1*(nu3*q2 + nu2*q3)*pow(pow(q1,2),0.5)) + P11(q1,p)*(q2*P11(q3,p)*((pow(q1,2) + pow(q3,2))*(nu1*nu3 + pow(1 - pow(nu1,2),0.5)*pow(1 - pow(nu3,2),0.5)) + f1*nu2*(nu3*q1 + nu1*q3)*pow(pow(q2,2),0.5)) + q3*P11(q2,p)*(pow(q1,2)*(nu1*nu2 + pow(1 - pow(nu1,2),0.5)*pow(1 - pow(nu2,2),0.5)) + f1*nu2*nu3*q1*pow(pow(q3,2),0.5) + q2*(nu1*nu2*q2 + q2*pow(1 - pow(nu1,2),0.5)*pow(1 - pow(nu2,2),0.5) + f1*nu1*nu3*pow(pow(q3,2),0.5))))) ;
			break ;

		case 2: // b2*
			return APfactor* (2*(P11(q2,p)*P11(q3,p)*(7 - 2*pow(nu3,2) + pow(nu2,2)*(-2 + 4*pow(nu3,2)) + 4*nu2*nu3*pow(1 - pow(nu2,2),0.5)*pow(1 - pow(nu3,2),0.5)) + P11(q1,p)*(P11(q2,p)*(7 - 2*pow(nu2,2) + pow(nu1,2)*(-2 + 4*pow(nu2,2)) + 4*nu1*nu2*pow(1 - pow(nu1,2),0.5)*pow(1 - pow(nu2,2),0.5)) + P11(q3,p)*(7 - 2*pow(nu3,2) + pow(nu1,2)*(-2 + 4*pow(nu3,2)) + 4*nu1*nu3*pow(1 - pow(nu1,2),0.5)*pow(1 - pow(nu3,2),0.5)))))/7. ;
			break ;

		case 3: // b4*
			return APfactor* 2*(P11(q2,p)*P11(q3,p) + P11(q1,p)*(P11(q2,p) + P11(q3,p))) ;
			break ;

		case 4: // b11*
			return APfactor* 105*(P11(q1,p) + P11(q2,p) + P11(q3,p)) ;
			break ;

		case 5: // b8*b8* // ANALYTICAL
			return APfactor * 11025. ;

	}
}



void APtransform (const double & aperp, const double & apar, 
	const double & k1, const double & k2, const double & k3, const double & mu1, const double & phi, 
	double & q1, double & q2, double & q3, double & nu1, double & nu2, double & nu3) {


	double F = apar/aperp ;

	double cosa5k = (k1*k1 + k2*k2 - k3*k3) / (2.*k1*k2) ;
	double mu2 = mu1 * cosa5k - sqrt(1.-mu1*mu1)*sqrt(1.-cosa5k*cosa5k)*cos(phi) ;
	double mu3 = - mu2 * cosa5k + sqrt(1.-mu2*mu2)*sqrt(1.-cosa5k*cosa5k) ;

	q1 = k1/aperp * sqrt(1.+ mu1*mu1 * (pow(F,-2)-1.)) ;
	q2 = k2/aperp * sqrt(1.+ mu2*mu2 * (pow(F,-2)-1.)) ;
	q3 = k3/aperp * sqrt(1.+ mu3*mu3 * (pow(F,-2)-1.)) ;

	nu1 = mu1/F * 1./sqrt(1.+ mu1*mu1 * (pow(F,-2)-1.)) ;
	nu2 = mu2/F * 1./sqrt(1.+ mu2*mu2 * (pow(F,-2)-1.)) ;
	nu3 = mu3/F * 1./sqrt(1.+ mu3*mu3 * (pow(F,-2)-1.)) ;
}