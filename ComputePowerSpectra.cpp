#include "ResumCLASS.h"


///////

void ComputePowerSpectraLinearNoResum (const ParamsP11 & params, PowerSpectraNoResum * Ps) {
	
	double f1 = params.f ;

	for (unsigned int i = 0 ; i < 3 ; i++) {
		for (unsigned int m = 0 ; m < Nk ; m++) {
			double P11k = P11(klist[m],params) ;

			(*Ps)[i][0][m] = m4[i] *f1*f1 *P11k ; 	// 1
			(*Ps)[i][1][m] = m2[i] *2.*f1 *P11k ; 	// b1
			(*Ps)[i][2][m] = m0[i] *P11k ; 			// b1*b1
		}
	}
}