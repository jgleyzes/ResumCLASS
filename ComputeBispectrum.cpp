#include "RedshiftBiasEFT.h"
#include <cuba.h>

#define NDIM 2
#define NCOMP 1
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 10000000
#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define KEY 0

#define NVEC 1
#define VERBOSE 0
#define LAST 4
#define SPIN NULL

#define EPSABS 0.1

void ComputeBispectrumMonopole (const double & EpsRel, const string & PathToFolder, const string & PathToTriangles, const ParamsP11 & params, const double & aperp, const double & apar) {

	// Import Triangles configuration
	ifstream triangle(PathToTriangles, ios::in) ;

	if (!triangle.is_open()) {
		cerr << "There was a problem opening the triangles configuration file:" << PathToTriangles << endl ;
		exit(EXIT_FAILURE) ;
	}

	else {
		// Export results in file
		ostringstream filename ;
		filename << PathToFolder << "/BispectrumTreeMonopole.dat" ;

		ofstream write (filename.str(), ios::out | ios::trunc) ;
		int pre = 13 ; // white spacing
		write << "#" << setw(pre-1) << "k1" << setw(pre) << "k2" << setw(pre) << "k3" << setw(pre) << "1" << setw(pre) << "b1" << setw(pre) << "b2" << setw(pre) << "b4" << setw(pre) << "b11" << setw(pre) << "b8*b8" << endl ;

		double k1, k2, k3 ;
		int neval, fail, nregions ; // Cuba integration thingy

		while(triangle >> k1 >> k2 >> k3) {

			write << setw(pre) << k1 << setw(pre) << k2 << setw(pre) << k3 ;

			ParamsIntegrBispectrumAP integr = { k1,k2,k3, params, aperp, apar, 0 } ;

			for (unsigned int i = 0 ; i < 5 ; i++) { // 4 integrations to perform

				integr.id = i ;
				double res, err, prob ; 

				// CUBA v4.2
				Cuhre (NDIM,NCOMP,IntegrandBispectrumAP_CUBA, &integr, NVEC, EpsRel, EPSABS, VERBOSE, MINEVAL, MAXEVAL, KEY, STATEFILE, SPIN, &nregions, &neval, &fail, &res, &err, &prob) ;
				
				write << setw(pre) << res ;
			}
			
			// One of the integrand does not depend on mu and phi so we just multiply it by 4*PI to get the result
			write << setw(pre) << 4.*M_PI*IntegrandBispectrumAP(5,k1,k2,k3,0,0,params,aperp,apar) << endl ;
		} ;
	}


}


static int IntegrandBispectrumAP_CUBA (const int *ndim, const double a[], const int *ncomp, double ff[], void *params) {

	ParamsIntegrBispectrumAP p = *(ParamsIntegrBispectrumAP *) params ;

	double mu = -1. + 2.*a[0] ;
	double phi = 0. + 2.*M_PI*a[1] ;

	ff[0] = 2.*2.*M_PI * IntegrandBispectrumAP (p.id, p.k1, p.k2, p.k3, mu, phi, p.p11, p.aperp, p.apar) ;
}


