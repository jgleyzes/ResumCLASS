#ifndef RESUMCLASS_H
#define RESUMCLASS_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring> 
#include <sstream>
#include <fstream>
#include <vector>
#include <complex>

#include <gsl/gsl_spline.h>

#include <cuba.h>

#include "klist.h"

using namespace std ;

// Math stuff
double Heaviside (const double & a) ;
const double Pi = M_PI ;
const double E = exp(1.) ;


// Struct declarations
struct InterpFunc {
	gsl_interp_accel * accel ;
	gsl_spline * interp ;
} ;

struct ParamsP11 {
	double f ;
	gsl_interp_accel * accel ;
	gsl_spline * interp ;
} ;

struct ParamsIntegr {
	double k ;
	ParamsP11 p11 ;
	int id ;
	bool UseCosmoRef ;
	ParamsP11 p11ref ;
	double RescaleFactor ;
} ;

typedef double redshift ;

// Precision of the evalution of the loop integrals: 2 for UseCosmoRef = Yes/No, 2 for EpsAbs,EpsRel 
typedef double PrecisionIntegr[2][2] ; 

/** Cosmology **/
const size_t Nc = 5 ; // Number of cosmological parameters 
typedef double ParametersCosmology[Nc] ; // cosmological parameters
const string ParametersCosmologyNames[Nc] = { "A_s", "n_s", "h", "omega_b", "omega_cdm" } ; // A_s refers to ln(10^10 A_s)

/* Reference cosmology: Planck2015 */
const ParametersCosmology Reference = { 3.094, 0.9645, 0.6727, 0.02225, 0.1198 } ;

// Linear Growth rate f: GrowthFunction.cpp
double LinearGrowthRate (const ParametersCosmology & p, const redshift & z);

// LinearPowerSpectrum.cpp
void LoadP11 (const string & LinearPowerSpectrumData, const ParametersCosmology & cosmo, const redshift & z, ParamsP11 & p) ;
double P11 (const double & q, const ParamsP11 & params) ;
void UnloadP11 (const ParamsP11 & params) ;


// Multipole expansion
const size_t Nl = 5 ;
typedef double MultipoleMoments[Nl] ; // l = 0, 2, 4, 6, 8
// We call Mi with i: power of mu
const MultipoleMoments m0 = { 1., 0., 0., 0., 0. } ;
const MultipoleMoments m2 = { 1./3., 2./3., 0., 0., 0. } ;
const MultipoleMoments m4 = { 1./5., 4./7., 8./35., 0., 0. } ;
const MultipoleMoments m6 = { 1./7., 10./21., 24./77., 16./231., 0. } ;
const MultipoleMoments m8 = { 1./9., 40./99., 48./148., 64./495., 128./6435. } ;

// ComputePowerSpectra.cpp
const size_t N0 = 3 ; 	// 3 linear terms
typedef double PowerSpectraNoResum[Nl][N0][Nk] ;

void ComputePowerSpectraLinearNoResum (const ParamsP11 & params, PowerSpectraNoResum * Ps) ;

// LoadConfigFile.cpp
typedef bool YesNo ;
void LoadConfigFile (char * ConfigFile, double & nbar, double & km, double & knl, redshift & z0, ParametersCosmology & cosmo, string & PathToFolder,string & PathToFolderRD,string & PathToFolderCosmoRef, string & PathToLinearPowerSpectrum, YesNo & ComputePowerSpectrum, YesNo & UseRef, YesNo & ImportM, YesNo & ExportM, PrecisionIntegr & Eps) ;

#include "kout.h"


struct ParamsResumIntegr {
	InterpFunc InterpP, InterpM ;
} ;



//////////////////////////////////////////////////////////////////
// Functions needed to create the M-matrices.
//////////////////////////////////////////////////////////////////

// ResumI.cpp
typedef double ResumI (const double & k, const double & X1q, const double & f1) ;

ResumI I00, I02, I04, I06, I08, I22, I24, I26, I28, I44, I46, I48, I66, I68, I88 ;
ResumI J00, J02, J04, J22, J24, J44 ;

ResumI * const ResumI0[5][5] = { 
	{ I00, I02, I04, I06, I08 } ,
	{ I02, I22, I24, I26, I28 } ,
	{ I04, I24, I44, I46, I48 } ,
	{ I06, I26, I46, I66, I68 } ,
	{ I08, I28, I48, I68, I88 } 
} ;

ResumI * const ResumI2[3][3] = {
	{ J00, J02, J04 } ,
	{ J02, J22, J24 } ,
	{ J04, J24, J44 } ,
} ; 

typedef complex<double> dcomplex ;

const double qInfinity = 10000. ;
const double alpha = 1./3. ;
const double beta = 1./3. ;


const size_t Nlout = 3 ; // l = 0,2,4

// M^o_l,lp (k,kp) : StoreM [new kp:0 | M:1][o][l][lp][k][kp]
typedef double StoreM [2][2][Nlout][Nl][Nout][Nkp] ;

// ResumX1Y1.cpp
double ResumX1 (const double & q, const ParamsP11 & InterpP11) ;
void GetX1Y1 (const ParamsP11 & InterpP11, InterpFunc & InterpX1, InterpFunc & InterpY1, double & X1Infinity) ;

// ResumQ.cpp
double ResumQ0 (const int & l, const int & lp, const double & k, const double & q, const InterpFunc & InterpX1, const InterpFunc & InterpY1, const double & X1Infinity, const double & f1) ;
double ResumQ1(const int & l, const int & lp, const double & k, const double & q, const InterpFunc & InterpX1, const InterpFunc & InterpY1, const double & X1Infinity, const double & f1) ;

// ResumM.cpp
void ResumM (const string & PathToFolderRD, const ParamsP11 & params, const bool & ExportM, StoreM * TableM) ;
double QInfinity (const unsigned & order, const double & k, const unsigned & l, const unsigned & lp, const double & X1Infinity, const double & f1) ;
double LoadQInfinity (const string & PathToFolderRD, const unsigned & Morder, const double & k, const unsigned & l, const unsigned & lp) ;


//////////////////////////////////////////////////////////////////
// IR-resummation
//////////////////////////////////////////////////////////////////

// LoadResummation.cpp
void LoadM (const string & PathToFolderRD, const unsigned int & order, const double & k, const unsigned int & l, const unsigned int & lp, InterpFunc & InterpM) ;
void UnloadInterp (InterpFunc & params) ;

void UnloadInterp (InterpFunc & params) ;

// Resummation.cpp
double P (const double & q, const InterpFunc & InterpP) ;
double M (const double & q, const InterpFunc & InterpM) ;
double X1 (const double & q, const InterpFunc & InterpX1) ;
double Y1 (const double & q, const InterpFunc & InterpX1) ;

void ResumPowerSpectra (const string & PathToFolder,const string & PathToFolderRD, const ParamsP11 & InterpP11, PowerSpectraNoResum * Linear, PowerSpectraNoResum * Loop, const YesNo & ImportM, const YesNo & ExportM, StoreM * TableM) ;
double Resum_Integrand_GSL (double q, void * params) ;




#endif
