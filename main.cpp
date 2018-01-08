#include "ResumCLASS.h"

#include <ctime>

int main(int argc, char *argv[]) {

	if (argc <= 1) {
		cerr << "Error: no configuration file specified." << endl ;
		exit(EXIT_FAILURE) ;
	}
	
	else {

		// Default values
		double nbar = 1./105. , km = 1. , knl = 1. ;
		redshift z0 = 0.67 ;
		ParametersCosmology cosmo ; for (unsigned int i = 0 ; i < Nc ; i++) cosmo[i] = Reference[i] ;
		YesNo UseRef = 0, ImportM = 1, ExportM = 0, ComputePowerSpectrum = 1 ;
		string PathToFolder = "./" ;
		string PathToFolderRD ;
		string PathToFolderCosmoRef ;
		string PathToLinearPowerSpectrum ;
		PrecisionIntegr Eps = { { 1.,1e-2 } , { 1.,5e-1 } } ;

		LoadConfigFile (argv[1], nbar, km, knl, z0, cosmo, PathToFolder,PathToFolderRD,PathToFolderCosmoRef, PathToLinearPowerSpectrum, ComputePowerSpectrum, UseRef, ImportM, ExportM,Eps) ;

		///////////////////////////
		//
		///////////////////////////

		ParamsP11 paramsP11 ;
		LoadP11 (PathToLinearPowerSpectrum, cosmo, z0, paramsP11) ;


		PowerSpectraNoResum PsLinear ;
		static StoreM TableM ;


		int start_s=clock() ;

		if (ImportM == false) ResumM (PathToFolderRD, paramsP11, ExportM, &TableM) ;

		if (ComputePowerSpectrum == true) {
			ComputePowerSpectraLinearNoResum  (paramsP11, &PsLinear) ;

			ResumPowerSpectra (PathToFolder,PathToFolderRD, paramsP11, &PsLinear, 0, ImportM, ExportM, &TableM) ;
		}

		int stop_s=clock() ;
		cout << "ResumCLASS ran in " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " seconds." << endl ;

		
		///////////////////////////
		///
		///////////////////////////
	}

	return 0 ;

}
