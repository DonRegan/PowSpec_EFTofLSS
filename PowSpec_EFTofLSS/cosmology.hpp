//
//  cosmology.hpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 22/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#ifndef cosmology_hpp
#define cosmology_hpp

#include "splining_routines.hpp"
#include "integration_routine.hpp"
#include "import_variables.hpp"

class zdepCosmol_FLRW;

class Cosmology_FLRW{//general cosmology variables
private:
    cosmolInfo cosmol;
    double h,ns,OmegaM,OmegaB,OmegaC,OmegaL,Omegar,OmegaK,TCMB,Neff,kpivot,H0;
    double deltaphi;
    float *kvals,*tf_full,*tf_baryon,*tf_cdm;int sizekvec;
    std::vector<double>kvec;
    Splining PowL;
public:
    Cosmology_FLRW(cosmolInfo cosmolIn=cosmolInfo(),double Neff_in=3.046,double kpivotIn=0.05);//setting defaults for the latter two values
    void set_othervals();
public://access routines
    double get_OmegaM(){return OmegaM;}
    cosmolInfo get_cosmolInfo(){return cosmol;}
public:
    //Next part is required for the no-wiggle power spectrum calculation
    double Ea(double a);//gives H(a)^2/H0^2
    double tk_eh(double k);
    std::vector<double> create_powspec_eh(int numpts,double kmin,double kmax,double z,std::vector<double>&kvals,std::vector<double>&Pkvals_eh,which_scale scale=logScale);
    double growth_function(double z);
    friend double growth_integrand(double aprime,void *p);
    friend class zdepCosmol_FLRW;
    
};
double growth_integrand(double aprime,void *p);

#endif /* cosmology_hpp */
