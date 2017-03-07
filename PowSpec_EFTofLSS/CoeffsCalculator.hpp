//
//  CoeffsCalculator.hpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 22/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#ifndef CoeffsCalculator_hpp
#define CoeffsCalculator_hpp

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "cosmology.hpp"
class Coeffs;
void solvesystem(Coeffs &coeffs_info);


class zdepCosmol_FLRW{
private:
    Cosmology_FLRW cosmol;
public:
    zdepCosmol_FLRW(Cosmology_FLRW cosmolIn=Cosmology_FLRW()):cosmol(cosmolIn){}
    ~zdepCosmol_FLRW(){}
    //specific functions for the calculation of the Coefficients
    double get_Omegam_z(double zval);
    double get_Omegam_z_prime(double zval);
    double get_epsilon(double zval);
    double  get_epsilon_prime(double zval);
    double get_c0(double zval);
    double get_c1(double zval);
    double get_c0prime(double zval);
    double get_c1prime(double zval);
    friend class Coeffs;
};

class Coeffs{
private:
    zdepCosmol_FLRW cosmoz;
    Cosmology_FLRW cosmol;
    std::vector<double> yinit,yfinal,yfinalEdS,zvals;
    int ode_size,numz;
    which_scale scale;
    double zinit,zfinal;
public:
    Coeffs(double zinitIn,double zfinalIn,int num_zsteps,Cosmology_FLRW &cosmolIn,which_scale scaleIn=which_scale::linearScale);
    Coeffs(){}
    ~Coeffs(){}
    double get_zfinal(){return zfinal;}
    void set_zvals();
    void initialise_y();
    zdepCosmol_FLRW get_cosmoz(){return cosmoz;}
    Cosmology_FLRW get_cosFLRW(){return cosmol;}
    int get_odesize(){return ode_size;}
    int get_numz(){return numz;}
    std::vector<double> get_yinit(){ return yinit;}
    double get_zval(int i){return zvals[i];}
    void set_yfinal(std::vector<double>ysolution);
    double get_yfinal(int i){return yfinal[i];}
    void set_yfinalEdS();
    double get_yfinalEdS(int i){return yfinalEdS[i];}
};


#endif /* CoeffsCalculator_hpp */
