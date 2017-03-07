//
//  PowDM_SPT_Resummation.hpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 23/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#ifndef PowDM_SPT_Resummation_hpp
#define PowDM_SPT_Resummation_hpp

#include <stdio.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_bessel.h>
#include "PowDM_SPT.hpp"

class IR_Resum_RSD;

class IR_Resum{
private:
    double g,f;
    int numevals;
    double kmax,kmin,q0PqInt;
    double Aval;// <<int dp/(6 pi^2) [P_L(p)(1-j0(pq))]>> , where P_L(p)=g^2 P*(p), and we average q where the wiggle part has support ~(10-300)Mpc/h
    std::vector<double>kvals;
    Coeffs coeffs;
    Splining P13,P22,P13_nw,P22_nw;
    Splining PowLinit,PowLinit_nw,PowL,PowL_nw,IR_lin,IR1loop,nonIR1loop;
    which_scale scale;
    finalData finaldata;
    calc_read calcopt;
    p13CalcDM p13,p13_nw;
    p22calcDM p22,p22_nw;
public:
    IR_Resum(import_class imports,Coeffs &coeffsIn,int numevalsIn,which_scale scaleIn=linearScale,EdSoption edsopt=EdSoption::Off);
public:
    void calc_IR_Resum_Linear();
    Splining set_IR_Resum_to1Loop_DM();
    Splining set_nonIR_Resum_to1Loop_DM();
    Splining get_IR_DM(){return IR1loop;}
    Splining get_nonIR_DM(){return nonIR1loop;}
    Splining get_PowL(){return PowL;}
    Splining get_PowL_IR(){return IR_lin;}
public:
    friend class IR_Resum_RSD;
};
struct params_Sigma{double q;Splining *PowL;};
double SigmaAverage(Splining &PowL,double qmin,double qmax);
double SigmaVal(double q,void *p);
double SigmaIntegrand(double pval,void *p);

#endif /* PowDM_SPT_Resummation_hpp */
