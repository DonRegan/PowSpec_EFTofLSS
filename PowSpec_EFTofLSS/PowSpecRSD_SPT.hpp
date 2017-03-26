//
//  PowSpecRSD_SPT.hpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 23/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#ifndef PowSpecRSD_SPT_hpp
#define PowSpecRSD_SPT_hpp

#include <stdio.h>
#include <fstream>
#include <sstream>
#include "PowDM_SPT.hpp"

//Much of the formulae here are in the Mathematica notebooks FabrikantImplemented.nb and RSD_Checks.nb

class p13CalcRSD{
private:
    double z,g,A,B,D,E,F,G,J;
    double f,Abarprime,Bbarprime,Dbarprime,Ebarprime,Fbarprime,Gbarprime,Jbarprime;
    double H,N,P,Q,R,S,K,L;
    double kmin,kmax,timedep_ct_dd,timedep_ct_dd_deriv,timedep_ct_tt;
    int numevals;
    Coeffs coeffs;
    Splining PowL;
    which_scale scale;
    std::vector<double>kvals;
    std::vector<Splining>J0,J1,J2,J3,J4; //Jn^m(k)=  1/(2 \[Pi]^2) int dq q^m F(n,q,k) P(q) , F(n,q,k)=\[Integral]dx x jl(k x) jl(q x)
    std::vector<Splining>P13dsds;
    Splining P13theta;
    NoWiggle_option wigOpt;
public:
    p13CalcRSD(p13CalcDM &p13DM,Splining &PowLIn);
public:
    void set_coeffs(EdSoption edsopt);
    void calc_P13_contribs(p13CalcDM &p13DM);
public://access routines
    Splining get_P13RSD(int mupow){assert(mupow==0||mupow==2||mupow==4||mupow==6);return P13dsds[mupow/2];}
    Splining get_P13theta(){return P13theta;}
};
class p22CalcRSD{
private:
    double z,g,A,B,K,L,f,Abarprime,Bbarprime;
    double kmin,kmax;
    int numevals;
    Coeffs coeffs;
    Splining PowL;
    which_scale scale;
    std::vector<double>kvals;
    std::vector<Splining>P22RSD_mu2n;
    Splining P22theta;
    NoWiggle_option wigOpt;
    finalData finaldata;
public:
    p22CalcRSD(p22calcDM &p22DM,Splining &PowLIn,calc_read calcoption=calc_read::calc);
public:
    void set_coeffs(EdSoption edsopt);
    double integrate_over_q(double k,int mupow,double lowerlim,double upperlim,double (*func)(double,void *),void *p);//for the double integral over q and r (with |k-q| <= r <= k+q )
    void calc_p22_RSD();
public://access routines
    Splining get_P22RSD(int mupow){assert(mupow==0||mupow==2||mupow==4||mupow==6||mupow==8);return P22RSD_mu2n[mupow/2];}
    Splining get_P22theta(){return P22theta;}
public://Input/Output routines
    void output_p22_RSD();
    void input_p22_RSD();
public://integrands for the 22 calculations
    friend double P22KKintegrand_r(double r,void *p);
    friend double P22vdvd_integrand(double r,void *p);
    friend double P22v2v2_integrand(double r,void *p);
    friend double P22Kvd_integrand(double r,void *p);
    friend double P22Kv2_integrand(double r,void *p);
    friend double P22v2vd_integrand(double r,void *p);
    friend double P22RSD_integrand(double r,void *p);
};
struct leg_paramsp22a{double &k;int &mupow;Splining &PowL;p22CalcRSD *p22class;};
struct leg_paramsp22b{double &k;double &q;int &mupow;Splining &PowL;double (*func)(double,void *);;p22CalcRSD *p22class;};
double integrateRSD_over_r(double q,void *p);

double P22KKintegrand_r(double r,void *p);
double P22vdvd_integrand(double r,void *p);
double P22v2v2_integrand(double r,void *p);
double P22Kvd_integrand(double r,void *p);
double P22Kv2_integrand(double r,void *p);
double P22v2vd_integrand(double r,void *p);
double P22RSD_integrand(double r,void *p);



#endif /* PowSpecRSD_SPT_hpp */
