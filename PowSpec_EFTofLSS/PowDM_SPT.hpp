//
//  PowDM_SPT.hpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 22/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#ifndef PowDM_SPT_hpp
#define PowDM_SPT_hpp

#include <stdio.h>
#include <boost/timer.hpp>
#include "CoeffsCalculator.hpp"
#include "import_powerspectra.hpp"
#include "Fabrikant.hpp"
#include "Wiggle_NoWiggle.hpp"

class p13CalcRSD;
class p22CalcRSD;

enum class calc_read{calc,read};
enum class EdSoption{On,Off};
enum class NoWiggle_option{On,Off};
//Much of the formulae here are in the Mathematica notebooks FabrikantImplemented.nb and RSD_Checks.nb


/*********************/
/***P13 calculation***/
/*********************/
class p13CalcDM{
private:
    double z,g,A,B,D,E,F,G,J;
    double kmin,kmax,timedep_ct_dd;//time dep counter term for delta delta
    int numevals;
    Coeffs coeffs;
    Splining PowL;
    which_scale scale;
    std::vector<double>kvals,P13DM;
    std::vector<Splining>J0,J1,J2,J3,J4; //Jn^m(k)=  1/(2 \[Pi]^2) int dq q^m F(n,q,k) P(q) , F(n,q,k)=\[Integral]dx x jl(k x) jl(q x)
    double q2PqInt,q0PqInt,q0PqIntIR; //   1/(2 \[Pi]^2) int dq q^2 P(q) &  1/(2 \[Pi]^2) int dq q^0 P(q)
    Splining p13_DM_spline;
    NoWiggle_option wigOpt;
    EdSoption edsopt;
public:
    p13CalcDM(import_class imports,Coeffs &coeffsIn,int numevalsIn,NoWiggle_option wigOptIn=NoWiggle_option::Off,which_scale scaleIn=linearScale,    EdSoption edsoptIn=EdSoption::Off);
    p13CalcDM(){}
    double get_timedep_ct_dd(){return timedep_ct_dd;}
    void calc_Jn_vectors();
    std::vector<double> calc_P13KKtypical(double A1,double B1,double C1,double D1,double E1,double G1);
public://access
    Splining get_p13_spline(){return p13_DM_spline;}
public://friends
    friend class p13CalcRSD;
};
struct leg_paramsp13{double &k;int &qindex;int & nindex;Splining &PowL;};
double Jintegrand(double qval,void *p);
double q2PqIntegrand(double qval,void *p);
double q0PqIntegrand(double qval,void *p);
double Fkernel(int n,double qval,double kval);

/*********************/
/***P22 calculation***/
/*********************/
class p22calcDM{
private:
    double z,g,A,B;
    int numpts;
    double kmax,kmin;
    std::vector<double>kvals;
    Splining PowL;
    std::vector<double>P22_DM;
    Splining P22_DM_spline;
    Coeffs coeffs;
    which_scale scale;
    std::string outputbase;
    finalData finaldata;
    NoWiggle_option wigOpt;
    EdSoption edsopt;
public:
    p22calcDM(import_class imports,Coeffs &coeffsIn,int numevalsIn,NoWiggle_option wigOptIn=NoWiggle_option::Off,which_scale scaleIn=linearScale,EdSoption edsoptIn=EdSoption::Off);
    p22calcDM(){}
    double P22KKtypical_integrand(double k,double q,double r,double A1,double B1,double A2,double B2);
public://access to the spline
    Splining get_p22_spline(){return P22_DM_spline;}
public://the integrands
    double integrate_over_q(double k,double lowerlim,double upperlim,double (*func)(double,void *),void *p);//for the double integral over q and r (with |k-q| <= r <= k+q )
    friend double P22DMintegrand_r(double r,void *p);
public://friends
    friend class p22CalcRSD;
};
struct paramsp22b{double &k;double &q;Splining &PowL;double (*func)(double,void *);p22calcDM *p22class;};
double integrate_over_r(double q,void *p);
double P22DMintegrand_r(double r,void *p);
double P22KKtypical_integrand(double k,double q,double r,double A1,double B1,double A2,double B2);

#endif /* PowDM_SPT_hpp */
