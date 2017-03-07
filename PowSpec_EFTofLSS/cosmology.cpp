//
//  cosmology.cpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 22/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#include "cosmology.hpp"
using namespace std;

Cosmology_FLRW::Cosmology_FLRW(cosmolInfo cosmolIn,double Neff_in,double kpivotIn):cosmol(cosmolIn),Neff(Neff_in),kpivot(kpivotIn){
    set_othervals();
}

void Cosmology_FLRW::set_othervals(){
    OmegaM=cosmol.OmegaM,OmegaB=cosmol.OmegaB,h=cosmol.h,ns=cosmol.ns,deltaphi=cosmol.deltaphi;
    OmegaC=OmegaM-OmegaB;
    double rho_critical=1.87847*h*h*1E5;
    double rho_gamma=4.64511*pow(TCMB/2.7255,4); 
    double rho_neutrino=Neff*(7./8.)*pow(4./11.,4./3.)*rho_gamma;
    double Omega_gamma=rho_gamma/rho_critical;
    double Omega_neutrino=rho_neutrino/rho_critical;
    Omegar=Omega_gamma+Omega_neutrino;
    OmegaL=1.0-OmegaM-Omegar;
    OmegaK=0.;
    H0=100.*h/3.0E5; //km/s/Mpc -> c /Mpc (we use c=1)
}

double Cosmology_FLRW::Ea(double a){//gives H(a)^2/H0^2
    assert(a<=1.);
    double res=OmegaM/pow(a,3)+Omegar/pow(a,4)+OmegaK/pow(a,2)+OmegaL;
    return res;
}
double Cosmology_FLRW::tk_eh(double k){//need k in h/Mpc //Martin White function
    double q,theta,ommh2,ombh2,a,s,gamma,L0,C0,tmp;
    ombh2=OmegaB*h*h;
    ommh2=OmegaM*h*h;
    theta = 2.728/2.7;
    s = 44.5 * log(9.83 / ommh2) / sqrt(1. + 10. * exp(0.75 * log(ombh2))) * h;
    a = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2
    + 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * (ombh2 / ommh2);
    gamma = a + (1. - a) / (1. + exp(4 * log(0.43 * k * s)));
    gamma *= OmegaM * h;
    q = k * theta * theta / gamma;
    L0 = log(2. * exp(1.) + 1.8 * q);
    C0 = 14.2 + 731. / (1. + 62.5 * q);
    tmp = L0 / (L0 + C0 * q * q);
    return (tmp);
}

vector<double> Cosmology_FLRW::create_powspec_eh(int numpts,double kmin,double kmax,double z,vector<double>&kvals,vector<double>&Pkvals_eh,which_scale scale){
    create_vec(kmin*0.999, kmax*1.001, numpts, kvals, scale);//CHANGE on 12/9/2016

    Pkvals_eh=vector<double>(numpts);
    for (int i=0; i<numpts; i++) {
        double k=kvals[i];
        double Mkz_eh=tk_eh(k)*(2./3.)*k*k*growth_function(z)*h*h/(OmegaM*H0*H0);
        Pkvals_eh[i]=deltaphi*pow(kvals[i]/kpivot,ns-1.)/pow(kvals[i],3);
        Pkvals_eh[i]*=Mkz_eh*Mkz_eh;
    }
    return Pkvals_eh;
}


double Cosmology_FLRW::growth_function(double z){ //this is only good for w=-1
    double a=1./(1+z);
    double ha_h0=sqrt(Ea(a));
    double growth_fn=integrate(0.,a,&growth_integrand,this)*(5.*OmegaM/2.)*ha_h0;
    return growth_fn;
}

double growth_integrand(double aprime,void *p){
    Cosmology_FLRW *p1=(Cosmology_FLRW *)p;
    if(aprime==0.){return 0.;}
    double ha_h0=sqrt(p1->Ea(aprime));
    return 1./pow(aprime*ha_h0,3.);
}
