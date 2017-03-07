//
//  Wiggle_NoWiggle.cpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 22/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#include "Wiggle_NoWiggle.hpp"

using namespace std;
Splining get_Power_nowiggle(double z,Splining &PowSpec,Cosmology_FLRW cosFLRW,double lambda){
    double kmin=PowSpec.get_xmin(),kmax=PowSpec.get_xmax();
    vector<double>kvals,Pkvals_eh;
    int numpts=1000;
    vector<double>Pk_ratio(numpts);
    cosFLRW.create_powspec_eh(numpts, kmin, kmax, z, kvals, Pkvals_eh);
    for (int i=0; i<numpts; i++) {
        Pk_ratio[i]=PowSpec.get_yval(kvals[i])/Pkvals_eh[i]-1.;
    }
    Splining Pkratio(kvals, Pk_ratio);
    
    vector<double>Pk_nw(numpts),Pk_w(numpts);
    for (int i=0; i<numpts; i++) {
        double lambdanew=lambda*(1-0.5*i/numpts);//(1-0.7*i/numpts);//adding a slight scale dependence and increasing by up to 10%
        double k=kvals[i];
        struct params_wigglepower params={Pkratio,lambdanew,log10(k)};
        double lowerlim=log10(kmin*1.0001),upperlim=log10(kmax*0.9999);
        double denomin=integrate(lowerlim,upperlim, &gaussianfiltered_power_denomin, &params,1E-8);
        Pk_nw[i]=integrate(lowerlim,upperlim, &gaussianfiltered_power, &params,1E-8)*Pkvals_eh[i]/denomin+Pkvals_eh[i];
    }
    Splining PowSpec_nw(kvals, Pk_nw);
    return PowSpec_nw;

}
Splining get_Power_wiggle(Splining PowSpec,Splining PowerNoWiggle){
    double kmin=PowSpec.get_xmin(),kmax=PowSpec.get_xmax();
    int numpts=1000;
    vector<double>kvals,Pkvals_w(numpts);
    create_vec(kmin, kmax, numpts, kvals, logScale);
    for (int i=0; i<numpts; i++) {
        Pkvals_w[i]=PowSpec.get_yval(kvals[i])-PowerNoWiggle.get_yval(kvals[i]);
    }
    Splining PowSpec_w(kvals, Pkvals_w);
    return PowSpec_w;
}
double gaussianfiltered_power(double log10q,void *p){
    struct params_wigglepower *params=(struct params_wigglepower *)p;
    double log10k=params->log10k;
    double lambda=params->lambda;
    double q=pow(10,log10q);
    return (1./(sqrt(2.*M_PI)*lambda))*exp(-0.5*(log10k-log10q)*(log10k-log10q)/(lambda*lambda))*params->PowL.get_yval(q) ;
}
double gaussianfiltered_power_denomin(double log10q,void *p){
    struct params_wigglepower *params=(struct params_wigglepower *)p;
    double log10k=params->log10k;
    double lambda=params->lambda;
    return (1./(sqrt(2.*M_PI)*lambda))*exp(-0.5*(log10k-log10q)*(log10k-log10q)/(lambda*lambda));
}
