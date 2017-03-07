//
//  PowDM_SPT_Resummation.cpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 23/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#include "PowDM_SPT_Resummation.hpp"
using namespace std;
IR_Resum::IR_Resum(import_class imports,Coeffs &coeffsIn,int numevalsIn,which_scale scaleIn,EdSoption edsopt):coeffs(coeffsIn),numevals(numevalsIn),scale(scaleIn){

    finaldata=imports.get_finaldata();
    p13=p13CalcDM(imports, coeffs, numevals,NoWiggle_option::Off,scale,edsopt);
    p13_nw=p13CalcDM(imports, coeffs, numevals,NoWiggle_option::On,scale,edsopt);
    calcopt=calc_read::calc;
    if(finaldata.calcP22==false){
        calcopt=calc_read::read;
    }
    p22=p22calcDM(imports, coeffs, numevals,NoWiggle_option::Off,scale,edsopt);
    p22_nw=p22calcDM(imports, coeffs, numevals,NoWiggle_option::On,scale,edsopt);
    P13=p13.get_p13_spline(),P13_nw=p13_nw.get_p13_spline();
    P22=p22.get_p22_spline(),P22_nw=p22_nw.get_p22_spline();

    kmin=finaldata.kminCalc;
    kmax=finaldata.kmaxCalc;
    create_vec(kmin, kmax, numevals, kvals, scale);

    PowL=PowSpec(finaldata.linearpow_DM).get_PowSpec();
    PowL_nw=get_Power_nowiggle(finaldata.zfinal, PowL,coeffs.get_cosFLRW(),0.25 );
    
    initData initdata=imports.get_initialdata();
    PowLinit=PowSpec(initdata.initpow_file).get_PowSpec();
    PowLinit_nw=get_Power_nowiggle(initdata.zinit, PowLinit,coeffs.get_cosFLRW(),0.25 );

    Aval=SigmaAverage(PowL_nw, 10, 300);
    g=coeffs.get_yfinal(0);
    f=-(1+finaldata.zfinal)*coeffs.get_yfinal(1)/g;
    calc_IR_Resum_Linear();
    IR1loop=set_IR_Resum_to1Loop_DM();
    nonIR1loop=set_nonIR_Resum_to1Loop_DM();
}
void IR_Resum::calc_IR_Resum_Linear(){
    vector<double>IR_Lin(kvals.size());
    for (int i=0; i<kvals.size(); i++) {
        double k=kvals[i];
        double PLk=PowL.get_yval(k);
        double PLk_nw=(PowL_nw.get_yval(k));
        double PLk_w=PLk-PLk_nw;
        IR_Lin[i]=PLk_nw+exp(-k*k*Aval)*(PLk_w*(1+k*k*Aval));//note the extra bit here which is needed at one loop order a la Matsubara
    }
    IR_lin=Splining(kvals, IR_Lin);
}
Splining IR_Resum::set_IR_Resum_to1Loop_DM(){
    vector<double>IR_p13(kvals.size()),IR_p22(kvals.size()),IR_1loop(kvals.size());
    for (int i=0; i<kvals.size(); i++) {
        double k=kvals[i];
        double P22k=P22.get_yval(k),P22k_nw=P22_nw.get_yval(k);
        double P13k=P13.get_yval(k),P13k_nw=P13_nw.get_yval(k);

        double P22k_w=P22k-P22k_nw,P13k_w=P13k-P13k_nw;
        IR_p22[i]=P22k_nw+exp(-k*k*Aval)*P22k_w;
        IR_p13[i]=P13k_nw+exp(-k*k*Aval)*P13k_w;
    }
    for (int i=0; i<kvals.size(); i++) {
        double k=kvals[i];
        IR_1loop[i]=IR_lin.get_yval(k)+IR_p13[i]+IR_p22[i];
    }
    return Splining(kvals,IR_1loop);
}

Splining IR_Resum::set_nonIR_Resum_to1Loop_DM(){
    vector<double>nonIR(kvals.size());
    for (int i=0; i<kvals.size(); i++) {
        double k=kvals[i];
        double P22k=P22.get_yval(k),P13k=P13.get_yval(k);
        nonIR[i]=PowL.get_yval(k)+P13k+P22k;
    }
    return Splining(kvals,nonIR);
}








double SigmaAverage(Splining &PowL,double qmin,double qmax){
    return integrate(qmin, qmax, &SigmaVal, &PowL)/(qmax-qmin);
}
double SigmaVal(double q,void *p){
    Splining *PowL=(Splining *)p;
    struct params_Sigma params={q,PowL};
    double lowerlim=PowL->get_xmin()*1.001,upperlim=0.07;
    return integrate(lowerlim, upperlim, &SigmaIntegrand, &params);
}
double SigmaIntegrand(double pval,void *p){
    struct params_Sigma *params=(struct params_Sigma *)p;
    double qval=params->q;
    return params->PowL->get_yval(pval)*(1-gsl_sf_bessel_j0(pval*qval))/(6.*M_PI*M_PI);
}
