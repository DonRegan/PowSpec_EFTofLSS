//
//  PowSpecRSD_SPT.cpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 23/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#include "PowSpecRSD_SPT.hpp"
using namespace std;

/*********************/
/***P13 calculation***/
/*********************/
p13CalcRSD::p13CalcRSD(p13CalcDM &p13DM,Splining &PowLIn):PowL(PowLIn){
    //just copy overstuff from p13 calculation
    numevals=p13DM.numevals;
    kmin=p13DM.kmin;
    kmax=p13DM.kmax;
    coeffs=p13DM.coeffs;
    scale=p13DM.scale;
    wigOpt=p13DM.wigOpt;
    z=p13DM.z;

    set_coeffs(p13DM.edsopt);
    J0=p13DM.J0,J1=p13DM.J1,J2=p13DM.J2,J3=p13DM.J3,J4=p13DM.J4;
    kvals=p13DM.kvals;
    P13dsds=vector<Splining>(4);
    boost::timer t;
    calc_P13_contribs(p13DM);
    std::cout<<"calculations of P13_RSD done in  "<<t.elapsed()<<" seconds" <<std::endl;
}

void p13CalcRSD::set_coeffs(EdSoption edsopt){
    vector<double>yfinal(16);
    for (int i=0; i<16; i++) {
        if(edsopt==EdSoption::On){
            yfinal[i]=coeffs.get_yfinalEdS(i);
        }else{
            yfinal[i]=coeffs.get_yfinal(i);
        }
    }
    g=yfinal[0],A=yfinal[2],B=yfinal[4],D=yfinal[6],E=yfinal[8],F=yfinal[10],G=yfinal[12],J=yfinal[14];
    f=-(1+z)*yfinal[1]/g,Abarprime=-(1+z)*yfinal[3],Bbarprime=-(1+z)*yfinal[5],Dbarprime=-(1+z)*yfinal[7],Ebarprime=-(1+z)*yfinal[9],Fbarprime=-(1+z)*yfinal[11],Gbarprime=-(1+z)*yfinal[13],Jbarprime=-(1+z)*yfinal[15];
    H=Dbarprime-Jbarprime,N=Ebarprime,P=Fbarprime+Jbarprime-f*g*A,Q=Gbarprime-f*g*B,R=f*g*A-g*Abarprime-2*Jbarprime+f*g*g*g,S=Jbarprime-g*Bbarprime+f*g*B;
    K=Bbarprime,L=Abarprime-f*g*g;
    timedep_ct_dd=-18*D-28.*E+7.*F+2.*G+13.*J;
    timedep_ct_dd_deriv=-18*Dbarprime-28.*Ebarprime+7.*Fbarprime+2.*Gbarprime+13.*Jbarprime;
    timedep_ct_tt=-18.*H-28.*N+7.*P+2.*Q+12.*R+12.*S;
}
void  p13CalcRSD::calc_P13_contribs(p13CalcDM &p13DM){
    double q0PqInt=p13DM.q0PqInt;
    double q2PqInt=p13DM.q2PqInt;
    //First for the KxK terms
    vector<double>P13KKvalsmu0(numevals),P13KKvalsmu2(numevals),P13KKvalsmu4(numevals);
    P13KKvalsmu0=p13DM.calc_P13KKtypical(D-J,E,F+J,G,-J,J/2.);
    P13KKvalsmu4=p13DM.calc_P13KKtypical(H,N,P,Q,R/2.,S/2.);//must multiply by f
    for (int i=0; i<numevals; i++) {//must divide by two for the following because I'd included the extra factor in calc_P13KKtypical
        P13KKvalsmu2[i]=(P13KKvalsmu0[i]*f+P13KKvalsmu4[i])/2;
        P13KKvalsmu4[i]*=f/2;
        P13KKvalsmu0[i]=P13KKvalsmu0[i]/2;
    }
    //Second for Kx[v delta]
    vector<double>P13Kvdmu2(numevals),P13Kvdmu4(numevals);
    for (int i=0; i<numevals; i++) {
        double k=kvals[i],k2=kvals[i]*kvals[i];
        P13Kvdmu2[i]=(2*g*g*((L + 4*K/3.)*(k2*J0[2].get_yval(k) - k*J1[3].get_yval(k)) + (L/2. + K)*((J0[4].get_yval(k) + 2.*J2[4].get_yval(k))/3. +(k2/3.)*(J0[2].get_yval(k) + 2*J2[2].get_yval(k)) - k*(J1[3].get_yval(k) + k2*J1[1].get_yval(k))) + (2*K/3.)*(k2*J2[2].get_yval(k) - (2/5.)*k*J1[3].get_yval(k) - 3*k*J3[3].get_yval(k)/5.) - (1/3.)*(f*A/2. + f*B)*(k2*q0PqInt + q2PqInt) ))*PowL.get_yval(k);
        P13Kvdmu4[i]=f*P13Kvdmu2[i];
    }
    //Third for K x[v v]
    vector<double>P13Kv2mu2(numevals),P13Kv2mu4(numevals),P13Kv2mu6(numevals);
    for (int i=0; i<numevals; i++) {
        double k=kvals[i],k2=kvals[i]*kvals[i];
        P13Kv2mu2[i]=(-2*f*g*g*k2*( (L + 4*K/3.)*(J0[2].get_yval(k) - J2[2].get_yval(k))/3. - (L/2. + K)*(k2*(J1[1].get_yval(k) - J3[1].get_yval(k)) + J1[3].get_yval(k) - J3[3].get_yval(k))/(5.*k) + (2*K/3.)*(-1/3.)*(-J2[2].get_yval(k) + J0[2].get_yval(k)/5. + 2*J2[2].get_yval(k)/7. + 18*J4[2].get_yval(k)/35.)))*PowL.get_yval(k);
        P13Kv2mu6[i]=(-2*f*f*g*g*k2*( (L + 4*K/3.)*(J2[2].get_yval(k) - k*J1[1].get_yval(k)) + (L/2. + K)*((J0[2].get_yval(k) + k2*J0[0].get_yval(k))/3. +2*(J2[2].get_yval(k) + k2*J2[0].get_yval(k))/3. - (2*(J1[3].get_yval(k) + k2*J1[1].get_yval(k)))/(5.*k) - (3*(J3[3].get_yval(k) + k2*J3[1].get_yval(k)))/(5.*k)) + (2*K/3.)*(J0[2].get_yval(k)/5. + 2*J2[2].get_yval(k)/7. +18*J4[2].get_yval(k)/35. - 2*k*J1[1].get_yval(k)/5. - 3*k*J3[1].get_yval(k)/5.)))*PowL.get_yval(k);
        P13Kv2mu4[i]=f*P13Kv2mu2[i]+P13Kv2mu6[i]/f;
    }
    //Fourth for K x[v^2 delta + v^3] = Kx[other]
    vector<double>P13Kothermu2(numevals),P13Kothermu4(numevals),P13Kothermu6(numevals);
    for (int i=0; i<numevals; i++) {
        P13Kothermu2[i]=(-f*g*g*f*g*g*kvals[i]*kvals[i]*q0PqInt/6.)*PowL.get_yval(kvals[i]);
        P13Kothermu4[i]=2*f*P13Kothermu2[i];
        P13Kothermu6[i]=f*f*P13Kothermu2[i];
    }
    
    //Finally let's gather everything up and create the splines
    vector<double>P13mu0(numevals),P13mu2(numevals),P13mu4(numevals),P13mu6(numevals),P13thet(numevals);
    for (int i=0; i<numevals; i++) {//Multiplying by factor of 2 so that P=PL+P13+P22
        P13mu0[i]=2*P13KKvalsmu0[i];
        P13mu2[i]=2*(P13KKvalsmu2[i]+P13Kvdmu2[i]+P13Kv2mu2[i]+P13Kothermu2[i]);
        P13mu4[i]=2*(P13KKvalsmu4[i]+P13Kvdmu4[i]+P13Kv2mu4[i]+P13Kothermu4[i]);
        P13mu6[i]=2*(P13Kv2mu6[i]+P13Kothermu6[i]);
        P13thet[i]=2.*P13KKvalsmu4[i];
    }
    P13dsds[0]=Splining(kvals, P13mu0);
    P13dsds[1]=Splining(kvals, P13mu2);
    P13dsds[2]=Splining(kvals, P13mu4);
    P13dsds[3]=Splining(kvals, P13mu6);
    P13theta=Splining(kvals,P13thet);
}






/*********************/
/***P22 calculation***/
/*********************/
p22CalcRSD::p22CalcRSD(p22calcDM &p22DM,Splining &PowLIn,calc_read calcoption):PowL(PowLIn){
    //just copy overstuff from p13 calculation
    numevals=p22DM.numpts;
    kmin=p22DM.kmin;
    kmax=p22DM.kmax;
    kvals=p22DM.kvals;
    coeffs=p22DM.coeffs;
    scale=p22DM.scale;
    wigOpt=p22DM.wigOpt;
    finaldata=p22DM.finaldata;
    z=p22DM.z;
    
    set_coeffs(p22DM.edsopt);
        
    if(calcoption==calc_read::calc){
        boost::timer t;
        calc_p22_RSD();
        std::cout<<"calculations of P22_RSD done in  "<<t.elapsed()<<" seconds" <<std::endl;
        output_p22_RSD();
    }else{
        input_p22_RSD();
    }
}
void p22CalcRSD::set_coeffs(EdSoption edsopt){
    vector<double>yfinal(16);
    for (int i=0; i<16; i++) {
        if(edsopt==EdSoption::On){
            yfinal[i]=coeffs.get_yfinalEdS(i);
        }else{
            yfinal[i]=coeffs.get_yfinal(i);
        }
    }
    g=yfinal[0],A=yfinal[2],B=yfinal[4];
    f=-(1+z)*yfinal[1]/g,Abarprime=-(1+z)*yfinal[3],Bbarprime=-(1+z)*yfinal[5];
    K=Bbarprime,L=Abarprime-f*g*g;
}
double p22CalcRSD::integrate_over_q(double k,int mupow,double lowerlim,double upperlim,double (*func)(double,void *),void *p){//int dq dr func(r,...) where ... has the q and k value. |k-q| <=r<=k+q
    struct leg_paramsp22b parsb{k,k,mupow,PowL,func,this};
    return integrate(lowerlim,upperlim,&integrateRSD_over_r,&parsb);
}

void p22CalcRSD::output_p22_RSD(){
    string filename=finaldata.baseoutputs+"P22RSD_z"+to_string(z);
    string NoWiggle="";
    if(wigOpt==NoWiggle_option::On){NoWiggle="_nowiggle";}
    string filenameRSD=filename+NoWiggle+".csv";
    std::ofstream myfile(filenameRSD);
    if (!myfile.is_open()) {
        std::cout<<"Problem opening file "<<filenameRSD<<std::endl;
        return;
    }
    myfile<<"kvalue"<<","<<"P22_ell0"<<","<<"P22_ell2"<<","<<"P22_ell4"<<","<<"P22_ell6"<<","<<"P22_ell8"<<"\n";
    for (int i=0; i<numevals; i++) {
        myfile<<kvals[i];
        for (int mupow=0; mupow<5; mupow++) {
            myfile<<","<<P22RSD_mu2n[mupow].get_yval(kvals[i])<<std::scientific;
        }
        myfile<<"\n";
    }
    myfile.close();
    //now output the theta file
    string filenameth=filename+NoWiggle+"theta.csv";
    std::ofstream myfile2(filenameth);
    if (!myfile2.is_open()) {
        std::cout<<"Problem opening file "<<filenameth<<std::endl;
        return;
    }
    myfile2<<"kvalue"<<","<<"Pkvalue"<<"\n";
    for (int i=0; i<numevals; i++) {
        myfile2<<kvals[i]<<","<<P22theta.get_yval(kvals[i])<<std::scientific<<"\n";
    }
    myfile2.close();
}
void p22CalcRSD::input_p22_RSD(){
    string filename=finaldata.baseoutputs+"P22RSD_z"+to_string(z);
    string NoWiggle="";
    if(wigOpt==NoWiggle_option::On){NoWiggle="_nowiggle";}
    string filenameRSD=filename+NoWiggle+".csv";
    string filenameth=filename+NoWiggle+"theta.csv";
    vector<double>kvalsIn,Pk_ell0,Pk_ell2,Pk_ell4,Pk_ell6,Pk_ell8,Pk_theta;
    ifstream infile(filenameRSD);
    if (!infile.is_open()) {
        std::cout<<"Problem opening file "<<filenameRSD<<std::endl;
        return;
    }
    std::string name1,name2;

    getline(infile,name2);
    double k1,Pk1;
    int numpoints=0;
    while(getline(infile,name1)){
        std::stringstream ss(name1);
        ss>>k1;
        ss.get();
        ss>>Pk1;
        kvalsIn.push_back(k1);
        Pk_ell0.push_back(Pk1);
        ss.get();
        ss>>Pk1;
        Pk_ell2.push_back(Pk1);
        ss.get();
        ss>>Pk1;
        Pk_ell4.push_back(Pk1);
        ss.get();
        ss>>Pk1;
        Pk_ell6.push_back(Pk1);
        ss.get();
        ss>>Pk1;
        Pk_ell8.push_back(Pk1);
        numpoints++;
    }
    infile.close();

    std::ifstream myfile2(filenameth);
    if (!myfile2.is_open()) {
        std::cout<<"Problem opening file "<<filenameth<<std::endl;
        return;
    }

    getline(myfile2,name1);
    numpoints=0;
    while(getline(myfile2,name1)){
        std::stringstream ss(name1);
        ss>>k1;
        ss.get();
        ss>>Pk1;
        Pk_theta.push_back(Pk1);
        numpoints++;
    }
    myfile2.close();
    //overwrite kvals and numevals to reflect these values
    kvals=kvalsIn;
    kmin=kvals[0],kmax=kvals[kvals.size()-1];
    numevals=numpoints;
    //Fit all of these to Splines
    P22RSD_mu2n=vector<Splining>(5);
    P22RSD_mu2n[0]=Splining(kvals, Pk_ell0);
    P22RSD_mu2n[1]=Splining(kvals, Pk_ell2);
    P22RSD_mu2n[2]=Splining(kvals, Pk_ell4);
    P22RSD_mu2n[3]=Splining(kvals, Pk_ell6);
    P22RSD_mu2n[4]=Splining(kvals, Pk_ell8);
    P22theta=Splining(kvals,Pk_theta);

    return;
}
void p22CalcRSD::calc_p22_RSD(){
    vector<double>P22RSD_mu0(numevals),P22RSD_mu2(numevals);
    vector<double>P22RSD_mu4(numevals),P22RSD_mu6(numevals);
    vector<double>P22RSD_mu8(numevals),P22thet(numevals);

    for (int i=0; i<numevals; i++) {
        P22RSD_mu0[i]=integrate_over_q(kvals[i], 0, kmin, kmax, &P22RSD_integrand, &f);
        P22RSD_mu2[i]=integrate_over_q(kvals[i], 2, kmin, kmax, &P22RSD_integrand, &f);
        P22RSD_mu4[i]=integrate_over_q(kvals[i], 4, kmin, kmax, &P22RSD_integrand, &f);
        P22RSD_mu6[i]=integrate_over_q(kvals[i], 6, kmin, kmax, &P22RSD_integrand, &f);
        P22RSD_mu8[i]=integrate_over_q(kvals[i], 8, kmin, kmax, &P22RSD_integrand, &f);
        P22thet[i]=integrate_over_q(kvals[i],4,kmin,kmax,&P22KKintegrand_r,&f);
    }
    P22RSD_mu2n=vector<Splining>(5);
    P22RSD_mu2n[0]=Splining(kvals, P22RSD_mu0);
    P22RSD_mu2n[1]=Splining(kvals, P22RSD_mu2);
    P22RSD_mu2n[2]=Splining(kvals, P22RSD_mu4);
    P22RSD_mu2n[3]=Splining(kvals, P22RSD_mu6);
    P22RSD_mu2n[4]=Splining(kvals, P22RSD_mu8);
    P22theta=Splining(kvals,P22thet);
}

double P22KKintegrand_r(double r,void *p){
    struct leg_paramsp22b *pars=(struct leg_paramsp22b *)p;
    double k=pars->k;
    int mupow=pars->mupow;

    double q=pars->q;
    double kmin=pars->PowL.get_xmin(),kmax=pars->PowL.get_xmax();
    double val;
    double A=pars->p22class->A ,B=pars->p22class->B,L=pars->p22class->L,K=pars->p22class->K;
    //P22KK terms are already symmetric in q and r
    switch (mupow) {
        case 0:
            val= P22KKtypical_integrand(k,q,r,A,B,A,B);
            break;
        case 2:
            val= 2.*P22KKtypical_integrand(k,q,r,A,B,L,K);
            break;
        case 4:
            val= P22KKtypical_integrand(k,q,r,L,K,L,K);
            break;
        default:
            val= 0.;
            break;
    }
    if(r<kmin||r>kmax){return 0.;}
    return val*pars->PowL.get_yval(r);
}

double P22vdvd_integrand(double r,void *p){
    struct leg_paramsp22b *pars=(struct leg_paramsp22b *)p;
    double k=pars->k;
    int mupow=pars->mupow;
    double q=pars->q;
    double kmin=pars->PowL.get_xmin(),kmax=pars->PowL.get_xmax();
    double val;
    double f=pars->p22class->f,g=pars->p22class->g;
    
    double q2=q*q,r2=r*r,k2=k*k,q3=q2*q,r3=r2*r,k3=k2*k,q4=q2*q2,r4=r2*r2,k4=k2*k2;
    switch (mupow) {
        case 2:
            val=-((M_PI* (k - q - r) *(q - r)*(q-r)* (k + q - r)* (k - q + r) *(q + r)*(q+r)* (k + q + r))/(64* k3 *q3 *r3));
            break;
        case 4:
            val=(M_PI* (3.* pow(q2 - r2,4) - 6.* k2* (q2 - r2)*(q2-r2)* (q2 + r2) +k4* (3.* q4 + 2.* q2* r2 + 3.* r4)))/(64.* k3* q3* r3);
            break;
        default:
            val= 0.;
            break;
    }
    if(r<kmin||r>kmax){return 0.;}
    return val*pars->PowL.get_yval(r)*pow(k*f*g*g,2);
}

double P22v2v2_integrand(double r,void *p){
    struct leg_paramsp22b *pars=(struct leg_paramsp22b *)p;
    double k=pars->k;
    int mupow=pars->mupow;
    double q=pars->q;
    double kmin=pars->PowL.get_xmin(),kmax=pars->PowL.get_xmax();
    double val;
    double f=pars->p22class->f,g=pars->p22class->g;
    
    double k2=k*k,k4=k2*k2,k6=k4*k2,k8=k4*k4,r2=r*r,r4=r2*r2,q2=q*q,q4=q2*q2;
    double r3=r2*r,q3=q2*q,k5=k4*k;
    switch (mupow) {//already symmetrised
        case 4:
            val=(M_PI* (3 *k8 - 12 *k6 *q2 + 18 *k4 *q4 - 12 *k6 *r2 +12 *k4 *q2 *r2 + 18 *k4 *r4 + 3 *pow(q2 - r2,4) -12 *k2 *(q2 - r2)*(q2 - r2)* (q2 + r2)))/(512 *k5 *q3* r3);
            break;
        case 6:
            val=(M_PI* (2* k8 + 8.* k6* q2 - 52* k4* q4 + 8 *k6* r2 + 8 *k4* q2* r2 -52 *k4* r4 - 30 *pow(q2 - r2,4) +72* k2* (q2 - r2)*(q2-r2)* (q2 + r2)))/(512.* k5* q3 *r3);
            break;
        case 8:
            val=(M_PI* (3. *k8 + 18.* k4* q4 + 12* k4 *q2* r2 + 18 *k4* r4 +35 *pow(q2 - r2,4) + 4* k6* (q2 + r2) -60* k2* (q2 - r2)*(q2-r2)* (q2 + r2)))/(512* k5 *q3* r3);
            break;
        default:
            val= 0.;
            break;
    }
    if(r<kmin||r>kmax){return 0.;}
    return val*pars->PowL.get_yval(r)*pow(k*f*g,4)/2.;
}
double P22Kvd_integrand(double r,void *p){
    struct leg_paramsp22b *pars=(struct leg_paramsp22b *)p;
    double k=pars->k;
    int mupow=pars->mupow;
    double q=pars->q;
    double kmin=pars->PowL.get_xmin(),kmax=pars->PowL.get_xmax();
    double val;
    double f=pars->p22class->f,g=pars->p22class->g,A=pars->p22class->A,B=pars->p22class->B;
    double fA=pars->p22class->Abarprime/A,fB=pars->p22class->Bbarprime/B;
    
    double k2=k*k,q2=q*q,r2=r*r,k4=k2*k2,r4=r2*r2,q4=q2*q2;
    switch (mupow) {
        case 2:
            val=(M_PI* (B* k4 +A* k2* q2 - (A + B)* q4 + (A *k2 + 2 *(A + B)* q2)* r2 - (A +B)* r4)* (-(q2 - r2)*(q2-r2) + k2* (q2 + r2)))/(64* k2* q4* r4);
            break;
        case 4:
            val=(M_PI* (-(q2 - r2)*(q2-r2) +k2* (q2 + r2))* (B* fB* (k4 - (q2 - r2)*(q2-r2)) + (A *fA -f* g*g)* (-(q2 - r2)*(q2-r2) + k2* (q2 + r2))))/(64* k2* q4* r4);
            break;
        default:
            val= 0.;
            break;
    }
    if(r<kmin||r>kmax){return 0.;}
    return val*pars->PowL.get_yval(r)*4.*f*g*g*k*r*q;
}
double P22Kv2_integrand(double r,void *p){
    struct leg_paramsp22b *pars=(struct leg_paramsp22b *)p;
    double k=pars->k;
    int mupow=pars->mupow;
    double q=pars->q;
    double kmin=pars->PowL.get_xmin(),kmax=pars->PowL.get_xmax();
    double val;
    double f=pars->p22class->f,g=pars->p22class->g,A=pars->p22class->A,B=pars->p22class->B;
    double fA=pars->p22class->Abarprime/A,fB=pars->p22class->Bbarprime/B;
    
    double q2=q*q,q4=q2*q2,r2=r*r,r4=r2*r2,k2=k*k,k4=k2*k2;
    
    switch (mupow) {//ensure symmetrised
        case 2:
            val=(M_PI *(k - q - r) *(k + q - r) *(k - q + r) *(k + q +r)* (B *k4 +A * k2 *q2 - (A + B) *q4 + (A *k2 + 2 *(A + B) *q2) *r2 - (A + B) *r4))/(64 *k *q2*q* r2*r);
            break;
        case 4:
            val=(M_PI/(64. *k *q2*q *r2* r))* (B* (k4 - (q2 - r2)*(q2 - r2)) *((1 + fB) *k4 + (-3 + fB)* (q2 - r2)*(q2 - r2) -2 *(-1 + fB) *k2 *(q2 + r2)) + (-(q2 - r2)*(q2 - r2) +k2 *(q2 + r2)) *(-f *g*g* (k - q - r) *(k + q - r)* (k - q + r) *(k + q + r) +A *((1 + fA)* k4 + (-3 + fA) *(q2 - r2)*(q2 - r2) -2 *(-1 + fA)* k2* (q2 + r2))));
            break;
        case 6:
            val=( M_PI *(k4 - 3 *(q2 - r2)*(q2 - r2) +2 *k2* (q2 + r2))* (B*fB *(k4 - (q2 - r2)*(q2 - r2)) + (A* fA -f *g*g) *(-(q2 - r2)*(q2 - r2) +k2* (q2 + r2))))/(64 *k *q2*q* r2*r);
            break;
        default:
            val= 0.;
            break;
    }
    if(r<kmin||r>kmax){return 0.;}
    return val*pars->PowL.get_yval(r)*(f*g)*(f*g);
}
double P22v2vd_integrand(double r,void *p){
    struct leg_paramsp22b *pars=(struct leg_paramsp22b *)p;
    double k=pars->k;
    int mupow=pars->mupow;
    
    double q=pars->q;
    double kmin=pars->PowL.get_xmin(),kmax=pars->PowL.get_xmax();
    double val;
    double f=pars->p22class->f,g=pars->p22class->g;
    
    double q2=q*q,r2=r*r,k2=k*k,q3=q2*q,r3=r2*r,q4=q2*q2,r4=r2*r2,k4=k2*k2,k6=k4*k2;
    switch (mupow) {//ensure symmetrised
        case 4:
            val=(M_PI*(k - q - r)*(k + q - r)*(k - q + r)*(k + q +r)*(-3*(q2 - r2)*(q2-r2) + k2* (q2 + r2)))/(128.*k4*q3*r3);
            break;
        case 6:
            val=(M_PI*(5.* pow(q2 - r2,4) + k6* (q2 + r2) -9* k2* (q2 - r2)*(q2-r2)* (q2 + r2) +k4* (3* q4 + 2* q2* r2 + 3 *r4)))/(128.* k4* q3* r3);
            break;
        default:
            val= 0.;
            break;
    }
    if(r<kmin||r>kmax){return 0.;}
    return val*pars->PowL.get_yval(r)*2.*pow(f*g*k,3)*g;
}
double P22RSD_integrand(double r,void *p){
    return P22KKintegrand_r(r, p)+P22v2vd_integrand(r, p)+P22v2v2_integrand(r, p)+P22Kvd_integrand(r, p)+P22Kv2_integrand(r, p)+(P22vdvd_integrand(r, p));
}
double integrateRSD_over_r(double q,void *p){//double k,double q,double lowerlim,double upperlim,double (*func)(double,void *),void *p){//int dq dr func(r,...) where ... has the q and k value. |k-q| <=r<=k+q
    struct leg_paramsp22b *parsb=(struct leg_paramsp22b *)p;
    double k=parsb->k;
    int mupow=parsb->mupow;

    struct leg_paramsp22b parsc={k,q,mupow,parsb->PowL,parsb->func,parsb->p22class};
    double upperlim=k+q,lowerlim=sqrt((k-q)*(k-q));
    return integrate(lowerlim,upperlim,parsb->func,&parsc)*parsb->PowL.get_yval(q)/pow(M_PI,3);
}

