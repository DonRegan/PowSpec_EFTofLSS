//
//  PowDM_SPT.cpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 22/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#include "PowDM_SPT.hpp"
using namespace std;


//Wiggle Option... not currently being used
/*********************/
/***P13 calculation***/
/*********************/
p13CalcDM::p13CalcDM(import_class imports,Coeffs &coeffsIn,int numevalsIn,NoWiggle_option wigOptIn,which_scale scaleIn,EdSoption edsoptIn):coeffs(coeffsIn),numevals(numevalsIn),scale(scaleIn),wigOpt(wigOptIn),edsopt(edsoptIn){
    initData initdata=imports.get_initialdata();
    finalData finaldata=imports.get_finaldata();
    z=finaldata.zfinal;
    kmin=finaldata.kminCalc;
    kmax=finaldata.kmaxCalc;
    PowL= PowSpec(initdata.initpow_file).get_PowSpec();
    //Change this to the nowiggle version if told to do so!!
    if(wigOpt==NoWiggle_option::On){
        Splining PowLnw=get_Power_nowiggle(initdata.zinit, PowL,coeffs.get_cosFLRW(),0.25 );
        PowL=PowLnw;
    }

    vector<double>yfinal(16);
    for (int i=0; i<16; i++) {
        if(edsopt==EdSoption::On){
            yfinal[i]=coeffs.get_yfinalEdS(i);
        }else{
            yfinal[i]=coeffs.get_yfinal(i);
        }
    }
    g=yfinal[0],A=yfinal[2],B=yfinal[4],D=yfinal[6],E=yfinal[8],F=yfinal[10],G=yfinal[12],J=yfinal[14];
    timedep_ct_dd=-18*D-28.*E+7.*F+2.*G+13.*J;
    J0=vector<Splining>(6),J1=vector<Splining>(6),J2=vector<Splining>(6),J3=vector<Splining>(6),J4=vector<Splining>(6);
    create_vec(kmin, kmax, numevals, kvals, scale);
    boost::timer t;
    calc_Jn_vectors();
    P13DM=calc_P13KKtypical(D-J,E,F+J,G,-J,J/2.);
    std::cout<<"calculations of P13_DM done in  "<<t.elapsed()<<" seconds" <<std::endl;
    p13_DM_spline=Splining(kvals, P13DM);
}

void p13CalcDM::calc_Jn_vectors(){
    vector<double>J0vec(numevals),J1vec(numevals),J2vec(numevals),J3vec(numevals),J4vec(numevals);
    double kmin4int=PowL.get_xmin(),kmax4int=kmax;
    double tol=1E-6;//need it to be more accurate than the P22 calculation as we'll multiply by PL(k) after
    for (int j=0; j<6; j++) {
        for (int i=0; i<numevals; i++) {
            int nind0=0,nind1=1,nind2=2,nind3=3,nind4=4;
            struct leg_paramsp13 params0={kvals[i],j,nind0,PowL};
            struct leg_paramsp13 params1={kvals[i],j,nind1,PowL};
            struct leg_paramsp13 params2={kvals[i],j,nind2,PowL};
            struct leg_paramsp13 params3={kvals[i],j,nind3,PowL};
            struct leg_paramsp13 params4={kvals[i],j,nind4,PowL};
            J0vec[i]=(integrate(kmin4int, kmax4int, &Jintegrand, &params0,tol))/(2*M_PI*M_PI);
            J1vec[i]=(integrate(kmin4int, kmax4int, &Jintegrand, &params1,tol))/(2*M_PI*M_PI);
            J2vec[i]=(integrate(kmin4int, kmax4int, &Jintegrand, &params2,tol))/(2*M_PI*M_PI);
            J3vec[i]=(integrate(kmin4int, kmax4int, &Jintegrand, &params3,tol))/(2*M_PI*M_PI);
            J4vec[i]=(integrate(kmin4int, kmax4int, &Jintegrand, &params4,tol))/(2*M_PI*M_PI);
        }
        J0[j]=Splining(kvals, J0vec);J1[j]=Splining(kvals, J1vec);J2[j]=Splining(kvals, J2vec);
        J3[j]=Splining(kvals, J3vec);J4[j]=Splining(kvals, J4vec);
    }
    q0PqInt=integrate(kmin4int,kmax4int,&q0PqIntegrand,&PowL);
    q0PqIntIR=integrate(kmin4int,0.066,&q0PqIntegrand,&PowL);
    q2PqInt=integrate(kmin4int,kmax4int,&q2PqIntegrand,&PowL);
}

vector<double> p13CalcDM::calc_P13KKtypical(double A1,double B1,double C1,double D1,double E1,double G1){
    vector<double>P13KKvals(numevals);
    double val1= (-A1/6. - B1/3. +5.*C1/12. +D1/2. +E1 +4.*G1/3.)*q2PqInt,val2= (-A1/6. - B1/3. - C1/12. - D1/6.)*q0PqInt;
    for (int i=0; i<numevals; i++) {
        double k=kvals[i];
        double k2=k*k,k3=k2*k;
        P13KKvals[i]=val1+k2*val2+k2*J0[2].get_yval(k)*(A1/2. + 13*B1/15. - C1/12. - D1/6. - E1/6. - G1/3.)+J0[4].get_yval(k)*(A1/6. + B1/3. - 7*C1/12. - 5*D1/6. - 7*E1/6. - 5*G1/3.) +k3*J1[1].get_yval(k)*(-3*A1/10. - 3*B1/5.) +k*J1[3].get_yval(k)*(-13*A1/10. - 11*B1/5. + 3*C1/4. + 13*D1/10. + 3*E1/2. + 13*G1/5.) +(J1[5].get_yval(k)/k)*(C1/4. + D1/2. + E1/2. + G1) +k2*J2[2].get_yval(k)*(A1 + 40*B1/21. - C1/6. - D1/3. - E1/3. - 2*G1/3.) +J2[4].get_yval(k)*(A1/3. + 2*B1/3. - C1/6. - 2*D1/3. - E1/3. - 4*G1/3.) +k3*J3[1].get_yval(k)*(-A1/5. - 2*B1/5.) + k*J3[3].get_yval(k)*(-A1/5. - 4*B1/5. + D1/5. + 2*G1/5.) +k2*J4[2].get_yval(k)*(8/35.)*B1;
        P13KKvals[i]*=8.*g*PowL.get_yval(k);//extra factor of 2 here P13->P13new=2P13
    }
    return P13KKvals;
}
double Fkernel(int n,double qval,double kval){//Int dx j_n(k x) j_n(q x)
    double argval=sqrt((1 + kval/qval)*(1 + kval/qval)/((1 - kval/qval) *(1 - kval/qval)));
    if(sqrt((1 - kval/qval) *(1 - kval/qval))<1.0E-4){argval=0.;}
    switch (n) {
        case 0:
            return log(argval)/(2.*kval*qval);
            break;
        case 1:
            return (-kval*qval + (kval*kval + qval*qval)*log(argval)/2.)/(2.*kval*kval*qval*qval);
            break;
        case 2:
            return (-3.*kval*qval*(kval*kval + qval*qval) + (3*pow(kval,4) + 2.*kval*kval*qval*qval + 3.*pow(qval,4))*log(argval)/2.)/(8.*pow(kval*qval,3.));
            break;
        case 3:
            return (-kval*qval*(15.*pow(kval,4) + 14. *kval*kval*qval*qval + 15.*pow(qval,4)) + 3.*(5.*pow(kval,6) + 3.*pow(kval*kval*qval,2) + 3.*pow(kval*qval*qval,2) + 5.*pow(qval,6))*log(argval)/2.)/(48.*pow(kval*qval,4));
            break;
        case 4:
            return (-5.*kval*qval*(kval*kval + qval*qval)*(21.*pow(kval,4) - 2.*kval*kval*qval*qval + 21.*pow(qval,4)) +
                    3.*(35.*pow(kval,8.) + 20.*pow(kval,6)*qval*qval + 18.*pow(kval*qval,4) + 20.*kval*kval*pow(qval,6) + 35.*pow(qval,8.))*log(argval)/2.)/(384.*pow(kval*qval,5));
            break;
        default:
            cout<<"Have a problem with the index of the Fkernal\n";
            return 0.;
            break;
    }
}
double Jintegrand(double qval,void *p){
    struct leg_paramsp13 *params =(struct leg_paramsp13 *)p;
    double kval=params->k;
    return params->PowL.get_yval(qval)*pow(qval,params->qindex)*Fkernel(params->nindex,qval,kval);
}
double q2PqIntegrand(double qval,void *p){
    Splining *PowL=(Splining *)p;
    return PowL->get_yval(qval)*qval*qval/(2.*M_PI*M_PI);
}
double q0PqIntegrand(double qval,void *p){
    Splining *PowL=(Splining *)p;
    return PowL->get_yval(qval)/(2.*M_PI*M_PI);
}

/*********************/
/***P22 calculation***/
/*********************/
p22calcDM::p22calcDM(import_class imports,Coeffs &coeffsIn,int numevalsIn,NoWiggle_option wigOptIn,which_scale scaleIn,EdSoption edsoptIn):coeffs(coeffsIn),numpts(numevalsIn),scale(scaleIn),wigOpt(wigOptIn),edsopt(edsoptIn){
    initData initdata=imports.get_initialdata();
    finaldata=imports.get_finaldata();
    outputbase=finaldata.baseoutputs;
    
    z=finaldata.zfinal;
    kmin=finaldata.kminCalc;
    kmax=finaldata.kmaxCalc;
    PowL= PowSpec(initdata.initpow_file).get_PowSpec();
    //Change this to the nowiggle version if told to do so!!
    if(wigOpt==NoWiggle_option::On){
        Splining PowLnw=get_Power_nowiggle(initdata.zinit, PowL,coeffs.get_cosFLRW(),0.25 );
        PowL=PowLnw;
    }
    
    vector<double>yfinal(16);
    for (int i=0; i<16; i++) {
        if(edsopt==EdSoption::On){
            yfinal[i]=coeffs.get_yfinalEdS(i);
        }else{
            yfinal[i]=coeffs.get_yfinal(i);
        }
    }
    g=yfinal[0],A=yfinal[2],B=yfinal[4];
    create_vec(kmin, kmax, numpts, kvals, scale);
    P22_DM=vector<double>(numpts);
    boost::timer t;
    for (int i=0; i<numpts; i++) {
        P22_DM[i]=integrate_over_q(kvals[i], kmin, kmax, &P22DMintegrand_r, &g);
    }
    std::cout<<"calculations of P22_DM done in  "<<t.elapsed()<<" seconds" <<std::endl;
    
    P22_DM_spline=Splining(kvals, P22_DM);
}

double p22calcDM::integrate_over_q(double k,double lowerlim,double upperlim,double (*func)(double,void *),void *p){//int dq dr func(r,...) where ... has the q and k value. |k-q| <=r<=k+q
    struct paramsp22b parsb{k,k,PowL,func,this};
    return integrate(lowerlim,upperlim,&integrate_over_r,&parsb);
}
double P22DMintegrand_r(double r,void *p){
    struct paramsp22b *pars=(struct paramsp22b *)p;
    double k=pars->k;
    
    double q=pars->q;
    double kmin=pars->PowL.get_xmin(),kmax=pars->PowL.get_xmax();
    double val;
    double A=pars->p22class->A ,B=pars->p22class->B;
    //P22KK terms are already symmetric in q and r
  
    val= P22KKtypical_integrand(k,q,r,A,B,A,B);
    if(r<kmin||r>kmax){return 0.;}
    return val*pars->PowL.get_yval(r);
}

double P22KKtypical_integrand(double k,double q,double r,double A1,double B1,double A2,double B2){
    double r2=r*r,q2=q*q;
    double I000=Fabrik000(k, q, r),I011=Fabrik011(k,q,r),I022=Fabrik022(k,q,r),I033=Fabrik033(k,q,r),I044=Fabrik044(k,q,r);
    
    double val=2*r2*q2*((7*A1*A2/6. + (5/3.)*(A1*B2+A2*B1) + (38/15.)*B1*B2)*I000 + (A1*A2/12. + (1/6.)*(A1*B2+A2*B1) + (1/3.)*B1*B2)*(q2/r2 + r2/q2)*I000 - (1/5.)*(5.*A1*A2+9.*(A1*B2+A2*B1)+16.*B1*B2)*(q/r + r/q)*I011 + (A1*A2/3. + (4./3.)*(A1*B2+A2*B1) + (68/21.)*B1*B2)*I022 + (A1*A2/6. + (1/3.)*(A1*B2+A2*B1) + (2/3.)*B1*B2)*(q2/r2 + r2/q2)*I022 - (1./5.)*(A1*B2+A2*B1+4.*B1*B2)*(q/r + r/q)*I033 + (8*B1*B2/35.)*I044);
    return val;
}

double integrate_over_r(double q,void *p){//double k,double q,double lowerlim,double upperlim,double (*func)(double,void *),void *p){//int dq dr func(r,...) where ... has the q and k value. |k-q| <=r<=k+q
    struct paramsp22b *parsb=(struct paramsp22b *)p;
    double k=parsb->k;
    
    struct paramsp22b parsc={k,q,parsb->PowL,parsb->func,parsb->p22class};
    double upperlim=k+q,lowerlim=sqrt((k-q)*(k-q));
    return integrate(lowerlim,upperlim,parsb->func,&parsc)*parsb->PowL.get_yval(q)/pow(M_PI,3);
}
