//
//  PowSpecRSD_SPT_resummation.cpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 23/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#include "PowSpecRSD_SPT_resummation.hpp"
using namespace std;

IR_Resum_RSD::IR_Resum_RSD(IR_Resum &irresum){
    
    kvals=irresum.kvals;
    p13CalcRSD p13RSD(irresum.p13,irresum.PowLinit); //must set the wiggle or no wiggle
    p13CalcRSD p13RSD_nw(irresum.p13_nw,irresum.PowLinit_nw);
    p22CalcRSD p22RSD(irresum.p22,irresum.PowLinit,irresum.calcopt);
    p22CalcRSD p22RSD_nw(irresum.p22_nw,irresum.PowLinit_nw,irresum.calcopt);
    
    P13dsds=vector<Splining>(4),P13dsds_nw=vector<Splining>(4),P13dsds_w=vector<Splining>(4);//for mu^0,mu^2,mu^4,mu^6
    P22dsds=vector<Splining>(5),P22dsds_nw=vector<Splining>(5),P22dsds_w=vector<Splining>(5);//for mu^0,mu^2,mu^4,mu^6,mu^8
    for (int n=0; n<4; n++) {
        P13dsds[n]=p13RSD.get_P13RSD(2*n);
        P22dsds[n]=p22RSD.get_P22RSD(2*n);
        P13dsds_nw[n]=p13RSD_nw.get_P13RSD(2*n);
        P22dsds_nw[n]=p22RSD_nw.get_P22RSD(2*n);

    }
    P22dsds[4]=p22RSD.get_P22RSD(8);
    P22dsds_nw[4]=p22RSD_nw.get_P22RSD(8);
    calc_wiggle_contribs();
    
    PowL=irresum.PowL;
    PowL_nw=irresum.PowL_nw;
    PowL_IR=irresum.IR_lin;
    Aval=irresum.Aval;
    g=irresum.g;
    f=irresum.f;
    calc_IR_Resum_Legendre0();
    calc_IR_Resum_Legendre2();
    calc_IR_Resum_Legendre4();
}
void IR_Resum_RSD::calc_wiggle_contribs(){
    vector<double>P13_w_vals(kvals.size()),P22_w_vals(kvals.size());
    for (int n=0; n<4; n++) {
        for (int i=0; i<kvals.size(); i++) {
            double k=kvals[i];
            P13_w_vals[i]=P13dsds[n].get_yval(k)-P13dsds_nw[n].get_yval(k);
            P22_w_vals[i]=P22dsds[n].get_yval(k)-P22dsds_nw[n].get_yval(k);
        }
        P13dsds_w[n]=Splining(kvals,P13_w_vals);
        P22dsds_w[n]=Splining(kvals,P22_w_vals);
    }
    for (int i=0; i<kvals.size(); i++) {
        double k=kvals[i];
        P22_w_vals[i]=P22dsds[4].get_yval(k)-P22dsds_nw[4].get_yval(k);
    }
    P22dsds_w[4]=Splining(kvals,P22_w_vals);
}
void IR_Resum_RSD::calc_IR_Resum_Legendre0(){
    vector<double>IR_Leg0(kvals.size());
    vector<double>Non_IR_Leg0(kvals.size());
    
    cout<<"getting the Legendre0 \n";
    for (int i=0; i<kvals.size(); i++) {
        double k=kvals[i];
        double Aksq=k*k*Aval;
        double Alpha=Aksq*f*(f+2);
        double F1=(erf(sqrt(Alpha))/sqrt(Alpha))*sqrt(M_PI)/2.;
        double F2=exp(-Alpha),F3=exp(-Aksq);
        double F1overF2=F1/F2;
        double PLk=PowL.get_yval(k);
        double PLk_nw=(PowL_nw.get_yval(k));
        double PLk_w=PLk-PLk_nw;
        
        double Pbar0=PLk_w+P13dsds_w[0].get_yval(k)+P22dsds_w[0].get_yval(k)+Aksq*PLk_w;
        double Pbar2=2*PLk_w*f+P13dsds_w[1].get_yval(k)+P22dsds_w[1].get_yval(k)+(Alpha+2*f*Aksq)*PLk_w;
        double Pbar4=PLk_w*f*f+P13dsds_w[2].get_yval(k)+P22dsds_w[2].get_yval(k)+(f*f*Aksq+2*f*Alpha)*PLk_w;
        double Pbar6=P13dsds_w[3].get_yval(k)+P22dsds_w[3].get_yval(k)+(f*f*Alpha)*PLk_w+f*f*Alpha*PLk_w;
        double Pbar8=P22dsds_w[4].get_yval(k);
        double Alpha2=Alpha*Alpha,Alpha3=Alpha2*Alpha,Alpha4=Alpha2*Alpha2;
        IR_Leg0[i]=F3*F2*(Pbar0*F1overF2  + Pbar2*(F1overF2-1.)*0.5/Alpha +  Pbar4*(3*(F1overF2-1.)-2*Alpha)/(4*Alpha2) + Pbar6*(15*(F1overF2-1.)+(-10*Alpha-4*Alpha2))/(8*Alpha3) + Pbar8*(105*(F1overF2-1.)+(-70*Alpha-28*Alpha2-8*Alpha3))/(16*Alpha4));
        
        double PnonIR_nw=(PLk_nw+P13dsds_nw[0].get_yval(k)+P22dsds_nw[0].get_yval(k))+(2*PLk*f+P13dsds_nw[1].get_yval(k)+P22dsds_nw[1].get_yval(k))/3. +(PLk_nw*f*f+P13dsds_nw[2].get_yval(k)+P22dsds_nw[2].get_yval(k))/5.+(P13dsds_nw[3].get_yval(k)+P22dsds_nw[3].get_yval(k))/7.+P22dsds_nw[4].get_yval(k)/9.;
        
        IR_Leg0[i]+=PnonIR_nw;
        double PnonIR=(PLk+P13dsds[0].get_yval(k)+P22dsds[0].get_yval(k))+(2*PLk*f+P13dsds[1].get_yval(k)+P22dsds[1].get_yval(k))/3. +(PLk*f*f+P13dsds[2].get_yval(k)+P22dsds[2].get_yval(k))/5.+(P13dsds[3].get_yval(k)+P22dsds[3].get_yval(k))/7.+P22dsds[4].get_yval(k)/9.;
        Non_IR_Leg0[i]=PnonIR;
        
        if(Alpha<4E-3){
            IR_Leg0[i]=Non_IR_Leg0[i];
        }
        //cout<<kvals[i]<<'\t'<<IR_Leg0[i]<<'\t'<<IR_Leg0[i]/Non_IR_Leg0[i]<<'\n';
    }
    IR_leg0=Splining(kvals, IR_Leg0);
    Non_IR_leg0=Splining(kvals, Non_IR_Leg0);
}

void IR_Resum_RSD::calc_IR_Resum_Legendre2(){
    vector<double>IR_Leg2(kvals.size());
    vector<double>Non_IR_Leg2(kvals.size());
    
    cout<<"getting the Legendre2 \n";
    for (int i=0; i<kvals.size(); i++) {
        double k=kvals[i];
        double Aksq=k*k*Aval;
        double Alpha=Aksq*f*(f+2);
        double F1=(erf(sqrt(Alpha))/sqrt(Alpha))*sqrt(M_PI)/2.;
        double F2=exp(-Alpha),F3=exp(-Aksq);
        double F1overF2=F1/F2,F1overF2minus1=F1/F2-1.;
        
        double PLk=PowL.get_yval(k);
        double PLk_nw=(PowL_nw.get_yval(k));
        double PLk_w=PLk-PLk_nw;
        
        double Pbar0=PLk_w+P13dsds_w[0].get_yval(k)+P22dsds_w[0].get_yval(k)+Aksq*PLk_w;
        double Pbar2=2*PLk_w*f+P13dsds_w[1].get_yval(k)+P22dsds_w[1].get_yval(k)+(Alpha+2*f*Aksq)*PLk_w;
        double Pbar4=PLk_w*f*f+P13dsds_w[2].get_yval(k)+P22dsds_w[2].get_yval(k)+(f*f*Aksq+2*f*Alpha)*PLk_w;
        double Pbar6=P13dsds_w[3].get_yval(k)+P22dsds_w[3].get_yval(k)+(f*f*Alpha)*PLk_w+f*f*Alpha*PLk_w;
        double Pbar8=P22dsds_w[4].get_yval(k);
        double Alpha2=Alpha*Alpha,Alpha3=Alpha2*Alpha,Alpha4=Alpha2*Alpha2,Alpha5=Alpha3*Alpha2;
        
        
        IR_Leg2[i]=5.*F3*F2*(Pbar0*(3*(F1overF2minus1)-2*Alpha*F1overF2)/(4*Alpha)  + Pbar2*F2*(9*(F1overF2minus1)-2*Alpha*F1overF2-4*Alpha)/(8*Alpha2) +  Pbar4*(45*(F1overF2minus1)-6*Alpha*F1overF2-24*Alpha-8*Alpha2)/(16*Alpha3) + Pbar6*((F1overF2minus1)*315 -30*Alpha*F1overF2+(-180*Alpha-64*Alpha2- 16*Alpha3))/(32*Alpha4) + Pbar8*((F1overF2minus1)*2835 -210*Alpha*F1overF2 +(-1680*Alpha - 616*Alpha2 - 160*Alpha3 - 32*Alpha4))/(64*Alpha5));
        
        if(Alpha<1E-3)IR_Leg2[i]=2.*Pbar2/3 + 4.*Pbar4/7. + 10.*Pbar6/21. + 40.*Pbar8/99.;
        
        double PnonIR_nw=2.*(2*PLk_nw*f+P13dsds_nw[1].get_yval(k)+P22dsds_nw[1].get_yval(k))/3. +4.*(PLk*f*f+P13dsds_nw[2].get_yval(k)+P22dsds_nw[2].get_yval(k))/7.+10.*(P13dsds_nw[3].get_yval(k)+P22dsds_nw[3].get_yval(k))/21.+40.*P22dsds_nw[4].get_yval(k)/99.;
        IR_Leg2[i]+=PnonIR_nw;
        
        double PnonIR=2.*(2*PLk*f+P13dsds[1].get_yval(k)+P22dsds[1].get_yval(k))/3. +4.*(PLk*f*f+P13dsds[2].get_yval(k)+P22dsds[2].get_yval(k))/7.+10.*(P13dsds[3].get_yval(k)+P22dsds[3].get_yval(k))/21.+40.*P22dsds[4].get_yval(k)/99.;
        Non_IR_Leg2[i]=PnonIR;
        if(Alpha<4E-3){
            IR_Leg2[i]=Non_IR_Leg2[i];
        }
    }
    IR_leg2=Splining(kvals, IR_Leg2);
    Non_IR_leg2=Splining(kvals, Non_IR_Leg2);
}
void IR_Resum_RSD::calc_IR_Resum_Legendre4(){
    vector<double>IR_Leg4(kvals.size());
    vector<double>Non_IR_Leg4(kvals.size());
    
    cout<<"getting the Legendre4 \n";
    for (int i=0; i<kvals.size(); i++) {
        double k=kvals[i];
        double Aksq=k*k*Aval;
        double Alpha=Aksq*f*(f+2);
        double F1=(erf(sqrt(Alpha))/sqrt(Alpha))*sqrt(M_PI)/2.;
        double F2=exp(-Alpha),F3=exp(-Aksq);
        double F1overF2=F1/F2;
        
        double PLk=PowL.get_yval(k);
        double PLk_nw=(PowL_nw.get_yval(k));
        double PLk_w=PLk-PLk_nw;
        
        double Pbar0=PLk_w+P13dsds_w[0].get_yval(k)+P22dsds_w[0].get_yval(k)+Aksq*PLk_w;
        double Pbar2=2*PLk_w*f+P13dsds_w[1].get_yval(k)+P22dsds_w[1].get_yval(k)+(Alpha+2*f*Aksq)*PLk_w;
        double Pbar4=PLk_w*f*f+P13dsds_w[2].get_yval(k)+P22dsds_w[2].get_yval(k)+(f*f*Aksq+2*f*Alpha)*PLk_w;
        double Pbar6=P13dsds_w[3].get_yval(k)+P22dsds_w[3].get_yval(k)+(f*f*Alpha)*PLk_w+f*f*Alpha*PLk_w;
        double Pbar8=P22dsds_w[4].get_yval(k);
        double Alpha2=Alpha*Alpha,Alpha3=Alpha2*Alpha,Alpha4=Alpha2*Alpha2,Alpha5=Alpha3*Alpha2,Alpha6=Alpha3*Alpha3;
        
        IR_Leg4[i]=9*F3*F2*(((12 *Alpha2 *F1overF2 + 105 *(-1. + F1overF2) -10 *Alpha*(1. + 6 *F1overF2)) /(32* Alpha2))*Pbar0 );
        IR_Leg4[i]+=9*F3*F2*( (525 *(-1.+ F1overF2) + 4 *Alpha2 *(-8 + 3 *F1overF2) -10 *Alpha *(17 + 18 *F1overF2))/(64*Alpha3)*Pbar2);
        IR_Leg4[i]+=-9*F3*F2*((64* Alpha3 + Alpha2 *(416-36* F1overF2) -3675 *(-1.+F1overF2) + 50*Alpha*(31 + 18* F1overF2))/(128*Alpha4)*Pbar4);
        IR_Leg4[i]+=-9*F3*F2*(((960* Alpha3 + 128* Alpha4) - 33075* (-1. + F1overF2) +3150 *Alpha *(5 + 2 *F1overF2) -60 *Alpha2 *(-80 + 3 *F1overF2))/(256*Alpha5)*Pbar6);
        IR_Leg4[i]+=-9*F3*F2*(1/(512*Alpha6)*((13440* Alpha3 + 2176 *Alpha4 + 256 *Alpha5) - 1260*Alpha2 *(-48 + F1overF2) -363825* (-1. + F1overF2) +3150 *Alpha* (59 + 18* F1overF2))*Pbar8 );
        
        if(Alpha<4E-3)IR_Leg4[i]=8.*Pbar4/35. + 24.*Pbar6/77. + 48.*Pbar8/143.;
        
        double PnonIR_nw=8.*(PLk_nw*f*f+P13dsds_nw[2].get_yval(k)+P22dsds_nw[2].get_yval(k))/35.+24.*(P13dsds_nw[3].get_yval(k)+P22dsds_nw[3].get_yval(k))/77.+48.*P22dsds_nw[4].get_yval(k)/143.;
        IR_Leg4[i]+=PnonIR_nw;
        
        double PnonIR=8.*(PLk*f*f+P13dsds[2].get_yval(k)+P22dsds[2].get_yval(k))/35.+24.*(P13dsds[3].get_yval(k)+P22dsds[3].get_yval(k))/77.+48.*P22dsds[4].get_yval(k)/143.;
        Non_IR_Leg4[i]=PnonIR;
        
        //cout<<kvals[i]<<'\t'<<IR_Leg4[i]<<'\t'<<IR_Leg4[i]/Non_IR_Leg4[i]<<'\n';
    }
    IR_leg4=Splining(kvals, IR_Leg4);
    Non_IR_leg4=Splining(kvals, Non_IR_Leg4);
}
Splining IR_Resum_RSD::get_RSD_1loop_IR(int legMode){
    Splining splineout;
    switch (legMode) {
        case 0:
            return IR_leg0;
            break;
        case 2:
            return IR_leg2;
            break;
        case 4:
            return IR_leg4;
            break;
        default:
            cout<<"Trying to get a Legendre mode that doesn't exist for the 1 loop IR RSD\n";
            break;
    }
    return splineout;
}
Splining IR_Resum_RSD::get_RSD_1loop_nonIR(int legMode){
    Splining splineout;
    switch (legMode) {
        case 0:
            return Non_IR_leg0;
            break;
        case 2:
            return Non_IR_leg2;
            break;
        case 4:
            return Non_IR_leg4;
            break;
        default:
            cout<<"Trying to get a Legendre mode that doesn't exist for the 1 loop non-IR RSD\n";
            break;
    }
    return splineout;
}
