//
//  output_routines.cpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 24/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#include "output_routines.hpp"
using namespace std;

void output_RSD_13_22_mu2n(int numpts,import_class &imports,double f,Splining &PowL,std::vector<Splining>&P13dsds, std::vector<Splining>&P22dsds){
    string filename=imports.get_finaldata().baseoutputs+"Pk_RSD_contribs_z"+to_string(imports.get_finaldata().zfinal)+".csv";
    
    double kmin=imports.get_finaldata().kminCalc,kmax=imports.get_finaldata().kmaxCalc;
    vector<double>kvals;
    create_vec(kmin, kmax, numpts, kvals, logScale);
    
    std::ofstream myfile(filename);
    if (!myfile.is_open()) {
        std::cout<<"Problem opening file "<<filename<<std::endl;
        return;
    }
    myfile<<"kvalue"<<","<<"Ptree_mu0"<<","<<"P13_mu0"<<","<<"P22_mu0"<<","<<"Ptree_mu2"<<","<<"P13_mu2"<<","<<"P22_mu2"<<","<<"Ptree_mu4"<<","<<"P13_mu4"<<","<<"P22_mu4"<<'\n';
    for (int i=0; i<kvals.size(); i++) {
        myfile<<kvals[i]<<","<<PowL.get_yval(kvals[i])<<","<<P13dsds[0].get_yval(kvals[i])<<","<<P22dsds[0].get_yval(kvals[i])<<",";
        myfile<<2*f*PowL.get_yval(kvals[i])<<","<<P13dsds[1].get_yval(kvals[i])<<","<<P22dsds[1].get_yval(kvals[i])<<",";
        myfile<<f*f*PowL.get_yval(kvals[i])<<","<<P13dsds[2].get_yval(kvals[i])<<","<<P22dsds[2].get_yval(kvals[i])<<std::scientific<<"\n";
        
    }
    myfile.close();
}

void output_DM_powspec(int numpts, import_class &imports,IR_Resum &irresum){
    string filename_IR=imports.get_finaldata().baseoutputs+"Pk_DM_z"+to_string(imports.get_finaldata().zfinal)+".csv";
    string filename_nonIR=imports.get_finaldata().baseoutputs+"Pk_DM_IR_z"+to_string(imports.get_finaldata().zfinal)+".csv";

    Splining PkDM_nonIR=irresum.get_nonIR_DM();
    Splining PkDM_IR=irresum.get_IR_DM();
    Splining PkL_nonIR=irresum.get_PowL();
    Splining PkL_IR=irresum.get_PowL_IR();
    
    //Also read in the NL values and put them to a spline
    PowSpec PowNLclass(imports.get_finaldata().nonlinearpow_DM);
    Splining PkNL=PowNLclass.get_PowSpec();
    
    double kmin=imports.get_finaldata().kminCalc,kmax=imports.get_finaldata().kmaxCalc;
    vector<double>kvals;
    create_vec(kmin, kmax, numpts, kvals, logScale);
    //Now let's output the nonIR resummed values
    output_DM_values(filename_nonIR,kvals,PkDM_nonIR,PkL_nonIR,PkNL);
    output_DM_values(filename_IR,kvals,PkDM_IR,PkL_IR,PkNL);

}
void output_DM_values(string filename,vector<double>kvals,Splining Pk1loop,Splining PkL,Splining PkNL){
    std::ofstream myfile(filename);
    if (!myfile.is_open()) {
        std::cout<<"Problem opening Pk DM output file "<<filename<<std::endl;
        return;
    }
    myfile<<"kvalue"<<","<<"PkNL"<<","<<"Pk1loop"<<","<<"PkL"<<'\n';
    for (int i=0; i<kvals.size(); i++) {
        myfile<<kvals[i]<<","<<PkNL.get_yval(kvals[i])<<","<<Pk1loop.get_yval(kvals[i])<<","<<PkL.get_yval(kvals[i])<<std::scientific<<"\n";
    }
    myfile.close();
}

void output_RSD_powspec(int numpts,import_class &imports,IR_Resum_RSD &irresum_rsd){
    string filename_IR=imports.get_finaldata().baseoutputs+"Pk_RSD_z"+to_string(imports.get_finaldata().zfinal)+".csv";
    string filename_nonIR=imports.get_finaldata().baseoutputs+"Pk_RSD_IR_z"+to_string(imports.get_finaldata().zfinal)+".csv";

    double kmin=imports.get_finaldata().kminCalc,kmax=imports.get_finaldata().kmaxCalc;
    vector<double>kvals;
    create_vec(kmin, kmax, numpts, kvals, logScale);
    Splining PowL_IR=irresum_rsd.get_PowL_IR();
    Splining PowL_nw=irresum_rsd.get_PowL_nw();
    Splining PowL_nonIR=irresum_rsd.get_PowL();
    vector<Splining>PowRSD_nonIR(3),PowRSD_IR(3);
    for (int i=0; i<3; i++) {
        PowRSD_nonIR[i]=irresum_rsd.get_RSD_1loop_nonIR(2*i);
        PowRSD_IR[i]=irresum_rsd.get_RSD_1loop_IR(2*i);
    }
    output_RSD_values(filename_IR, kvals, PowL_IR,PowL_nw, PowRSD_IR);
    output_RSD_values(filename_nonIR, kvals, PowL_nonIR,PowL_nw, PowRSD_nonIR);

}
void output_RSD_values(string filename,vector<double>kvals,Splining PowL,Splining PowL_nw,vector<Splining> PowRSD){
    std::ofstream myfile(filename);
    if (!myfile.is_open()) {
        std::cout<<"Problem opening file "<<filename<<std::endl;
        return;
    }
    myfile<<"kvalue"<<","<<"PkL"<<","<<"PkL_nw"<<","<<"PkIR_ell_0"<<","<<"PkIR_ell_2"<<","<<"PkIR_ell_4"<<'\n';
    for (int i=0; i<kvals.size(); i++) {
        myfile<<kvals[i]<<","<<PowL.get_yval(kvals[i])<<","<<PowL_nw.get_yval(kvals[i])<<","<<PowRSD[0].get_yval(kvals[i])<<","<<PowRSD[1].get_yval(kvals[i])<<","<<PowRSD[2].get_yval(kvals[i])<<std::scientific<<"\n";
    }
    myfile.close();
}



void output_timedep_cts(import_class &imports,Coeffs &coeffs){
    double z=imports.get_finaldata().zfinal;
    double g,A,B,D,E,F,G,J,f,Abarprime,Bbarprime,Dbarprime,Ebarprime,Fbarprime,Gbarprime,Jbarprime;
    vector<double>yfinal(16);
    for (int i=0; i<16; i++) {
        yfinal[i]=coeffs.get_yfinal(i);
    }
    g=yfinal[0],A=yfinal[2],B=yfinal[4],D=yfinal[6],E=yfinal[8],F=yfinal[10],G=yfinal[12],J=yfinal[14];
    f=-(1+z)*yfinal[1]/g,Abarprime=-(1+z)*yfinal[3],Bbarprime=-(1+z)*yfinal[5];
    //see RSD_Checks.nb for the derivation of these time dependencies
    double T1,T1bar,g2=g*g,g3=g2*g,f2=f*f,f3=f2*f;
    T1=-18*D-28*E+7*F+2*G+13*J;
    T1bar=-18*Dbarprime-28*Ebarprime+7*Fbarprime+2*Gbarprime+13*Jbarprime;
    double ct_mu0=-2*T1/g;
    double ct_mu2=-2/g *(f*T1 + T1bar  - 8 *f*g*(Abarprime+Bbarprime) +5.5*f2*g3);
    double ct_mu4=-2/g *(      f*T1bar - 8 *f2*g*(Abarprime+Bbarprime) -f*g*(Abarprime+6*Bbarprime) +f2*g3*(1+3*f) );
    double ct_mu6=-2/g *(                                              -f2*g*(Abarprime+6*Bbarprime) +f3*g3*(1-2.5*f));

    string filename=imports.get_finaldata().baseoutputs +"cts_z"+to_string(z)+".csv";
    std::ofstream myfile(filename);
    if (!myfile.is_open()) {
        std::cout<<"Problem opening time dependence of the counterterms file "<<filename<<std::endl;
        return;
    }
    myfile<<"ct_mu0"<<","<<"ct_mu2"<<","<<"ct_mu4"<<","<<"ct_mu6"<<","<<"g"<<","<<"f"<<'\n';
    myfile<<ct_mu0<<","<<ct_mu2<<","<<ct_mu4<<","<<ct_mu6<<","<<g<<","<<f<<'\n';
    myfile.close();
}
