//
//  fitting_routine.cpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 24/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#include "fitting_routine.hpp"
using namespace std;

bool polynomialfit(int obs, int deg_start,int degree,std::vector<double>dx, std::vector<double>dy, std::vector<double>&store) /* n, p */
{ //fitting sum_i store[i]*dx[j]^i =y[j] for j labelling the observations.
    //i runs from deg_start to degree-1 (inclusive)
    gsl_multifit_linear_workspace *ws;
    gsl_matrix *cov, *X;
    gsl_vector *y, *c;
    double chisq;
    
    int i, j;
    
    X = gsl_matrix_alloc(obs, degree);
    y = gsl_vector_alloc(obs);
    c = gsl_vector_alloc(degree);
    
    cov = gsl_matrix_alloc(degree, degree);
    
    for(i=0; i < obs; i++) {
        if(deg_start==0){
            gsl_matrix_set(X, i, 0, 1.0); //change here ...was 1.0
        }
        for(j=deg_start; j < degree; j++) {//change omit j=0 case
            gsl_matrix_set(X, i, j, pow(dx[i], j));
        }
        gsl_vector_set(y, i, dy[i]);
    }
    
    ws = gsl_multifit_linear_alloc(obs, degree);
    gsl_multifit_linear(X, y, c, cov, &chisq, ws);
    
    /* store result ... */
    for(i=0; i < degree; i++)
    {
        //store[i] = gsl_vector_get(c, i);
        store.push_back(gsl_vector_get(c, i));
    }
    
    gsl_multifit_linear_free(ws);
    gsl_matrix_free(X);
    gsl_matrix_free(cov);
    gsl_vector_free(y);
    gsl_vector_free(c);
    return true; /* we do not "analyse" the result (cov matrix mainly)
                  to know if the fit is "good" */
}

double fit_cssquared(int numpts,double kmin,double kmax,import_class imports,IR_Resum irresum,bool printvals){
    Splining Pow1loop_IR=irresum.get_IR_DM();
    Splining Pow1loop_nonIR=irresum.get_nonIR_DM();
    Splining PowL_IR=irresum.get_PowL_IR();
    Splining PowL_nonIR=irresum.get_PowL();
    Splining PowNL=PowSpec(imports.get_finaldata().nonlinearpow_DM).get_PowSpec();
    vector<double>P_tofit_IR(numpts),P_tofit_nonIR(numpts),kvals;
    create_vec(kmin, kmax, numpts, kvals, linearScale);
    for (int i=0; i<numpts; i++) {
        P_tofit_IR[i]=PowNL.get_yval(kvals[i])-Pow1loop_IR.get_yval(kvals[i]) ;
        P_tofit_nonIR[i]=PowNL.get_yval(kvals[i])-Pow1loop_nonIR.get_yval(kvals[i]) ;

        P_tofit_IR[i]/=-2.*kvals[i]*kvals[i]*PowL_IR.get_yval(kvals[i]);
        P_tofit_nonIR[i]/=-2.*kvals[i]*kvals[i]*PowL_IR.get_yval(kvals[i]);
    }
    int deg(1),degstart(0);
    std::vector<double>alpha_IRfit,alpha_nonIRfit;
    polynomialfit(numpts, degstart,deg, kvals, P_tofit_IR, alpha_IRfit);
    polynomialfit(numpts, degstart,deg, kvals, P_tofit_nonIR, alpha_nonIRfit);
    if(printvals){
        cout<<"Finding the fit and the error using a polynomial fitting routine (error is deviation of curve to the fit)"<<'\n';
        double error_csquared_IR=error_estim(numpts,alpha_IRfit[0],P_tofit_IR);
        cout<<"csquared (IR) "<<alpha_IRfit[0]<<" ± "<<error_csquared_IR<<endl;
        double error_csquared_nonIR=error_estim(numpts,alpha_nonIRfit[0],P_tofit_nonIR);
        cout<<"csquared (nonIR) "<<alpha_nonIRfit[0]<<" ± "<<error_csquared_nonIR<<endl;
    }
    return alpha_IRfit[0];
}
void fit_cssquared_2(int numpts,double kmin,double kmax,import_class imports,IR_Resum irresum){
    
    cout<<"Alternatively we get the error by putting a 1% error on PNL and see how cs^2 changes across many realisations\n";
    Splining PowNL=PowSpec(imports.get_finaldata().nonlinearpow_DM).get_PowSpec();
    vector<double>kvals;
    create_vec(kmin, kmax, numpts, kvals, linearScale);
    
    Splining Pow1loop_IR=irresum.get_IR_DM(),PowL_IR=irresum.get_PowL_IR();
    
    double cs_sq=fit_cssquared(numpts, kmin, kmax, imports, irresum,false);
    std::default_random_engine gen;
    vector<std::normal_distribution<double>>distribution(numpts);
    for (int i=0; i<numpts; i++) {
        distribution[i]=std::normal_distribution<double>(PowNL.get_yval(kvals[i]),PowNL.get_yval(kvals[i])*0.01);
    }
    int num_realisations=100;
    std::vector<double>PkNL_realisation(numpts),P_tofit(numpts);
    
    double errval=0.0;
    for(int n=0;n<num_realisations;n++){
        int deg(1),degstart(0);
        std::vector<double>alpha;
        for (int i=0; i<numpts; i++) {
            PkNL_realisation[i]=distribution[i](gen);
            P_tofit[i]=-(PkNL_realisation[i]-Pow1loop_IR.get_yval(kvals[i]))/(2.*kvals[i]*kvals[i]*PowL_IR.get_yval(kvals[i]));
        }
        polynomialfit(numpts, degstart,deg, kvals, P_tofit, alpha);
        double cs_sq_realisation=alpha[0];
        errval+=(cs_sq-cs_sq_realisation)*(cs_sq-cs_sq_realisation);
    }
    errval=sqrt(errval/(num_realisations-1));
    cout<<"csq ± error "<<cs_sq<<'\t'<<errval<<'\n';
}

double error_estim(int numpts,double mean,std::vector<double>yvals){
    double sum=0.0;
    for (int i=0; i<numpts; i++) {
        sum+=(yvals[i]-mean)*(yvals[i]-mean);
    }
    return sqrt(sum/(numpts-1));
}


