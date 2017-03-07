//
//  integration_routine.cpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 22/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#include "integration_routine.hpp"
double integrate(double lowerlim,double upperlim,double (*func)(double,void *),void *p,double tol){
    gsl_function F={func,p};
    gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc (100);
    double result,error;
    size_t nevals;
    gsl_integration_cquad (&F, lowerlim, upperlim, 0, tol,w, &result, &error,&nevals);
    gsl_integration_cquad_workspace_free (w);
    return result;
}
