//
//  fitting_routine.hpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 24/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#ifndef fitting_routine_hpp
#define fitting_routine_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <gsl/gsl_multifit.h>
#include <stdbool.h>
#include <math.h>
#include "PowSpecRSD_SPT_resummation.hpp"
//https://www.gnu.org/software/gsl/manual/html_node/Multi_002dparameter-fitting.html
//http://rosettacode.org/wiki/Polynomial_regression#C
bool polynomialfit(int obs, int deg_start,int degree,std::vector<double>dx, std::vector<double>dy, std::vector<double>&store);

double fit_cssquared(int numpts,double kmin,double kmax,import_class imports,IR_Resum irresum,bool printvals=true);
double error_estim(int numpts,double mean,std::vector<double>yvals);

void fit_cssquared_2(int numpts,double kmin,double kmax,import_class imports,IR_Resum irresum);

#endif /* fitting_routine_hpp */
