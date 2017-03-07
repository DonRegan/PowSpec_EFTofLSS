//
//  integration_routine.hpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 22/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#ifndef integration_routine_hpp
#define integration_routine_hpp

#include <stdio.h>
#include <gsl/gsl_integration.h>

double integrate(double lowerlim,double upperlim,double (*func)(double,void *),void *p,double tol=1E-4);

#endif /* integration_routine_hpp */
