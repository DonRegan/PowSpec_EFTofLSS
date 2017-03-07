//
//  Wiggle_NoWiggle.hpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 22/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//
/************
 The goal here is to find the no-wiggle linear power spectrum.
 As indicated in the appendix of 1509.02120 I divide the linear power spectrum by the BBKS approximation.
 Then use a Gaussian filter (in log space) on [Plin/PBBKS] to get the wiggle and no-wiggle parts (after multiplying back by PBBKS.
 ***********/

#ifndef Wiggle_NoWiggle_hpp
#define Wiggle_NoWiggle_hpp

#include "cosmology.hpp"

Splining get_Power_nowiggle(double z,Splining &PowSpec,Cosmology_FLRW cosFLRW,double lambda=0.25);//eq A.4 of 1509.02120
Splining get_Power_wiggle(Splining PowSpec,Splining PowerNoWiggle);
double gaussianfiltered_power(double log10q,void *p);
double gaussianfiltered_power_denomin(double log10q,void *p);
struct params_wigglepower{Splining PowL;double lambda;double log10k;};


#endif /* Wiggle_NoWiggle_hpp */
