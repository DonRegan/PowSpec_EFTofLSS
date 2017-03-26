//
//  output_routines.hpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 24/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#ifndef output_routines_hpp
#define output_routines_hpp

#include <stdio.h>
#include "PowSpecRSD_SPT_resummation.hpp"

void output_DM_powspec(int numpts,import_class &imports,IR_Resum &irresum);
void output_DM_values(std::string filename,std::vector<double>kvals,Splining Pk1loop,Splining PkL,Splining PkNL);

void output_RSD_powspec(int numpts,import_class &imports,IR_Resum_RSD &irresum_rsd);
void output_RSD_values(std::string filename,std::vector<double>kvals,Splining PowL,Splining PowL_nw,std::vector<Splining> PowRSD);

void output_timedep_cts(import_class &imports,Coeffs &coeff);
#endif /* output_routines_hpp */
