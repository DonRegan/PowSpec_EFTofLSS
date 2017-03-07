//
//  PowSpecRSD_SPT_resummation.hpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 23/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#ifndef PowSpecRSD_SPT_resummation_hpp
#define PowSpecRSD_SPT_resummation_hpp

#include <stdio.h>
#include "PowDM_SPT_Resummation.hpp"
#include "PowSpecRSD_SPT.hpp"

class IR_Resum_RSD{
private:
    double Aval,f,g;
    std::vector<double>kvals;
    std::vector<Splining>P13dsds,P22dsds,P13dsds_nw,P22dsds_nw,P13dsds_w,P22dsds_w;
    Splining PowL,PowL_nw,PowL_IR;
    Splining IR_leg0,IR_leg2,IR_leg4;
    Splining Non_IR_leg0,Non_IR_leg2,Non_IR_leg4;
public:
    IR_Resum_RSD(IR_Resum &irresum);
public:
    void calc_wiggle_contribs();
    void calc_IR_Resum_Legendre0();
    void calc_IR_Resum_Legendre2();
    void calc_IR_Resum_Legendre4();
public://access routines
    Splining get_RSD_1loop_IR(int legMode);
    Splining get_RSD_1loop_nonIR(int legMode);
    Splining get_PowL(){return PowL;}
    Splining get_PowL_IR(){return PowL_IR;}
};
#endif /* PowSpecRSD_SPT_resummation_hpp */
