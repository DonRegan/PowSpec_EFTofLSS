//
//  import_powerspectra.hpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 22/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#ifndef import_powerspectra_hpp
#define import_powerspectra_hpp

#include <stdio.h>
#include "splining_routines.hpp"


class PowSpec{
private:
    std::vector<double> kvals;
    std::vector<double> Pkvals;
    Splining powerspec;
    int sizePk;
    std::string filename;
public:
    PowSpec(std::string filenameIn):filename(filenameIn){
        read_camb_pk();
        powerspec=Splining(kvals, Pkvals);
    }
    PowSpec();
    ~PowSpec(){}
    void read_camb_pk();
    std::vector<double> get_kvals(){return kvals;}
    std::vector<double> get_Pkvals(){return Pkvals;}
    int get_sizePk(){return sizePk;}
    double get_kmin(){ return kvals[0];}
    double get_kmax(){ return kvals[sizePk-1];}
    Splining get_PowSpec(){return powerspec;}
};

#endif /* import_powerspectra_hpp */
