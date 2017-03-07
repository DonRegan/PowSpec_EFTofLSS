//
//  import_variables.cpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 22/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#include "import_variables.hpp"

import_class::import_class(std::string inifile){
    reader=INIReader(inifile);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load "<<inifile<<'\n';
    }
    defaultString="UNKNOWN";
    defaultInt=-1;
    defaultDouble=-1;
    //naming the groups
    initialdata="initial_data";
    finaldata="final_data";
    
    set_initial_data();
    set_final_data();
}
//Now let's get the initial data
void import_class::set_initial_data(){
    initdata.zinit=reader.GetReal(initialdata, "zinit", defaultDouble);
    initdata.initpow_file=reader.Get(initialdata, "powspecInit", defaultString);
    set_cosmology_data();
    initdata.cosmol=cosmol;
}
void import_class::set_cosmology_data(){
    
    cosmol.deltaphi= reader.GetReal(initialdata, "deltaphi", defaultDouble);
    cosmol.h= reader.GetReal(initialdata, "h", defaultDouble);
    cosmol.OmegaM= reader.GetReal(initialdata, "OmegaM", defaultDouble);
    cosmol.ns= reader.GetReal(initialdata, "ns", defaultDouble);
    cosmol.OmegaB= reader.GetReal(initialdata, "OmegaB", defaultDouble);
    cosmol.TCMB= reader.GetReal(initialdata, "TCMB", defaultDouble);
}
//Now let's get the final data
void import_class::set_final_data(){
    final_data.zfinal=reader.GetReal(finaldata, "zfinal", defaultDouble);
    final_data.linearpow_DM=reader.Get(finaldata, "pklinear", defaultString);
    final_data.nonlinearpow_DM=reader.Get(finaldata, "pknonlinear", defaultString);
    final_data.calcP22=reader.GetBoolean(finaldata, "calcP22", true);
    final_data.baseoutputs=reader.Get(finaldata, "baseoutputs", defaultString);
    final_data.kminCalc=reader.GetReal(finaldata, "kminCalc", defaultDouble);
    final_data.kmaxCalc=reader.GetReal(finaldata, "kmaxCalc", defaultDouble);
}
