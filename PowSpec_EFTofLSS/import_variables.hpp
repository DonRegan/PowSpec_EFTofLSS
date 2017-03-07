//
//  import_variables.hpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 22/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#ifndef import_variables_hpp
#define import_variables_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include "INIReader.hpp"
struct cosmolInfo{
    double deltaphi;
    double h;
    double OmegaM;
    double ns;
    double OmegaB;
    double TCMB;
    cosmolInfo(){//default values from Table 4. Col 2. Planck 2015 parameters
        deltaphi=1.49228E-8;h=0.6704;OmegaM=0.3168;ns=0.9619;OmegaB=0.049;TCMB=2.7255;
    }
};
struct initData{
    double zinit;
    std::string initpow_file;
    cosmolInfo cosmol;
};
struct finalData{
    double zfinal;
    std::string linearpow_DM;
    std::string nonlinearpow_DM;
    bool calcP22; //this is true if we need to calculate P22 for the DM or RSDs
    std::string baseoutputs;
    double kminCalc;
    double kmaxCalc;
};

class import_class
{
public:
    import_class(std::string inifile);
public: //initial data
    void set_initial_data();
    void set_cosmology_data();
    cosmolInfo get_cosmolvals(){return cosmol;}
    initData get_initialdata(){return initdata;}
private:
    cosmolInfo cosmol;
    initData initdata;
public:
    void set_final_data();
    finalData get_finaldata(){return final_data;}
private:
    finalData final_data;
private:
    INIReader reader;
    std::string initialdata,finaldata;
    std::string defaultString;
    int defaultInt;
    double defaultDouble;
};
#endif /* import_variables_hpp */
