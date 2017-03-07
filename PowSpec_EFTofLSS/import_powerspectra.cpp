//
//  import_powerspectra.cpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 22/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#include "import_powerspectra.hpp"

using namespace std;

void PowSpec::read_camb_pk(){
    ifstream infile(filename);
    if(!infile.is_open()){
        cout<<"Problem opening the CAMB file \n";
        return;
    }
    double k1,Pk1;
    sizePk=0;
    std::string name1,name2;
    getline(infile,name1);
    getline(infile,name2);
    
    while(infile>>k1>>Pk1){
        kvals.push_back(k1);
        Pkvals.push_back(Pk1);
        sizePk++;
    }
}
