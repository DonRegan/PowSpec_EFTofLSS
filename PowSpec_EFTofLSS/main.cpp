//
//  main.cpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 22/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#include <iostream>
#include "output_routines.hpp"
#include "fitting_routine.hpp"

using namespace std;
int main(int argc, const char * argv[]) {
    

    //Read in the variables
    import_class imports("/Users/donoughregan/Desktop/Projects_2017/EFTofLSS_Code4Paper/EFTofLSS_C++/PowSpec_EFTofLSS/paramFiles.ini");
    initData initdata= imports.get_initialdata();
    finalData finaldata=imports.get_finaldata();
    
    /**** Create the cosmology and calculate the necessary coefficients ****/
    Cosmology_FLRW cosFLRW(initdata.cosmol);
    Coeffs coeff(initdata.zinit,finaldata.zfinal,5,cosFLRW);
    solvesystem(coeff);
    
    /****** Calculation + Resummation of the Dark Matter Power Spectrum ******/
    int numevals=200;
    IR_Resum irresum(imports,coeff,numevals,logScale);
    
    /****** Calculation + Resummation of the RSD Power Spectra ******/
    IR_Resum_RSD irresum_rsd(irresum);
    irresum_rsd.output_PowSpec_contribs(100, imports);
    
    /***** Outputs here ******/
    output_DM_powspec(1000,imports,irresum);
    output_RSD_powspec(1000,imports,irresum_rsd);
    output_timedep_cts(imports,coeff);
    
    /***** On the fly post-processing here ******/
    fit_cssquared(100, 0.2, 0.4, imports, irresum);//cs^2 calculation for the DM power spectrum
    fit_cssquared_2(100, 0.2, 0.4, imports, irresum);
    return 0;
}
