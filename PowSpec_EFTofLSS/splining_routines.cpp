//
//  splining_routines.cpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 22/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#include "splining_routines.hpp"
void create_vec(double minval,double maxval,int size,std::vector<double>&xvals,which_scale scale){
    if (scale==linearScale) {
        double step=(maxval-minval)/(double)(size-1);
        for (int i=0; i<size; i++) {
            double val=minval+step*i;
            xvals.push_back(val);
        }
    }else{
        double step=pow(maxval/minval,1./(size-1));
        for (int i=0; i<size; i++) {
            xvals.push_back(minval*pow(step,i));
        }
    }
    return;
}

void Splining::make_spline(){
    
    spline=(gsl_spline *)gsl_spline_alloc(gsl_interp_cspline,size_vec);
    acc=(gsl_interp_accel *)gsl_interp_accel_alloc();
    gsl_spline_init(spline,&xvals[0],&yvals[0],size_vec);
}


Spline_To_Vec::Spline_To_Vec(Splining &splineIn,double xminIn,double xmaxIn,int numpointsIn,which_scale scaleIn):xmin(xminIn),xmax(xmaxIn),spline(splineIn),numpoints(numpointsIn),scale(scaleIn){
    set_xvals();
    set_yvals();
}
void Spline_To_Vec::set_yvals(){
    for (int i=0; i<numpoints; i++) {
        yvals.push_back(spline.get_yval(xvals[i]));
    }
}
void Spline_To_Vec::set_xvals(){
    xvals.push_back(xmin);
    if(scale==linearScale){
        double step=(xmax-xmin)/(double)(numpoints-1);
        for (int i=1; i<numpoints; i++) {
            xvals.push_back(xmin+i*step);
        }
    }else{//log spacing
        double step=pow(xmax/xmin,1.0/(double)(numpoints-1));
        for (int i=1; i<numpoints; i++) {
            xvals.push_back(xmin*pow(step,i));
        }
    }
    
    return;
}
