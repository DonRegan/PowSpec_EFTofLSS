//
//  splining_routines.hpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 22/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#ifndef splining_routines_hpp
#define splining_routines_hpp

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <assert.h>
#include <fstream>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>


//using namespace std;
enum which_scale{logScale,linearScale};
void create_vec(double minval,double maxval,int size,std::vector<double>&xvals,which_scale scale);

class Splining{
private:
    std::vector<double> xvals;
    std::vector<double> yvals;
    int size_vec;
    gsl_spline *spline;
    gsl_interp_accel *acc;
    double xmin,xmax;
public:
    Splining(){
        
    }
    Splining(std::vector<double> xvalsIn,std::vector<double> yvalsIn):xvals(xvalsIn),yvals(yvalsIn){
        size_vec=static_cast<int>(xvals.size());
        xmin=xvals[0];
        xmax=xvals[size_vec-1];
        make_spline();
    }
    ~Splining(){
    }
    void make_spline();
    double get_yval(double xval){
        // assert(xval<=xvals[size_vec-1]);
        if(xval>xvals[size_vec-1])return 0.0;//TEMPORARY CHANGE DMR
        
        if(xval<xvals[0])return 0.0;//TEMPORARY CHANGE DMR
        return gsl_spline_eval(spline, xval, acc);
    }
    double get_xmin(){ return xmin;}
    double get_xmax(){ return xmax;}
    void free_gsl_structures(){
        gsl_spline_free (spline);
        gsl_interp_accel_free (acc);
    }
};

class Spline_To_Vec{
private:
    Splining spline;
    which_scale scale;
    std::vector<double>xvals,yvals,yvalsfac;
    double xmin,xmax;
    int numpoints;
public:
    Spline_To_Vec(Splining &splineIn,double xminIn,double xmaxIn,int numpoints,which_scale scaleIn=linearScale);
    ~Spline_To_Vec(){}
    void set_xvals();
    void set_yvals();
    std::vector<double>get_xvals(){return xvals;}
    std::vector<double>get_yvals(){return yvals;}
    std::vector<double>get_yvals_fac(double fac){
        for (int i=0; i<numpoints; i++) {
            yvalsfac.push_back(fac*yvals[i]);
        }
        return yvalsfac;
    }
    
};

#endif /* splining_routines_hpp */
