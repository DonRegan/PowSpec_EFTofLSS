//
//  CoeffsCalculator.cpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 22/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#include "CoeffsCalculator.hpp"
using namespace std;

double zdepCosmol_FLRW::get_Omegam_z(double z){
    double fac(1+z),fac2(fac*fac),fac3(fac2*fac),fac4(fac2*fac2);
    double denomin=cosmol.OmegaM*fac3+cosmol.Omegar*fac4+cosmol.OmegaL;
    return ((cosmol.OmegaM*fac3) /denomin);
}
double zdepCosmol_FLRW::get_Omegam_z_prime(double z){
    double fac(1+z),fac2(fac*fac),fac3(fac2*fac),fac4(fac2*fac2);
    double denomin=cosmol.OmegaM*fac3+cosmol.Omegar*fac4+cosmol.OmegaL;
    double Omegam_z=get_Omegam_z(z);
    return (3.0*fac2*Omegam_z/denomin - Omegam_z*(3.0*fac2*cosmol.OmegaM+4.0*fac3*cosmol.Omegar)/(denomin));
}
double zdepCosmol_FLRW::get_epsilon(double z){
    double fac(1+z),fac2(fac*fac),fac3(fac2*fac),fac4(fac2*fac2);
    double denomin=cosmol.OmegaM*fac3+cosmol.Omegar*fac4+cosmol.OmegaL;
    return (1.5*get_Omegam_z(z) +2.0*cosmol.Omegar*fac4/denomin);
}
double  zdepCosmol_FLRW::get_epsilon_prime(double z){
    double fac(1+z),fac2(fac*fac),fac3(fac2*fac),fac4(fac2*fac2);
    double denomin=cosmol.OmegaM*fac3+cosmol.Omegar*fac4+cosmol.OmegaL;
    
    return (1.5*get_Omegam_z_prime(z) + 8.0*cosmol.Omegar*fac3/denomin - (2.0*cosmol.Omegar*fac4/denomin)*(3.0*fac2*cosmol.OmegaM+4.0*fac3*cosmol.Omegar)/(denomin));
}
double zdepCosmol_FLRW::get_c0(double zval){
    return ((1-get_epsilon(zval))/(1+zval));
}
double zdepCosmol_FLRW::get_c1(double zval){
    return (1.5*get_Omegam_z(zval)/((1+zval)*(1+zval)));
}
double zdepCosmol_FLRW::get_c0prime(double z){
    return (-(1-get_epsilon(z))/((1+z)*(1+z)) - get_epsilon_prime(z)/(1+z));
}
double zdepCosmol_FLRW::get_c1prime(double z){
    double fac(1+z),fac2(fac*fac),fac3(fac2*fac);
    return (-3.0*get_Omegam_z(z)/fac3 + 1.5*get_Omegam_z_prime(z)/fac2);
}
/******Coeffs class *********/
/****ODES to solve*******
 g''-c0 g' -c1 g =0
 A''-c0 A' -c1 A = c1 g^2
 ...
 G''-c0 G' -c1 G = c1 g B
 H''-c0 H' -c1 H = (g')^2 g
 y[0]=g,y[1]=g'
 y[2]=A,y[3]=A'
 ...
 y[12]=G,y[13]=G'
 y[14]=J,y[15]=J' ... new line
 HERE c0 = (1-epsilon)/(1+z), c1 = (3/2) Omegam(z)/(1+z)^2
 *************************/
struct params2{Coeffs coeffs;};
Coeffs::Coeffs(double zinitIn,double zfinalIn,int num_zsteps,Cosmology_FLRW &cosmolIn,which_scale scaleIn):zinit(zinitIn),zfinal(zfinalIn),numz(num_zsteps),scale(scaleIn),cosmol(cosmolIn),cosmoz(zdepCosmol_FLRW(cosmolIn)){
    ode_size=16;
    initialise_y();
    set_zvals();
}
void Coeffs::set_zvals(){
    zvals.push_back(zinit);
    if(scale==which_scale::linearScale){
        double step=(zfinal-zinit)/(double)(numz);
        for (int i=1; i<=numz; i++) {
            zvals.push_back(zinit+i*step);
        }
    }else{//log spacing
        double step=pow(zfinal/zinit,1.0/(double)(numz));
        for (int i=1; i<=numz; i++) {
            zvals.push_back(zinit*pow(step,i));
        }
    }
    
    return;
}



void Coeffs::initialise_y(){
    //For matching to input spectrum we need g(z*)=1
    //gLCDM[z] = Exp[Integrate[OmegaM[zp]^(gamma-1)/(1+zp),{zp,z,1100}]
    
    yinit.push_back(1.0);//corresponding to g(z*)=1
    double glcdmPrime=-pow(cosmoz.get_Omegam_z(zinit),0.55-1)/(1+zinit);
    //yinit.push_back(-1.0/(1.0+zinit));//g'(z)=a(z)/a(init)=-1/(1+zinit)
    yinit.push_back(glcdmPrime);
    
    for (int i=2; i<ode_size; i++) {
        yinit.push_back(0.0);
    }
    yinit[2]=-1.;//CHANGE DMR
    yinit[4]=-(2./3.)*glcdmPrime*glcdmPrime*(1+zinit)*(1+zinit)/cosmoz.get_Omegam_z(zinit);
    yfinal.resize(ode_size);
}
void Coeffs::set_yfinal(std::vector<double>ysolution){
    yfinal=ysolution;
}
void Coeffs::set_yfinalEdS(){
    //just need g and f to get these...
    double g=yfinal[0];
    yfinalEdS=vector<double>(yfinal.size());
    double g2=g*g,g3=g2*g;
    double g2deriv=2.*g*yfinal[1];
    double g3deriv=3.*g2*yfinal[1];
    
    yfinalEdS[0]=yfinal[0];yfinalEdS[1]=yfinal[1];
    //A,B,C,D,F,G,J and derivs
    yfinalEdS[2]=3.*g2/7.;yfinalEdS[3]=3.*g2deriv/7.; //A and Aprime
    yfinalEdS[4]=2.*g2/7.;yfinalEdS[5]=2.*g2deriv/7.; //B and Bprime
    yfinalEdS[6]=2.*g3/21.;yfinalEdS[7]=2.*g3deriv/21.; //D and Dprime
    yfinalEdS[8]=4.*g3/63.;yfinalEdS[9]=4.*g3deriv/63.; //E and Eprime
    yfinalEdS[10]=1.*g3/14.;yfinalEdS[11]=1.*g3deriv/14.; //F and Fprime
    yfinalEdS[12]=1.*g3/21.;yfinalEdS[13]=1.*g3deriv/21.; //G and Gprime
    yfinalEdS[14]=1.*g3/9.;yfinalEdS[15]=1.*g3deriv/9.; //J and Jprime
}

int odesystem_dydz(double z,const double y[],double f[],void *p){//the Coeffs stucture is passed in via params
    struct params2 *pars=(struct params2 *)p;
    Coeffs coeffs=pars->coeffs;
    
    double c0(coeffs.get_cosmoz().get_c0(z)),c1(coeffs.get_cosmoz().get_c1(z));
    
    for (int i=0; i<coeffs.get_odesize()/2; i++) {
        f[i*2]=y[i*2+1];//e.g. f[0]=y[1]\equiv dy[0]/dz
        f[i*2+1]=c0*y[i*2+1]+c1*y[i*2];//e.g. f[1]=c0*y[1]+c1*y[0]
    }
    //extra bits
    f[3]+=c1*y[0]*y[0];
    f[5]+=y[1]*y[1];
    f[7]+=y[1]*y[3];
    f[9]+=y[1]*y[5];
    f[11]+=c1*y[0]*y[2];
    f[13]+=c1*y[0]*y[4];
    f[15]+=y[1]*y[1]*y[0];//possible issue with J...
    
    return GSL_SUCCESS;
}
int odesystem_jacobian(double z,const double y[],double *dfdy,double dfdz[],void *p){
    struct params2 *pars=(struct params2 *)p;
    Coeffs coeffs=pars->coeffs;
    double c0(coeffs.get_cosmoz().get_c0(z)),c1(coeffs.get_cosmoz().get_c1(z));
    //double c0p(coeffs.get_cosmoz().get_c0prime()),c1p(coeffs.get_cosmoz().get_c1prime());
    int dim=coeffs.get_odesize();
    
    gsl_matrix_view dfdy_mat=gsl_matrix_view_array(dfdy,dim,dim);
    gsl_matrix *jacmat=&dfdy_mat.matrix;
    
    for (int i=0; i<dim/2; i++) {
        gsl_matrix_set(jacmat,i*2,i*2+1,1.0);//e.g. df0/dy1 = dy1/dy1 = dg'/dg' =1
        gsl_matrix_set(jacmat,i*2+1,i*2+1,c0);//df[i*2+1]/dy[i*2+1]
        gsl_matrix_set(jacmat,i*2+1,i*2,c1);
        
    }
    //extrabits
    gsl_matrix_set(jacmat,3,0, 2*c1*y[0] );
    gsl_matrix_set(jacmat,5,1, 2*y[1] );//change from 2*y[0]
    gsl_matrix_set(jacmat,7,1, y[3] );     gsl_matrix_set(jacmat,7,3, y[1] );
    gsl_matrix_set(jacmat,9,1, y[5] );     gsl_matrix_set(jacmat,9,5, y[1] );
    gsl_matrix_set(jacmat,11,0, c1*y[2] ); gsl_matrix_set(jacmat,11,2, c1*y[0] );
    gsl_matrix_set(jacmat,13,0, c1*y[4] ); gsl_matrix_set(jacmat,13,4, c1*y[0] );
    gsl_matrix_set(jacmat,15,0, y[1]*y[1] ); gsl_matrix_set(jacmat,15,1, 2*y[0]*y[1] );
    
    
    for (int i=0; i<dim; i++) {
        dfdz[i]=0.0;
    }
    
    return 0;
}

void solvesystem(Coeffs &coeffs_info){
    //I'll put the updated yvalues into the Coeffs structure
    int odesize=coeffs_info.get_odesize();
    int numtimesteps=coeffs_info.get_numz();
    struct params2 myparams={coeffs_info};
    
    gsl_odeiv2_system sys_a = {odesystem_dydz, odesystem_jacobian,(size_t) odesize,&myparams};
    gsl_odeiv2_driver * d =gsl_odeiv2_driver_alloc_y_new (&sys_a, gsl_odeiv2_step_rk8pd,-1e-10, 1e-10, 0.0);
    double *ysolution=new double[odesize];
    std::vector<double>ysol;
    std::copy(coeffs_info.get_yinit().begin(),coeffs_info.get_yinit().end(),ysolution);//copy vector into array
    for(int i=0;i<odesize;i++){
       // cout<<"ysolution "<<i<<'\t'<<ysolution[i]<<endl;
    }
    double zval=coeffs_info.get_zval(0);
    
    for (int i=1; i<=numtimesteps; i++) { //don't actually need all the z values. Just got them for checking. I'll save the values at the final z only
        double zi=coeffs_info.get_zval(i);
        int status=gsl_odeiv2_driver_apply(d,&zval,zi,ysolution );//using ICs ysolution at zval previous to get solution at zi//z will run up to tstep_i from zval at end of previous loop
        assert(status==GSL_SUCCESS);
    }
    for (int i=0; i<odesize; i++) {
        cout<<"ysol "<<i<<'\t'<<ysolution[i]<<'\n';
        ysol.push_back(ysolution[i]);
    }
    coeffs_info.set_yfinal(ysol);
    coeffs_info.set_yfinalEdS();
    gsl_odeiv2_driver_free(d);
    free(ysolution);
    return;
    
}



