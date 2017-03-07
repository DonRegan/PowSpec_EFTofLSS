//
//  Fabrikant.cpp
//  PowSpec_EFTofLSS
//
//  Created by Donough Regan on 22/02/2017.
//  Copyright © 2017 Donough Regan. All rights reserved.
//

#include "Fabrikant.hpp"
//Fabrikant are Fabrikant_mnp=int x^2 jm(kx)jn(qx)jp(rx)
double Fabrik000(double k,double q,double r){
    return M_PI/(4.*k*q*r);
}
double Fabrik011(double k,double q,double r){
    return (M_PI*(-k*k + q*q + r*r))/(8.*k*q*q*r*r);
}
double Fabrik112(double k,double q,double r){
    double ksq=k*k,rsq=r*r,qsq=q*q;
    return (M_PI*(-3.*(ksq - qsq)*(ksq - qsq) + 2.*(ksq + qsq)* rsq + rsq*rsq))/(32.*ksq*qsq*rsq*r);
}
double Fabrik022(double k,double q,double r){
    double ksq=k*k,rsq=r*r,qsq=q*q;
    return (M_PI*(3.*(ksq - qsq)*(ksq - qsq) + 2.*(-3.*ksq + qsq)*rsq + 3.*rsq*rsq))/(32.*k*qsq*q*rsq*r);
}
double Fabrik222(double k,double q,double r){
    double k2=k*k,r2=r*r,q2=q*q,k4=k2*k2,q4=q2*q2,r4=r2*r2;
    return (M_PI*(-3.*(k2 - q2)*(k2 - q2)*(k2 + q2) + (3.*k4 + 2.*k2*q2 +3.*q4)* r2 + 3.*(k2 + q2)* r4 - 3.*r4*r2))/(64.*k*q*r*k2*q2*r2);
}
double Fabrik123(double k,double q,double r){
    double k2=k*k,r2=r*r,q2=q*q,k4=k2*k2,q4=q2*q2,r4=r2*r2;
    
    return (M_PI*(5.*pow((k2 - q2),3) +3.*(-3.*k4 + 2.*k2*q2 + q4)* r2 + (3.*k2 + q2)*r4 +r2*r4))/(64*k2*q2*q*r4);
}
double Fabrik242(double k,double q,double r){
    double k2=k*k,r2=r*r,q2=q*q,k3=k2*k,q3=q2*q,r3=r2*r,k4=k2*k2,q4=q2*q2,r4=r2*r2,q5=q3*q2,k6=k4*k2,q6=q4*q2,r6=r4*r2,k8=k4*k4;
    
    return (1./(512.*k3*q5*r3))*M_PI*(35.*k8 - 20.*k6*(3.*q2 + 7.*r2) +6.*k4 *(3. *q4 + 10.* q2* r2 + 35.* r4) + (q2 - r2)*(q2 - r2)*(3.*q4 +10.*q2*r2 + 35.*r4) + 4.*k2 *(q6 + 3.*q4*r2 + 15.*q2*r4 - 35.*r6));
}
double Fabrik033(double k,double q,double r){
    double k2=k*k,r2=r*r,q2=q*q,k4=k2*k2,q4=q2*q2,r4=r2*r2,r6=r4*r2;
    return (M_PI*(-5.*pow((k2 - q2),3) + 3.* (5.*k4 - 6.* k2* q2 + q4)* r2 +3.*(-5.*k2 + q2)* r4 + 5.* r6))/(64.* k *q4 *r4);
}
double Fabrik044(double k,double q,double r){
    double k2=k*k,r2=r*r,q2=q*q,q3=q2*q,r3=r2*r,k4=k2*k2,q4=q2*q2,r4=r2*r2,q5=q3*q2,r5=r3*r2,r6=r4*r2,r8=r4*r4;
    
    return (M_PI*(35.* pow((k2 - q2),4) - 20.*(k2 - q2)*(k2 - q2)*(7.*k2 - q2)*r2 +6.* (35.* k4 - 30.* k2 *q2 + 3.* q4)* r4 + 20.* (-7.* k2 + q2)* r6 +35.* r8))/(512.* k *q5 *r5);
}
