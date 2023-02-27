#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <cylinder2024.h>


// Floating-point numbers
const double PI = 3.14159265358979323846;

// Cylinder parameters
const double rcyl        = 0.5;
const double Atheta      = PI;
const double Aa          = 1.5;
const double Ag          = 0.15;
const double alpha_max   = 0.3;
const double omega_alpha = 6.0;
const double tref        = 1.0;




double complex c_atan2(complex double y, complex double x){
    double a, b, c, d;
    double complex theta;

    a = creal(y);
    b = cimag(y);
    c = creal(x);
    d = cimag(x);

    if ((pow(a,2.0) + pow(c,2.0)) < 1.0e-13){
    //if ((pow(a,2.0) + pow(c,2.0)) == 0.0){
        theta = atan2(a,c) + 0.0*I;
    } else {
        theta = atan2(a,c) + (c*b-a*d)/(cpow(a,2.0) + cpow(c,2.0))*I;
    }
    return theta;
}


double complex alpha1(double complex t){
    double complex res;
    res = t*t*t*(8.0-3.0*t)/16.0;
    return res;
}



double complex alpha2(double complex t){
    double complex res;
    res = alpha_max*csin(omega_alpha*t)*(cpow(t,6.0)/(cpow(t,6.0)+cpow(tref,6.0)));
    return res;
}

double complex alpha_generic(int alpha, double complex t){
    double complex res;

    if (alpha == ALPHA_SHORT){
        res = alpha1(t);
    } else if (alpha == ALPHA_LONG){
        res = alpha2(t);
    } else{
        printf("Invalid selection of alpha time-activation function.");
        exit(1);
    }
    return res;
}


double complex psi(double complex t, int alpha){
    double complex res;
    res = 1.0 + (Aa - 1.0)*alpha_generic(alpha,t);
    return res;
}


double complex eta(double complex lam, double omega, double tau){
    double complex res;
    res = csin(omega*lam + tau*(1.0 - ccos(omega*lam)));
    return res; 
}


double complex fg(double complex t, double complex r0, double complex theta0){
    double complex res;
    res = (16.0*cpow(r0,4.0) + (cpow(t,6.0)/(cpow(t,6.0) + 0.01)) * \
          eta(t,10.0,0.7)*(ccos(32.0*PI*cpow(r0,4.0)) - 1.0))*eta(theta0,1.0,0.7);
    return res; 
}



double complex thetag(double complex t, double complex r0, double complex theta0){
    double complex res;
    res = theta0 + Ag*fg(t,r0,theta0);
    return res;
}


    
void coords(int alpha,double complex xref,double complex yref, double complex t, \
            double complex *x,double complex *y){
    double complex r0, theta0, xdef0, ydef0;

    r0     = csqrt(xref*xref + yref*yref);
    theta0 = c_atan2(yref,xref);

    xdef0 = r0*ccos(thetag(t,r0,theta0));
    ydef0 = r0*csin(thetag(t,r0,theta0));

    *x = ccos(Atheta*alpha_generic(alpha,t))*psi(t,alpha)*xdef0 - \
        (csin(Atheta*alpha_generic(alpha,t))/psi(t,alpha))*ydef0;
    *y = csin(Atheta*alpha_generic(alpha,t))*psi(t,alpha)*xdef0 + \
        (ccos(Atheta*alpha_generic(alpha,t))/psi(t,alpha))*ydef0 + alpha_generic(alpha,t);

}




void coords_dot(int alpha,double complex xref,double complex yref,double complex t, \
                double *dxdx,double *dxdy,double *dxdt,double *dydx,double *dydy,double *dydt){
    double complex  eps, my_psi, my_alpha, \
                    dxdx_c, dxdy_c, dxdt_c,\
                    dydx_c, dydy_c, dydt_c;

    double h; 

    h   = 1.e-30;
    eps = 0.0 + h*I;

    // Issues related to the use of atan2 in coords function make the derivative 
    // at 0,0 undefined. The form of the motion simplifies a bit at (0,0) and we 
    // hard-code the spatial derivatives here
    if ( cabs(creal(xref)) < 1.0e-13 && \
         cabs(creal(yref)) < 1.0e-13 ){
        my_psi = psi(t,alpha);
        my_alpha = alpha_generic(alpha,t);
        *dxdx = creal(my_psi)*cos(Atheta*creal(my_alpha));
        *dxdy = -(1.0/creal(my_psi))*sin(Atheta*creal(my_alpha));
        *dydx = creal(my_psi)*sin(Atheta*creal(my_alpha));
        *dydy = (1.0/creal(my_psi))*cos(Atheta*creal(my_alpha));

    // Away from (0,0) we evaluate derivatives via complex-step
    } 
    else {
        // Evaluate function with complex perturbations in x,y
        coords(alpha,xref+eps,yref    ,t, &dxdx_c, &dydx_c);
        coords(alpha,xref    ,yref+eps,t, &dxdy_c, &dydy_c);

        // Divide imaginary part by step size
        *dxdx = cimag(dxdx_c)/h;
        *dxdy = cimag(dxdy_c)/h;
        *dydx = cimag(dydx_c)/h;
        *dydy = cimag(dydy_c)/h;
    }

    // There are no complicating issues for the time derivative, evaluate everywhere with complex-step
    coords(alpha,xref,    yref,    t+eps,&dxdt_c, &dydt_c);
    *dxdt = cimag(dxdt_c)/h;
    *dydt = cimag(dydt_c)/h;

}

