#include <cylinder2024.hpp>
#include <cmath>
#include <cstdio>

const complex<double> I(0.0,1.0);

// Cylinder parameters
const double rcyl        = 0.5;
const double Atheta      = M_PI;
const double Aa          = 1.5;
const double Ag          = 0.15;
const double alpha_max   = 0.3;
const double omega_alpha = 6.0;
const double tref        = 1.0;

complex<double> c_atan2(complex<double> y, complex<double> x) {
    double a, b, c, d;
    complex<double> theta;

    a = real(y);
    b = imag(y);
    c = real(x);
    d = imag(x);

    if ((pow(a,2.0) + pow(c,2.0)) < 1.0e-13) {
    //if ((pow(a,2.0) + pow(c,2.0)) == 0.0) {
        theta = atan2(a,c) + 0.0*I;
    } else {
        theta = atan2(a,c) + (c*b-a*d)/(pow(a,2.0) + pow(c,2.0))*I;
    }
    return theta;
}


complex<double> alpha1(complex<double> t) {
    complex<double> res;
    res = t*t*t*(8.0-3.0*t)/16.0;
    return res;
}



complex<double> alpha2(complex<double> t) {
    complex<double> res;
    res = alpha_max*sin(omega_alpha*t)*(pow(t,6.0)/(pow(t,6.0)+pow(tref,6.0)));
    return res;
}

complex<double> alpha_generic(AlphaFunction alpha, complex<double> t) {
    complex<double> res;
    switch (alpha) {
    case ALPHA_SHORT:
        res = alpha1(t);
        break;
    case ALPHA_LONG:
        res = alpha2(t);
        break;
    }
    return res;
}


complex<double> psi(complex<double> t, AlphaFunction alpha) {
    complex<double> res;
    res = 1.0 + (Aa - 1.0)*alpha_generic(alpha,t);
    return res;
}


complex<double> eta(complex<double> lam, double omega, double tau) {
    complex<double> res;
    res = sin(omega*lam + tau*(1.0 - cos(omega*lam)));
    return res; 
}


complex<double> fg(complex<double> t, complex<double> r0, complex<double> theta0) {
    complex<double> res;
    res = (16.0*pow(r0,4.0) + (pow(t,6.0)/(pow(t,6.0) + 0.01)) *
          eta(t,10.0,0.7)*(cos(32.0*M_PI*pow(r0,4.0)) - 1.0))*eta(theta0,1.0,0.7);
    return res; 
}



complex<double> thetag(complex<double> t, complex<double> r0, complex<double> theta0) {
    complex<double> res;
    res = theta0 + Ag*fg(t,r0,theta0);
    return res;
}

void coords(AlphaFunction alpha,complex<double> xref,complex<double> yref, complex<double> t,
            complex<double> &x,complex<double> &y) {
    complex<double> r0, theta0, xdef0, ydef0;

    r0     = sqrt(xref*xref + yref*yref);
    theta0 = c_atan2(yref,xref);

    xdef0 = r0*cos(thetag(t,r0,theta0));
    ydef0 = r0*sin(thetag(t,r0,theta0));

    x = cos(Atheta*alpha_generic(alpha,t))*psi(t,alpha)*xdef0 -
        (sin(Atheta*alpha_generic(alpha,t))/psi(t,alpha))*ydef0;
    y = sin(Atheta*alpha_generic(alpha,t))*psi(t,alpha)*xdef0 +
        (cos(Atheta*alpha_generic(alpha,t))/psi(t,alpha))*ydef0 + alpha_generic(alpha,t);

}

// Real inputs / outputs
void coords(AlphaFunction alpha, double xref, double yref, double t, double &x, double &y) {
    complex<double> xref_c(xref), yref_c(yref), t_c(t), x_c, y_c;
    coords(alpha, xref_c, yref_c, t_c, x_c, y_c);
    x = real(x_c);
    y = real(y_c);
}
    

void coords_dot(AlphaFunction alpha,complex<double> xref,complex<double> yref,complex<double> t,
                double &dxdx,double &dxdy,double &dxdt,double &dydx,double &dydy,double &dydt) {
    complex<double> eps, my_psi, my_alpha,
                    dxdx_c, dxdy_c, dxdt_c,
                    dydx_c, dydy_c, dydt_c;

    double h; 

    h   = 1.e-30;
    eps = 0.0 + h*I;

    // Issues related to the use of atan2 in coords function make the derivative 
    // at 0,0 undefined. The form of the motion simplifies a bit at (0,0) and we 
    // hard-code the spatial derivatives here
    if ( abs(real(xref)) < 1.0e-13 &&
         abs(real(yref)) < 1.0e-13 ) {
        my_psi = psi(t,alpha);
        my_alpha = alpha_generic(alpha,t);
        dxdx = real(my_psi)*cos(Atheta*real(my_alpha));
        dxdy = -(1.0/real(my_psi))*sin(Atheta*real(my_alpha));
        dydx = real(my_psi)*sin(Atheta*real(my_alpha));
        dydy = (1.0/real(my_psi))*cos(Atheta*real(my_alpha));

    // Away from (0,0) we evaluate derivatives via complex-step
    } 
    else {
        // Evaluate function with complex perturbations in x,y
        coords(alpha,xref+eps,yref    ,t, dxdx_c, dydx_c);
        coords(alpha,xref    ,yref+eps,t, dxdy_c, dydy_c);

        // Divide imaginary part by step size
        dxdx = imag(dxdx_c)/h;
        dxdy = imag(dxdy_c)/h;
        dydx = imag(dydx_c)/h;
        dydy = imag(dydy_c)/h;
    }

    // There are no complicating issues for the time derivative, evaluate everywhere with complex-step
    coords(alpha,xref,    yref,    t+eps, dxdt_c, dydt_c);
    dxdt = imag(dxdt_c)/h;
    dydt = imag(dydt_c)/h;
}

// Real inputs / outputs
void coords_dot(AlphaFunction alpha, double xref, double yref, double t,
                double &dxdx, double &dxdy, double &dxdt,
                double &dydx, double &dydy, double &dydt) {
    complex<double> xref_c(xref), yref_c(yref), t_c(t);
    coords_dot(alpha, xref_c, yref_c, t_c, dxdx, dxdy, dxdt, dydx, dydy, dydt);
}
