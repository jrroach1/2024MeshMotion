#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>

#include <cylinder2024.h>

int main(){

    double fd_eps, \
           x0, y0, t, \
           dxdx_fd, dydx_fd, dxdy_fd,  \
           dydy_fd, dxdt_fd, dydt_fd,  \
           dxdx_cs, dydx_cs, dxdy_cs,  \
           dydy_cs, dxdt_cs, dydt_cs;

    double complex x0_c, y0_c, t_c, xa, ya,       \
                   x_xdot, x_ydot, x_tdot,  \
                   y_xdot, y_ydot, y_tdot;

    for (x0 = -0.5; x0<0.5; x0=x0+0.1 ){
    for (y0 = -0.5; y0<0.5; y0=y0+0.1 ){
    for (t  =  0.0; t <1.0; t =t +0.1 ){

    x0_c   = x0     + 0.0*I;
    y0_c   = y0     + 0.0*I;
    t_c    = t      + 0.0*I;
    fd_eps = 1.e-10 + 0.0*I;

    // Compute complex-step derivatives
    coords_dot(ALPHA_SHORT,x0,y0,t,&dxdx_cs,&dxdy_cs,&dxdt_cs,&dydx_cs,&dydy_cs,&dydt_cs);

    // Compute finite-difference derivatives
    coords(ALPHA_SHORT,x0,       y0,       t       , &xa,     &ya    );
    coords(ALPHA_SHORT,x0+fd_eps,y0       ,t       , &x_xdot, &y_xdot);
    coords(ALPHA_SHORT,x0       ,y0+fd_eps,t       , &x_ydot, &y_ydot);
    coords(ALPHA_SHORT,x0       ,y0       ,t+fd_eps, &x_tdot, &y_tdot);


    dxdx_fd = (creal(x_xdot) - creal(xa))/fd_eps;
    dydx_fd = (creal(y_xdot) - creal(ya))/fd_eps;

    dxdy_fd = (creal(x_ydot) - creal(xa))/fd_eps;
    dydy_fd = (creal(y_ydot) - creal(ya))/fd_eps;

    dxdt_fd = (creal(x_tdot) - creal(xa))/fd_eps;
    dydt_fd = (creal(y_tdot) - creal(ya))/fd_eps;

    if ( abs(dxdx_fd - dxdx_cs) > 0.0001 || \
         abs(dxdy_fd - dxdy_cs) > 0.0001 || \
         abs(dxdt_fd - dxdt_cs) > 0.0001 || \
         abs(dydx_fd - dydx_cs) > 0.0001 || \
         abs(dydy_fd - dydy_cs) > 0.0001 || \
         abs(dydt_fd - dydt_cs) ){
        printf("Derivative error detected...");
        printf("x,y,t: %.15f %.15f %.15f\n", x0, y0, t);
        printf("dxdx %.15f %.15f %.15f\n", dxdx_fd, dxdx_cs, dxdx_fd - dxdx_cs);
        printf("dxdy %.15f %.15f %.15f\n", dxdy_fd, dxdy_cs, dxdy_fd - dxdy_cs);
        printf("dxdt %.15f %.15f %.15f\n", dxdt_fd, dxdt_cs, dxdt_fd - dxdt_cs);
        printf("dydx %.15f %.15f %.15f\n", dydx_fd, dydx_cs, dydx_fd - dydx_cs);
        printf("dydy %.15f %.15f %.15f\n", dydy_fd, dydy_cs, dydy_fd - dydy_cs);
        printf("dydt %.15f %.15f %.15f\n", dydt_fd, dydt_cs, dydt_fd - dydt_cs);}

    }}}

}


