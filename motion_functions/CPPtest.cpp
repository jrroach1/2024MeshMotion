#include <cstdio>
#include <cylinder2024.hpp>

const double fd_eps = 1.e-10;
const double fd_tol = 1.e-5;
AlphaFunction alpha = ALPHA_SHORT;

int main() {
    for (double x0 = -0.5; x0<0.5; x0=x0+0.1 ) {
        for (double y0 = -0.5; y0<0.5; y0=y0+0.1 ) {
            for (double t  =  0.0; t <1.0; t =t +0.1 ) {
                double
                    x, y,
                    x_xdot, x_ydot, x_tdot,
                    y_xdot, y_ydot, y_tdot,
                    dxdx_fd, dydx_fd, dxdy_fd,
                    dydy_fd, dxdt_fd, dydt_fd,
                    dxdx_cs, dydx_cs, dxdy_cs,
                    dydy_cs, dxdt_cs, dydt_cs;

                // Compute complex-step derivatives
                coords_dot(alpha, x0, y0, t,
                           dxdx_cs, dxdy_cs, dxdt_cs,
                           dydx_cs, dydy_cs, dydt_cs);

                // Compute finite-difference derivatives
                coords(alpha, x0,        y0,        t,        x,      y     );
                coords(alpha, x0+fd_eps, y0,        t,        x_xdot, y_xdot);
                coords(alpha, x0,        y0+fd_eps, t,        x_ydot, y_ydot);
                coords(alpha, x0,        y0,        t+fd_eps, x_tdot, y_tdot);

                dxdx_fd = (x_xdot - x)/fd_eps;
                dydx_fd = (y_xdot - y)/fd_eps;
                
                dxdy_fd = (x_ydot - x)/fd_eps;
                dydy_fd = (y_ydot - y)/fd_eps;
                
                dxdt_fd = (x_tdot - x)/fd_eps;
                dydt_fd = (y_tdot - y)/fd_eps;

                if (abs(dxdx_fd - dxdx_cs) > fd_tol ||
                    abs(dxdy_fd - dxdy_cs) > fd_tol ||
                    abs(dxdt_fd - dxdt_cs) > fd_tol ||
                    abs(dydx_fd - dydx_cs) > fd_tol ||
                    abs(dydy_fd - dydy_cs) > fd_tol ||
                    abs(dydt_fd - dydt_cs) > fd_tol ) {
                    printf("Derivative error detected...");
                    printf("x,y,t: %.15f %.15f %.15f\n", x0, y0, t);
                    printf("dxdx %.15f %.15f %.15f\n", dxdx_fd, dxdx_cs, dxdx_fd - dxdx_cs);
                    printf("dxdy %.15f %.15f %.15f\n", dxdy_fd, dxdy_cs, dxdy_fd - dxdy_cs);
                    printf("dxdt %.15f %.15f %.15f\n", dxdt_fd, dxdt_cs, dxdt_fd - dxdt_cs);
                    printf("dydx %.15f %.15f %.15f\n", dydx_fd, dydx_cs, dydx_fd - dydx_cs);
                    printf("dydy %.15f %.15f %.15f\n", dydy_fd, dydy_cs, dydy_fd - dydy_cs);
                    printf("dydt %.15f %.15f %.15f\n", dydt_fd, dydt_cs, dydt_fd - dydt_cs);
                }
            }
        }
    }
    return 0;
}
