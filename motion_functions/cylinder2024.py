import math
import cmath


rcyl        = 0.5
Atheta      = math.pi
Aa          = 1.5
Ag          = 0.15
alpha_max   = 0.3
omega_alpha = 6.
tref        = 1.


def c_atan2(y,x):
    """
    Approach detailed in:
    https://mdolab.engin.umich.edu/misc/complex-step-guide-python

    Note, convention for x,y order in atan2 utilized in the above reference as of 
    7Feb2023 is atan2(x,y), which is not consistent with general practice. Function 
    updated here to follow standard atan2(y,x) interface.
    """
    a=y.real
    b=y.imag
    c=x.real
    d=x.imag
    return complex(math.atan2(a,c),(c*b-a*d)/(a**2+c**2))


def alpha1(t):
    return t*t*t*(8.-3.*t)/16.

def alpha2(t):
    return alpha_max*cmath.sin(omega_alpha*t)*(t**6./(t**6.+tref**6.))

def psi(t,alpha):
    return 1. + (Aa - 1.)*alpha(t)

def eta(lam,omega,tau):
    return cmath.sin(omega*lam + tau*(1. - cmath.cos(omega*lam)))

def fg(t,r0,theta0):
    return (16.*(r0**4.) + (t**6./(t**6. + 0.01))*eta(t,10.,0.7)*(cmath.cos(32.*cmath.pi*r0**4.) - 1.))*eta(theta0,1.,0.7)


def thetag(t,r0,theta0):
    return theta0 + Ag*fg(t,r0,theta0)


def coords(alpha,x0,y0,t):

    r0     = cmath.sqrt(x0*x0 + y0*y0)
    theta0 = c_atan2(y0,x0)

    xdef0 = r0*cmath.cos(thetag(t,r0,theta0))
    ydef0 = r0*cmath.sin(thetag(t,r0,theta0))

    x = cmath.cos(Atheta*alpha(t))*psi(t,alpha)*xdef0 - (cmath.sin(Atheta*alpha(t))/psi(t,alpha))*ydef0
    y = cmath.sin(Atheta*alpha(t))*psi(t,alpha)*xdef0 + (cmath.cos(Atheta*alpha(t))/psi(t,alpha))*ydef0 + alpha(t)

    # cmath call is always converting things into complex, so we turn things back
    # into floats if inputs are floats.
    if isinstance(x0,complex) or isinstance(y0,complex) or isinstance(t,complex):
        x_return = x
        y_return = y
    else:
        x_return = x.real
        y_return = y.real

    return x_return, y_return


def coords_dot(alpha,x0,y0,t):

    h = 1.e-30
    eps = complex(0.,h)

    # Evaluate function with complex perturbations in x,y,t
    dxdx,dydx = coords(alpha,x0+eps,y0    ,t    )
    dxdy,dydy = coords(alpha,x0    ,y0+eps,t    )
    dxdt,dydt = coords(alpha,x0    ,y0    ,t+eps)

    # Divide imaginary part by step size
    dxdx = dxdx.imag/h
    dxdy = dxdy.imag/h
    dxdt = dxdt.imag/h
    dydx = dydx.imag/h
    dydy = dydy.imag/h
    dydt = dydt.imag/h

    return dxdx, dxdy, dxdt, dydx, dydy, dydt


def test_derivatives():

    x0 = 1.
    y0 = 0.3
    t = 1.
    fd_eps = 1.e-9

    # Compute complex-step derivatives
    dxdx_cs, dxdy_cs, dxdt_cs, dydx_cs, dydy_cs, dydt_cs = coords_dot(alpha1,x0,y0,t)

    # Compute finite-difference derivatives
    xa,ya         = coords(alpha1,x0,       y0,       t       )
    x_xdot,y_xdot = coords(alpha1,x0+fd_eps,y0       ,t       )
    x_ydot,y_ydot = coords(alpha1,x0       ,y0+fd_eps,t       )
    x_tdot,y_tdot = coords(alpha1,x0       ,y0       ,t+fd_eps)

    dxdx_fd = (x_xdot - xa)/fd_eps
    dydx_fd = (y_xdot - ya)/fd_eps

    dxdy_fd = (x_ydot - xa)/fd_eps
    dydy_fd = (y_ydot - ya)/fd_eps

    dxdt_fd = (x_tdot - xa)/fd_eps
    dydt_fd = (y_tdot - ya)/fd_eps


    print(f"Cylinder2024: Case1 short-time composite motion with deformation extension (x0,y0,t) = ({x0},{y0},{t})")
    print("{: >30} {: >30} {: >30} {: >30}".format(*[' ', f'Finite-difference(step={fd_eps})', 'Complex-step', 'Error']))
    print("{: >30} {: >30} {: >30} {: >30}".format(*['dxdx',dxdx_fd, dxdx_cs, dxdx_fd - dxdx_cs]))
    print("{: >30} {: >30} {: >30} {: >30}".format(*['dxdy',dxdy_fd, dxdy_cs, dxdy_fd - dxdy_cs]))
    print("{: >30} {: >30} {: >30} {: >30}".format(*['dxdt',dxdt_fd, dxdt_cs, dxdt_fd - dxdt_cs]))
    print("{: >30} {: >30} {: >30} {: >30}".format(*['dydx',dydx_fd, dydx_cs, dydx_fd - dydx_cs]))
    print("{: >30} {: >30} {: >30} {: >30}".format(*['dydy',dydy_fd, dydy_cs, dydy_fd - dydy_cs]))
    print("{: >30} {: >30} {: >30} {: >30}".format(*['dydt',dydt_fd, dydt_cs, dydt_fd - dydt_cs]))



    # Compute complex-step derivatives
    dxdx_cs, dxdy_cs, dxdt_cs, dydx_cs, dydy_cs, dydt_cs = coords_dot(alpha2,x0,y0,t)

    # Compute finite-difference derivatives
    xa,ya         = coords(alpha2,x0,       y0,       t       )
    x_xdot,y_xdot = coords(alpha2,x0+fd_eps,y0       ,t       )
    x_ydot,y_ydot = coords(alpha2,x0       ,y0+fd_eps,t       )
    x_tdot,y_tdot = coords(alpha2,x0       ,y0       ,t+fd_eps)

    dxdx_fd = (x_xdot - xa)/fd_eps
    dydx_fd = (y_xdot - ya)/fd_eps

    dxdy_fd = (x_ydot - xa)/fd_eps
    dydy_fd = (y_ydot - ya)/fd_eps

    dxdt_fd = (x_tdot - xa)/fd_eps
    dydt_fd = (y_tdot - ya)/fd_eps

    print(" ")
    print(f"Cylinder2024: Case2 long-time study of conservation. (x0,y0,t) = ({x0},{y0},{t})")
    print("{: >30} {: >30} {: >30} {: >30}".format(*[' ', f'Finite-difference(step={fd_eps})', 'Complex-step', 'Error']))
    print("{: >30} {: >30} {: >30} {: >30}".format(*['dxdx',dxdx_fd, dxdx_cs, dxdx_fd - dxdx_cs]))
    print("{: >30} {: >30} {: >30} {: >30}".format(*['dxdy',dxdy_fd, dxdy_cs, dxdy_fd - dxdy_cs]))
    print("{: >30} {: >30} {: >30} {: >30}".format(*['dxdt',dxdt_fd, dxdt_cs, dxdt_fd - dxdt_cs]))
    print("{: >30} {: >30} {: >30} {: >30}".format(*['dydx',dydx_fd, dydx_cs, dydx_fd - dydx_cs]))
    print("{: >30} {: >30} {: >30} {: >30}".format(*['dydy',dydy_fd, dydy_cs, dydy_fd - dydy_cs]))
    print("{: >30} {: >30} {: >30} {: >30}".format(*['dydt',dydt_fd, dydt_cs, dydt_fd - dydt_cs]))

if __name__ == "__main__":
    test_derivatives()






