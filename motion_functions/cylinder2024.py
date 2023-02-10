import math
import cmath

""" Guide for using mesh motion functions

Usage:
--------------------
    import cylinder2024
    x,y,z                         = cylinder2024.coords(    cylinder2024.alpha1,xref,yref,t)
    dxdx,dxdy,dxdt,dydx,dydy,dydt = cylinder2024.coords_dot(cylinder2024.alpha1,xref,yref,t)


Top-level functions:
--------------------
    coords
    coords_dot


Low-level functions:
--------------------
    c_atan2
    alpha1
    alpha2
    psi
    eta
    fg
    thetag


See test_derivatives function for examples of usage.

"""

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


    if (a**2+c**2) == 0:
        x_sign = math.copysign(1,x.real)
        y_sign = math.copysign(1,y.real)

        a = a + 1.e-12
        #c = c 

    #    #if x_sign > 0:
    #    #    if y.imag != 0:
    #    #        result = complex(math.atan2(a,c),0.)
    #    #    else:
    #    #        result = complex(math.atan2(a,c),0.)
    #    #else:
    #    result = complex(math.atan2(a,c),-.0000000000001)
    #    result = complex(math.atan2(a,c),(c*b-a*d)/(a**2+c**2))

    result = complex(math.atan2(a,c),(c*b-a*d)/(a**2+c**2))

    return result


def alpha1(t):
    """ Time activation function for cylinder case1

    short-time composite motion with deformation extension 

    """
    return t*t*t*(8.-3.*t)/16.

def alpha2(t):
    """ Time activation function for cylinder case2

    long-time composite motion with deformation extension 

    """
    return alpha_max*cmath.sin(omega_alpha*t)*(t**6./(t**6.+tref**6.))

def psi(t,alpha):
    return 1. + (Aa - 1.)*alpha(t)

def eta(lam,omega,tau):

    return cmath.sin(omega*lam + tau*(1. - cmath.cos(omega*lam)))

def fg(t,r0,theta0):
    """ Deformation function for non-unit geometry mapping Jacobian 

    Parameters
    ----------
    t: float or complex
        Time
    r0: float or complex
        Reference radius
    theta0: float or complex
        Reference theta


    """
    return (16.*(r0**4.) + (t**6./(t**6. + 0.01))*eta(t,10.,0.7)*(cmath.cos(32.*cmath.pi*r0**4.) - 1.))*eta(theta0,1.,0.7)


def thetag(t,r0,theta0):
    return theta0 + Ag*fg(t,r0,theta0)


def coords(alpha,x0,y0,t):
    """ Return the physical coordinates

    Parameters
    ----------
    alpha: function
        Time-activation function. Definded above (either alpha1 or alpha2)
    x0: float or complex
        Reference x-coordinate
    y0: float or complex
        Reference y-coordinate
    t: float or complex
        time

    """

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
    """ Return the derivatives of physical coordinates

    Parameters
    ----------
    alpha: function
        Time-activation function. Definded above (either alpha1 or alpha2)
    x0: float or complex
        Reference x-coordinate
    y0: float or complex
        Reference y-coordinate
    t: float or complex
        time

    """

    h = 1.e-30
    eps = complex(0.,h)

    # Issues related to the use of atan2 in coords function make the derivative 
    # at 0,0 undefined. The form of the motion simplifies a bit at (0,0) and we 
    # hard-code the spatial derivatives here
    if (x0 == 0) and (y0 == 0):
        dxdx = psi(t,alpha).real*math.cos(Atheta*alpha(t).real)
        dxdy = -(1./psi(t,alpha).real)*math.sin(Atheta*alpha(t).real)
        dydx = psi(t,alpha).real*math.sin(Atheta*alpha(t).real)
        dydy = (1./psi(t,alpha).real)*math.cos(Atheta*alpha(t).real)

    # Away from (0,0) we evaluate derivatives via complex-step
    else:
        # Evaluate function with complex perturbations in x,y
        dxdx,dydx = coords(alpha,x0+eps,y0    ,t    )
        dxdy,dydy = coords(alpha,x0    ,y0+eps,t    )

        # Divide imaginary part by step size
        dxdx = dxdx.imag/h
        dxdy = dxdy.imag/h
        dydx = dydx.imag/h
        dydy = dydy.imag/h

    # There are no complicating issues for the time derivative, evaluate everywhere with complex-step
    dxdt,dydt = coords(alpha,x0    ,y0    ,t+eps)
    dxdt = dxdt.imag/h
    dydt = dydt.imag/h

    return dxdx, dxdy, dxdt, dydx, dydy, dydt


def test_derivatives():

    x0 = -0.0
    y0 = 0.0
    t = 0.
    fd_eps = 1.e-10

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






