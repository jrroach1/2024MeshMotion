program test_derivatives
    use cylinder2024,   only: rk, ZERO, coords, coords_dot, ALPHA_SHORT, ALPHA_LONG
    implicit none

    real(rk) :: fd_eps, &
                dxdx_fd, dydx_fd, dxdy_fd,  &
                dydy_fd, dxdt_fd, dydt_fd,  &
                dxdx_cs, dydx_cs, dxdy_cs,  &
                dydy_cs, dxdt_cs, dydt_cs

    complex(rk) :: x0, y0, t, xa, ya,       &
                   x_xdot, x_ydot, x_tdot,  &
                   y_xdot, y_ydot, y_tdot

    x0     = complex(0.0_rk, ZERO)
    y0     = complex(0.0_rk, ZERO)
    t      = complex(1._rk, ZERO)
    fd_eps = complex(1.e-10_rk, ZERO)

    ! Compute complex-step derivatives
    call coords_dot(ALPHA_SHORT,x0,y0,t,dxdx_cs,dxdy_cs,dxdt_cs,dydx_cs,dydy_cs,dydt_cs)

    ! Compute finite-difference derivatives
    call coords(ALPHA_SHORT,x0,       y0,       t       , xa,     ya    )
    call coords(ALPHA_SHORT,x0+fd_eps,y0       ,t       , x_xdot, y_xdot)
    call coords(ALPHA_SHORT,x0       ,y0+fd_eps,t       , x_ydot, y_ydot)
    call coords(ALPHA_SHORT,x0       ,y0       ,t+fd_eps, x_tdot, y_tdot)

    dxdx_fd = (x_xdot%re - xa%re)/fd_eps
    dydx_fd = (y_xdot%re - ya%re)/fd_eps

    dxdy_fd = (x_ydot%re - xa%re)/fd_eps
    dydy_fd = (y_ydot%re - ya%re)/fd_eps

    dxdt_fd = (x_tdot%re - xa%re)/fd_eps
    dydt_fd = (y_tdot%re - ya%re)/fd_eps

    print*, "Cylinder2024: Case1 short-time composite motion with deformation extension." 
    print*, " ", "Finite-difference(step=1e-10)", "Complex-step", "Error"
    print*, "dxdx", dxdx_fd, dxdx_cs, dxdx_fd - dxdx_cs
    print*, "dxdy", dxdy_fd, dxdy_cs, dxdy_fd - dxdy_cs
    print*, "dxdt", dxdt_fd, dxdt_cs, dxdt_fd - dxdt_cs
    print*, "dydx", dydx_fd, dydx_cs, dydx_fd - dydx_cs
    print*, "dydy", dydy_fd, dydy_cs, dydy_fd - dydy_cs
    print*, "dydt", dydt_fd, dydt_cs, dydt_fd - dydt_cs


    ! Compute complex-step derivatives
    call coords_dot(ALPHA_LONG,x0,y0,t,dxdx_cs,dxdy_cs,dxdt_cs,dydx_cs,dydy_cs,dydt_cs)

    ! Compute finite-difference derivatives
    call coords(ALPHA_LONG,x0,       y0,       t       , xa,     ya    )
    call coords(ALPHA_LONG,x0+fd_eps,y0       ,t       , x_xdot, y_xdot)
    call coords(ALPHA_LONG,x0       ,y0+fd_eps,t       , x_ydot, y_ydot)
    call coords(ALPHA_LONG,x0       ,y0       ,t+fd_eps, x_tdot, y_tdot)

    dxdx_fd = (x_xdot%re - xa%re)/fd_eps
    dydx_fd = (y_xdot%re - ya%re)/fd_eps

    dxdy_fd = (x_ydot%re - xa%re)/fd_eps
    dydy_fd = (y_ydot%re - ya%re)/fd_eps

    dxdt_fd = (x_tdot%re - xa%re)/fd_eps
    dydt_fd = (y_tdot%re - ya%re)/fd_eps

    print*, "Cylinder2024: Case2 long-time study of conservation."
    print*, " ", "Finite-difference(step=1e-10)", "Complex-step", "Error"
    print*, "dxdx", dxdx_fd, dxdx_cs, dxdx_fd - dxdx_cs
    print*, "dxdy", dxdy_fd, dxdy_cs, dxdy_fd - dxdy_cs
    print*, "dxdt", dxdt_fd, dxdt_cs, dxdt_fd - dxdt_cs
    print*, "dydx", dydx_fd, dydx_cs, dydx_fd - dydx_cs
    print*, "dydy", dydy_fd, dydy_cs, dydy_fd - dydy_cs
    print*, "dydt", dydt_fd, dydt_cs, dydt_fd - dydt_cs



end program test_derivatives

