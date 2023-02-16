!>
! COMPILE:
!   compile:  gfortran cylinder2024.f90
!   produces: cylinder2024.mod, a.out
!  
!   cylinder2024.mod can be linked with your fortran application 
!   and "use"ed as follows
!
! USAGE EXAMPLE:
!   use cyliner2024,    only: coords, coords_dot, ALPHA_SHORT
!
!   call coords(ALPHA_SHORT,xref,yref,t,x,y)
!   call coords_dot(ALPHA_SHORT,xref,yref,t,dxdx,dxdy,dxdt,dydx,dydy,dydt)
!
!
! FURTHER EXAMPLE:
! See the program "test_derivatives" at the bottom of this file
! as an example of usage.
!
! NOTE:
! Since this implementation utilizes the complex-step method
! to obtain derivatives of the prescribed motion functions,
! the input reference coordinates and time must be declared
! and input as complex variables with 0 imaginary component.
!
!-----------------------------------------------------------------
module cylinder2024
    implicit none

    integer, parameter :: rk = selected_real_kind(15, 307 )

    integer, parameter :: ALPHA_SHORT = 1
    integer, parameter :: ALPHA_LONG  = 2

    ! Floating-point numbers
    real(rk), parameter     :: ZERO     = 0._rk
    real(rk), parameter     :: ONE      = 1._rk
    real(rk), parameter     :: TWO      = 2._rk
    real(rk), parameter     :: THREE    = 3._rk
    real(rk), parameter     :: FOUR     = 4._rk
    real(rk), parameter     :: FIVE     = 5._rk
    real(rk), parameter     :: SIX      = 6._rk
    real(rk), parameter     :: SEVEN    = 7._rk
    real(rk), parameter     :: EIGHT    = 8._rk
    real(rk), parameter     :: NINE     = 9._rk
    real(rk), parameter     :: TEN      = 10._rk
    real(rk), parameter     :: PI       = 3.14159265358979323846_rk

    ! Cylinder parameters
    real(rk) :: rcyl        = 0.5_rk
    real(rk) :: Atheta      = PI
    real(rk) :: Aa          = 1.5_rk
    real(rk) :: Ag          = 0.15_rk
    real(rk) :: alpha_max   = 0.3_rk
    real(rk) :: omega_alpha = 6._rk
    real(rk) :: tref        = 1._rk


contains

    !>
    !!
    !!
    !----------------------------------------------
    function c_atan2(y,x) result(theta)
        complex(rk),   intent(in)  :: y
        complex(rk),   intent(in)  :: x
        
        real(rk) :: a, b, c, d
        complex(rk) :: theta

        a = y%re
        b = y%im
        c = x%re
        d = x%im

        if ((a**TWO + c**TWO) == 0.) then
            theta = cmplx(atan2(a,c),ZERO, rk)
        else
            theta = cmplx(atan2(a,c),(c*b-a*d)/(a**TWO + c**TWO), rk)
        end if

    end function c_atan2
    !***********************************************




    !>
    !!
    !!
    !-----------------------------------------------
    function alpha_generic(alpha,t) result(res)
        integer,     intent(in) :: alpha
        complex(rk), intent(in) :: t

        complex(rk) :: res

        if (alpha == ALPHA_SHORT) then
            res = alpha1(t)
        else if (alpha == ALPHA_LONG) then
            res = alpha2(t)
        else
            print*, "Invalid selection of alpha time-activation function."
            stop 
        end if

    end function alpha_generic
    !***********************************************


    !>
    !!
    !!
    !-----------------------------------------------
    function alpha1(t) result(res)
        complex(rk),    intent(in)  :: t

        complex(rk) :: res

        res = t*t*t*(EIGHT-THREE*t)/16._rk

    end function alpha1
    !***********************************************


    !>
    !!
    !!
    !-----------------------------------------------
    function alpha2(t) result(res)
        complex(rk),    intent(in)  :: t

        complex(rk) :: res

        res = alpha_max*sin(omega_alpha*t)*(t**SIX/(t**SIX+tref**SIX))

    end function alpha2
    !***********************************************


    !>
    !!
    !--------------------------------------------------
    function psi(t,alpha) result(res)
        complex(rk),    intent(in)  :: t
        integer,        intent(in)  :: alpha

        complex(rk) :: res

        res = ONE + (Aa - ONE)*alpha_generic(alpha,t)

    end function psi
    !***********************************************
    

    !>
    !!
    !--------------------------------------------------
    function eta(lam,omega,tau) result(res)
        complex(rk),    intent(in)  :: lam
        real(rk),       intent(in)  :: omega
        real(rk),       intent(in)  :: tau

        complex(rk) :: res
    
        res = sin(omega*lam + tau*(ONE - cos(omega*lam)))

    end function eta
    !***********************************************
    

    !> Deformation function for non-unit geometry mapping Jacobian 
    !!
    !!
    !--------------------------------------------------
    function fg(t,r0,theta0) result(res)
        complex(rk),    intent(in)  :: t
        complex(rk),    intent(in)  :: r0
        complex(rk),    intent(in)  :: theta0

        complex(rk) :: res

        res = (16._rk*(r0**FOUR) + (t**SIX/(t**SIX + 0.01_rk)) * &
              eta(t,TEN,0.7_rk)*(cos(32._rk*PI*r0**FOUR) - ONE))*eta(theta0,ONE,0.7_rk)
    
    end function fg
    !**************************************************
    

    !>
    !!
    !--------------------------------------------------
    function thetag(t,r0,theta0) result(res)
        complex(rk),    intent(in)  :: t
        complex(rk),    intent(in)  :: r0
        complex(rk),    intent(in)  :: theta0

        complex(rk) :: res

        res = theta0 + Ag*fg(t,r0,theta0)
    
    end function thetag
    !**************************************************
    


    !>
    !!  """ Return the physical coordinates
    !!
    !!  Parameters
    !!  ----------
    !!  alpha: function
    !!      Time-activation function. Definded above (either alpha1 or alpha2)
    !!  x0: float or complex
    !!      Reference x-coordinate
    !!  y0: float or complex
    !!      Reference y-coordinate
    !!  t: float or complex
    !!      time
    !!
    !--------------------------------------------------
    subroutine coords(alpha,xref,yref,t,x,y) 
        integer,        intent(in)      :: alpha
        complex(rk),    intent(in)      :: xref
        complex(rk),    intent(in)      :: yref
        complex(rk),    intent(in)      :: t
        complex(rk),    intent(inout)   :: x
        complex(rk),    intent(inout)   :: y

        complex(rk) :: r0, theta0, xdef0, ydef0


        r0     = sqrt(xref*xref + yref*yref)
        theta0 = c_atan2(yref,xref)

        xdef0 = r0*cos(thetag(t,r0,theta0))
        ydef0 = r0*sin(thetag(t,r0,theta0))

        x = cos(Atheta*alpha_generic(alpha,t))*psi(t,alpha)*xdef0 - &
           (sin(Atheta*alpha_generic(alpha,t))/psi(t,alpha))*ydef0
        y = sin(Atheta*alpha_generic(alpha,t))*psi(t,alpha)*xdef0 + &
           (cos(Atheta*alpha_generic(alpha,t))/psi(t,alpha))*ydef0 + alpha_generic(alpha,t)

    end subroutine coords
    !**************************************************



    !>
    !!
    !!
    !!
    !--------------------------------------------------
    subroutine coords_dot(alpha,xref,yref,t,dxdx,dxdy,dxdt,dydx,dydy,dydt)
        integer,        intent(in)      :: alpha
        complex(rk),    intent(in)      :: xref
        complex(rk),    intent(in)      :: yref
        complex(rk),    intent(in)      :: t
        real(rk),       intent(inout)   :: dxdx
        real(rk),       intent(inout)   :: dxdy
        real(rk),       intent(inout)   :: dxdt
        real(rk),       intent(inout)   :: dydx
        real(rk),       intent(inout)   :: dydy
        real(rk),       intent(inout)   :: dydt

        complex(rk) :: eps, my_psi, my_alpha, &
                       dxdx_c, dxdy_c, dxdt_c, &
                       dydx_c, dydy_c, dydt_c
        real(rk) :: h
    
        h = 1.e-30_rk
        eps = complex(ZERO,h)
    
        ! Issues related to the use of atan2 in coords function make the derivative 
        ! at 0,0 undefined. The form of the motion simplifies a bit at (0,0) and we 
        ! hard-code the spatial derivatives here
        if ( (abs(xref%re) < 1.e-13_rk) .and. &
             (abs(yref%re) < 1.e-13_rk) ) then 
            my_psi = psi(t,alpha)
            my_alpha = alpha_generic(alpha,t)
            dxdx = my_psi%re*cos(Atheta*my_alpha%re)
            dxdy = -(ONE/my_psi%re)*sin(Atheta*my_alpha%re)
            dydx = my_psi%re*sin(Atheta*my_alpha%re)
            dydy = (ONE/my_psi%re)*cos(Atheta*my_alpha%re)
    
        ! Away from (0,0) we evaluate derivatives via complex-step
        else
            ! Evaluate function with complex perturbations in x,y
            call coords(alpha,xref+eps,yref    ,t, dxdx_c, dydx_c)
            call coords(alpha,xref    ,yref+eps,t, dxdy_c, dydy_c)
    
            ! Divide imaginary part by step size
            dxdx = dxdx_c%im/h
            dxdy = dxdy_c%im/h
            dydx = dydx_c%im/h
            dydy = dydy_c%im/h
        end if
    
        ! There are no complicating issues for the time derivative, evaluate everywhere with complex-step
        call coords(alpha,xref,yref,t+eps, dxdt_c, dydt_c)
        dxdt = dxdt_c%im/h
        dydt = dydt_c%im/h

    end subroutine coords_dot
    !**************************************************


end module cylinder2024



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

    x0     = complex(0.5_rk, ZERO)
    y0     = complex(0.5_rk, ZERO)
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
