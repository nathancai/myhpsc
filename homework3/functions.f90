! $MYHPSC/homework3/functions.f90
! Nathan Juswanto

module functions

    implicit none
    real(kind=8) :: pi
    save

contains

real(kind=8) function f_sqrt(x)
    implicit none
    real(kind=8), intent(in) :: x

    f_sqrt = x**2 - 4.d0

end function f_sqrt


real(kind=8) function fprime_sqrt(x)
    implicit none
    real(kind=8), intent(in) :: x
    
    fprime_sqrt = 2.d0 * x

end function fprime_sqrt

real(kind=8) function f_sqrt2(x)
    implicit none
    real(kind=8), intent(in) :: x

    f_sqrt2 = x*cos(pi*x)-1.d0+.6d0*x**2

    ! x*np.cos(np.pi*x)-1+.6*x**2
end function f_sqrt2


real(kind=8) function fprime_sqrt2(x)
    implicit none
    real(kind=8), intent(in) :: x
    
    fprime_sqrt2 = cos(pi*x)+1.2d0*x-x*pi*sin(pi*x)

    ! x*np.pi*-np.sin(np.pi*x)+np.cos(np.pi*x)+1.2*x
end function fprime_sqrt2

end module functions
