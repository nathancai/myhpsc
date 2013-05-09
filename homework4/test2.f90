! $MYHPSC/homework4/test2.f90
! Nathan Juswanto

program test2

    use quadrature, only: trapezoid, error_table

    implicit none
    real(kind=8) :: a2, b2, int_true2, k
    integer :: nvals(12), i

    a2 = 0.d0
    b2 = 2.d0
    k  = 1000.d0
    int_true2 = (b2 - a2) + (b2**4 - a2**4)/4.d0 - (1.d0/k)*(cos(k*b2) - cos(k*a2))

    print 10, int_true2
 10 format("true integral: ", es22.14)
    print *, " "  ! blank line

    ! values of n to test:
    do i = 1, size(nvals)
        nvals(i) = 5 * 2**(i-1)
        enddo

    call error_table(f2, a2, b2, nvals, int_true2)

contains

    real(kind=8) function f2(x)
        implicit none
        real(kind=8), intent(in) :: x 
        
        f2 = 1.d0 + x**3 + sin(k*x)
    end function f2

end program test2
 
