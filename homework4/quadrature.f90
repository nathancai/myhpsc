! $MYHPSC/homework4/quadrature.f90
! Nathan Juswanto

module quadrature

contains

real(kind=8) function trapezoid(f, a, b, n)

    implicit none
    integer, intent(in) :: n
    real(kind = 8), external :: f
    real(kind = 8), intent(in) :: a, b
    
    integer :: i
    real(kind = 8) :: fj(n), h, xj(n)    
    
    h  = (b - a)/(n - 1.d0)
    
    xj(1) = a
    fj(1) = .5d0*f(xj(1))
    do i = 2, n - 1
        xj(i) = xj(i - 1) + h
        fj(i) = f(xj(i))
        end do
    xj(n) = b
    fj(n) = .5d0*f(xj(n))
    
    trapezoid = h*sum(fj)
    
end function trapezoid

subroutine error_table(f, a, b, nvals, int_true)

    implicit none
    real(kind=8), external :: f
    real(kind=8), intent(in) :: a, b, int_true
    integer, dimension(:), intent(in) :: nvals
    
    integer :: n
    real(kind=8) :: error, int_trap, last_error, ratio
    
    print *, "    n         trapezoid            error       ratio"

    last_error = 0.d0  ! need something for first ratio
    do n = 1, size(nvals)
        int_trap = trapezoid(f, a, b, nvals(n))
        error = abs(int_trap - int_true)
        ratio = last_error/error
        last_error = error ! for next n
        print 11, nvals(n), int_trap, error, ratio
11      format(i8, es22.14, es13.3, es13.3)
        
        enddo
    
end subroutine error_table

end module quadrature
 
