! $MYHPSC/homework3/am583/newton.f90
! Nathan Juswanto

module newton

    ! module parameters:
    implicit none
    integer, parameter :: maxiter = 40
    real(kind=8) :: tol 
    save

contains

subroutine solve(f, fp, x0, x, iters, debug)

    ! Estimate the zero of f(x) using Newton's method. 
    ! Input:
    !   f:  the function to find a root of
    !   fp: function returning the derivative f'
    !   x0: the initial guess
    !   debug: logical, prints iterations if debug=.true.
    ! Returns:
    !   the estimate x satisfying f(x)=0 (assumes Newton converged!) 
    !   the number of iterations iters
     
    implicit none
    real(kind=8), intent(in) :: x0
    real(kind=8), external :: f, fp
    logical, intent(in) :: debug
    real(kind=8), intent(out) :: x
    integer, intent(out) :: iters

    ! Declare any local variables:
    real(kind=8) :: deltax, fx, fxprime
    integer :: k

    ! initial guess
    x = x0
    k = 1

    if (debug) then
        print 11, x
 11     format('Initial guess: x = ', e22.15)
        endif

    ! Newton iteration to find a zero of f(x) 

    fx = f(x)
    do while (abs(fx) > tol) 

        ! evaluate function and its derivative:
        fx = f(x)
        fxprime = fp(x)

        if ((abs(fx) < tol) .or. (k > maxiter)) then
            exit  ! jump out of do loop
            endif

        ! compute Newton increment x:
        deltax = fx/fxprime

        ! update x:
        x = x - deltax

        if (debug) then
            print 12, k, x
 12         format('After', i3, ' iterations, x = ', e22.15)
            endif
            
        k = k + 1

        enddo

    ! number of iterations taken:
    iters = k - 1

end subroutine solve

end module newton
