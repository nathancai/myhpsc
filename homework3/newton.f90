! $MYHPSC/homework3/newton.f90
! Nathan Juswanto

module newton

    ! module parameters:
    implicit none
    integer, parameter :: maxiter = 20
    real(kind=8), parameter :: tol = 1.d-14

contains

subroutine initialize()

    ! Set the value of pi used elsewhere.
    use functions, only: pi
    pi = acos(-1.d0)

end subroutine initialize

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

    if (debug) then
        print 11, x
 11     format('Initial guess: x = ', e22.15)
        endif

    ! Newton iteration to find a zero of f(x) 

    do k=1,maxiter

        ! evaluate function and its derivative:
        fx = f(x)
        fxprime = fp(x)

        if (abs(fx) < tol) then
            exit  ! jump out of do loop
            endif

        ! compute Newton increment x:
        deltax = fx/fxprime

        ! update x:
        x = x - deltax

        if (debug) then
            print 12, k,x
 12         format('After', i3, ' iterations, x = ', e22.15)
            endif

        enddo

    if (k > maxiter) then
        ! might not have converged

        fx = f(x)
        if (abs(fx) > tol) then
            print *, '*** Warning: has not yet converged'
            endif
        endif 

    ! number of iterations taken:
    iters = k-1

end subroutine solve

subroutine solve2(f, fp, x0, x, iters, debug)

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
        print 21, x
 21     format('Initial guess: x = ', es22.15)
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
            print 22, k,x
 22         format('After', i3, ' iterations, x = ', es22.15)
            endif
        
        k = k + 1

        enddo

    if (k > maxiter) then
        ! might not have converged

        fx = f(x)
        if (abs(fx) > tol) then
            print *, '*** Warning: has not yet converged'
            endif
        endif 

    ! number of iterations taken:
    iters = k-1

end subroutine solve2

end module newton
