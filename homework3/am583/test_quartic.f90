! $MYHPSC/homework3/am583/test1.f90
! Nathan Juswanto

program test_quartic

    use newton, only: solve, tol
    use functions, only: f_quartic, fprime_quartic, epsilon

    implicit none
    real(kind=8) :: x, x0, fx, x0vals, xstar, epsilonvals(3), tolvals(4)
    integer :: iters, itest, jtest
    logical :: debug         ! set to .true. or .false.

    debug = .false.

    ! values to test as x0:
    x0vals      = .4d1
    epsilonvals = (/ .1d-3, .1d-7, .1d-11 /)
    tolvals     = (/ .1d-4, .1d-9, .1d-13, .1d-18/)

    print 11, x0vals
11  format('Starting with initial guess', d24.15)
    print *, ' '  ! blank line
    print *, '    epsilon        tol    iters          x                 f(x)        x-xstar'
    
    do itest=1,3
        
        ! xstar
        epsilon = epsilonvals(itest)
        tol     = tolvals(4)
        
        call solve(f_quartic, fprime_quartic, x0, x, iters, debug)
        
        xstar = x
           
        do jtest=1,3
            x0      = x0vals
            tol     = tolvals(jtest)
            
            call solve(f_quartic, fprime_quartic, x0, x, iters, debug)
            fx = f_quartic(x)
                        
            print 12, epsilon, tol, iters, x, fx, x-xstar
12          format(2d13.3, i4, d24.15, 2d13.3)

            enddo       

            print *, ' '  ! blank line
            
            enddo

end program test_quartic
