! $MYHPSC/homework3/intersections.f90
! Nathan Juswanto

program intersections

    use newton, only: initialize, solve2
    use functions, only: f_sqrt2, fprime_sqrt2

    implicit none
    real(kind=8) :: x, x0, fx
    real(kind=8) :: x0vals(4)
    integer :: iters, itest
	logical :: debug         ! set to .true. or .false.

    call initialize()   ! sets pi
    
    print *, "Test routine for computing zero of f"
    debug = .true.

    ! values to test as x0:
    x0vals = (/-2.d0, -1.5d0, -1.d0, 1.5d0/)

    do itest=1,4
        x0 = x0vals(itest)
		print *, ' '  ! blank line
        call solve2(f_sqrt2, fprime_sqrt2, x0, x, iters, debug)

        print 11, x, iters
11      format('solve returns x = ', es22.15, ' after', i3, ' iterations')

        fx = 1d0-.6d0*x**2
        print 12, fx
12      format('the value of f(x) is ', es22.15)

        enddo

end program intersections