
program test

    use omp_lib

    use quadrature, only: trapezoid, error_table
    use functions, only: f, g, fevals, gevals, k

    implicit none
    real(kind=8) :: a,b,c,d,int_true
    integer :: nxvals(10), nyvals(10)
    integer :: i, nthreads

    real(kind=8) :: t1, t2, elapsed_time
    integer(kind=8) :: tclock1, tclock2, clock_rate

    nthreads = 1      ! for serial mode
    !$ nthreads = 4   ! for openmp
    !$ call omp_set_num_threads(nthreads)
    print 100, nthreads
100 format("Using ",i2," threads")

    fevals = 0

    k = 1.d3   ! functions module variable for function f2
    a = 0.d0
    b = 2.d0
    c = 1.d0
    d = 4.d0
    int_true = -sin(b + 4) + sin(a + 4) + sin(b + 1) - sin(a + 1)
    

    print 10, int_true
 10 format("true integral: ", es22.14)
    print *, " "  ! blank line

    ! values of n to test:   (larger values than before)
    do i=1,10
        nxvals(i) = 5 * 2**(i-1)
        nyvals(i) = 1000
        enddo
  
    ! time the call to error_table:
    call system_clock(tclock1)  
    call cpu_time(t1)
    call error_table(f, a, b, nxvals, int_true, trapezoid)
    call cpu_time(t2)   
    call system_clock(tclock2, clock_rate)

    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, " "
    print 11, elapsed_time
 11 format("Elapsed time = ",f12.8, " seconds")

    print 12, t2-t1
 12 format("CPU time = ",f12.8, " seconds")

    
    ! print the number of function evaluations by each thread:
    do i=0,nthreads-1
        print 101,  i, fevals(i)
101     format("fevals by thread ",i2,": ",i13)
        enddo

    print 102, sum(fevals)
102 format("Total number of fevals: ",i10)

    do i=0,nthreads-1
        print 103,  i, gevals(i)
103     format("gevals by thread ",i2,": ",i13)
        enddo

    print 104, sum(gevals)
104 format("Total number of gevals: ",i10)

end program test
