
program laplace_mc

    use problem_description, only: ax, bx, ay, by, nx, ny, dx, dy, utrue, uboundary
    use mc_walk, only: nwalks, random_walk, many_walks
    use random_util, only: init_random_seed

    implicit none
    real(kind=8) :: error, x0, y0, u_mc, u_mc_total, u_true, u_sum_new, u_sum_old
    integer :: i, i0, j0, max_steps, n_mc, n_success, n_total, seed1 = 12345

    nwalks = 0
    open(unit=25, file='mc_laplace_error.txt', status='unknown')
    
    ! Try it out from a specific (x0,y0):
    x0 = 0.9d0
    y0 = 0.6d0
    
    i0 = idnint((x0-ax)/dx)
    j0 = idnint((y0-ay)/dy)

    ! shift (x0,y0) to a grid point if it wasn't already:
    x0 = ax + i0*dx
    y0 = ay + j0*dy

    u_true = utrue(x0,y0)

    call init_random_seed(seed1)  

    ! maximum number of step in each before giving up:
    max_steps = 100*max(nx, ny)

    ! initial number of Monte-Carlo walks to take:
    n_mc = 10

    call many_walks(i0, j0, max_steps, n_mc, u_mc, n_success)

    error = abs((u_mc - u_true) / u_true)

    print '("After ", i8, " random walks, u = ", es15.9, " rel. error = ", es15.6)', n_success, u_mc, error
    write(25,'(i10,e23.15,e15.6)') n_success, u_mc, error

    ! start accumulating totals:
    u_mc_total = u_mc
    n_total = n_success
    
    do i=1,12
        
        u_sum_old = u_mc_total * n_total
        call many_walks(i0, j0, max_steps, n_mc, u_mc, n_success)
        u_sum_new = u_mc * n_success
        n_total = n_total + n_success
        u_mc_total = (u_sum_old + u_sum_new) / n_total
        error = abs((u_mc_total - u_true) / u_true)

        print '("After ", i8, " random walks, u = ", es15.9, " rel. error = ", es15.6)', n_total, u_mc_total, error
        write(25,'(i10,e23.15,e15.6)') n_success, u_mc, error
        n_mc = 2*n_mc   ! double number of trials for next iteration
        
        enddo
        
    print '("Final approximation to u(x0,y0) : ", es22.14)', u_mc_total
    print *, "Total number of random walks: ", nwalks
    
end program laplace_mc