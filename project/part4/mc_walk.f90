
module mc_walk

    use mpi
    use problem_description, only: ax, bx, ay, by, nx, ny, dx, dy, uboundary

    implicit none
    integer :: nwalks
    save

contains

    subroutine random_walk(i0, j0, max_steps, ub, iabort)
    
      !"""
      !Take one random walk starting at (i0,j0) until we reach the boundary or
      !exceed max_steps steps.
      !Return the value at the boundary point reached, or nan if we failed.
      !"""
    
      implicit none
      integer, intent(in) :: i0, j0, max_steps
      real(kind=8), intent(out) :: ub
      integer, intent(out) :: iabort
      
      integer :: i, j, iold, jold, istep
      real(kind=8) :: r(max_steps), xb, yb

      nwalks = nwalks + 1
      iabort = 0
      
      ! starting point:
      i = i0
      j = j0

      ! generate as many random numbers as we could possibly need
      ! for this walk, since this is much faster than generating one at a time:
      call random_number(r)     

      do istep=1, max_steps
	  iold = i
	  jold = j  ! needed for plotting only
      
	  ! Take the next random step with equal probability in each
	  ! direction:

	  if (r(istep) < 0.25d0) then
	      i = i-1   ! step left
	  else if (r(istep) < 0.5d0) then
	      i = i+1   ! step right
	  else if (r(istep) < 0.75d0) then
	      j = j-1   ! step down
	  else   
	      j = j+1   ! step up
	      endif



	  ! check if we hit the boundary:
	  if (i*j*(nx+1-i)*(ny+1-j) == 0) then
	      xb = ax + i*dx
	      yb = ay + j*dy
	      ub = uboundary(xb, yb)
		  
	      exit  ! end the walk
	      
	      endif

	  if (istep==max_steps) then
	    
	      iabort = 1
	      endif
	      
	      
	  enddo
	      
    end subroutine random_walk

    subroutine many_walks(i0, j0, max_steps, n_mc, u_mc, n_success)
    
      implicit none
      integer, intent(in) :: i0, j0, max_steps, n_mc
      integer, intent(out) :: n_success
      real(kind=8), intent(out) :: u_mc
      
      integer :: i, j, k, ierr, iabort, jj, num_procs, proc_num, numsent, sender
      real(kind=8) :: ub, ub_sum
      integer, dimension(MPI_STATUS_SIZE) :: status

      call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, proc_num, ierr)
   
      ub_sum = 0.d0   ! to accumulate boundary values reached from all walks
      n_success = 0    ! to keep track of how many walks reached boundary

      numsent = 0
      
  
      if (proc_num==0) then	  

	  do j=1,min(num_procs-1, n_mc)
	    call MPI_SEND(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION,&
				j, j, MPI_COMM_WORLD, ierr)
	    numsent = numsent + 1
	    enddo

	  do j=1,n_mc
	      call MPI_RECV(ub, 1, MPI_DOUBLE_PRECISION, &
			    MPI_ANY_SOURCE, MPI_ANY_TAG, &
			    MPI_COMM_WORLD, status, ierr)
	      sender = status(MPI_SOURCE)
	      
	      jj = status(MPI_TAG)
	      
	      if (jj==0) then
	          ub_sum = ub_sum + ub
	          n_success = n_success + 1
	      endif
	      
	      if (numsent < n_mc) then
		! still more work to do, the next column will be sent and
		! this index also used as the tag:
		i = numsent + 1
		
		call MPI_SEND(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION,&
				sender, i, MPI_COMM_WORLD, ierr)
				
		numsent = numsent + 1
	      else
		! send an empty message with tag=0 to indicate this worker
		! is done:
		call MPI_SEND(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION,&
				sender, 0, MPI_COMM_WORLD, ierr)
		endif

	      enddo
      endif

    ! -----------------------------------------
    ! code for Workers (Processors 1, 2, ...):
    ! -----------------------------------------
    if (proc_num /= 0) then
    
        if (proc_num > n_mc) go to 99   ! no work expected
        
        do while (.true.)

	    call MPI_RECV(MPI_ANY_TAG, 0, MPI_DOUBLE_PRECISION, &
			  0, MPI_ANY_TAG, &
			  MPI_COMM_WORLD, status, ierr)

	    j = status(MPI_TAG)   ! this is the subinterval number

	    if (j==0) go to 99
	    
	    call random_walk(i0, j0, max_steps, ub, iabort)

	    call MPI_SEND(ub, 1, MPI_DOUBLE_PRECISION, &
			0, iabort, MPI_COMM_WORLD, ierr)
	    
	   
	    enddo

        endif

99    continue
        
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      
      if (proc_num == 0) then

          u_mc = ub_sum / n_success   ! average over successful walks
          endif
      
      call MPI_BCAST(u_mc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(n_success, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    end subroutine many_walks
    
end module mc_walk