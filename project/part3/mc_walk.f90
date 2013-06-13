
module mc_walk

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
      
      integer :: i, j, k, iabort
      real(kind=8) :: ub, ub_sum

   
      ub_sum = 0.d0   ! to accumulate boundary values reached from all walks
      n_success = 0    ! to keep track of how many walks reached boundary

      do k=1,n_mc
	  i = i0
	  j = j0
	  call random_walk(i0, j0, max_steps, ub, iabort)
	  if (iabort==0) then
	      ! use this result unless walk didn't reach boundary
	      ub_sum = ub_sum + ub
	      n_success = n_success + 1
	      endif
	  enddo

      u_mc = ub_sum / n_success   ! average over successful walks

    end subroutine many_walks
    
end module mc_walk