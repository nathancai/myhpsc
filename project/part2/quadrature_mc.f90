
module quadrature_mc

contains

    function quad_mc(g, a, b, ndim, npoints)

    implicit none
    real(kind=8) :: quad_mc, f, f2, x(ndim), volume
    real(kind=8), external :: g
    real(kind=8), intent(in), dimension(ndim) :: a, b
    integer, intent(in) :: ndim, npoints
    integer :: i
     
    quad_mc = 0.d0

    f= 0.d0
    
    volume = product(b-a)  ! =  product of b(i)-a(i) of ndim dimensions
    
    do i=1, npoints
        
        call random_number(x)
        
        x = x*volume/npoints

        f = f + g(x, ndim)
        

        end do
    
    
    
    quad_mc = f/npoints
    
    end function quad_mc

end module quadrature_mc
