
module functions

    use omp_lib
    implicit none
    integer :: fevals(0:7), gevals(0:7)
    real(kind=8) :: k
    save

contains

    real(kind=8) function f(x)
        implicit none
        real(kind=8), intent(in) :: x 
        integer thread_num

        ! keep track of number of function evaluations by
        ! each thread:
        thread_num = 0   ! serial mode
        !$ thread_num = omp_get_thread_num()
        fevals(thread_num) = fevals(thread_num) + 1
        
        f = -cos(x + 4.d0) + cos(x + 1.d0)
        
    end function f
    
    real(kind=8) function g(x, y)
        implicit none
        real(kind=8), intent(in) :: x, y
        integer thread_num

        ! keep track of number of function evaluations by
        ! each thread:
        thread_num = 0   ! serial mode
        !$ thread_num = omp_get_thread_num()
        gevals(thread_num) = gevals(thread_num) + 1
        
        g = sin(x + y)
        
    end function g    

end module functions
