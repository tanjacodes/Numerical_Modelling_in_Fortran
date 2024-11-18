module derivatives
contains
    subroutine second_derivative (array, N, h, primeprime)
        integer, intent(in) :: N
        integer :: i
        real, intent(in) :: h
        real, intent(in) :: array(N)
        real, intent(out) :: primeprime(N)

        primeprime(1) = 0.0
        do i = 2,N-1
            primeprime(i) = (array(i+1) - 2* array(i) + array(i-1))/(h*h) 
        end do
        primeprime(N) = 0.0
    
    end subroutine second_derivative
end module derivatives
