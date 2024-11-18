program main
    
    !Test 1 with x**2
    real :: quadrat(11) = [-25.0, -16.0, -9.0, -4.0, -1.0, 0.0, 1.0, 4.0, 9.0, 16.0, 25.0]
    real :: calc_sec_deriv_quadr(11)
    real :: sec_deriv_quadr(11) = [2.0,2.0, 2.0, 2.0,2.0, 2.0, 2.0, 2.0, 2.0,2.0, 2.0]
    !Test 2 with cos
    real :: sine(5) = [0.0,1.0, 0.0, -1.0, 0.0]
    real :: calc_secderiv_sine(5)
    real :: secderiv_sine(5) = [0.0, -1.0, 0.0, 1.0, 0.0]
    real :: h = 3.1415926/2.0

    !Test 1
    call second_derivative(quadrat,11,1.0 ,calc_sec_deriv_quadr)
    print*, "Calculated second derivatives of quadratic function", calc_sec_deriv_quadr
    print*, "Expected second derivative of quadratic function", sec_deriv_quadr 
    print*, "Error", sec_deriv_quadr - calc_sec_deriv_quadr

    !Test 2
    call second_derivative(sine, 5, h,calc_secderiv_sine )
    print*, "Calculated second dericatives of cosine", calc_secderiv_sine
    print*, "Expected second derivative of cosine", secderiv_sine
    print*, "Error", secderiv_sine - calc_secderiv_sine
end program main


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
