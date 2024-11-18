module mean_and_SD
contains    
    subroutine statistical_moments(array1D,N, mean, SD )
   
        implicit none
        integer, intent(in) :: N
        real, intent(in):: array1D(N)
        real, intent(out):: mean, SD
        real :: number
        integer :: i
        !If the input is valid: calculating mean and standard deviation SD
        if (N > 0) then
            ! will store sum of input values
            mean = 0.0
            ! will store sum of squares
            SD = 0.0
            do i = 1,N
                number = array1D(i)
                mean = mean + number
                SD = SD + number*number 
            end do
            ! the mean is the sum of input values divided by the number of input values
            mean = mean/N
            ! the SD is the sum of squares divided by N - mean*mean and the square root of the whole expression
            SD = (SD/N - mean*mean)**(0.5) 
        !invalid input
        else 
            print*, "Error:The input number is 0 or negative"
        end if
    end subroutine statistical_moments
end module mean_and_SD
