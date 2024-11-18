program statistical_moments
   
     implicit none
    integer :: N,i
    real :: mean, SD, input_number

    ! Asking for and storing the number of input numbers
    print*, "Please enter the number of real numbers in the series you would like to enter:)."
    read* , N

    !If the input is valid: calculating mean and standard deviation SD
    if (N > 0) then
        ! will store sum of input values
        mean = 0.0
        ! will store sum of squares
        SD = 0.0
        do i = 1,N
            print*, "Please enter the next number in the series"
            read*, input_number
            mean = mean + input_number
            SD = SD + input_number*input_number 
        end do
        ! the mean is the sum of input values divided by the number of input values
        mean = mean/N
        ! the SD is the sum of squares divided by N - mean*mean and the square root of the whole expression
        SD = (SD/N - mean*mean)**(0.5)

        !printing the results
        print*, "The mean is:", mean, "The standard deviation is:", SD 
    !invalid input
    else 
        print*, "Error:The input number is 0 or negative"
    end if


end program statistical_moments
