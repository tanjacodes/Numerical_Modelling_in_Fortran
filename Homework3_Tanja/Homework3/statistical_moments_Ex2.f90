module mean_and_SD
contains    
    subroutine statistical_moments(array1D,N, mean, SD )
   
        implicit none
        integer, intent(in) :: N
        real(8), intent(in):: array1D(N)
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
            SD =  sqrt(SD/N - mean*mean) 
        !invalid input
        else 
            print*, "Error:The input number is 0 or negative"
        end if
    end subroutine statistical_moments
end module mean_and_SD


program main_program
    use mean_and_SD
    implicit none
    integer :: N
    real(8), allocatable :: array(:)
    real :: mean, SD, b
    integer :: i

    

    open(1, file = 'MeanStdInput.txt', status = 'old')

    N = 0
    do 
        read(1,*,iostat = i) b
        if (i<0) exit
        if (i /= 0) stop 'error reading data'
        N = N + 1
    end do
    print*, 'found', N, 'values'
    
    allocate(array(N))

    rewind(1) ! moves file pointer back to start
    
    do i = 1,N
        read(1,*) array(i)
        print*, array(i)
    end do
    
    call statistical_moments(array,N,mean,SD)
    print*, "The mean is:", mean, "The standard deviation is:", SD
    deallocate(array)
end program main_program

