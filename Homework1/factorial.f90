program factorial
    implicit none
    integer :: N, i, result

    print*, "Calculating N!; Please enter the non-negative integer N:"
    read*, N

    if (N < 0) then
        print*, "Error: You typed in a non-negative number."
    else if (N == 0) then
        print*, "0! = 1"
    else
        result = 1
        Do i = 1,N
            result = result * i
        end Do
        print*, N, '! = ', result
    end if

end program factorial
