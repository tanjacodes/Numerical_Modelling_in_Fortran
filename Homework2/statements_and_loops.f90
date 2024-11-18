program statements
    implicit none
    
    !declarations for STATEMENTS
    !Declaring a string of length 15 char
    character(len=15) :: string
    !Declaring an integer called parameter with value 5
    integer :: parameter = 5
    !Declaring an 1D array with indices running from -1 to 10
    real :: array(-1:10)
    !Declaring 4 dim allocatable array
    real, allocatable:: alloc(:,:,:,:)

    !Declaring a &b for remainder calculation 
    real :: a,b

    !declarations  for LOOP
    !Declaring an array arr(1:00) for the last subexercise
    real :: arr(1:100)
    !Declaring sum to sum up the sum from 12 to 124
    integer :: sum
    !Declaring "Laufvariable" i, for the loops
    integer :: i

    !STATEMENT
    !Converting a real number to the nearest integer
    real :: number
    print*, "Type in a real number that should be rounded to the next int"
    read*, number
    print*, "Nearest integer is", Nint(number)

    !STATEMENT
    !Reading two numbers and calculating the remainder
    print*, "Please insert two number to calculate the remainder of a/b, where a is the first number and b the second"
    read*, a,b
    print*, "The remainder is given by", mod(a,b)

    !LOOPs
    !adding all numbers from 12 to 124
    do i = 12,124,2
        sum = sum + i
    end do
    print*, "Adding all even number from 12 to 124 results in:", sum

    !Test each element of array a(1:100) 
    do i = 1, 100
        if (arr(i) > 0) then
            print*, "A message is printed, because the array a(1:100) contains a positive element"
            stop
        end if
    end do
end program statements
