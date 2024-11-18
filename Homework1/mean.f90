program mean
    implicit none
    real :: first, second, third, arithmetic, harmonic, geometric 

    print*, "Please type in the first number."
    read*, first

    print*, "Please type in the second number."
    read*, second

    print*, "Please type in the third number."
    read*, third

    arithmetic = (first + second + third) / 3.0 
    harmonic = 1.0 / ((1.0/first + 1.0/second + 1.0/third)/ 3.0)
    geometric = (first*second*third)**(1.0/3.0)
    print*, "Arithmetic mean:", arithmetic ,"Harmonic mean:", harmonic, "Geometric mean:", geometric 
end program mean
