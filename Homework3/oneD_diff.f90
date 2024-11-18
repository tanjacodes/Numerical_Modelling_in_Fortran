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


program oneDdiffusion
    use derivatives
    implicit none
    integer:: nx, nsteps,i,j
    real:: dx, k = 1.0, total_time = 0.1, a, dt
    character(len =10):: Tinit
    real, allocatable :: T(:), DD(:)
    
    
    namelist /inputs/ nx,k,total_time,a,Tinit
    
    open(1,file = 'oneD_diffusion_input.txt', status= 'old')
    read(1,inputs)
    close(1)

    
    
!divide by nx because, there are actually have nx+1 gridpoints
!used this convention from slide "Finite Difference grid in 1D"
    dx = 1.0/nx 
    dt = a*dx*dx/k
    nsteps = total_time/dt
    allocate(T(0:nx))
    allocate(DD(0:nx))
    if (Tinit == 'spike') then
    	do i = 0,nx	
    		T(i) = 0.0
    	end do
    	T(nx/2) = 1.0
    else if (Tinit == 'random') then
    	T(0) = 0.0
    	T(nx) = 0.0
    	do i = 1,nx-1	
    		call random_number(T(i))
    	end do
    else
    	print*, "Error, Tinit must be spike or random"
    end if
    
    
    do i = 1,nsteps
    	call second_derivative(T,nx+1,dx,DD)
    	do j = 1,nx-1 !makes sure boundaries at 0 and nx are constant 0
    		T(j) = T(j) + dt*k*DD(j)
    	end do
    end do
    
    deallocate(T)
    deallocate(DD)
end program oneDdiffusion
