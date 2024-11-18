module derivatives
contains
    subroutine second_derivative (array, Nx,Ny, dx, dy, primeprime)
        integer, intent(in) :: Nx,Ny
        integer :: i,j
        real, intent(in) :: dx, dy
        real, intent(in) :: array(Nx,Ny)
        real, intent(out) :: primeprime(Nx,Ny)

    	forall(i=1:Ny) primeprime(1,i) = 0.0
        forall(i=1:Nx) primeprime(i,1) = 0.0
    	forall(i=1:Nx) primeprime(i,ny) = 0.0
        forall(i=1:Ny) primeprime(nx,i) = 0.0
    
        do i = 2,Nx-1
            do j = 2, Ny-1
                primeprime(i,j) = (array(i-1,j) + array(i+1,j) - 2*array(i,j))/(dx*dx)
                primeprime(i,j) = primeprime(i,j) + (array(i,j-1) + array(i, j+1)-2*array(i,j))/(dy*dy) 
            end do    
        end do
    
    end subroutine second_derivative
end module derivatives


program twoDdiffusion
    use derivatives
    implicit none
    integer:: nx,ny, nsteps,i,j,h
    real:: dx,dy, k = 1.0, total_time = 0.1, a, dt
    character(len =10):: Tinit
    real, allocatable :: T(:,:), DD(:,:)
    
    
    namelist /inputs/ nx,ny,k,total_time,a,Tinit
    
    open(1,file = 'twoD_diffusion_input.txt', status= 'old')
    read(1,inputs)
    close(1)
	print*, "Reading inputs done"
    dx = 1.0/ny
    dy= dx
    dt = a*dx*dx/k
    nsteps = total_time/dt
    allocate(T(0:nx,0:ny))
    allocate(DD(0:nx, 0:ny))
    if (Tinit == 'spike') then
    	T = 0.0
        T(nx/2, ny/2) = 1.0
    else if (Tinit == 'random') then
    	forall(i=0:Ny) T(0,i) = 0.0
        forall(i=0:Nx) T(i,0) = 0.0
    	forall(i=0:Nx) T(i,ny) = 0.0
        forall(i=0:Ny) T(nx,i) = 0.0
    	do i = 1,nx-1	
            do j= 1,ny-1
    		    call random_number(T(i,j))
            end do
    	end do
    else
    	print*, "Error, Tinit must be spike or random"
    end if
    
    print*, "Initialisation done."
    
    do h  = 1,nsteps
    	call second_derivative(T,nx+1, ny+1 ,dx,dy, DD)
    	do i  = 1,nx-1 !makes sure boundaries at 0 and nx are constant 0
    		do j = 1,ny-1
                T(i,j) = T(i,j) + dt*k*DD(i,j)
            end do
    	end do
    end do
    if (Tinit == 'spike') then
        open(1,file='twoD_results_spike.dat', status='replace')
        do j=0,ny
            write(1,'(*(es13.5))') T(:,j)
        end do
        close(1)
    endif

    if (Tinit == 'random') then
        open(1,file='twoD_results_random.dat', status='replace')
        do j=0,ny
            write(1,'(*(es13.5))') T(:,j)
        end do
        close(1)
    endif
    deallocate(T)
    deallocate(DD)
end program twoDdiffusion
