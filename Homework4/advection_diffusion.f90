module derivatives
 !   implicit none
!    real, parameter :: PI  = 3.1415926536
contains
    !calculates velocity given stream function
    subroutine vel(S,vx,vy,nx,ny)
        implicit none
        real, dimension(:,:),intent(in)::S
        integer ::nx,ny,i,j
        real, intent(out)::vx(nx,ny),vy(nx,ny)
        vx = 0
        vy = 0
        do i= 2, nx-1
            do j = 1, ny
                !vy = -dS/dx using centered finite difference approx
                vy(i,j) = (S(i-1,j) - S(i+1,j))*(nx-1)/2
            end do
        end do

        do i= 1, nx
            do j = 2, ny-1
                !vx = dS/dy using centered finite difference approx
                vx(i,j) = (S(i,j+1) - S(i, j-1))*(ny-1)/2
            end do
        end do        
    end subroutine vel
    !calculates vgradT
    subroutine vgrad(T, vx, vy, dx, dy,nx,ny, vgradT)
        implicit none
        real, dimension(:,:), intent(in):: T,vx,vy
        real, intent(in):: dx, dy 
        real :: vgradTdx(nx,ny), vgradTdy(nx,ny)
        real,intent(out):: vgradT(nx,ny)
        real, parameter :: PI  = 3.1415926536
        integer, intent(in):: nx, ny
        integer :: i,j
         
        do i=1,nx
            do j = 1,ny
                !first calculating vgradTdx
                if(vx(i,j) > 0 .AND. i>1) then
                    vgradTdx(i,j) = vx(i,j)*(T(i,j)-T(i-1,j))/dx
                else if (vx(i,j) < 0 .AND. i<nx) then
                    vgradTdx(i,j) = vx(i,j)*(T(i+1,j)-T(i,j))/dx
                else 
                    vgradTdx(i,j) = 0 
                end if
                !then vgradTdy
                if(vy(i,j) > 0 .AND. j > 1) then
                    vgradTdy(i,j)= vy(i,j)*(T(i,j)-T(i,j-1))/dy
                else if (vy(i,j) < 0 .AND. j< ny) then
                    vgradTdy(i,j)= vy(i,j)*(T(i,j+1)-T(i,j))/dy
                else 
                    vgradTdy(i,j) = 0 
                end if
            end do
        end do

        !final result is adding the two subresults
        vgradT = vgradTdx + vgradTdy
    end subroutine vgrad
    
    
    subroutine second_derivative (array, Nx,Ny, dx, dy, primeprime)
        implicit none
        integer, intent(in) :: Nx,Ny
        integer :: i,j
        real, intent(in) :: dx, dy
        real, intent(in) :: array(Nx,Ny)
        real, intent(out) :: primeprime(Nx,Ny)

    	primeprime(1,:) = 0.0
        primeprime(:,1) = 0.0
    	primeprime(nx,:) = 0.0
        primeprime(:,ny) = 0.0
    
        do i = 2,Nx-1
            do j = 2, Ny-1
                primeprime(i,j) = (array(i-1,j) + array(i+1,j) - 2*array(i,j))/(dx*dx)
                primeprime(i,j) = primeprime(i,j) + (array(i,j-1) + array(i, j+1)-2*array(i,j))/(dy*dy) 
            end do    
        end do
    
    end subroutine second_derivative
end module derivatives


program advection
    use derivatives
    implicit none
    integer:: nx,ny, nsteps,i,j,h
    real:: dx,dy,B, k, total_time, a_adv, a_dif, dt, dt_dif, dt_adv, x_max, y_max, vx_max, vy_max
    real, allocatable :: T(:,:), DD(:,:),S(:,:),vgradT(:,:),vx(:,:), vy(:,:) 
    real, parameter :: PI  = 3.1415926536
    !Declaring namelist  
    namelist /inputs/ nx,ny,B,k,a_adv,a_dif,total_time
    
    open(1,file = 'AdvDif_parameters.dat', status= 'old')
    read(1,inputs)
    close(1)
	print*, "Reading inputs done"
    
    dx = 1.0/(nx-1.0)
    dy = 1.0 /(ny-1.0)

    print*, "calculated dx, dy"


    print*, "calculated nsteps"
    allocate(T(nx,ny))
    T = 0.0
    allocate(DD(nx, ny))
    DD = 0.0
    allocate(S(nx,ny))
    S = 0.0
    allocate(vgradT(nx,ny))
    vgradT = 0.0
    allocate(vx(nx,ny))
    vx = 0.0
    allocate(vy(nx,ny))
    vy = 0.0
    
    print*, "Allocated all arrays"
    !random initialisation of T (because it shouldn't matter for the result)
    call random_number(T)
    T(:,ny) = 0.0
    T(:,1) = 1    
    !We are using a 1x1 grid, therefore max values for x and y are 1
    x_max = 1.0
    y_max = 1.0
    !initialise S using the convention that S =  B*sin(pi*x/x_max)*sin(pi*y/y_max)
    do i = 1,nx
        do j = 1,ny
            S(i,j) = B*sin(PI*dx*(i-1)/x_max)*sin(PI*dy*(j-1)/y_max)
        end do
    end do
    !calculating velocities
    call vel(S,vx,vy,nx,ny)

    !get max value of vx, and vy
    vx_max = MAXVAL(vx)
    vy_max = MAXVAL(vy)

    !calculate timestep from diffusion and from advection timestep
    dt_dif = a_dif*MIN(dx,dy)**2/k
    dt_adv = a_adv*MIN(dx/vx_max, dy/vy_max)
    dt = MIN(dt_dif, dt_adv)
    !total number of steps
    nsteps = total_time/dt
    print*, "Initialisation done."
    !print*, "vxmax", vx_max
    !print*, "vymax", vy_max
    !print*, "adif", a_dif
    !print*, "aadv", a_adv
    !print*, "t_dif", dt_dif
    !print*, "t_adv", dt_adv
    !print*, "nsteps", nsteps

    !main loop with timesteps
    do h  = 1,nsteps
        !boundary condition for top and bottom
        T(:,ny) = 0.0
        T(:,1) = 1    
        !Derivative condition for sides
        T(1,:) = T(2,:)
        T(nx,:) = T(nx-1,:) 
    	call second_derivative(T,nx, ny ,dx,dy, DD)
        call vgrad(T, vx, vy, dx, dy,nx,ny,vgradT)
        T = T + dt*(k*DD-vgradT)
    end do
    
    print*, "Main loop done."

    !writing results and stream function to output files
    open(1,file='adv_diff_result.dat', status='replace')
    do j=1,ny
        write(1,'(*(es13.5))') T(:,j)
    end do
    close(1)

    open(1,file='stream_function.dat', status='replace')
    do j=1,ny
        write(1,'(*(es13.5))') S(:,j)
    end do
    close(1)
    print*, "writing into file done"
    
    !Deallocating all arrays
    deallocate(T)
    deallocate(DD)
    deallocate(S)
    deallocate(vgradT)
    deallocate(vx)
    deallocate(vy)
end program advection
