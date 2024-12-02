module derivatives
contains
	!calculates dT/dx
	function partialdx(T,dx, nx,ny)
	    implicit none
        integer, intent(in) :: nx, ny
        real, intent(in) :: T(nx,ny), dx
        
        integer :: i
        
        real, dimension(nx,ny):: partial
        real, dimension(nx,ny):: partialdx

        do i = 2,nx-1
            partial(i,:) = (T(i+1,:) - T(i-1,:))/dx/2.0
        end do      
        partial(1,:) = 0.0
        partial(nx,:) = 0.0
        !partial(:,1) = 0.0
        !partial(:,ny) =0.0
        partialdx = partial
        
    end function partialdx
    !calculates velocity given stream function
    function  vel(S,vy,nx,ny)
        implicit none
        real, dimension(:,:),intent(in)::S
        integer ::nx,ny,i,j
        real, intent(out)::vy(nx,ny)
        real, dimension(nx,ny) :: vel
        vel = 0
        vy = 0
        do i= 2, nx-1
            do j = 1, ny
                !vy = -dS/dx using centered finite difference approx
                vy(i,j) = (S(i-1,j) - S(i+1,j))*(ny-1)/2
            end do
        end do

        do i= 1, nx
            do j = 2, ny-1
                !vx = dS/dy using centered finite difference approx
                vel(i,j) = (S(i,j+1) - S(i, j-1))*(ny-1)/2
            end do
        end do        
    end function vel
    !calculates vgradT
    function vgrad(T, vx, vy, dx, dy,nx,ny)
        implicit none
        real, dimension(:,:), intent(in):: T,vx,vy
        real, intent(in):: dx, dy 
        real :: vgradTdx(nx,ny), vgradTdy(nx,ny)
        real, dimension(nx,ny):: vgrad
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
        vgrad = vgradTdx + vgradTdy
    end function vgrad
    
    function second_derivative (array, Nx,Ny, dx, dy)
        implicit none
        integer, intent(in) :: Nx,Ny
        integer :: i,j
        real, intent(in) :: dx, dy
        real, intent(in) :: array(Nx,Ny)
        real, dimension(Nx,Ny) :: second_derivative
        real,dimension(Nx,Ny) :: delsq
        

    	second_derivative(1,:) = 0.0
        second_derivative(:,1) = 0.0
    	second_derivative(nx,:) = 0.0
        second_derivative(:,ny) = 0.0
    
        do i = 2,Nx-1
            do j = 2, Ny-1
                second_derivative(i,j) = (array(i-1,j) + array(i+1,j) - 2*array(i,j))/(dx*dx) &
                 + (array(i,j-1) + array(i, j+1)-2*array(i,j))/(dy*dy) 
            end do    
        end do
    end function second_derivative
end module derivatives

module Poisson_solver
    contains
        
    real function iteration_2DPoisson(u,f,h,alpha) result(rms_residue)
        implicit none
        real ::  res
        real,intent(inout)::u(:,:)
        real,intent(in)::f(:,:),h,alpha
        integer :: i,j,nx,ny
        real :: dx,dy
        
        nx = size(u,1)
        ny = size(u,2)
        
        res = 0.0
        rms_residue = 0.0
        do i = 2,nx-1
            do j = 2,ny-1
                !Gauss-Seidel Slide 37
                res =  (u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4.0*u(i,j))/(h**2)-f(i,j)
                rms_residue = rms_residue + res*res
                u(i,j) = u(i,j) + alpha*res*h*h/4.0
            end do
        end do
        rms_residue = sqrt(rms_residue/nx/ny)

    end function iteration_2DPoisson
    
    subroutine residue_2DPoisson(u,f,h,res)
        implicit none
        real, intent(in) :: u(:,:), f(:,:),h
        real, intent(out) ::res(:,:)
        real, allocatable :: DD(:,:)
        integer :: nx,ny,i,j
        
        nx = size(u,1)
        ny = size(u,2)

        !Calculating residue as del^2 u - f
        do i=2,nx-1
            do j=2,ny-1
                res(i,j)= (u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j))/(h**2)-f(i,j)
            end do
        end do
        res(1,:)=0.0
        res(nx,:)=0.0
        res(:,1)=0.0
        res(:,ny)=0.0
    end subroutine residue_2DPoisson
    
    subroutine restrict(fine,coarse)
        implicit none
        real, intent(in)::fine(:,:)
        real, intent(out)::coarse(:,:)
        integer :: nxf, nyf, nxc, nyc, hx, hy,i,j
        !copies every other point in fine into coarse
        nxf = size(fine,1)
        nyf = size(fine,2)
        nxc = size(coarse,1)
        nyc = size(coarse,2)

        do i = 1,nxc
            do j = 1,nyc
                coarse(i,j) = fine(i*2-1,j*2-1)
            end do
        end do
    end subroutine restrict

    subroutine prolongate(coarse,fine)
        implicit none
        real, intent(in)::coarse(:,:)
        real, intent(out)::fine(:,:)
        integer :: nxf, nyf, nxc,nyc, i,j,hx,hy
        !copies coarse into every other point in fine
        !linear interpolation
        
        nxf = size(fine,1)
        nyf = size(fine,2)
        nxc = size(coarse,1)
        nyc = size(coarse,2)
        hx = (nxf-1)/(nxc-1) !hx = 2 in our case
        hy = (nyf-1)/(nyc-1) !hy = 2 in our case

        !assigning to every second entry (also called odd column, odd row case)
        do i=1,nxc
            do j=1,nyc
                fine(i*2-1,j*2-1)=coarse(i,j)
            end do
        end do
    
        !interpolating the entries inbetween
        ! even column index, odd row index
        do j=1,nyc*2-1,2
            do i=2,2*nxc-2,2
                fine(i,j) = (fine(i-1,j)+fine(i+1,j))/2
            end do
        end do
    
        ! odd column index, even row index
        do j=2,nyc*2-2,2
            do i=1,2*nxc-1,2
                fine(i,j) = (fine(i,j-1)+fine(i,j+1))/2
                !fine(i*2,2*ny-1) = (fine(i*2-1,2*ny-1)+fine(i*2+1,2*ny-1))/2
            end do
        end do
        !even column and row index
        do i=2,2*nxc-1,2
            do j=2,2*nyc-1,2
                fine(i,j)=(fine(i-1,j)+fine(i+1,j)+fine(i,j+1)+fine(i,j-1))/4
            end do
        end do
        
    end subroutine prolongate
    
    recursive function Vcycle_2DPoisson(u_f,rhs,h) result (resV)
        implicit none
        real resV
        real,intent(inout):: u_f(:,:)  ! arguments
        real,intent(in)   :: rhs(:,:),h
        integer         :: nx,ny,nxc,nyc, i  ! local variables
        real,allocatable:: res_c(:,:),corr_c(:,:),res_f(:,:),corr_f(:,:)
        real            :: alpha=0.7, res_rms

        nx=size(u_f,1); ny=size(u_f,2)  ! must be power of 2 plus 1
        if( nx-1/=2*((nx-1)/2) .or. ny-1/=2*((ny-1)/2) ) &
            stop 'ERROR:not a power of 2'
        nxc=1+(nx-1)/2; nyc=1+(ny-1)/2  ! coarse grid size

        if (min(nx,ny)>5) then  ! not the coarsest level

           allocate(res_f(nx,ny),corr_f(nx,ny), &
                corr_c(nxc,nyc),res_c(nxc,nyc))

           !---------- take 2 iterations on the fine grid--------------
           res_rms = iteration_2DPoisson(u_f,rhs,h,alpha) 
           res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)

           !---------- restrict the residue to the coarse grid --------
           call residue_2DPoisson(u_f,rhs,h,res_f) 
           call restrict(res_f,res_c)

           !---------- solve for the coarse grid correction -----------
           corr_c = 0.  
           res_rms = Vcycle_2DPoisson(corr_c,res_c,h*2) ! *RECURSIVE CALL*

           !---- prolongate (interpolate) the correction to the fine grid 
           call prolongate(corr_c,corr_f)

           !---------- correct the fine-grid solution -----------------
           u_f = u_f - corr_f  

           !---------- two more smoothing iterations on the fine grid---
           res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
           res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)

           deallocate(res_f,corr_f,res_c,corr_c)

        else  

           !----- coarsest level (ny=5): iterate to get 'exact' solution

           do i = 1,100
              res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
           end do

        end if

        resV = res_rms   ! returns the rms. residue

  end function Vcycle_2DPoisson

end module Poisson_solver

module grid_mod
    implicit none
    type grid
        integer::nx,ny
        real::h
        real, allocatable :: S(:,:), W(:,:), T(:,:), vx(:,:), vy(:,:)
    end type grid

end module grid_mod

!main program
program convection_simulation
    use derivatives
    use poisson_solver
    use grid_mod

    implicit none
    !real, allocatable :: T(:,:), S(:,:), W(:,:), vx(:,:), vy(:,:), rhs(:,:)
    type(grid):: grids
    real :: Pr = 0.1, Ra = 1e6, err = 1.e-3 ,a_adv = 0.4, a_dif = 0.15, total_time = 0.1, dt_dif, dt_adv, dt 
    real :: time = 0.0, res_rms,vmax = 0.0,vxmax = 0.0, vymax = 0.0, W_rms = 0.0, EPS = 1.e-4
    integer :: nx = 257, ny = 65, i,j
    character(len=10) :: Tinit = 'cosine'
    character(len=70) :: filename = 'lowPrandt_parameters.txt'   
    real :: PI = 3.14159265359
    !read input parameter
    namelist /inputs/ Pr,nx, ny, Ra, total_time, err, a_adv, a_dif, Tinit 

    if(command_argument_count()>0) &
        call get_command_argument(1,filename)
    print*, filename
    
    !read input paramter from file
    open(1, file = filename, status= 'old')
    read(1,inputs)
    close(1)     

    !initialise some variables
    !initialise grid
    grids%nx = nx
    grids%ny = ny
    call Initialise_grid(grids)
    !grids%h = 1.0/(grids%ny-1)
    dt_dif = a_dif*(grids%h**2)/MAX(Pr,1.0)
    print*, "Initialised grid"
    
    !open file to write time vs vxmax 
    open(unit=10, file='time_vmax.dat', status="replace")
    write(10,*) time, vmax
    !mainloop - timesteps
    print*, "Start do loop"
    do while (time <=  total_time) 
        
        
        !boundary conditions for W    
        grids%W(:,1)=1.
        grids%W(:,ny)=0.
        grids%W(1,:)=grids%T(2,:)
        grids%W(nx,:)=grids%T(nx-1,:)        
        !calculate rms of W
        W_rms = 0.0
        do concurrent (i = 1:grids%nx, j = 1:grids%ny)
            W_rms = W_rms + grids%W(i,j)**2
        end do
        W_rms = sqrt(W_rms/grids%nx/grids%ny)
        !if W_rms = 0, then W ==0, then it is not necessary to calculate S using Vcycle_2DPoisson
        !we then know that S = 0
        if(W_rms < EPS) then
            grids%S = 0
        else
            !poisson solve to get S from W
            res_rms = Vcycle_2DPoisson(grids%S,grids%W,grids%h)
            do while(res_rms >  err*W_rms)
                res_rms = Vcycle_2DPoisson(grids%S,grids%W,grids%h)

            !boundary conditions for S
            grids%S(1,:) = 0.0
            grids%S(nx,:) = 0.0
            grids%S(:,1) = 0.0
            grids%S(:,ny) = 0.0
            end do
        end if
         

        !calculating vecolicites from S
        grids%vx =  vel(grids%S,grids%vy,grids%nx,grids%ny)
        vxmax = MAXVAL(ABS(grids%vx)) 
        vymax = MAXVAL(ABS(grids%vy))
        vmax = MAX(vxmax,vymax)
        !calculate timestep from diffusion and from advection timestep
        !if both velocities are = 0, then vxmax = vymax = 0
        !then the timestep of the diffusion should be taken
        if(vxmax < EPS .or. vymax < EPS)then
            dt = dt_dif
        !if only vxmax == 0 then avoid dividing by vxmax
        !else if (vxmax < EPS) then
        !    dt_adv = a_adv*grids%h/vymax
        !    dt = MIN(dt_dif,dt_adv)
        !else if (vymax < EPS) then
        !    dt_adv = a_adv*grids%h/vxmax
        !    dt = MIN(dt_dif, dt_adv)                        
        else 
            dt_adv = a_adv*MIN(grids%h/vxmax, grids%h/vymax)
            dt = MIN(dt_dif, dt_adv)                        
        end if
        grids%T = grids%T + dt*(second_derivative(grids%T,grids%nx,grids%ny,grids%h,grids%h) &
                 - vgrad(grids%T,grids%vx,grids%vy,grids%h, grids%h,grids%nx,grids%ny))
        grids%W = grids%W + dt*(Pr*second_derivative(grids%W, grids%nx, grids%ny,grids%h,grids%h) &
                -vgrad(grids%W, grids%vx, grids%vy, grids%h, grids%h, grids%nx, grids%ny) &
                -Pr*Ra*partialdx(grids%T,grids%h,grids%nx,grids%ny))     
        
        time = time + dt
        write(10,*) time, vmax
        !boundary condition at bottom and top
        grids%T(:,1) = 1.0
        grids%T(:,ny) = 0.0
        !boundary condition at sides
        grids%T(nx,:) = grids%T(nx-1,:)
        grids%T(1,:) = grids%T(2,:)
        !boundary conditions for W    
        grids%W(:,1)=0.0
        grids%W(:,ny)=0.0
        grids%W(1,:)=0.0
        grids%W(nx,:)=0.0       
       end do
    close(10)
    
    print*, "End of do loop"
    
    open(unit= 10, file = 'T_output.dat', status = "replace")
    do i = 1, grids%ny
        write(10,*)(grids%T(j,i), j = 1,grids%nx)
    end do
    close(10)

    open(unit= 10, file = 'W_output.dat', status = "replace")
    do i = 1, grids%ny
        write(10,*)(grids%W(j,i), j = 1,grids%nx)
    end do
    close(10)
    open(unit= 10, file = 'S_output.dat', status = "replace")
    do i = 1, grids%ny
        write(10,*)(grids%S(j,i), j = 1,grids%nx)
    end do
    close(10)
    call deallocate_grid(grids)

    contains

    

    subroutine Initialise_grid(a)
        implicit none
        integer :: i,j
        type(grid),intent(inout)::a

        a%h=1./(a%ny-1.0)

        ! Allocate space for all grids
        allocate(a%T(a%nx,a%ny), a%W(a%nx,a%ny), a%S(a%nx,a%ny), a%vx(a%nx,a%ny), a%vy(a%nx,a%ny))
        ! Initialize all grids
        a%W=0.0
        a%S=0.0
        a%vx=0.0
        a%vy=0.0

        ! Initialize T random
        if (Tinit == 'random') then
            call random_number(a%T)

        ! Initialize T according to slide 20
        else if(Tinit == 'cosine') then
            do i=1,a%nx
                do j=1,a%ny
                    a%T(i,j)=0.5*(1+cos(3*pi*(i-1)/(a%nx-1)))
                    !print*, a%T(i,j)
                end do
            end do
        else
            print*, "Error, init has to be either random or cosine"
        end if


        ! Boundary Conditions for T
        a%T(:,1)=1.0
        a%T(:,a%ny)=0.0
        a%T(1,:)=a%T(2,:)
        a%T(a%nx,:)=a%T(a%nx-1,:)

    end subroutine Initialise_grid

    subroutine deallocate_grid(a)
        implicit none
        type(grid),intent(inout)::a
        deallocate(a%T)
        deallocate(a%W)
        deallocate(a%S)
        deallocate(a%vx)
        deallocate(a%vy)
    end subroutine deallocate_grid
end program convection_simulation
