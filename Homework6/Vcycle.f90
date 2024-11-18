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

program main
    use Poisson_solver
    implicit none
    integer :: nx,ny,i,j,n_iterations = 1
    character(len = 10)::init
    logical :: multigrid
    real :: alpha, res_rms, f_rms, TOL,h
    real, allocatable :: f(:,:), u(:,:)
    
    !reading in input
    namelist /in/ nx, ny, init,alpha,multigrid

    open(1, file = 'relaxIN.txt', status = 'old')
    read(1,in)
    close(1)
    allocate(f(nx,ny))
    allocate(u(nx,ny))
    
    h = 1.0/(ny-1)
    u = 0.0
    TOL = 0.00001
    f_rms = 0.0
    n_iterations = 0

    if (Init == 'spike') then
        f  = 0.0
        f(nx/2+1, ny/2+1) = 1.0/(h**2)
    else if (Init == 'random') then
        call random_number(f)
    else 
        print*, "Init has to be spike or random"
    end if

    ! Boundary Conditions for u and f
    u(1,:)=0
    u(nx,:)=0
    u(:,1)=0
    u(:,ny)=0
    f(1,:)=0
    f(nx,:)=0
    f(:,1)=0
    f(:,ny)=0


    f_rms = 0.0
    do i=1,nx
        do  j=1,ny
            f_rms = f_rms +  f(i,j)**2
        end do    
    end do
    f_rms = sqrt(f_rms/nx/ny)

    if (Multigrid .eqv. .false.) then
        res_rms = iteration_2DPoisson(u,f,h,alpha)
        do while(res_rms > TOL*f_rms)
            res_rms = iteration_2DPoisson(u,f,h,alpha)
            n_iterations = n_iterations + 1
        end do
   
    else
        res_rms = Vcycle_2DPoisson(u,f,h)        
        do while(res_rms > TOL*f_rms)
            res_rms = Vcycle_2DPoisson(u,f,h)
            n_iterations = n_iterations + 1
        end do
    end if

    print*, "number of iterations:", n_iterations
      ! Write results to two files
    open(1,file='u_vcycle_random.dat', status='replace')
        do j=1,ny
            write(1,'(*(es13.5))') u(:,j)
        end do
        close(1)

        open(1,file='f_vcycle_random.dat', status='replace')
        do j=1,ny
            write(1,'(*(es13.5))') f(:,j)
        end do
        close(1)

    deallocate(u)
    deallocate(f)

    print*, "Done with loop."
end program main
