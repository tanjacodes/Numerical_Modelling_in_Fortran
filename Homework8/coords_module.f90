

module coords
  use iso_fortran_env
  implicit none

  private
  public :: point, pointT
  
  type point
    real :: x,y,z
   contains
     procedure,private :: pointplus, pointminus, pointseparation, dotproduct, crossproduct,&
                         rightmult32, rightmult64, divide32, divide64
     procedure,private,pass(this) :: absvec32,absvec64,leftmult32,leftmult64
     generic,public :: operator(+) => pointplus
     generic,public :: operator(-) => pointminus
     generic,public :: operator(.distance.) => pointseparation
     generic,public :: assignment(=) => absvec32, absvec64
     generic,public :: operator(.dot.) => dotproduct
     generic,public :: operator(.cross.) => crossproduct
     generic,public :: operator(*) => leftmult32, leftmult64, rightmult32, rightmult64
     generic,public :: operator(/) => divide32, divide64
  end type point

  type, extends(point) :: pointT
    real :: temp
   contains 
    procedure,public :: Tplus, TMinus
  end type pointT


contains

  type(point) function pointplus(this,b)   ! for +
    class(point),intent(in):: this,b
    pointplus%x = this%x + b%x
    pointplus%y = this%y + b%y
    pointplus%z = this%z + b%z
  end function pointplus
  
  type(point) function pointminus(this,b)  !  for -
    class(point),intent(in):: this,b
    pointminus%x = this%x - b%x
    pointminus%y = this%y - b%y
    pointminus%z = this%z - b%z
  end function pointminus

  real function pointseparation(this,b) ! for .distance.
    class(point),intent(in):: this,b
    pointseparation = sqrt(  &
         (this%x-b%x)**2+(this%y-b%y)**2+(this%z-b%z)**2)
  end function pointseparation

  subroutine absvec32(a,this)  ! for = (distance
    implicit none
    real(REAL32),intent(out):: a     !  from origin)
    class(point),intent(in):: this
    a = sqrt(this%x**2+this%y**2+this%z**2)
  end subroutine absvec32
  
  subroutine absvec64(a,this)  ! for = (distance
    implicit none
    real(REAL64),intent(out):: a     !  from origin)
    class(point),intent(in):: this
    a = sqrt(this%x**2+this%y**2+this%z**2)
  end subroutine absvec64

  real function dotproduct(this,b)
    class(point), intent(in)::this,b
    dotproduct = this%x*b%x + this%y*b%y + this%z*b%z
  end function dotproduct

  type(point) function crossproduct(this,b)
    class(point),intent(in):: this,b
    !crossproduct = this
    crossproduct%x = this%y*b%z - this%z*b%y
    crossproduct%y = this%z*b%x - this%x*b%z
    crossproduct%z = this%x*b%y - this%y*b%x
  end function crossproduct

  type(point) function leftmult32(a,this)
    real(real32), intent(in) :: a
    class(point), intent(in):: this
    !leftmult = this
    leftmult32%x = a*this%x
    leftmult32%y = a*this%y
    leftmult32%z = a*this%z
  end function leftmult32
  type(point) function leftmult64(a,this)
    real(real64), intent(in) :: a
    class(point), intent(in):: this
    !leftmult = this
    leftmult64%x = a*this%x
    leftmult64%y = a*this%y
    leftmult64%z = a*this%z
  end function leftmult64

  type(point) function rightmult32(this,b)
    class(point), intent(in):: this
    real(real32), intent(in) :: b
    !rightmult = this
    rightmult32%x = this%x*b
    rightmult32%y = this%y*b
    rightmult32%z = this%z*b
  end function rightmult32

  type(point) function rightmult64(this,b)
    class(point), intent(in):: this
    real(real64), intent(in) :: b
    !rightmult = this
    rightmult64%x = this%x*b
    rightmult64%y = this%y*b
    rightmult64%z = this%z*b
  end function rightmult64

  type(point) function divide32(this,b)
    class(point), intent(in) :: this
    real(real32), intent(in) :: b
    !divide = this
    divide32%x = this%x/b
    divide32%y = this%y/b
    divide32%z = this%z/b
     
  end function divide32

  type(point) function divide64(this,b)
    class(point), intent(in) :: this
    real(real64), intent(in) :: b
    !divide = this
    divide64%x = this%x/b
    divide64%y = this%y/b
    divide64%z = this%z/b
     
  end function divide64
  
  real function Tplus(this,b)
    class(pointT), intent(in) :: this,b
    Tplus = this%temp + b%temp
  end function Tplus

  real function Tminus(this,b)
    class(pointT), intent(in) :: this,b
    Tminus = this%temp - b%temp
  end function Tminus

  !real function crossproduct(b,this)
  !end function crossproduct
end module coords

program main
  use coords
  use iso_fortran_env

  implicit none
  real(real32):: r4 
  real(real64):: r8
  type(point)::p1,p2,p3
  type(pointT):: pT1, pT2

  p1%x = 2.0
  p1%y = 3.0
  p1%z = 4.0
  p2%x = 1.0
  p2%y = 1.0
  p2%z = 0.0

  pT1%x = 0.0
  pT1%y = 0.0
  pT1%z = 0.0
  pT1%temp = 10.0
  pT2%x = 1.0
  pT2%y = 1.0
  pT2%z = 1.0
  pT2%temp = 5.0

  r4 = p1
  r8 = p1
  r4 = p1.dot.p2
  !print*, r4
  
  p3 = p1.cross.p2
  !print*, p3
  
  p3 = p1*r4
  !print*, p3
  
  p3 = p1*r8
  !print*, p3
  
  p3 = r4*p1
  !print*, p3
  
  p3 = r8*p1
  !print*, p3
  
  p3 = p1/r4
  !print*, p3
  
  p3 = p1/r8
  !print*, p3
  
  r4 = pT1%Tplus(pT2)
  !print*, r4
  
  r4 = pT1%Tminus(pT2)
  !print*, r4
end program main
