
module polar_coords
implicit none

private
public :: polar
real :: PI = 3.14159265359

type polar
    real :: z, theta
   contains
     procedure, private, pass(this) :: polar_to_ord, ord_to_polar, mult, divide
     generic,public :: assignment(=) => polar_to_ord, ord_to_polar    
     generic,public :: operator(*) => mult
     generic,public :: operator(/) => divide
end type polar

contains

    subroutine ord_to_polar(this,x)
        complex, intent(in) :: x
        class(polar), intent(out) :: this
        real :: PI = 3.14159265359
        this%z = sqrt(x%re**2 + x%im**2)
        this%theta = atan(x%im, x%re)/2./PI*360.0
    end subroutine ord_to_polar

    subroutine polar_to_ord(x,this)
        class(polar), intent(in) ::this
        complex, intent(out) :: x
        x = cmplx(this%z*cos(this%theta*2*PI/360),this%z*sin(this%theta*2*PI/360))
    end subroutine polar_to_ord

    type(polar) function mult(this,p2)
        class(polar), intent(in) :: this,p2
        mult%z = this%z*p2%z
        !using mod to make sure mult%theta is a value between -180 and 180
        mult%theta = mod(this%theta + p2%theta + 180, 360.0) - 180.0
    end function mult

    type(polar) function divide(this,p2)
        class(polar), intent(in) :: this, p2
        divide%z = this%z/p2%z
        !using mod to make sure mult%theta is a value between -180 and 180
        divide%theta = mod(this%theta - p2%theta + 180, 360.0)- 180.0
    end function divide

end module polar_coords

program main
    use polar_coords
    implicit none
    type(polar):: p1,p2,p3
    complex :: c1 = cmplx(-1.0, 2.0), c2 = cmplx(2.0,-0.5)

    !test type conversion
    p1 = c1
    c2 = p1
    print*, p1
    print*, c2/c1

    !test multiply and divide
    p2 = c2
    p3 = p1*p2
    print*, p3
    p3 = p2/p1
    print*, p3
end program main
