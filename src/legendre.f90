subroutine legendre(x,n,p)
!----------------------------------------------------------------
! 
!	Compute the Legendre polynomials
!   
!
!	Given:
!		x		d			x=cos(fai)
!		n		i			order of Legendre Polynomials
!
!   Returned:
!        p		d(0:n)		Legendre Polynomials
!
!----------------------------------------------------------------
  implicit none
  real*8, intent(in) :: x
  integer,intent(in) :: n
  real*8, dimension(0:n),intent(out) :: p
  integer :: ind
  
  p = 0.0d0
  p(0) = 1.0d0
  p(1) = x
  p(2) = 1.5d0*x*x - 0.5d0
  p(3) = 2.5d0*x*x*x - 1.5d0*x

  do ind = 4,n
  p(ind) = ((2.0d0*ind-1.0d0)*x*p(ind-1) - (ind-1.0d0)*p(ind-2))/ind
  end do

  end subroutine
