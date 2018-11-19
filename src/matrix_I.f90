subroutine matrix_I(x,n,matrixx)
!---------------------------------------------------
!	matrixx = x * Id(n,n)
!	form a diagonal matrix with x in the diagonal
!
!	Given:
!		x		d
!		n		i
!
!	Returned:
!		matrixx		d(n,n)
!
!----------------------------------------------------

implicit none
integer, intent(in) :: n
real*8, intent(in) :: x
real*8, dimension(n,n), intent(out) :: matrixx
integer :: ind1, ind2

do ind1 = 1,n
do ind2 = 1,n
if(ind1.eq.ind2) then
matrixx(ind1,ind2) = x
else
matrixx(ind1,ind2) = 0.0d0
end if
end do
end do


end subroutine