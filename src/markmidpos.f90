subroutine markmidpos(setupinfo,midpos, flag)
!----------------------------------------------------------------
!
!	Mark the  midpoint within the interested area
!
!
!   Given:
!		midpos		d(n_obs, 2)		B,L coordinates of the midpoint of the twin satellites in CTS
!
!   Returned:
!		flag		i(n_obs)			flag(i) = 0 / 1
!									if flag(i) = 1, at epoch i, the midpoint overfly
!                                   the interested area, otherwise the vice
!
!	Notes:
!		The boundary of the interested area is
!
!			ooooooooooooooooooooo*(lat2+border, lon2+border)
!			oooooooooooooooooooooo		
!			oo!!!!!!!!!!!!!!!!!!oo
!			oo!!!!!!!!!!!!!!!!!!oo
!			oo!!!!!!!!!!!!!!!!!!oo
!			oo!!!!!!!!!!!!!!!!!!oo
!			oo!!!!!!!!!!!!!!!!!!oo
!			oo!!!!!!!!!!!!!!!!!!oo
!			oooooooooooooooooooooo
!			*ooooooooooooooooooooo	
!(lat1-border, lon1-border)
!
!	The declarations of lat1,lat2,lon1,lon2 and border are in module localdata(commondata.f90)
!
!-----------------------------------------------------------------
implicit none
include 'commondata.h'
include 'rbfconstypes.h'

type(setup) :: setupinfo
real*8, dimension(setupinfo%n_obs, 2), intent(in) :: midpos
integer, dimension(setupinfo%n_obs), intent(inout) :: flag
integer :: ind
real*8 :: lon1, lon2, lat1, lat2,border

lon1 = setupinfo%regions(1)%lon1
lon2 = setupinfo%regions(1)%lon2
lat1 = setupinfo%regions(1)%lat1
lat2 = setupinfo%regions(1)%lat2
border = setupinfo%border

flag = 0
do ind = 1, setupinfo%n_obs
if((midpos(ind, 2) >= lon1 - setupinfo%border).and.(midpos(ind, 2) <= lon2 + border) &
.and.(midpos(ind, 1) >= lat1 - setupinfo%border).and.(midpos(ind, 1) <= lat2 + border)) then

flag(ind) = 1

end if
end do

end subroutine