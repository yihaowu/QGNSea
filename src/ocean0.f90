subroutine ocean0(x1,x,blocknum,setupinfo)
!------------------------------------------------------------------------------------
!
!	Setting the EWH of blocks in the ocean to be zero 
!
!	Considering the severe signal leakage near seas and the EWH of water variation
!	in the ocean ought to approximate zero, we simply set the EWH of blocks in the ocean to be zero 
!
!	Given:
!		x1	d	EWH recovered from orbit data with the spatial constraint SCF_han
!
!	Returned:
!		x	d	EWH, after setting the value of blocks in the ocean to be zero 
!
!	Inputfile: 
!		'In/ocean.txt'	1 stands for continent, 0 stands for ocean, 
!						I write this file manually
!
!-----------------------------------------------------------------------------------------------------
implicit none
include 'commondata.h'
include 'rbfconstypes.h'

integer, intent(in) :: blocknum
real*8, dimension(blocknum,1), intent(in) :: x1
real*8, dimension(blocknum,1), intent(out) :: x
type(setup), intent(in) :: setupinfo
integer, dimension(blocknum) :: ocean
integer :: status, ind

!---------------------------------------------------------------------
!	read file, and store the data in ocean
!	ocean	d(blocknum)		0: oceans	1: continents
!---------------------------------------------------------------------
open(unit = 18, file = setupinfo%oceanfile, status = 'old', action = 'read', iostat = status)
if(status == 0) then
   read(18, *) ocean
   rewind(18)
   close(18)
 end if

 do ind = 1,blocknum

 if(ocean(ind).eq.1) then
 x(ind,1) = x1(ind,1)
 else
 x(ind,1) = 0.0d0
 end if

 end do
end subroutine