subroutine relativeerror(x,setupinfo)
implicit none
include 'commondata.h'
include 'rbfconstypes.h'

type(setup), intent(in) :: setupinfo
real*8, dimension(setupinfo%blocknum,1), intent(in) :: x
real*8, allocatable, dimension(:,:) :: synthmodel
real*8 :: cerr1,cerr2,commissionerr
integer, dimension(setupinfo%blocknum) :: ocean
integer :: status, ind1,ind2


!-----------------------------------------------------------------------------------------------------------
!	Read the true EWH in interested area
!	synthmodel		d(3,blocknum)	the data in each rows are L(in degree), B(in degree), synthsis EWH(in mm)
!-----------------------------------------------------------------------------------------------------------
open(unit = 10, file = setupinfo%synthsetupfile,iostat = status, action = 'read', status = 'old')
if(status.eq.0) then
allocate(synthmodel(3,setupinfo%blocknum),stat = status)
if(status.eq.0) then
read(10,*) synthmodel
end if
end if
rewind(10)
close(10)


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

!-----------------------------------------------------------------------------------------------
!	ceer1		RMS of error	(error: the difference between signal and recoverd EWH)
!	ceer2		RMS of signal	(signal: the EWH of synthsis model)
!	ind2		number of blocks on land
!-----------------------------------------------------------------------------------------------
 cerr1 = 0.0d0
 cerr2 = 0.0d0
ind2 = 0
do ind1 = 1,setupinfo%blocknum
if(ocean(ind1).eq.1) then
ind2 = ind2 + 1

 cerr1 = cerr1 + (x(ind1,1)*1000 - synthmodel(3,ind1))**2
 cerr2 = cerr2 +  (synthmodel(3,ind1))**2
end if
end do

 cerr1 = sqrt(cerr1 / ind2)
 cerr2 = sqrt(cerr2 / ind2)
 commissionerr = cerr1/cerr2


!----------------------------------------------------------------------------------------------------
!	The RMS of error, RMS of signal and commission error are stored in file 'Out/commissionerr.txt'
!----------------------------------------------------------------------------------------------------
open(unit = 19, file = trim(setupinfo%outputdir)//'commissionerr.txt', iostat = status, action = 'write', status = 'replace')
if(status.eq.0) then
write(19,'(e25.15)') cerr1,cerr2,commissionerr
endif
rewind(19)
close(19)

end subroutine
