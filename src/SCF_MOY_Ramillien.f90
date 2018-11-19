subroutine SCF_MOY_Ramillien(setupinfo,ctc)
!--------------------------------------------------------------------------------------
!    
!	Calculate the uniform weighting(MOY) spatial constraint functions
!
!	References: 
!		Ramillien,2012
!
!	Given: 
!		setupinfo				type(setup)
!		setupinfo%regpar_alpha(alpha)		d	regularization parameter, in most cases between 1.0e-7 to 1.0e-2
!		setupinfo%collength(fai0)		d	maximum radius of correlation, in km, in most cases between 200 to 1800
!		setupinfo%masconsblockfile		c*30	filename of the L, B coordinates of block midpoints in interested area
!
!	Returned:
!		ctc			d(blocknum,blocknum)		regularizer matrix
!
!	Called:
!		matrix_I	form a diagonal matrix
!
!	Inputfile:
!		setupinfo%masconsblockfile
!
!
!
!	Notes:
!	1) the regularizer matrix assists the solutions on reducing striping error,
!		without the regularizer matrix, x = N^(-1)*W
!		with the matrix, x = (N + ctc)^(-1)*W
!
!	2) alpha and fai0 should be selected manually
!
!
!-------------------------------------------------------------------------------------
implicit none
include 'commondata.h'
include 'rbfconstypes.h'

type(setup), intent(in) :: setupinfo
real*8, dimension(setupinfo%blocknum,setupinfo%blocknum), intent(out) :: ctc
!------------------------------------------------------------------------
!	blockxyz	d(blocknum,3)	
!	the XYZ coordinates of the midpoint of blocks in interested area
!
!	diag_alpha = alpha * Id
!	We use this matrix just for the convenience of computation
!
!	fai_ij		d(blocknum,blocknum)
!	spherical distance between blocks, in km
!
!	Bij			d(blocknum, blocknum)
!	x = Bij * x
!
!	C			d(blocknum, blocknum)	
!	coefficient matrix of spatial contraint functions, 0 = C * x
!-------------------------------------------------------------------------

real*8 :: alpha,fai0
real*8, allocatable, dimension(:,:) :: temp, blockxyz,diag_alpha
real*8, allocatable, dimension(:,:) :: fai_ij, Bij, C
real*8, allocatable, dimension(:) :: l_fai0
real*8 :: dij, theta
integer :: status, blocknum, ind1, ind2
real*8, allocatable, dimension(:,:) :: bc, bc1
real*8, allocatable, dimension(:) :: wc

blocknum = setupinfo%blocknum
alpha = setupinfo%regpar_alpha
fai0 = setupinfo%collength
!----------------------------------------------------------------------------
!	read the L,B,r coordinates of the midpoint of blocks in interested area
!	and calculate the XYZ coordinates
!----------------------------------------------------------------------------
open(unit = 17, file = setupinfo%masconsblockfile, status = 'old', action = 'read', iostat = status)
if(status == 0) then
 allocate(temp(5, blocknum), stat = status)
 if(status == 0) then
   read(17, *) temp
   rewind(17)
   close(17)
 end if
end if
allocate(blockxyz(blocknum,3))
do ind1 = 1,blocknum
call BLR2XYZ(temp(3,ind1),temp(2,ind1),temp(4,ind1),blockxyz(ind1,1),blockxyz(ind1,2),blockxyz(ind1,3))
end do
deallocate(temp)

!---------------------------------------------------------------------
!	Calculate the spherical distance between any pair of blocks
!	the results are stored in matrix fai_ij(blocknum,blocknum)
!	the diagonal values of fai_ij are zeros
!---------------------------------------------------------------------
allocate(fai_ij(blocknum,blocknum),stat = status)
if(status.eq.0) then
fai_ij = 0.0d0
do ind1 = 1,blocknum-1
do ind2 = ind1+1,blocknum

!---------------------------------------------------------------------
!	space distance between block ind1 and ind2, in meters
!---------------------------------------------------------------------
dij = sqrt((blockxyz(ind1,1)-blockxyz(ind2,1))**2 + &
           (blockxyz(ind1,2)-blockxyz(ind2,2))**2 + &
		   (blockxyz(ind1,3)-blockxyz(ind2,3))**2)

!---------------------------------------------------------------------
!	angular distance between block ind1 and ind2, in radius
!---------------------------------------------------------------------
theta = 2.0d0 * asin(dij/(2.0d0*a))

!---------------------------------------------------------------------
!	spherical distance between block ind1 and ind2, in km
!---------------------------------------------------------------------
fai_ij(ind1,ind2) = a * theta * 1.0d-3
fai_ij(ind2,ind1) = fai_ij(ind1,ind2)

end do
end do
end if

deallocate(blockxyz)

!----------------------------------------------------------------------
!	Calculate l_fai0
!	l_fai0(ind1) is the number of blocks inside the geographical disc 
!	of radius fai0 and of which the block ind1 is the center
!----------------------------------------------------------------------
allocate(l_fai0(blocknum),stat = status)
if(status.eq.0) then
l_fai0 = 0          !initialization
do ind1 = 1,blocknum
do ind2 = 1,blocknum
if(ind2.ne.ind1) then
if(fai_ij(ind1,ind2) < fai0) then
l_fai0(ind1) = l_fai0(ind1) + 1
end if
end if
end do
end do
end if

!--------------------------------------------------------------------------
!	calculate Bij
!	references: Ramillien, 2012, sec.3.4
!--------------------------------------------------------------------------
allocate(Bij(blocknum,blocknum),stat = status)
if(status.eq.0) then
Bij = 0               !initialization
do ind1 = 1,blocknum
do ind2 = 1,blocknum

if(ind2 .ne. ind1) then
if(fai_ij(ind1,ind2) < fai0) then
Bij(ind1,ind2) = 1.0d0/l_fai0(ind1)
else
Bij(ind1,ind2) = 0
end if
end if

end do
end do
endif

!----------------------------------------------------------------
!	Calculate matrix C, C = Id - Bij
!	C is a diagonally dominant matrix
!----------------------------------------------------------------

allocate(c(blocknum,blocknum),stat = status)
if(status.eq.0) then
do ind1 = 1,blocknum
do ind2 = 1,blocknum
if(ind1.eq.ind2) then
c(ind1,ind2) = 1.0d0
else
c(ind1,ind2) = -bij(ind1,ind2)
end if
end do
end do
endif

!----------------------------------------------------------------
!	diag_alpha = alpha * Id
!	We use this matrix just for the convenience of computation
!	Id denotes identity matrix
!----------------------------------------------------------------
allocate(diag_alpha(blocknum,blocknum))
call matrix_I(alpha,blocknum,diag_alpha)

deallocate(fai_ij)
deallocate(l_fai0)
deallocate(bij)

!-----------------------------------------
!	ctc = alpha * transpose(c)*c
!	diag_alpha = alpha * Id
!	so ctc = transpose(c) * c * diag_alpha
!-----------------------------------------
ctc = matmul(transpose(c),c)
ctc = matmul(ctc,diag_alpha)

deallocate(diag_alpha)
deallocate(c)

end subroutine
