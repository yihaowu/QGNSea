subroutine SCF_Rowlands(setupinfo,bctwbc)
implicit none
include 'commondata.h'
include 'rbfconstypes.h'
!--------------------------------------------------------------------------------------
!    
!	Calculate the spatial constraint functions
!
!	Given: 
!		setupinfo			type(setup)
!		setupinfo%regpar_alpha(s)	d		scale factor, in most cases between 1.0e-7 to 1.0e-2
!		setupinfo%collength(d)		d		correction distance, in km, in most cases between 100 to 600
!		setupinfo%masconsblockfile	c*30		filename of the L, B coordinates of block midpoints in interested area
!
!	Returned:
!		bctwbc		d(blocknum,blocknum)		regularizer matrix
!
!	Inputfile:
!		setupinfo%masconsblockfile
!
!
!	Notes:
!	1) the regularizer matrix assists the solutions on reducing striping error,
!		without the regularizer matrix, x = N^(-1)*W
!		with the matrix, x = (N + bctwbc)^(-1)*W
!
!	2) s and d should be selected manually
!-------------------------------------------------------------------------------------

type(setup), intent(in) :: setupinfo
real*8, dimension(setupinfo%blocknum,setupinfo%blocknum), intent(out) :: bctwbc
real*8 :: s,d
real*8 :: dij,theta,sij
integer :: status,blocknum,i,j,i1,j1
!---------------------------------------------------------------------------------
!	bc			d(blocknum*(blocknum - 1)/2,blocknum)	
!	coefficient matrix of spatial contraint functions
!
!	bc1			d(blocknum*(blocknum - 1)/2,blocknum)
!	we use this matrix just for computation efficiency
!
!	blockxyz	d(blocknum,3)	
!	the XYZ coordinates of the midpoint of blocks in interested area
!
!	wc			d(blocknum*(blocknum - 1)/2)			
!	the weight matrix of spatial contraint functions is a diagonal matrix
!	with the diagonal values of wc
!---------------------------------------------------------------------------------
real*8, allocatable, dimension(:,:) :: temp, bc, bc1, blockxyz
real*8, allocatable, dimension(:) :: wc

blocknum = setupinfo%blocknum
s = setupinfo%regpar_alpha
d = setupinfo%collength
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
do i1 = 1,blocknum
call BLR2XYZ(temp(3,i1),temp(2,i1),temp(4,i1),blockxyz(i1,1),blockxyz(i1,2),blockxyz(i1,3))
end do
deallocate(temp)

!--------------------------------------------------------------------------------
!	Constraint equation: for any pair of blocks i and i1, 0 = x(i1) - x(i)
!	x stands for the unknown parameters
!	There are up to blocknum*(blocknum - 1)/2 pairs(constraint equations).
!	the partical derivative of the jth constraint equation is:
!	bc(j,i1) = 1, bc(j,i) = -1, 0 for others
!	The constraint equations can be sumed as: 0 = bc * x 	
!	With fixed s and d, weight of the jth constraint equation wc(j) is determined by 
!	the spherical distance between block i and i1.
!--------------------------------------------------------------------------------
allocate(bc(blocknum*(blocknum - 1)/2,blocknum))
allocate(wc(blocknum*(blocknum - 1)/2))
j = 0
do i = 1, blocknum - 1
do j1 = 1, blocknum - i
j = j + 1
if(j > (blocknum*(blocknum - 1)/2)) then
write(*,*) 'exceed error! i = ',i
exit
endif
i1 = i + j1
bc(j,i) = -1
bc(j,i1) = 1

!---------------------------------------------------
!	space distance between i1 and i, in meters
!---------------------------------------------------
dij = sqrt((blockxyz(i1,1)-blockxyz(i,1))**2 + &
           (blockxyz(i1,2)-blockxyz(i,2))**2 + &
		   (blockxyz(i1,3)-blockxyz(i,3))**2)

!----------------------------------------------------
!	angular distance between i1 and i, in radius
!----------------------------------------------------
theta = 2.0d0 * asin(dij/(2.0d0*a))

!----------------------------------------------------
!	spherical distance between i1 and i, in km
!----------------------------------------------------
sij = a * theta * 1.0e-3

!----------------------------------------------------------
!	weight of the jth constraint funtion: 0 = H(i1) - H(i)
!----------------------------------------------------------
wc(j) = s * exp(2.0d0 - sij/d)
end do
end do


!-----------------------------------------------------------
!	bctwbc = transpose(bc) * p * bc
!	p is a diagonal matrix with the diagonal values of wc
!	bc1 = p * bc
!	so bctwbc = transpose(bc1) * bc
!-----------------------------------------------------------
allocate(bc1(blocknum*(blocknum - 1)/2,blocknum))

do j = 1,blocknum*(blocknum-1)/2
bc1(j,:) = bc(j,:) * wc(j)
end do

bctwbc = matmul(transpose(bc1),bc)
!call MATMULT ( transpose(bc1),bc, blocknum,blocknum*(blocknum - 1)/2,blocknum,bctwbc )

deallocate(bc)
deallocate(wc)
deallocate(bc1)
deallocate(blockxyz)

end subroutine
