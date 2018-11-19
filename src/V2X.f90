subroutine V2X(setupinfo)
implicit none
include 'commondata.h'
include 'rbfconstypes.h'

type(setup), intent(in) :: setupinfo
integer*4 :: status, ind1, ind2, ind3, flagnum,num,blocknum
!----------------------------------------------------------------------------------------
!	block_lbrarea(blocknum, 5) 
!	NO.	L	B	r	area
!----------------------------------------------------------------------------------------
real*8, allocatable, dimension(:, :) :: block_lbrarea
!--------------------------------------------------------------------
!	mark the epoch at which the midpoint overfly the interested area
!--------------------------------------------------------------------
integer*4, allocatable, dimension(:) :: flag

real*8, allocatable, dimension(:,:) :: temp, midposbl,pab, blra,t,y,t1,y1
real*8, allocatable, dimension(:,:) :: btb, btl, x, x1, ctc, matrix_d, priorimodel, truemodel,invN
real*8, allocatable, dimension(:) :: temp1, cos_fai_a


PI = 4.0D0*ATAN(1.0D0)
num = setupinfo%n_obs		!number of the observations
blocknum = setupinfo%blocknum	!number of rbf

open(unit = 13, file = setupinfo%obssetupfile, status = 'old', action = 'read', iostat = status)
if(status == 0) then
 allocate(pab(7, num), stat = status)
 if(status == 0) then
   read(13, *) pab
   rewind(13)
   close(13)
 end if
end if

allocate(blra(3,num))

do ind1 = 1, num
    call XYZ2BLR(pab(1,ind1),pab(2,ind1),pab(3,ind1),blra(1,ind1),blra(2,ind1),blra(3,ind1))
end do
blra(1:2,:) = blra(1:2,:)/180.0d0*pi

allocate(midposbl(num,2))
midposbl = transpose(temp(4:5,:))
deallocate(temp)

allocate(flag(num))
call markmidpos(setupinfo,midposbl,flag)
flagnum = sum(flag)

!------------------------------------------------------------------------------
!	read the L, B coordinates of the midpoints of blocks in interested area
!	store the data in matrix block_lbrarea(blocknum,5) temporarily
!------------------------------------------------------------------------------
open(unit = 18, file = setupinfo%masconsblockfile, status = 'old', action = 'read', iostat = status)
if(status == 0) then
 allocate(temp(5, blocknum))
   read(18, *) temp
   rewind(18)
   close(18)
end if
allocate(block_lbrarea(blocknum,5))
allocate(cos_fai_a(blocknum))
block_lbrarea = transpose(temp)
deallocate(temp)
block_lbrarea(:,2:3) = block_lbrarea(:,2:3) / 180.0d0 * pi

write(*,'(A)') 'forming design matrix...'
allocate(t1(flagnum, blocknum))
allocate(y1(flagnum, 1))

ind1 = 0
do ind2 = 1, num
if(flag(ind2).eq.1) then
ind1 = ind1 + 1

        cos_fai_a = dot_product(blra(:,ind2), block_lbrarea(:,1:3))
!----------------------------------------
!	form t1 and y1
!----------------------------------------
call buildline(cos_fai_a,setupinfo,block_lbrarea,t1(ind1,:))
y1(ind1,1) = pab(7,ind2)
end if
end do

deallocate(pab)
deallocate(blra)
deallocate(blrb)
deallocate(midposbl)
deallocate(flag)

!------------------------------------------------------------------------
!	t * x = y
!	double differenciate t1 and y1 to get rid of the energy constant
!------------------------------------------------------------------------
allocate(t(flagnum-1, blocknum))
allocate(y(flagnum-1, 1))
do ind1 = 1,flagnum - 1
t(ind1,:) = t1(ind1+1,:) - t1(ind1,:)
y(ind1,1) = y1(ind1+1,1) - y1(ind1,1)
end do

deallocate(t1)
deallocate(y1)

!-------------------------------------------------------------------
!	btb = transpose(t) * t
!	btl = transpose(t) * y
!-------------------------------------------------------------------
allocate(btb(blocknum, blocknum))
allocate(btl(blocknum, 1))
btb = matmul(transpose(t),t)
btl = matmul(transpose(t),y)

deallocate(t)
deallocate(y)

!-------------------------------------------------------------------------------------------------------------------------
!	Read the a priori information  in interested area
!-------------------------------------------------------------------------------------------------------------------------
open(unit = 10, file = setupinfo%priorisetupfile,iostat = status, action = 'read', status = 'old')
if(status.eq.0) then
allocate(priorimodel(3,blocknum))
read(10,*) priorimodel
end if
rewind(10)
close(10)

!-------------------------------------------------------------------------------------------------------------------------
!	Read the synthsis values in interested area
!-------------------------------------------------------------------------------------------------------------------------
open(unit = 100, file = setupinfo%synthsetupfile,iostat = status, action = 'read', status = 'old')
if(status.eq.0) then
allocate(truemodel(3,blocknum))
read(100,*) truemodel
end if
rewind(100)
close(100)

!-----------------------------------------------------------
!	form constraint matrix
!-----------------------------------------------------------
write(*,'(A)') 'forming constraint matrix...'
allocate(ctc(blocknum,blocknum))
allocate(matrix_d(blocknum, 1))

if(setupinfo%typeconstmat .eq. 1) then
!-------------------------------------------------------------
! adding spatial constraint based on rowlands, 2005
!-------------------------------------------------------------
call scf_rowlands(setupinfo,ctc)
matrix_d = 0.0d0

elseif(setupinfo%typeconstmat .eq. 2) then
!------------------------------------------------------
! adding spatial constraint based on Ramillien, 2012
!------------------------------------------------------
call scf_moy_ramillien(setupinfo,ctc)
matrix_d = 0.0d0

endif

!---------------------------------------------
!	Calculate x(EWH)
!	x = (btb+ctc)^(-1)*(btl+matrix_d)
!---------------------------------------------
allocate(x(blocknum,1))
allocate(invN(blocknum,blocknum))

if(setupinfo%const_solution) then
invN = btb + ctc
btl = btl + matrix_d
else
invN = btb
endif

call matinv(invN,blocknum,ind1)
x = matmul(invN, btl)

deallocate(ctc)
deallocate(btb)
deallocate(btl)
deallocate(matrix_d)

!--------------------------------------------------------------------------------------
! calculate the posterior correlation coefficients of a particular mascon with others
! modification will be done later
!--------------------------------------------------------------------------------------
call correlation(setupinfo,invN, 110, block_lbrarea)
deallocate(block_lbrarea)
deallocate(invN)

allocate(x1(blocknum,1))
x1 = x
deallocate(x1)

!----------------------------------
!	calculate relative error
!----------------------------------
call relativeerror(x,setupinfo)


open(unit = 19, file = setupinfo%recoveredEWH, iostat = status, action = 'write', status = 'replace')
if(status.eq.0) then
do ind1 = 1, blocknum
write(19,20) priorimodel(1,ind1),priorimodel(2,ind1),x(ind1,1)
end do
endif
rewind(19)
close(19)

!--------------------------------------------------------------------------------------------------------
!	The differences between recovered data the apriori information are stored in file 'Out/filenameherr'
!--------------------------------------------------------------------------------------------------------
open(unit = 22, file = setupinfo%recoveredERR, iostat = status, action = 'write', status = 'replace')
if(status.eq.0) then
do ind1 = 1, blocknum
write(22,20) priorimodel(1,ind1),priorimodel(2,ind1),x(ind1,1)-truemodel(3,ind1)
end do
endif
rewind(22)
close(22)

deallocate(x)
deallocate(priorimodel)
deallocate(truemodel)

20 format(2e15.5, e25.15)

end subroutine
