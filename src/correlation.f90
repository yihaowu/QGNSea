subroutine correlation(setupinfo, m_covariance, ind, block_lbrarea)
implicit none
include 'rbfconstypes.h'
include 'commondata.h'
type(setup), intent(in) :: setupinfo
real*8, dimension(setupinfo%blocknum, setupinfo%blocknum), intent(in) :: m_covariance
integer, intent(in) :: ind
real*8, dimension(setupinfo%blocknum,5), intent(in) :: block_lbrarea
real*8, dimension(setupinfo%blocknum) :: m_correlation

real*8 :: sqrtvarind
integer :: n, status

pi = 4.0d0 * atan(1.0d0)
sqrtvarind = sqrt(m_covariance(ind, ind))

!-------------------------------------------------------------------
!	m_correlation	correlation coefficient of mascon n with others
!-------------------------------------------------------------------
do n = 1, setupinfo%blocknum
m_correlation(n) = m_covariance(n,ind) / sqrt(m_covariance(n, n)) / sqrtvarind
end do

!---------------------------------------------------------------------------
!	the input covariance matrix is stored in file 'Out/covariance.txt'
!---------------------------------------------------------------------------
open(unit = 23, file = trim(setupinfo%outputdir)//'covariance.txt', iostat = status, action = 'write', status = 'replace')
if(status.eq.0) then
write(23,'(e25.15)') m_covariance
endif
rewind(23)
close(23)

!------------------------------------------------------------------------------------------------------
!	the correlation coefficient of mascon ind with other mascons are stored in file 'Out/plotcor.txt'
!------------------------------------------------------------------------------------------------------
open(unit = 18, file = trim(setupinfo%outputdir)//'plotcor.txt', iostat = status, action = 'write', status = 'replace')
if(status.eq.0) then
do n = 1, setupinfo%blocknum
write(18, 20) block_lbrarea(n,2) * 180.0d0 / pi, block_lbrarea(n,3) * 180.0d0 / pi, m_correlation(n)
end do
endif
rewind(18)
close(18)

20 format(2e15.5, e25.15)

end subroutine