subroutine buildline(cos_fai_a,numval,setupinfo,blockinfo,linet)
!-------------------------------------------------------------------------------------------------------------
!	form a single line of coefficient matrix t1
!	t connects observations with local mascon parameters(method = 1) or     
! 	radial basis function parameters(method = 2,...)
!
!
!	Notes:
!	1) typeRBFs		1		mascon parameters
!					2		......
!------------------------------------------------------------------------------------------------------------
implicit none
include 'rbfconstypes.h'
include 'commondata.h'

type(setup), intent(in) :: setupinfo
real*8, dimension(setupinfo%blocknum,5) :: blockinfo
real*8, dimension(1,setupinfo%blocknum), intent(out) :: linet

integer :: j, ind, status, typeRBFs, numval
integer, parameter :: type_operator = 2
real*8 :: cos_fai_a(numval)
!-------------------------------------------------------------------------
!	build linet which connect observations with local rbf parameters
!-------------------------------------------------------------------------
if (typeRBFs .eq. 2) then
  
call SBF_analytical(wavevec,cos_fai_a,blockinfo(:,1:3),type_operator)
!-----------------
!	form linet
!-----------------
linet(1,:) = wavevec

end if

end subroutine


! =========================================================================================
subroutine SBF_analytical(sbf,cost,PLRrbf,type_operator)
! =========================================================================================
!
!     Defines a sbf (i.e.,Poisson wavelets (PW)) in analytical form and  applies some differential operator to the sbf;
!
!     Input
!     =====
!     type_operator  type of the differential operator [integer]
!                   =1       gravity disturbance operator in spherical approx
!                   =2       gravity anomaly operator in spherical approx
!                   =3 or 4  Identity (i.e. kernel is computed)
!     cost          cosine of the spherical distance between computation point and
!                   SBF [real*8]
!     par           order of the SBF for type_sbf = 1,2 [integer]
!
!     Output
!     ======
!     sbf           differential operator applied to a normalized SBF.
!
! ====================================================================================
    implicit none
    include 'rbfconstants.h'
    integer*4 :: l,i,numval,idx
    real*8 :: r,rb,sbf,lambda,distance
    real*8 :: b0,b1,b_nm1,b_nm2
    real*8 :: db0,db1,db_nm1,db_nm2
    real*8 :: term1,term2,term3
    real*8 :: r_SBF
    real*8 :: SBF_analytical_kernel
    real*8 :: PLRrbf(numval,3)
    real*8 :: cost(numval)
    integer,parameter :: par = 3

    real*8, allocatable :: P(:),c(:)
    real*8, allocatable :: b_vector(:),db_vector(:)
    real*8, allocatable :: beta(:),rbr(:)

    interface
        RECURSIVE REAL*8 FUNCTION binom(par,k) RESULT (b)
        INTEGER*4 :: par,k
        END FUNCTION binom
    end interface

    allocate(b_vector(0:par+1),db_vector(0:par+1))
    allocate(beta(1:par)
        
        beta(1) = 3.0d0
        beta(2) = 17.0d0
        beta(3) = 13.0d0

        do idx=1,numval
        
            r_SBF = PLRrbf(idx,3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            distance = dsqrt(r*r+r_SBF*r_SBF-2.0d0*r*r_SBF*cost(idx))
            b0 = 1.0d0/distance
            b1 = -((r_SBF - r*cost(idx))/(distance**3.0d0))
            b_nm1 = b1
            b_nm2 = b0

            ! recursive computation of Poisson wavelet of order n (par)

            b_vector(1) = b1
            b_vector(0) = b0

            do i=2,par+1
                b_vector(i) = (dble(2*i-1))*distance*b1*b_nm1 - ((dble(i-1))**2.0d0)*b_nm2*(b0**2.0d0)
                b_nm2 = b_nm1
                b_nm1 = b_vector(i)
            end do

            ! recursive computation of radial derivative of Poisson wavelet of order n (par)

            db0 = -((r - r_SBF*cost(idx))/(distance**3.0d0))
            db1 = (cost(idx)/(distance**3.0d0)) - (3.0d0/(distance**5.0d0))*((r*cost(idx) - r_SBF)*(r - r_SBF*cost(idx)))
            db_vector(1) = db1
            db_vector(0) = db0
            
            b_nm2 = b0
            b_nm1 = b1
            db_nm2 = db0
            db_nm1 = db1

            do i=2,par+1
                term1 = dble(2*i-1)*((r - r_SBF*cost(idx))/distance)*b1*b_vector(i-1)
                term2 = dble(2*i-1)*distance*(db1*b_nm1 + b1*db_nm1)
                term3 = -((dble(i-1))**2.0d0)*(2.0d0*b0*db0*b_nm2 + b0*b0*db_nm2)
                db_vector(i) = term1 + term2 + term3
                db_nm2 = db_nm1
                db_nm1 = db_vector(i)
                b_nm2 = b_nm1
                b_nm1 = b_vector(i)
            end do
            
            SBF_analytical_kernel = 0.0d0

                if (type_operator.eq.1) then
                    SBF_analytical_kernel = -2.0d0*(r_SBF**(dble(par+1)))*db_vector(par+1)
                    do i=1,par
                        SBF_analytical_kernel = SBF_analytical_kernel - beta(i)*(r_SBF**(i))*db_vector(i)
                    end do
                else if (type_operator.eq.2) then
                    SBF_analytical_kernel = -2.0d0*(r_SBF**(dble(par+1)))*(db_vector(par+1) + (2.0d0/r)*b_vector(par+1))
                    do i=1,par
                        SBF_analytical_kernel = SBF_analytical_kernel - beta(i)*(r_SBF**(i))*(db_vector(i) + (2.0d0/r)*b_vector(i))
                    end do
                else if ((type_operator.eq.3).or.(type_operator.eq.4)) then
                    SBF_analytical_kernel = 2.0d0*(r_SBF**(dble(par+1)))*b_vector(par+1)
                    do i=1,par
                        SBF_analytical_kernel = SBF_analytical_kernel + beta(i)*(r_SBF**(i))*b_vector(i)
                    end do
                end if
            SBF_analytical_kernel = SBF_analytical_kernel
        end do        
    deallocate(P,c,b_vector,db_vector,beta,rbr)
end subroutine SBF_analytical
