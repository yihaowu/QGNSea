! ------------------------------------------------------------------------------
!
!                           SUBROUTINE MATMULT
!
!  this subroutine multiplies two matricies together.
!
!  Inputs          Description                    
!    Mat1        - Matrix number 1
!    Mat2        - Matrix number 2
!    Mat1r       - Matrix number 1 rows
!    Mat1c       - Matrix number 1 columns
!    Mat2c       - Matrix number 2 columns
!
!  OutPuts       :
!    Mat3        - Matrix result of Mat1 * Mat2 of size mat1r x mat2c
!
!  Locals        :
!    Row         - Row Index
!    Col         - Column Index
!    ktr         - Index
!
!  Coupling      :
!
! ------------------------------------------------------------------------------  

      SUBROUTINE MATMULT ( Mat1,Mat2, Mat1r,Mat1c,Mat2c,Mat3 )
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c,Mat2c
        REAL*8 Mat1(Mat1r,Mat1c),Mat2(Mat1c,Mat2c),Mat3(Mat1r,Mat2c)
! -----------------------------  Locals  ------------------------------
        INTEGER Row,Col,ktr

!         ! --------------------  Implementation   ----------------------
        DO Row=1, Mat1r
            DO Col= 1, Mat2c
                Mat3(Row,Col) = 0.0D0
                DO ktr= 1, Mat1c
                    Mat3(Row,Col)= Mat3(Row,Col)+ Mat1(Row,ktr)* Mat2(ktr,Col)
                  ENDDO
              ENDDO

          ENDDO
      RETURN
      END
!.
! ------------------------------------------------------------------------------
!
!                           SUBROUTINE MATADD
!
!  this subroutine adds two matricies together.
!
!  Inputs          Description                    
!    Mat1        - Matrix number 1
!    Mat2        - Matrix number 2
!    Mat1r       - Matrix number 1 rows
!    Mat1c       - Matrix number 1 columns
!
!  OutPuts       :
!    Mat3        - Matrix result of Mat1 + Mat2
!                    of size mat1r x mat1c
!
!  Locals        :
!    Row         - Row Index
!    Col         - Column Index
!
!  Coupling      :
!
! ------------------------------------------------------------------------------  

      SUBROUTINE MATADD ( Mat1,Mat2,Mat1r,Mat1c,Mat3 )
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c
        REAL*8 Mat1(Mat1r,Mat1c),Mat2(Mat1r,Mat1c),Mat3(Mat1r,Mat1c)
! -----------------------------  Locals  ------------------------------
        INTEGER Row,Col

!         ! --------------------  Implementation   ----------------------
        DO Row = 1, Mat1r

            DO Col =1, Mat1c

                  Mat3(Row,Col)= Mat1(Row,Col) + Mat2(Row,Col)
              ENDDO

          ENDDO

      RETURN
      END
!.
! ------------------------------------------------------------------------------
!
!                           SUBROUTINE MATSUB
!
!  this subroutine subtracts two matricies together.
!
!  Inputs          Description                    
!    Mat1        - Matrix number 1
!    Mat2        - Matrix number 2
!    Mat1r       - Matrix number 1 rows
!    Mat1c       - Matrix number 1 columns
!
!  OutPuts       :
!    Mat3        - Matrix result of Mat1 + Mat2
!                    of size mat1r x mat1c
!
!  Locals        :
!    Row         - Row Index
!    Col         - Column Index
!
!  Coupling      :
!
! ------------------------------------------------------------------------------

      SUBROUTINE MATSUB ( Mat1,Mat2,Mat1r,Mat1c,Mat3 )
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c
        REAL*8 Mat1(Mat1r,Mat1c),Mat2(Mat1r,Mat1c),Mat3(Mat1r,Mat1c)
! -----------------------------  Locals  ------------------------------
        INTEGER Row,Col

!         ! --------------------  Implementation   ----------------------
        DO Row= 1, Mat1r
            DO Col= 1, Mat1c
                  Mat3(Row,Col)= Mat1(Row,Col) - Mat2(Row,Col)
              ENDDO

          ENDDO

      RETURN
      END
!.
! ------------------------------------------------------------------------------
!
!                           SUBROUTINE MATTRANS
!
!  this subroutine finds the transpose of a matrix.
!
!  Inputs          Description                    
!    Mat1        - Matrix number 1
!    Mat1r       - Matrix number 1 rows
!    Mat1c       - Matrix number 1 columns
!
!  OutPuts       :
!    Mat2        - Matrix result of transpose Mat2
!
!  Locals        :
!    Row         - Row Index
!    Col         - Column Index
!
!  Coupling      :
!
! ------------------------------------------------------------------------------  

      SUBROUTINE MATTRANS ( Mat1,Mat1r,Mat1c,Mat2 )
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c
        REAL*8 Mat1(Mat1r,Mat1c),Mat2(Mat1c,Mat1r)
! -----------------------------  Locals  ------------------------------
        INTEGER Row,Col

!         ! --------------------  Implementation   ----------------------
        DO Row =1, Mat1r
            DO Col = 1, Mat1c
                  Mat2(Col,Row)= Mat1(Row,Col)
              ENDDO
          ENDDO
      RETURN
      END
!.
!--------------------------------------------------------------------------------
!
!                           SUBROUTINE MatVecMult
! ------------------------------------------------------------------------------
      SUBROUTINE MatVecMult (Mat1, Vec2, Mat1r, Mat1c, Vec3 )
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c
        REAL*8 Mat1(Mat1r,Mat1c), Vec2(Mat1c), Vec3(Mat1r)
!         ! -- local vars ---
        INTEGER Row, ktr

        DO Row=1, mat1r
            Vec3(Row) = 0.0D0
            DO ktr= 1, mat1c
                Vec3(Row) = Vec3(Row) + Mat1(Row,ktr) * Vec2(ktr)
              ENDDO
          ENDDO
      RETURN
      END   ! SUBROUTINE MatVecMult


! ------------------------------------------------------------------------------

      subroutine matinv(a,n,l)
! C
! C     gauss-jordan method to compute the inverse of matrix a
! C     a -- is the n*n matrix ,it's input and output(inverse of a)
! C     n--- is the n ,input
! C     l--- if the matrix a is dsingular ,l=0 ,else l !=0 
! C     is,js-----is a matrix is(n),and js(n)
! C
	dimension a(n,n),is(n),js(n)
	real*8 a,t,d,is,js
	l=1
	do 100 k=1,n
	  d=0.0
	  do 10 i=k,n
	  do 10 j=k,n
	    if (abs(a(i,j)).gt.d) then
	      d=abs(a(i,j))
	      is(k)=i
	      js(k)=j
	    end if
   10	  continue
	  if (d+1.0.eq.1.0) then
	    l=0
	    write(*,20)
	    return
	  end if
   20	  format(1x,'err**not inv')
	  do 30 j=1,n
	    t=a(k,j)
	    a(k,j)=a(is(k),j)
	    a(is(k),j)=t
   30	  continue
	  do 40 i=1,n
	    t=a(i,k)
	    a(i,k)=a(i,js(k))
	    a(i,js(k))=t
   40	  continue
	  a(k,k)=1/a(k,k)
	  do 50 j=1,n
	    if (j.ne.k) then
	      a(k,j)=a(k,j)*a(k,k)
	    end if
   50	  continue
	  do 70 i=1,n
	    if (i.ne.k) then
	      do 60 j=1,n
	        if (j.ne.k) then
	          a(i,j)=a(i,j)-a(i,k)*a(k,j)
	        end if
  60	      continue
	    end if
   70	  continue
	  do 80 i=1,n
	    if (i.ne.k) then
	      a(i,k)=-a(i,k)*a(k,k)
	    end if
   80	  continue
  100	continue
	do 130 k=n,1,-1
	  do 110 j=1,n
	    t=a(k,j)
	    a(k,j)=a(js(k),j)
	    a(js(k),j)=t
  110	  continue
	  do 120 i=1,n
	    t=a(i,k)
	    a(i,k)=a(i,is(k))
	    a(i,is(k))=t
  120	  continue
  130	continue
	return
	end

      subroutine BLR2XYZ(theta, lambda, r, x, y, z)
!-------------------------------------------------------------------------
!
!	From spherical coordinates theta,lambda,r to Cartesian coordinates X,Y,Z 
!
!
!   Given: 
!		theta		d		geocentric latitude, in degree
!		lambda		d		longitude, in degree
!		r		d		geocentric distance, in m 
!
!   Returned: 
!		x		d		Space rectangular coordinate X, in m
!		y		d		Space rectangular coordinate Y, in m
!		z		d		Space rectangular coordinate Z, in m
!
!
!-------------------------------------------------------------------------- 
      
      implicit none
      include 'commondata.h'
      
      real*8, intent(in) :: theta, lambda, r
      real*8, intent(out) :: x, y, z     
    
      pi = 4.0d0 * atan(1.0d0)

	x = r * cos(theta/180.0d0*pi) * cos(lambda/180.0d0*pi)
	y = r * cos(theta/180.0d0*pi) * sin(lambda/180.0d0*pi)
	z = r * sin(theta/180.0d0*pi)

      end subroutine

      subroutine XYZ2BLR(x,y,z,theta,lambda,r)
!-------------------------------------------------------------------------
!
!	From Cartesian coordinates X,Y,Z to spherical coordinates theta,lambda,r
!
!
!   Given: 
!		x		d		Space rectangular coordinate X, in m
!		y		d		Space rectangular coordinate Y, in m
!		z		d		Space rectangular coordinate Z, in m
!
!   Returned: 
!		theta	d		geocentric latitude, in degree
!		lambda	d		longitude, in degree
!		r		d		geocentric distance, in m
!
!
!--------------------------------------------------------------------------
      
      implicit none
      include 'commondata.h'

      real*8, intent(in) :: x, y, z     
      real*8, intent(out) :: theta, lambda, r
      pi = 4.0d0 * atan(1.0d0)
    

	!-----------------------------------------------------------------------
	! We don't take the x = y = 0 case into account
	!-----------------------------------------------------------------------
      r = sqrt(x*x + y*y + z*z)
      theta = atan(z /(sqrt(x*x + y*y)))*180.0d0/pi
      lambda = datan2(y,x)*180.0d0/pi;

      end subroutine