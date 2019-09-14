C *************************************************
C *************************************************
C *          UTILITY SUBROUTINES                  *
C *************************************************
C *************************************************



C *************************************************
C *       GENERALISED MATRIX MULTIPLICATION       *
C *************************************************
      
    ! A, B, OUT
      SUBROUTINE KMLTM(A,B,OUT,NCOL,M,N)
      INCLUDE 'aba_param.inc'
      DIMENSION A(M,NCOL),B(NCOL,N),OUT(M,N)

      DO I=1,M
        DO J=1,N
           OUT(i,j) = 0
           DO K=1,NCOL
              OUT(i,j) = OUT(i,j) + A(i, k)*B(k, j)
           END DO
        END DO
      END DO
      
      RETURN
      END
      
      
      ! multiply arbitrary Mx9 matrix by 9x1 vector
      SUBROUTINE KMLTMX991(DMIN,DVIN,DVOUT,M)
      INCLUDE 'aba_param.inc'
      PARAMETER (N=9)
      DIMENSION DMIN(M,N),DVIN(N),DVOUT(M)

      DO i=1,M
        x = 0.
        DO j=1,9
         x = x + DMIN(i,j)*DVIN(j)
        END DO
        
        DVOUT(i) = x
      END DO

      RETURN
      END
      

      SUBROUTINE KMLT3331(DMIN,DVIN,DVOUT)
      INCLUDE 'aba_param.inc'
      PARAMETER (M=3,N=3)
      DIMENSION DMIN(M,N),DVIN(N),DVOUT(M)

      DO i=1,3
        x = 0.
        DO j=1,3
         x = x + DMIN(i,j)*DVIN(j)
        END DO
        
        DVOUT(i) = x
      END DO

      RETURN
      END
           
      
      SUBROUTINE KMLT360601(DMIN,DVIN,DVOUT)
      INCLUDE 'aba_param.inc'
      PARAMETER (M=3,N=60)
      DIMENSION DMIN(M,N),DVIN(N),DVOUT(M)

      DO i=1,3
        x = 0.
        DO j=1,60
         x = x + DMIN(i,j)*DVIN(j)
        END DO
        
        DVOUT(i) = x
      END DO

      RETURN
      END

      SUBROUTINE KMLT60661(DMIN,DVIN,DVOUT)
      INCLUDE 'aba_param.inc'
      PARAMETER (M=60,N=6)
      DIMENSION DMIN(M,N),DVIN(N),DVOUT(M)

      DO i=1,60
        x = 0.
        DO j=1,6
         x = x + DMIN(i,j)*DVIN(j)
        END DO
        
        DVOUT(i) = x
      END DO

      RETURN
      END
      
      
      SUBROUTINE KMLT36991(DMIN,DVIN,DVOUT)
      INCLUDE 'aba_param.inc'
      PARAMETER (M=36,N=9)
      DIMENSION DMIN(M,N),DVIN(N),DVOUT(M)

      DO i=1,36
        x = 0.
        DO j=1,9
         x = x + DMIN(i,j)*DVIN(j)
        END DO
        
        DVOUT(i) = x
      END DO

      RETURN
      END
      
      
      SUBROUTINE KMLT6661(DMIN,DVIN,DVOUT)
      INCLUDE 'aba_param.inc'
      PARAMETER (M=6,N=6)
      DIMENSION DMIN(M,N),DVIN(N),DVOUT(M)

      DO i=1,6
        x = 0.
        DO j=1,6
         x = x + DMIN(i,j)*DVIN(j)
        END DO
        
        DVOUT(i) = x
      END DO

      RETURN
      END
      
      SUBROUTINE KMLT8881(DMIN,DVIN,DVOUT)
      INCLUDE 'aba_param.inc'
      PARAMETER (M=8,N=8)
      DIMENSION DMIN(M,N),DVIN(N),DVOUT(M)

      DO i=1,8
        x = 0.
        DO j=1,8
         x = x + DMIN(i,j)*DVIN(j)
        END DO
        
        DVOUT(i) = x
      END DO

      RETURN
      END



C *************************************************
C *              VECTOR PRODUCT                   *
C *************************************************

      SUBROUTINE CrossProd(A,B,C)
      INCLUDE 'aba_param.inc'
      PARAMETER (M=3)
      DIMENSION A(M),B(M),C(M)

      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)

      RETURN
      END

C *************************************************
C *              DYADIC PRODUCT                   *
C *************************************************

      SUBROUTINE DyadicProd(A,B,C)
      INCLUDE 'aba_param.inc'
      PARAMETER (M=3)
      DIMENSION A(M),B(M),C(M,M)
      
      DO i=1,3
        DO j=1,3
            C(i,j) = A(i)*B(j)
        END DO
      END DO
      
      RETURN
      END

C *************************************************
C *        CONVERT 2D MATRIX TO 1D COLUMN         *
C *************************************************

      SUBROUTINE Convert2Col(A,B)
      INCLUDE 'aba_param.inc'
      PARAMETER (M=3, N=9)
      DIMENSION A(M,M),B(N)
      
      B(1) = A(1,1)
      B(2) = A(2,1)
      B(3) = A(3,1)

      B(4) = A(1,2)
      B(5) = A(2,2)
      B(6) = A(3,2)

      B(7) = A(1,3)
      B(8) = A(2,3)
      B(9) = A(3,3)
      
      RETURN
      END

C *************************************************
C *           TRACE OF A 3X3 MATRIX               *
C *************************************************
      real*8 function trace(A)
      real*8, intent(in):: A(3,3)      
      trace = A(1,1)+A(2,2)+A(3,3)      
      return
      end function
      
      !Trace for any size of square matrix
!      real*8 function trace(A)
!      real*8, intent(in):: A(:,:)
!      real*8 :: X
!      integer :: N
!      
!      X = 0.
!      do i=1,ubound(A,1)      
!      X = X+ A(i,i)
!      end do
!      trace = X   
!      return
!      end function

C *************************************************
C *      TRANSFER 3X3 MATRIX TO 6X1 COLUMN VECTOR *         
C *************************************************
      SUBROUTINE KMATVEC6(DMIN,DVOUT)
       INCLUDE 'aba_param.inc'
        PARAMETER (M=3,N=3,K=6)
        DIMENSION DMIN(M,N),DVOUT(K)
        DO I=1,M
            DVOUT(I)=DMIN(I,I)
        END DO
            DVOUT(4) = DMIN(1,2)
            DVOUT(5) = DMIN(1,3)
            DVOUT(6) = DMIN(2,3)
        RETURN
      END
C *************************************************
C * TRANSFER 6X1 COLUMN VECTOR TO 3X3 MATRIX  *
C *************************************************

      SUBROUTINE KVECMAT6(DVIN,DMOUT)
       INCLUDE 'aba_param.inc'

        PARAMETER (M=3,N=3,K=6)

        DIMENSION DVIN(K),DMOUT(M,N)
        DO I=1,M
            DMOUT(I,I) = DVIN(I)
        END DO
            DMOUT(1,2) = DVIN(4)
            DMOUT(2,1) = DVIN(4)
            DMOUT(1,3) = DVIN(5)
            DMOUT(3,1) = DVIN(5)
            DMOUT(3,2) = DVIN(6)
            DMOUT(2,3) = DVIN(6)            
        RETURN
      END
C *************************************************
C *   VECTOR PRODUCT OF 3X1 WITH 3X1 GIVING 3X1   *
C *************************************************
      SUBROUTINE KVECPROD(DVIN1,DVIN2,DVOUT)
        INCLUDE 'ABA_PARAM.INC'
        PARAMETER(M=3)
        DIMENSION DVIN1(M),DVIN2(M),DVOUT(M)
        DVOUT(1)=DVIN1(2)*DVIN2(3)-DVIN2(2)*DVIN1(3)
        DVOUT(2)=DVIN2(1)*DVIN1(3)-DVIN1(1)*DVIN2(3)
        DVOUT(3)=DVIN1(1)*DVIN2(2)-DVIN2(1)*DVIN1(2)
        RETURN
      END 
C *************************************************
C *TRANSFER GENERAL 3X3 MATRIX TO 6X1 COLUMN VECTOR *         
C *************************************************
      SUBROUTINE KGMATVEC6(DMIN,DVOUT)
        INCLUDE 'ABA_PARAM.INC'
        PARAMETER (M=3,N=3,K=6)
        DIMENSION DMIN(M,N),DVOUT(K)
        DO I=1,M
            DVOUT(I)=DMIN(I,I)
        END DO
            DVOUT(4) = DMIN(1,2)+DMIN(2,1)
            DVOUT(5) = DMIN(3,1)+DMIN(1,3)
            DVOUT(6) = DMIN(2,3)+DMIN(3,2)            
        RETURN
      END
C *************************************************
C *        INVERSE OF A MATRIX  WITH LAPACK       *
C *************************************************
      !Checked inverse against python's numpy.linalg.inv
      subroutine lapinverse(xmatin,m,info,xmatout)
      implicit none
      
      integer,intent(in) :: m
      real*8,intent(in):: xmatin(m,m) !abaqus won't allow xmatin(:,:)
      
      integer,parameter :: LWORK = 8000
      real*8,parameter  :: zero=1.0e-6
      integer :: LDA,ialloc,i,n,j
      
      integer,intent(out) :: INFO
      real*8,intent(out):: xmatout(m,m)
      integer, allocatable ::  IPIV(:)
      real*8, allocatable :: A(:,:)
      real*8 :: WORK(LWORK)
      EXTERNAL DGETRI, DGETRF
      
      LDA = max(1,m); n = m
      xmatout = 0.
       
      allocate(A(LDA,n),IPIV(min(m,n)),stat=ialloc)      
      
      A = xmatin !Don't input array xmatin
          
      call DGETRF( M, N, A, LDA, IPIV, INFO )
      
      if(info < 0) then
!      write (6,*)"DGETRF, illegal value at = ",-info,". No inverse"
      xmatout = 0.
      
      else if(info > 0) then
!      write (6,*)"DGETRF, U(i,i) zero at i = ",info,". No inverse"
      xmatout = 0.
      
      else       
      call DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
      
      do i=1,m; do j=1,m; if (abs(A(i,j))<= zero) A(i,j) = 0. 
      end do; end do      
      xmatout = A
      
      end if
      
      deallocate(A,IPIV)
      return
      end subroutine        

C *************************************************
C *          SUBROUTINE MATRIX SQURAE ROOT        *
C *************************************************

      SUBROUTINE KMSQRT(a,b)

         INCLUDE 'ABA_PARAM.INC'

      PARAMETER (n=3,nx=3)

      DIMENSION a(n,n),diag(n,n),q(n,n),d(n),b(n,n),
     + qtrans(n,n),result(n,n)

             diag=0.; b=0.; q=0.;result=0.;d=0.
             
      call Kjacobi(a,n,d,q) 
      call Keigsrt(d,q,n)
C
      do i=1,n
         if (d(i) .ge. 0) then
            diag(i,i)=sqrt(d(i))
         else
            write (6,*) 'the matrix is not positive definite'
         end if
      end do

      result = matmul(q,diag)
      b = matmul(result,transpose(q))
      
      return
      end  
C
C
C *************************************************
C * SUBROUTINE MATRIX EIGENVALUES AND EIGENVECTORS*
C *************************************************
C
      SUBROUTINE KJACOBI(a,n,d,v)

        INCLUDE 'ABA_PARAM.INC'


      PARAMETER (NMAX=500) 


      DIMENSION a(n,n),d(n),v(n,n),b(NMAX),z(NMAX),dial(n,n)

      do ip=1,n
         do iq=1,n
            v(ip,iq)=0.
         end do
         v(ip,ip)=1.
      end do
C
      do ip=1,n
         b(ip)=a(ip,ip)
         d(ip)=b(ip)
         z(ip)=0.
      end do
C
      nrot=0
      do i=1,50
         sm=0.
         do ip=1,n-1
            do iq=ip+1,n
               sm=sm+abs(a(ip,iq))
            end do
         end do
C
         if (sm .eq. 0.) return
         if (i .lt. 4) then
           tresh=0.2*sm/n**2
         else
           tresh=0.
         end if
C
         do ip=1,n-1
           do iq=ip+1,n
              g=100.*abs(a(ip,iq))
              if ((i .gt. 4) .and. (abs(d(ip))+g .eq. abs(d(ip)))
     + .and. (abs(d(ip))+g .eq. abs(d(iq)))) then
                a(ip,iq)=0.
              else if (abs(a(ip,iq)) .gt. tresh) then
                h=d(iq)-d(ip)
                if (abs(h)+g .eq. abs(h)) then
                   t=a(ip,iq)/h 
                else
                   theta=0.5*h/a(ip,iq)
                   t=1./(abs(theta)+sqrt(1.+theta**2))
                   if (theta .lt. 0) t=-t
                end if
                c=1./sqrt(1+t**2)
                s=t*c
                ta=s/(1.+c)
                h=t*a(ip,iq)
                z(ip)=z(ip)-h
                z(iq)=z(iq)+h
                d(ip)=d(ip)-h
                d(iq)=d(iq)+h
                a(ip,iq)=0.
C
                do j=1,ip-1
                  g=a(j,ip)
                  h=a(j,iq)
                  a(j,ip)=g-s*(h+g*ta)
                  a(j,iq)=h+s*(g-h*ta)
                end do
C
                do j=ip+1,iq-1
                 g=a(ip,j)
                 h=a(j,iq)
                 a(ip,j)=g-s*(h+g*ta)
                 a(j,iq)=h+s*(g-h*ta)
                end do
C
                do j=iq+1,n
                 g=a(ip,j)
                 h=a(iq,j)
                 a(ip,j)=g-s*(h+g*ta)
                 a(iq,j)=h+s*(g-h*ta)
                end do
C
                do j=1,n
                  g=v(j,ip)
                  h=v(j,iq)
                  v(j,ip)=g-s*(h+g*ta)
                  v(j,iq)=h+s*(g-h*ta)
                end do
C
                nrot=nrot+1
              end if
           end do
         end do
C
         do ip=1,n
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip)
            z(ip)=0.
         end do
      end do
C      

C
      return
C
      END
C *************************************************
C *          SUBROUTINE SORT EIGENVALUES          *
C *************************************************
      SUBROUTINE kEIGSRT(D,V,N)
        INCLUDE 'ABA_PARAM.INC'
      DIMENSION D(N),V(N,N)
      DO I=1,N-1
        K=I
        P=D(I)
        DO J=I+1,N
          IF(D(J).GE.P)THEN
            K=J
            P=D(J)
          ENDIF
        END DO
        IF(K.NE.I)THEN
          D(K)=D(I)
          D(I)=P
          DO J=1,N
            P=V(J,I)
            V(J,I)=V(J,K)
            V(J,K)=P
          END DO
        ENDIF
      END DO
      RETURN
      END
C *************************************************
C *      THE DETERMINANT OF A 3X3 MATRIX          *
C *************************************************
      SUBROUTINE KDETER(DMIN,D)

        INCLUDE 'ABA_PARAM.INC'
        PARAMETER(M=3,N=3)
        DIMENSION DMIN(M,N)
        D=0.
        D=DMIN(1,1)*DMIN(2,2)*DMIN(3,3)+DMIN(1,2)*
     + DMIN(2,3)*DMIN(3,1)+DMIN(2,1)*DMIN(3,2)*
     + DMIN(1,3)-DMIN(1,3)*DMIN(2,2)*DMIN(3,1)-
     + DMIN(1,1)*DMIN(2,3)*DMIN(3,2)-DMIN(1,2)*
     + DMIN(2,1)*DMIN(3,3)
        RETURN
      END 

C *************************************************
C *      Build 4th order rotation matrix (tsigma) *
C *************************************************
      subroutine rotord4sig(xrot,tSigma) 
      
      real*8, intent(in) :: xrot(3,3)
      real*8, intent(out) :: tSigma(6,6)
      
      tSigma(1,1) = xRot(1,1)*xRot(1,1)
      tSigma(1,2) = xRot(2,1)*xRot(2,1)
      tSigma(1,3) = xRot(3,1)*xRot(3,1)
      tSigma(1,4) = 2.0*xRot(1,1)*xRot(2,1)
      tSigma(1,5) = 2.0*xRot(2,1)*xRot(3,1)
      tSigma(1,6) = 2.0*xRot(3,1)*xRot(1,1)

      tSigma(2,1) = xRot(1,2)*xRot(1,2)
      tSigma(2,2) = xRot(2,2)*xRot(2,2)
      tSigma(2,3) = xRot(3,2)*xRot(3,2)
      tSigma(2,4) = 2.0*xRot(1,2)*xRot(2,2)
      tSigma(2,5) = 2.0*xRot(2,2)*xRot(3,2)
      tSigma(2,6) = 2.0*xRot(3,2)*xRot(1,2)

      tSigma(3,1) = xRot(1,3)*xRot(1,3)
      tSigma(3,2) = xRot(2,3)*xRot(2,3)
      tSigma(3,3) = xRot(3,3)*xRot(3,3)
      tSigma(3,4) = 2.0*xRot(1,3)*xRot(2,3)
      tSigma(3,5) = 2.0*xRot(2,3)*xRot(3,3)
      tSigma(3,6) = 2.0*xRot(3,3)*xRot(1,3)

      tSigma(4,1) = xRot(1,1)*xRot(1,2)
      tSigma(4,2) = xRot(2,1)*xRot(2,2)
      tSigma(4,3) = xRot(3,1)*xRot(3,2)
      tSigma(4,4) = xRot(1,1)*xRot(2,2) + xRot(1,2)*xRot(2,1)
      tSigma(4,5) = xRot(2,1)*xRot(3,2) + xRot(2,2)*xRot(3,1) 
      tSigma(4,6) = xRot(3,1)*xRot(1,2) + xRot(3,2)*xRot(1,1)

      tSigma(5,1) = xRot(1,2)*xRot(1,3)
      tSigma(5,2) = xRot(2,2)*xRot(2,3)
      tSigma(5,3) = xRot(3,2)*xRot(3,3)
      tSigma(5,4) = xRot(1,2)*xRot(2,3) + xRot(1,3)*xRot(2,2)
      tSigma(5,5) = xRot(2,2)*xRot(3,3) + xRot(2,3)*xRot(3,2) 
      tSigma(5,6) = xRot(3,2)*xRot(1,3) + xRot(3,3)*xRot(1,2)

      tSigma(6,1) = xRot(1,3)*xRot(1,1)
      tSigma(6,2) = xRot(2,3)*xRot(2,1)
      tSigma(6,3) = xRot(3,3)*xRot(3,1)
      tSigma(6,4) = xRot(1,3)*xRot(2,1) + xRot(1,1)*xRot(2,3)
      tSigma(6,5) = xRot(2,3)*xRot(3,1) + xRot(2,1)*xRot(3,3) 
      tSigma(6,6) = xRot(3,3)*xRot(1,1) + xRot(3,1)*xRot(1,3)
      
      return
      end subroutine

C *************************************************
C *      Build 4th order rotation matrix (tstran) *
C *************************************************
      subroutine rotord4str(xrot,tStran)
      
      real*8, intent(in) :: xrot(3,3)
      real*8, intent(out) :: tStran(6,6)
      
      tStran(1,1) = xRot(1,1)*xRot(1,1)
      tStran(1,2) = xRot(2,1)*xRot(2,1)
      tStran(1,3) = xRot(3,1)*xRot(3,1)
      tStran(1,4) = xRot(1,1)*xRot(2,1)
      tStran(1,5) = xRot(2,1)*xRot(3,1)
      tStran(1,6) = xRot(3,1)*xRot(1,1)

      tStran(2,1) = xRot(1,2)*xRot(1,2)
      tStran(2,2) = xRot(2,2)*xRot(2,2)
      tStran(2,3) = xRot(3,2)*xRot(3,2)
      tStran(2,4) = xRot(1,2)*xRot(2,2)
      tStran(2,5) = xRot(2,2)*xRot(3,2)
      tStran(2,6) = xRot(3,2)*xRot(1,2)

      tStran(3,1) = xRot(1,3)*xRot(1,3)
      tStran(3,2) = xRot(2,3)*xRot(2,3)
      tStran(3,3) = xRot(3,3)*xRot(3,3)
      tStran(3,4) = xRot(1,3)*xRot(2,3)
      tStran(3,5) = xRot(2,3)*xRot(3,3)
      tStran(3,6) = xRot(3,3)*xRot(1,3)

      tStran(4,1) = 2.0*xRot(1,1)*xRot(1,2)
      tStran(4,2) = 2.0*xRot(2,1)*xRot(2,2)
      tStran(4,3) = 2.0*xRot(3,1)*xRot(3,2)
      tStran(4,4) = xRot(1,1)*xRot(2,2) + xRot(1,2)*xRot(2,1)
      tStran(4,5) = xRot(2,1)*xRot(3,2) + xRot(2,2)*xRot(3,1) 
      tStran(4,6) = xRot(3,1)*xRot(1,2) + xRot(3,2)*xRot(1,1)

      tStran(5,1) = 2.0*xRot(1,2)*xRot(1,3)
      tStran(5,2) = 2.0*xRot(2,2)*xRot(2,3)
      tStran(5,3) = 2.0*xRot(3,2)*xRot(3,3)
      tStran(5,4) = xRot(1,2)*xRot(2,3) + xRot(1,3)*xRot(2,2)
      tStran(5,5) = xRot(2,2)*xRot(3,3) + xRot(2,3)*xRot(3,2) 
      tStran(5,6) = xRot(3,2)*xRot(1,3) + xRot(3,3)*xRot(1,2)

      tStran(6,1) = 2.0*xRot(1,3)*xRot(1,1)
      tStran(6,2) = 2.0*xRot(2,3)*xRot(2,1)
      tStran(6,3) = 2.0*xRot(3,3)*xRot(3,1)
      tStran(6,4) = xRot(1,3)*xRot(2,1) + xRot(1,1)*xRot(2,3)
      tStran(6,5) = xRot(2,3)*xRot(3,1) + xRot(2,1)*xRot(3,3) 
      tStran(6,6) = xRot(3,3)*xRot(1,1) + xRot(3,1)*xRot(1,3)
      
      return
      end subroutine

C *************************************************************
C *      Build 4th order lower (and upper) tensor product     *
C * NB: When contracted from 4th to the format used for       *
C * C-matrix, it turns out that the upper and lower products  *
C * are the same!                                             *
C *************************************************************
      subroutine ltprod(A,B,C)
      
      real*8, intent(in) :: A(3,3),B(3,3)
      real*8, intent(out) :: C(6,6)
      
      C = 0.
      C(1,1) = 2*A(1,1)*B(1,1)
      C(1,2) = 2*A(1,2)*B(1,2)
      C(1,3) = 2*A(1,3)*B(1,3)
      C(1,4) = A(1,1)*B(1,2)+A(1,2)*B(1,1)
      C(1,5) = A(1,1)*B(1,3)+A(1,3)*B(1,1)
      C(1,6) = A(1,2)*B(1,3)+A(1,3)*B(1,2)
      
      C(2,2) = 2*A(2,2)*B(2,2)
      C(2,3) = 2*A(2,3)*B(2,3)
      C(2,4) = A(2,1)*B(2,2)+A(2,2)*B(2,1)
      C(2,5) = A(2,1)*B(2,3)+A(2,3)*B(2,1)
      C(2,6) = A(2,2)*B(2,3)+A(2,3)*B(2,2)
      
      C(3,3) = 2*A(3,3)*B(3,3)
      C(3,4) = A(3,1)*B(3,2)+A(3,2)*B(3,1)
      C(3,5) = A(3,1)*B(3,3)+A(3,3)*B(3,1)
      C(3,6) = A(3,2)*B(3,3)+A(3,3)*B(3,2)
      
      C(4,4) = A(1,1)*B(2,2)+A(1,2)*B(2,1)
      C(4,5) = A(1,1)*B(2,3)+A(1,3)*B(2,1)
      C(4,6) = A(1,2)*B(2,3)+A(1,3)*B(2,2)
      
      C(5,5) = A(1,1)*B(3,3)+A(1,3)*B(3,1)
      C(5,6) = A(1,2)*B(3,3)+A(1,3)*B(3,2)
      
      C(6,6) = A(2,2)*B(3,3)+A(2,3)*B(3,2)
      
      do i=2,6
         do j=1,i-1
            C(i,j)=C(j,i)
         end do
      end do
      
      C = 0.5*C
      
      return
      end subroutine      

C **************************************************
C *      Build 4th order tensor product (kronecker)*
C **************************************************
      subroutine tprod(A,B,C)
      
      real*8, intent(in) :: A(3,3),B(3,3)
      real*8, intent(out) :: C(6,6)
      
      C = 0.
      C(1,1) = 2*A(1,1)*B(1,1)
      C(1,2) = 2*A(1,1)*B(2,2)
      C(1,3) = 2*A(1,1)*B(3,3)
      C(1,4) = A(1,1)*B(1,2)+A(1,1)*B(2,1)
      C(1,5) = A(1,1)*B(1,3)+A(1,1)*B(3,1)
      C(1,6) = A(1,1)*B(2,3)+A(1,1)*B(3,2)
      
      C(2,2) = 2*A(2,2)*B(2,2)
      C(2,3) = 2*A(2,2)*B(3,3)
      C(2,4) = A(2,2)*B(1,2)+A(2,2)*B(2,1)
      C(2,5) = A(2,2)*B(1,3)+A(2,2)*B(3,1)
      C(2,6) = A(2,2)*B(2,3)+A(2,2)*B(3,2)
      
      C(3,3) = 2*A(3,3)*B(3,3)
      C(3,4) = A(3,3)*B(1,2)+A(3,3)*B(2,1)
      C(3,5) = A(3,3)*B(1,3)+A(3,3)*B(3,1)
      C(3,6) = A(3,3)*B(2,3)+A(3,3)*B(3,2)
      
      C(4,4) = A(1,2)*B(1,2)+A(1,2)*B(2,1)
      C(4,5) = A(1,2)*B(1,3)+A(1,2)*B(3,1)
      C(4,6) = A(1,2)*B(2,3)+A(1,2)*B(3,2)
      
      C(5,5) = A(1,3)*B(1,3)+A(1,3)*B(3,1)
      C(5,6) = A(1,3)*B(2,3)+A(1,3)*B(3,2)
      
      C(6,6) = A(2,3)*B(2,3)+A(2,3)*B(3,2)
      
      do i=2,6
         do j=1,i-1
            C(i,j)=C(j,i)
         end do
      end do
      
      C = 0.5*C
      
      return
      end subroutine      

