       SUBROUTINE PerformCurl(NOEL,NPT) 

      implicit none      
      INTEGER, PARAMETER  :: K=3,M=3,N=8,nnodes=8,knsdv=18
      REAL*8, PARAMETER  :: zero=1.0e-6,xgauss = 0.577350269189626
      INTEGER, INTENT(IN) :: NOEL, NPT
      
      REAL*8, PARAMETER :: xnat8(nnodes,3) =reshape([
     1 -1., -1., -1., -1.,  1.,  1.,  1.,  1., 
     2  1., -1.,  1., -1.,  1., -1.,  1., -1., 
     3  1.,  1., -1., -1.,  1.,  1., -1., -1. 
     1  ], [nnodes,3])
     
       REAL*8, PARAMETER :: gauss(nnodes,3) =reshape([
     1 -xgauss, -xgauss, -xgauss, -xgauss,  
     +  xgauss,  xgauss,  xgauss,  xgauss, 
     2  xgauss, -xgauss,  xgauss, -xgauss,
     +  xgauss, -xgauss,  xgauss, -xgauss, 
     3  xgauss,  xgauss, -xgauss, -xgauss,  
     +  xgauss,  xgauss, -xgauss, -xgauss 
     1  ], [nnodes,3])

c ------     
      real*8,intent(inout):: svars(144)
      INTEGER :: I,J,K, kintB, kintB2
      REAL*8 :: gausscoords(3,nnodes), yp, yq, yr

      real*8 :: xj(3,3),xjinv(3,3),dndloc(nnodes,3),dndx(nnodes,3),
     + fnode(3,nnodes), fmat1(3,3),dmout(3,3),xn(nnodes)
     
      !Extrapolation arrays
      real*8 :: z11i(nnodes),z11n(nnodes),z12i(nnodes),z12n(nnodes),
     + z13i(nnodes),z13n(nnodes),z21i(nnodes),z21n(nnodes),z22i(nnodes),
     + z22n(nnodes),z23i(nnodes),z23n(nnodes),z31i(nnodes),z31n(nnodes),
     + z32i(nnodes),z32n(nnodes),z33i(nnodes),z33n(nnodes),
     + xnmat(nnodes,nnodes),xnmatI(nnodes,nnodes)
     
      character(len=*),parameter :: fmt20 = "(' ',20(F6.4,1X))",
     + fmt3="(3(' ',(ES11.3,1X)))",fmt6 = "(' ',6(F12.3,1X))",
     + fmt5 = "(' ',5(F12.8,1X))",fmt8="(8(' ',(ES11.3,1X)))"

      include 'UserParameters.f'       
      REAL*8 :: kgausscoords, kFp, kcurlFp,      
      COMMON/UMPS/kgausscoords(TOTALELEMENTNUM,8,3),
     1 kFp(TOTALELEMENTNUM,8, 9),
     1 kcurlFp(TOTALELEMENTNUM, 8, 9)     

c =========================================================     
      fnode = 0.;fmat1 = 0.;xnmat=0.;xnmatI=0.
      z11i=0.;z11n=0.;z12i=0.;z12n=0.;z13i=0.;z13n=0.
      z21i=0.;z21n=0.;z22i=0.;z22n=0.;z23i=0.;z23n=0.
      z31i=0.;z31n=0.;z32i=0.;z32n=0.;z33i=0.;z33n=0.

c -----------------------------------------------
      DO J =1,8    
         DO i=1,3         
             gausscoords(i,J) = kgausscoords(noel,J,i)                          
         END DO
     
         DO i=1,9          
             svars(i + knsdv*(J-1)) = kFp(noel,J,i)         
         END DO
         END DO	  
	  
   
c ===================================================     
      do kintB2 = 1,8

C    SPECIFY yp,yq,yr - INTEGRATION POINT
      yp = gauss(kintB2,1)
      yq = gauss(kintB2,2)
      yr = gauss(kintB2,3)

C    SHAPE FUNCTIONS AND DERIVATIVES
      call kshapes8(yp,yq,yr,xnat8,xn,dndloc)      
!      write(6,*)"xn",kintB2; write(6,fmt8)xn
      
      xnmat(kintB2,:) = xn
!      write(6,*)"xnmat",kintB2; write(6,fmt8)(xnmat(i,:),i=1,8)      
C   
C     SET UP JACOBIAN           
C  
      xj = matmul(gausscoords,dndloc)
!      write(6,*)"xj"
!      write(6,*)xj
!      write(6,*)"xj"
!      write(6,*)xj
C
C    AND ITS INVERSE
C
      call KDETER(xj,det)
      
!      write(6,*)"kintB2,det",kintB2,det
      
      if (abs(det) <= 1.0e-6 .or. det /= det) then !last part true if det=NaN
               z11i(kintB2) = 0.0
               z21i(kintB2) = 0.0
               z31i(kintB2) = 0.0
               z12i(kintB2) = 0.0
               z22i(kintB2) = 0.0
               z32i(kintB2) = 0.0              
               z13i(kintB2) = 0.0
               z23i(kintB2) = 0.0
               z33i(kintB2) = 0.0               
      else  
         call lapinverse(xj,3,info,xjinv)
!         if(info /= 0) write(6,*) "inverse failure: xj in kcurl"
C
      dndx = matmul(dndloc,xjinv) 
C
C    DETERMINE first column of curlf: Read row, to determine column.
C
      do kintB = 1,8 !integration points are our nodes now
      fnode(1,kintB) = svars((kintB-1)*knsdv+1)
      fnode(2,kintB) = svars((kintB-1)*knsdv+4)
      fnode(3,kintB) = svars((kintB-1)*knsdv+7)
      end do
C
      fmat1 = matmul(fnode,dndx)
!      dmout(1,1) = fmat1(3,2) - fmat1(2,3)
!      dmout(2,1) = fmat1(1,3) - fmat1(3,1)
!      dmout(3,1) = fmat1(2,1) - fmat1(1,2)
      
      !Curlfp at integeration points
      z11i(kintB2) = fmat1(3,2) - fmat1(2,3)
      z21i(kintB2) = fmat1(1,3) - fmat1(3,1)
      z31i(kintB2) = fmat1(2,1) - fmat1(1,2)
      
!      write(6,*)"z11i",kintB2; write(6,fmt8)z11i
            
C
C    DETERMINE second column of curlf
C
      do kintB = 1,8
      fnode(1,kintB) = svars((kintB-1)*knsdv+2)
      fnode(2,kintB) = svars((kintB-1)*knsdv+5)
      fnode(3,kintB) = svars((kintB-1)*knsdv+8)
      end do
C
      fmat1 = matmul(fnode,dndx)
!      dmout(1,2) = fmat1(3,2) - fmat1(2,3)
!      dmout(2,2) = fmat1(1,3) - fmat1(3,1)
!      dmout(3,2) = fmat1(2,1) - fmat1(1,2)
      
      !Curlfp at integeration points
      z12i(kintB2) = fmat1(3,2) - fmat1(2,3)
      z22i(kintB2) = fmat1(1,3) - fmat1(3,1)
      z32i(kintB2) = fmat1(2,1) - fmat1(1,2)
      
!      write(6,*)"z12i",kintB2; write(6,fmt8)z12i
             
C
C    DETERMINE third column of curlf
C
      do kintB = 1,8
      fnode(1,kintB) = svars((kintB-1)*knsdv+3)
      fnode(2,kintB) = svars((kintB-1)*knsdv+6)
      fnode(3,kintB) = svars((kintB-1)*knsdv+9)
      end do
C
      fmat1 = matmul(fnode,dndx)
!      dmout(1,3) = fmat1(3,2) - fmat1(2,3)
!      dmout(2,3) = fmat1(1,3) - fmat1(3,1)
!      dmout(3,3) = fmat1(2,1) - fmat1(1,2)
      
      !Curlfp at integeration points
      z13i(kintB2) = fmat1(3,2) - fmat1(2,3)
      z23i(kintB2) = fmat1(1,3) - fmat1(3,1)
      z33i(kintB2) = fmat1(2,1) - fmat1(1,2)         
C
      end if

!  write curlfp into svars (for each integration point!)
! Assigning the values of the inner element's gauss point to the corresponding outer element's gauss point! This can be debated later!
!       do i=1,3
!        do j=1,3 
!        svars((kintB2-1)*knsdv+37+j+(i-1)*3) = dmout(i,j)
!        end do
!       end do
      
      end do 
C ******************************************************      
      !All integration points done. Extrapolation begins
!      write(6,*)"z11i",kintB2; write(6,fmt8)z11i
!      write(6,*)"z12i",kintB2; write(6,fmt8)z12i
!      write(6,*)"xnmat",kintB2; write(6,fmt8)(xnmat(i,:),i=1,8)
      
c      call lapinverse(xnmat,nnodes,info2,xnmatI)
!      if(info2 /= 0) write(6,*) "inverse failure: xnmat in kcurl"
!      write(6,*)"xnmatI",kintB2; write(6,fmt8)(xnmatI(i,:),i=1,8)
      
c      z11n = matmul(xnmatI,z11i)
c      z21n = matmul(xnmatI,z21i)
c      z31n = matmul(xnmatI,z31i)
      
c      z12n = matmul(xnmatI,z12i)
c      z22n = matmul(xnmatI,z22i)
c      z32n = matmul(xnmatI,z32i)
      
c      z13n = matmul(xnmatI,z13i)
c      z23n = matmul(xnmatI,z23i)
c      z33n = matmul(xnmatI,z33i)
 
      !The storage is done by row.
      do kintB=1,nnodes
         svars((kintB-1)*knsdv+10) = z11i(kintB)
         svars((kintB-1)*knsdv+13) = z12i(kintB)
         svars((kintB-1)*knsdv+16) = z13i(kintB)
         
         svars((kintB-1)*knsdv+11) = z21i(kintB)
         svars((kintB-1)*knsdv+14) = z22i(kintB)
         svars((kintB-1)*knsdv+17) = z23i(kintB)
         
         svars((kintB-1)*knsdv+12) = z31i(kintB)
         svars((kintB-1)*knsdv+15) = z32i(kintB)
         svars((kintB-1)*knsdv+18) = z33i(kintB)      
      end do !kintB 
	  
c     ==============================================
   
      call MutexLock( 4 )      
      DO J =1, 8
          DO i=1, 9
              kcurlFp(noel,J,i) = svars(9+i + 18*(J-1))
          END DO
      END DO
      call MutexUnlock( 4 )    
      
      RETURN
      END
C
      include 'kshapes.f'
