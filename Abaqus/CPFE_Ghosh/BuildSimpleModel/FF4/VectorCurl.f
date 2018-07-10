C     STRAIN GRADIENTS
C *******************************************************************************
C *                     Compute Curl of a 2nd order tensor                      *
C *              As a vector, curl of row written in corresponding column       *
C *                          Full integration                                   *
C *******************************************************************************
       SUBROUTINE VectorCurl(svars,xnat8,gauss,gausscoords) 

       INCLUDE 'ABA_PARAM.INC'

      integer, parameter:: K=3,M=3,N=8,nnodes=8
      real*8,parameter  :: zero=1.0e-6,xgauss = 0.577350269189626

      !scalars
      integer,parameter :: knsdv=6, nintpts = 8
      integer,parameter :: nsvars=nintpts * knsdv
      !arrays
      real*8,intent(inout):: svars(nsvars)
      real*8,intent(in) :: xnat8(nnodes,3),gauss(nnodes,3),
     + gausscoords(3,nnodes)

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

      fnode = 0.;fmat1 = 0.;xnmat=0.;xnmatI=0.
      z11i=0.;z11n=0.;z12i=0.;z12n=0.;z13i=0.;z13n=0.
      z21i=0.;z21n=0.;z22i=0.;z22n=0.;z23i=0.;z23n=0.
      z31i=0.;z31n=0.;z32i=0.;z32n=0.;z33i=0.;z33n=0.
C
C   USE EIGHT GAUSS POINT FOR GRADIENTS
C
C    LOOP OVER EIGHT INTEGRATION POINTS

!    Evaluate curlT at the integration points, and simultaneously populate xnmat which will be later inverted for extrapolation purposes.   

      do kint2 = 1,8

C    SPECIFY yp,yq,yr - INTEGRATION POINT
      yp = gauss(kint2,1)
      yq = gauss(kint2,2)
      yr = gauss(kint2,3)

C    SHAPE FUNCTIONS AND DERIVATIVES
      call kshapes8(yp,yq,yr,xnat8,xn,dndloc)      
!      write(6,*)"xn",kint2; write(6,fmt8)xn
      
      xnmat(kint2,:) = xn
!      write(6,*)"xnmat",kint2; write(6,fmt8)(xnmat(i,:),i=1,8)      
C   
C     SET UP JACOBIAN           
C  
      xj = matmul(gausscoords,dndloc)
c      write(6,*)"---------------------------"
c      write(6,*)"gausscoords"
c      write(6,*)gausscoords
c      write(6,*)"xj"
c      write(6,*)xj

C
C    AND ITS INVERSE
C
      call KDETER(xj,det)
c      write(6,*)"DET0=",abs(det)      
!      write(6,*)"kint2,det",kint2,det
      
      if (abs(det) <= zero .or. det /= det) then !last part true if det=NaN
         dmout = 0.0
         
      else  
         call lapinverse(xj,3,info,xjinv)
!         if(info /= 0) write(6,*) "inverse failure: xj in kcurl"
C
      dndx = matmul(dndloc,xjinv) 
C
C    DETERMINE first column of curlf: Read row, to determine column.
C
      do kint = 1,8 !integration points are our nodes now
      fnode(1,kint) = svars((kint-1)*knsdv+1)
      fnode(2,kint) = svars((kint-1)*knsdv+2)
      fnode(3,kint) = svars((kint-1)*knsdv+3)
      end do
C
      fmat1 = matmul(fnode,dndx)
!      dmout(1,1) = fmat1(3,2) - fmat1(2,3)
!      dmout(2,1) = fmat1(1,3) - fmat1(3,1)
!      dmout(3,1) = fmat1(2,1) - fmat1(1,2)
      
      !Curlfp at integeration points
      z11i(kint2) = fmat1(3,2) - fmat1(2,3)
      z21i(kint2) = fmat1(1,3) - fmat1(3,1)
      z31i(kint2) = fmat1(2,1) - fmat1(1,2)
      
!      write(6,*)"z11i",kint2; write(6,fmt8)z11i
            
C     
C
      end if

!  write curlfp into svars (for each integration point!)
! Assigning the values of the inner element's gauss point to the corresponding outer element's gauss point! This can be debated later!
!       do i=1,3
!        do j=1,3 
!        svars((kint2-1)*knsdv+37+j+(i-1)*3) = dmout(i,j)
!        end do
!       end do
      
      end do !kint2
      
      !All integration points done. Extrapolation begins
!      write(6,*)"z11i",kint2; write(6,fmt8)z11i
!      write(6,*)"z12i",kint2; write(6,fmt8)z12i
!      write(6,*)"xnmat",kint2; write(6,fmt8)(xnmat(i,:),i=1,8)
      
      call lapinverse(xnmat,nnodes,info2,xnmatI)
!      if(info2 /= 0) write(6,*) "inverse failure: xnmat in kcurl"
!      write(6,*)"xnmatI",kint2; write(6,fmt8)(xnmatI(i,:),i=1,8)
      
      z11n = matmul(xnmatI,z11i)
      z21n = matmul(xnmatI,z21i)
      z31n = matmul(xnmatI,z31i)
      
      
!      write(6,*)"z11i",kint2; write(6,fmt8)z11i
!      write(6,*)"z11n"; write(6,fmt8)z11n
!      write(6,*)"z12i",kint2; write(6,fmt8)z12i
!      write(6,*)"z12n"; write(6,fmt8)z12n

      !The storage is done by row.
      do kint=1,nnodes
c         svars((kint-1)*knsdv+4) = z11n(kint)
c         svars((kint-1)*knsdv+5) = z12n(kint)
c         svars((kint-1)*knsdv+6) = z13n(kint)
         svars((kint-1)*knsdv+4) = z11i(kint)
         svars((kint-1)*knsdv+5) = z12i(kint)
         svars((kint-1)*knsdv+6) = z13i(kint)          
      end do !kint 
C     
      RETURN
      END
C
C
