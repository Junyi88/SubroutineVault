C     STRAIN GRADIENTS
C *******************************************************************************
C *                     Compute Curl of a 2nd order tensor                      *
C *              As a vector, curl of row written in corresponding column       *
C *                          Full integration                                   *
C *******************************************************************************
       SUBROUTINE kcurl(svars,xnat8,gauss,gausscoords) 

       INCLUDE 'ABA_PARAM.INC'

      integer, parameter:: K=3,M=3,N=8,nnodes=8, knsdv=18
      real*8,parameter  :: zero=1.0e-6,xgauss = 0.577350269189626
      integer, parameter :: TOTALELEMENTNUM=100
      !scalars
      integer:: GG=0
      
      !arrays
      real*8,intent(inout):: svars(144)
      real*8,intent(in) :: xnat8(nnodes,3),gauss(nnodes,3),
     + gausscoords(3,nnodes)

      real*8 :: xj(3,3),xjinv(3,3),dndloc(nnodes,3),dndx(nnodes,3),
     + fnode(3,nnodes), fmat1(3,3),dmout(3,3),xn(nnodes)
     
      real*8 :: TEMP1(3,3)
      INTEGER :: NNN, NI,NJ,NG
     
      !Extrapolation arrays
      real*8 :: z11i(nnodes),z11n(nnodes),z12i(nnodes),z12n(nnodes),
     + z13i(nnodes),z13n(nnodes),z21i(nnodes),z21n(nnodes),z22i(nnodes),
     + z22n(nnodes),z23i(nnodes),z23n(nnodes),z31i(nnodes),z31n(nnodes),
     + z32i(nnodes),z32n(nnodes),z33i(nnodes),z33n(nnodes),
     + xnmat(nnodes,nnodes),xnmatI(nnodes,nnodes)
     
      character(len=*),parameter :: fmt20 = "(' ',20(F6.4,1X))",
     + fmt3="(3(' ',(ES11.3,1X)))",fmt6 = "(' ',6(F12.3,1X))",
     + fmt5 = "(' ',5(F12.8,1X))",fmt8="(8(' ',(ES11.3,1X)))"

      real*8 :: kgausscoords, kFp, kcurlFp, kDGA, kX
      COMMON/UMPS/kgausscoords(TOTALELEMENTNUM,8,3),
     1 kFp(TOTALELEMENTNUM,8, 9),
     1 kcurlFp(TOTALELEMENTNUM, 8, 9), 
     1 kDGA(TOTALELEMENTNUM, 8, 9),
     1 kX(TOTALELEMENTNUM, 8, 9)     
     
      fnode = 0.;fmat1 = 0.;xnmat=0.;xnmatI=0.
      z11i=0.;z11n=0.;z12i=0.;z12n=0.;z13i=0.;z13n=0.
      z21i=0.;z21n=0.;z22i=0.;z22n=0.;z23i=0.;z23n=0.
      z31i=0.;z31n=0.;z32i=0.;z32n=0.;z33i=0.;z33n=0.
C
C   USE EIGHT GAUSS POINT FOR GRADIENTS
C
C    LOOP OVER EIGHT INTEGRATION POINTS

!    Evaluate curlT at the integration points, and simultaneously populate xnmat which will be later inverted for extrapolation purposes.   
c      write(6,*) "kfp ------------------"
c      do kintB = 1,8 !integration points are our nodes now
c      write(6,*) "kintB = " , kintB
c      write(6,*) svars((kintB-1)*knsdv+1), svars((kintB-1)*knsdv+2),
c     1   svars((kintB-1)*knsdv+3)
c      write(6,*) svars((kintB-1)*knsdv+4), svars((kintB-1)*knsdv+5),
c     1   svars((kintB-1)*knsdv+6)	 
c      write(6,*) svars((kintB-1)*knsdv+7), svars((kintB-1)*knsdv+8),
c     1   svars((kintB-1)*knsdv+9)
	 
c      end do
      write(6,*) "XX+++++++++++++++++++++++++++++++++++++++++"	  
      write(6,*) "xNAT (8,3) "
      write(6,*) xnat8(1,:)
      write(6,*) xnat8(2,:)
      write(6,*) xnat8(3,:)
      write(6,*) xnat8(4,:)    
      write(6,*) xnat8(5,:)
      write(6,*) xnat8(6,:)
      write(6,*) xnat8(7,:)
      write(6,*) xnat8(8,:)        
      write(6,*) "   "      
      
      write(6,*) "Gauss (3,8) "
      write(6,*) gauss(1,:)
      write(6,*) gauss(2,:)
      write(6,*) gauss(3,:)
      write(6,*) gauss(4,:)
      write(6,*) gauss(5,:)
      write(6,*) gauss(6,:)
      write(6,*) gauss(7,:)
      write(6,*) gauss(8,:)
      write(6,*) "   "        
      
      write(6,*) "gausscoords(3,8) "
      write(6,*) gausscoords(:,1)
      write(6,*) gausscoords(:,2)
      write(6,*) gausscoords(:,3)
      write(6,*) gausscoords(:,4)
      write(6,*) gausscoords(:,5)
      write(6,*) gausscoords(:,6)
      write(6,*) gausscoords(:,7)
      write(6,*) gausscoords(:,8)
      write(6,*) "   "    

      do NG = 1,8
        do NI = 1,3              
        do NJ = 1,3  
           NNN = NI + (NJ-1)*3
           NNN = NNN + 18*(NG-1)
           TEMP1(NI,NJ) = svars(NNN)
        end do       
        end do       
        
      write(6,*) "kFP = ", NG      
      write(6,*) TEMP1(:,1)
      write(6,*) TEMP1(:,2)
      write(6,*) TEMP1(:,3)
      write(6,*) "   "    
      
        
      end do      
      
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
         dmout = 0.0
         write(6,*)  "DET is ZERO"
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
      fnode(2,kintB) = svars((kintB-1)*knsdv+2)
      fnode(3,kintB) = svars((kintB-1)*knsdv+3)
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
      fnode(1,kintB) = svars((kintB-1)*knsdv+4)
      fnode(2,kintB) = svars((kintB-1)*knsdv+5)
      fnode(3,kintB) = svars((kintB-1)*knsdv+6)
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
      fnode(1,kintB) = svars((kintB-1)*knsdv+7)
      fnode(2,kintB) = svars((kintB-1)*knsdv+8)
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
c ====   
        
C       write(6,*) "curlFP = ", kintB2     
C       write(6,*) z11i(kintB2), z12i(kintB2), z13i(kintB2)
C       write(6,*) z21i(kintB2), z22i(kintB2), z23i(kintB2)
C       write(6,*) z31i(kintB2), z32i(kintB2), z33i(kintB2)
C       write(6,*) "   "   
      
!  write curlfp into svars (for each integration point!)
! Assigning the values of the inner element's gauss point to the corresponding outer element's gauss point! This can be debated later!
!       do i=1,3
!        do j=1,3 
!        svars((kintB2-1)*knsdv+37+j+(i-1)*3) = dmout(i,j)
!        end do
!       end do
      
      end do !kintB2
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
      
!      write(6,*)"z11i",kintB2; write(6,fmt8)z11i
!      write(6,*)"z11n"; write(6,fmt8)z11n
!      write(6,*)"z12i",kintB2; write(6,fmt8)z12i
!      write(6,*)"z12n"; write(6,fmt8)z12n

      !The storage is done by row.
      do kintB=1,nnodes
         GG = (kintB-1)*knsdv+9
         svars(GG+1) = z11i(kintB)
         svars(GG+2) = z12i(kintB)
         svars(GG+3) = z13i(kintB)
         
         svars(GG+4) = z21i(kintB)
         svars(GG+5) = z22i(kintB)
         svars(GG+6) = z23i(kintB)
         
         svars(GG+7) = z31i(kintB)
         svars(GG+8) = z32i(kintB)
         svars(GG+9) = z33i(kintB)      
         
         
         write(6,*) "curlFP2 = ", kintB    
         write(6,*) z11i(kintB), z12i(kintB), z13i(kintB)
         write(6,*) z21i(kintB), z22i(kintB), z23i(kintB)
         write(6,*) z31i(kintB), z32i(kintB), z33i(kintB)


        
         write(6,*) "--"     
         write(6,*) GG+1, GG+2, GG+3
         write(6,*) GG+4, GG+5, GG+6
         write(6,*) GG+7, GG+8, GG+9
         write(6,*) "   "            
         
      end do !kintB 
      
      call MutexLock( 1 )      ! lock Mutex #1 

      do kintB=1,nnodes
         GG = (kintB-1)*knsdv+9
         KX(1,kintB,1) = z11i(kintB)
         KX(1,kintB,2) = z12i(kintB)
         KX(1,kintB,3) = z13i(kintB)
         
         KX(1,kintB,4) = z21i(kintB)
         KX(1,kintB,5) = z22i(kintB)
         KX(1,kintB,6) = z23i(kintB)
         
         KX(1,kintB,7) = z31i(kintB)
         KX(1,kintB,8) = z32i(kintB)
         KX(1,kintB,9) = z33i(kintB)      
         
         
      end do !kintB 
      
      call MutexUnlock( 1 )      ! lock Mutex #1 	  
c      do kintB=1,nnodes
c         svars((kintB-1)*knsdv+10) = svars((kintB-1)*knsdv+1)
c         svars((kintB-1)*knsdv+11) = svars((kintB-1)*knsdv+2)
c         svars((kintB-1)*knsdv+12) = svars((kintB-1)*knsdv+3)
         
c         svars((kintB-1)*knsdv+13) = svars((kintB-1)*knsdv+4)
c         svars((kintB-1)*knsdv+14) = svars((kintB-1)*knsdv+5)
c         svars((kintB-1)*knsdv+15) = svars((kintB-1)*knsdv+6)
         
c         svars((kintB-1)*knsdv+16) = svars((kintB-1)*knsdv+7)
c         svars((kintB-1)*knsdv+17) = svars((kintB-1)*knsdv+8)
c         svars((kintB-1)*knsdv+18) = svars((kintB-1)*knsdv+9)
c      end do !kintB 	  
C     
c      write(6,*) "kcurlfp ------------------"
c      do kintB = 1,8 !integration points are our nodes now
c      write(6,*) "kintB = " , kintB
c      write(6,*) svars((kintB-1)*knsdv+10), svars((kintB-1)*knsdv+11),
c     1   svars((kintB-1)*knsdv+12)
c      write(6,*) svars((kintB-1)*knsdv+13), svars((kintB-1)*knsdv+14),
c     1   svars((kintB-1)*knsdv+15)	 
c      write(6,*) svars((kintB-1)*knsdv+16), svars((kintB-1)*knsdv+17),
c     1   svars((kintB-1)*knsdv+18)
	 
c      end do

      RETURN
      END
C
C
