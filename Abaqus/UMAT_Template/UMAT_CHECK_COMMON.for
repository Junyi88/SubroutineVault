      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C


      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
C     ELASTIC USER SUBROUTINE
      PARAMETER (ONE=1.0D0, TWO=2.0D0)
c    -------------------- 	  
      integer, parameter :: TOTALELEMENTNUM = 4
      integer, parameter :: TOTALINT = 8
      integer, parameter :: TOTALSPACE = 32

      Integer :: I,J, NTest
      real*8 :: kgausscoords, kpoints
      COMMON/UMPS/kgausscoords(TOTALELEMENTNUM,TOTALINT,3),
     1  kpoints(TOTALELEMENTNUM,TOTALINT,3)

c     -------------------------------	  
      E=PROPS(1)
      ANU=PROPS(2)
      ALAMBDA=E/(ONE+ANU)/(ONE-TWO*ANU)
      BLAMBDA=(ONE-ANU)
         CLAMBDA=(ONE-TWO*ANU)
	     DO I=1,NTENS
	      DO J=1,NTENS
	      DDSDDE(I,J)=0.0D0
	      ENDDO
	     ENDDO
      DDSDDE(1,1)=(ALAMBDA*BLAMBDA)
      DDSDDE(2,2)=(ALAMBDA*BLAMBDA)
      DDSDDE(3,3)=(ALAMBDA*BLAMBDA)
      DDSDDE(4,4)=(ALAMBDA*CLAMBDA)
      DDSDDE(5,5)=(ALAMBDA*CLAMBDA)
      DDSDDE(6,6)=(ALAMBDA*CLAMBDA)
      DDSDDE(1,2)=(ALAMBDA*ANU)
      DDSDDE(1,3)=(ALAMBDA*ANU)
      DDSDDE(2,3)=(ALAMBDA*ANU)
      DDSDDE(2,1)=(ALAMBDA*ANU)
      DDSDDE(3,1)=(ALAMBDA*ANU)
      DDSDDE(3,2)=(ALAMBDA*ANU)
         DO I=1,NTENS
	      DO J=1,NTENS
	      STRESS(I)=STRESS(I)+DDSDDE(I,J)*DSTRAN(J)
	      ENDDO
	     ENDDO
		 
c      -----------
c      statev(1:32) = kgausscoordsX
c      statev(33:64) = kgausscoordsY
c      statev(65:96) = kgausscoordsZ
c      statev(97:128) = NEL
c      statev(129:160) = NPT
c      statev(161:192) = 1

c      STATEV(193:195) = coordinates
c      statev(196) = NEL
c      statev(197) = NPT

      IF ((KINC.LE.1)) THEN
      call MutexLock( 1 )      ! lock Mutex #2      
      ! use original co-ordinates X     
      do i =1,3
          kgausscoords(noel,npt,i) = coords(i)
      end do
          kpoints(noel,npt,1) = noel
          kpoints(noel,npt,2) = npt
		  kpoints(noel,npt,3) = 1.0
      call MutexUnlock( 1 )   ! unlock Mutex #2
      ENDIF
c----------------------------------------------------------
      IF (npt == 8 ) THEN    
c	   call MutexLock( 2 ) 
	     DO I=1,TOTALELEMENTNUM
	      DO J=1,TOTALINT
	        kpoints(I,J,3) = 2.0
	      ENDDO
	     ENDDO    	  
c       call MutexUnlock( 2 ) 
      ENDIF
c ---------------------------------------------------------
	   call MutexLock( 3 ) 
	     DO I=1,TOTALELEMENTNUM
	      DO J=1,TOTALINT
		     NTEST = ((I-1)*8) + J
			 statev(NTEST) = kgausscoords(I,J,1)
			 statev(NTEST+TOTALSPACE) = kgausscoords(I,J,2)
			 statev(NTEST+2*TOTALSPACE) = kgausscoords(I,J,3)
			 statev(NTEST+3*TOTALSPACE) = kpoints(I,J,1)
			 statev(NTEST+4*TOTALSPACE) = kpoints(I,J,2)
			 statev(NTEST+5*TOTALSPACE) = kpoints(I,J,3)
			 
	      ENDDO
	     ENDDO       
       call MutexUnlock( 3 ) 
		 
	     statev(1+6*TOTALSPACE) = coords(1)
		 statev(2+6*TOTALSPACE) = coords(2)
		 statev(3+6*TOTALSPACE) = coords(3)
		 statev(4+6*TOTALSPACE) = NOEL	 
		 statev(5+6*TOTALSPACE) = NPT		 
      RETURN
      END
	  
      include 'uexternaldb.f'