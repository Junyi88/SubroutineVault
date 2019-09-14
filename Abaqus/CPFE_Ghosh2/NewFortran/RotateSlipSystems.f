      subroutine RotateSlipSystems(GammaDot,dTIME,DSTRAN, SPIN_TENSOR,
     1 SLIP_S,SLIP_N,
     2 SLIP_SPE,SLIP_NPE,
     2 SLIP_SSE,SLIP_NSE,
     2 SLIP_SCB,SLIP_NCB)

C Subroutine to calculate forest parallel and mobile dislocations
      
      implicit none

      include 'ComplianceTensorMaps.f' 	 
	 
      real*8,intent(in) :: dTIME 
      real*8,intent(in) :: GammaDot(18),DStran(6), SPIN_TENSOR(3,3)
      real*8,intent(inout) :: SLIP_S(54),SLIP_N(54)
      real*8,intent(inout) :: SLIP_SPE(36),SLIP_NPE(36)
      real*8,intent(inout) :: SLIP_SSE(36),SLIP_NSE(36)
      real*8,intent(inout) :: SLIP_SCB(36),SLIP_NCB(36)	  
	  
      real*8:: DGA(18)
      real*8:: FLIPS(3,3,18)
      real*8:: DS(54),DN(54)
	  
      integer ISLIPS, IVAL, ISYS, I, J, K, L, ICOR

      Integer, PARAMETER::  FULL2VOIGT(3,3) =
     1 reshape([
     1  1, 4, 5, 
     1  4, 2, 6,
     1  5, 6, 3	 
     3  ], [3,3]
     4 )

      Integer, PARAMETER::  VOIGT2FULL(2,6) =
     1 reshape([
     1  1, 1,
     1  2, 2,
     1  3, 3,
     1  1, 2,
     1  1, 3,
     1  2, 3	 
     3  ], [2,6]
     4 )	 

      Integer, PARAMETER::  VOIGT2List(6,6) =
     1 reshape([
     1  1, 7, 8, 9 ,10, 11,
     2  7, 2, 12, 13 ,14, 15,
     3  8, 12, 3, 16 ,17, 18,
     4  9, 13, 16, 4 ,19, 20,
     5  10, 14, 17, 19 ,5, 21,
     6  11, 15, 18, 20 ,21, 6	 
     3  ], [6,6]
     4 )
	   
      DO ISLIPS=1,18
	        DGA(ISLIPS)=GammaDot(ISLIPS)*DTIME
      END DO

c-----------------------------------------
      DO ISLIPS=1,18
        DO I=1,3
        DO J=1,3
	        FLIPS(I,J,ISLIPS)=DSTRAN(FULL2VOIGT(I,J))+SPIN_TENSOR(I,J)
        END DO
        END DO
      END DO
	  
      DO ISLIPS=1,18
      DO ICOR=1,18
        IVAL=(ICOR-1)*3	  
        DO I=1,3
        DO J=1,3
	        FLIPS(I,J,ISLIPS)=FLIPS(I,J,ISLIPS)-
     1      SLIP_S(IVAL+I)*SLIP_N(IVAL+J)*DGA(ICOR)
        END DO
        END DO
      END DO	
      END DO	  
	  
C ----------------------------------------
c Calculate and update DS
      DO ISLIPS=1,54
	        DS(ISLIPS)=0.0
	        DN(ISLIPS)=0.0
      END DO	
      DO ISLIPS=1,18

        IVAL=(ISLIPS-1)*3	  
        DO I=1,3
        DO J=1,3
	        DS(IVAL+I)=DS(IVAL+I)+FLIPS(I,J,ISLIPS)*SLIP_S(IVAL+J)
	        DN(IVAL+J)=DN(IVAL+J)-FLIPS(J,I,ISLIPS)*SLIP_N(IVAL+J)
        END DO
        END DO
	
      END DO	
      DO ISLIPS=1,54
        SLIP_S(ISLIPS)=SLIP_S(ISLIPS)+DS(ISLIPS)
        SLIP_N(ISLIPS)=SLIP_N(ISLIPS)+DN(ISLIPS)       
	
      END DO	  

C ----------------------------------------
c Calculate and update PE
      DO ISLIPS=1,54
	        DS(ISLIPS)=0.0
	        DN(ISLIPS)=0.0
      END DO	
      DO ISLIPS=1,12

        IVAL=(ISLIPS-1)*3	  
        DO I=1,3
        DO J=1,3
	        DS(IVAL+I)=DS(IVAL+I)+FLIPS(I,J,ISLIPS)*SLIP_SPE(IVAL+J)
	        DN(IVAL+J)=DN(IVAL+J)-FLIPS(J,I,ISLIPS)*SLIP_NPE(IVAL+J)
        END DO
        END DO
	
      END DO	
      DO ISLIPS=1,36
        SLIP_SPE(ISLIPS)=SLIP_SPE(ISLIPS)+DS(ISLIPS)
        SLIP_NPE(ISLIPS)=SLIP_NPE(ISLIPS)+DN(ISLIPS)       
	
      END DO		  

C ----------------------------------------
c Calculate and update SE
      DO ISLIPS=1,54
	        DS(ISLIPS)=0.0
	        DN(ISLIPS)=0.0
      END DO	
      DO ISLIPS=1,12

        IVAL=(ISLIPS-1)*3	  
        DO I=1,3
        DO J=1,3
	        DS(IVAL+I)=DS(IVAL+I)+FLIPS(I,J,ISLIPS)*SLIP_SSE(IVAL+J)
	        DN(IVAL+J)=DN(IVAL+J)-FLIPS(J,I,ISLIPS)*SLIP_NSE(IVAL+J)
        END DO
        END DO
	
      END DO	
      DO ISLIPS=1,36
        SLIP_SSE(ISLIPS)=SLIP_SSE(ISLIPS)+DS(ISLIPS)
        SLIP_NSE(ISLIPS)=SLIP_NSE(ISLIPS)+DN(ISLIPS)       
	
      END DO	
	  
C ----------------------------------------
c Calculate and update CB
      DO ISLIPS=1,54
	        DS(ISLIPS)=0.0
	        DN(ISLIPS)=0.0
      END DO	
      DO ISLIPS=1,12

        IVAL=(ISLIPS-1)*3	  
        DO I=1,3
        DO J=1,3
	        DS(IVAL+I)=DS(IVAL+I)+FLIPS(I,J,ISLIPS)*SLIP_SCB(IVAL+J)
	        DN(IVAL+J)=DN(IVAL+J)-FLIPS(J,I,ISLIPS)*SLIP_NCB(IVAL+J)
        END DO
        END DO
	
      END DO	
      DO ISLIPS=1,36
        SLIP_SCB(ISLIPS)=SLIP_SCB(ISLIPS)+DS(ISLIPS)
        SLIP_NCB(ISLIPS)=SLIP_NCB(ISLIPS)+DN(ISLIPS)       
	
      END DO	
	  
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      return
      end subroutine RotateSlipSystems