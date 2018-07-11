      subroutine GetDSTRESSFP(DStress,GammaDot,DStran,Stress,dTIME, 
     1 SLIP_S,SLIP_N, dFP, FP,	   
     2 CinS)

C Subroutine to calculate forest parallel and mobile dislocations
      
      implicit none

      include 'ComplianceTensorMaps.f' 	 
	 
      real*8,intent(in) :: dTIME 
      real*8,intent(in) :: GammaDot(18),DStran(6),Stress(6)
      real*8,intent(in) :: SLIP_S(54),SLIP_N(54),FP(9)

      real*8,intent(out) :: DStress(6),dFP(9)

      real*8,intent(in) :: CinS(21)
      
      real*8:: HYDROSTRAIN, DGA(18), DUM1, MIU,OHM1,OHM2, LP(9)
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
	 
      HYDROSTRAIN=DSTRAN(1)+DSTRAN(2)+DSTRAN(3)
	  
	  
      DO ISLIPS=1,18
	        DGA(ISLIPS)=GammaDot(ISLIPS)*DTIME
      END DO

      DO ISLIPS=1,9
	        dFP(ISLIPS)=0.0
			LP(ISLIPS)=0.0
      END DO	  

C XXX	  

c      DO ISYS=1,6
c        DStress(ISYS)=0.0
c        I = VOIGT2FULL(1,ISYS)
c        J = VOIGT2FULL(2,ISYS)
c        DO K=1,3
c        DO L=1,3
c        IF (K.EQ.L) THEN
c          IVAL=FULL2VOIGT(K,L)
c          DStress(ISYS)=DStress(ISYS)+
c     1	     CinS(MFULL2LIST(I,J,K,L))*
c     1	     dSTRAN(IVAL)	 
c        else
c          IVAL=FULL2VOIGT(K,L)
c          DStress(ISYS)=DStress(ISYS)+
c     1	     0.5*CinS(MFULL2LIST(I,J,K,L))*
c     1	     dSTRAN(IVAL)			
c        END IF
c        END DO
c        END DO		
c      END DO

      DO ISYS=1,6
        DStress(ISYS)=0.0
      END DO	
	  
      DO ISYS=1,6
c        DStress(ISYS)=-STRESS(ISYS)*HYDROSTRAIN
        DO IVAL=1,6
           DStress(ISYS)=DStress(ISYS)+
     1     CINS(VOIGT2List(ISYS,IVAL))*
     1     DSTRAN(IVAL)	 
        END DO		
      END DO	
C ----------------------------------------
      DO ISYS=1,6
        I = VOIGT2FULL(1,ISYS)
        J = VOIGT2FULL(2,ISYS)
        DO ISLIPS=1,18
        ICOR=(ISLIPS-1)*3
        DO K=1,3
        DO L=1,3
          IVAL=FULL2VOIGT(K,L)
          MIU=0.5*(SLIP_S(K+ICOR)*SLIP_N(L+ICOR)+
     1	     SLIP_S(L+ICOR)*SLIP_N(K+ICOR))

          DStress(ISYS)=DStress(ISYS)-
     1	     CinS(MFULL2LIST(I,J,K,L))*
     1       (MIU)*
     1	     DGA(ISLIPS)	 

c -------------------------------------------------	 
c        IF (K.NE.L) THEN
c          DStress(ISYS)=DStress(ISYS)-
c     1	     CinS(MFULL2LIST(I,J,K,L))*
c     1       (MIU)*
c     1	     DGA(ISLIPS)	 
c        ELSE
c          DStress(ISYS)=DStress(ISYS)-
c     1	     CinS(MFULL2LIST(I,J,K,L))*
c     1       (MIU)*
c     1	     DGA(ISLIPS)/2.0	 		
		
c        END IF
        END DO
        END DO		
        END DO		
      END DO
c--------------------------------------
      DO ISLIPS=1,18
        ICOR=(ISLIPS-1)*3
        DO I=1,3
        DO J=1,3
        IVAL=3*(J-1)+I
        LP(IVAL)=LP(IVAL)+SLIP_S(I+ICOR)*SLIP_N(J+ICOR)*DGA(ISLIPS)
        END DO
        END DO		  
	        
      END DO	  


        DO I=1,3
        DO J=1,3
        DO K=1,3
		ISYS=3*(J-1)+I
        IVAL=3*(K-1)+I
		ICOR=3*(J-1)+K
        DFP(ISYS)=DFP(ISYS)+LP(IVAL)*FP(ICOR)
        END DO
        END DO		  
        END DO		         
	  
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      return
      end subroutine GetDSTRESSFP