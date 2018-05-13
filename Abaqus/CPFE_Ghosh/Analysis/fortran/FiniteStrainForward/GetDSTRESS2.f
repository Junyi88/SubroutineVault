      subroutine GetDSTRESS(DStress,GammaDot,DStran,Stress,dTIME, 
     1 SLIP_S,SLIP_N,	   
     2 CinS)

C Subroutine to calculate forest parallel and mobile dislocations
      
      implicit none

      include 'ComplianceTensorMaps.f' 	 
	 
      real*8,intent(in) :: dTIME 
      real*8,intent(in) :: GammaDot(18),DStran(6),Stress(6)
      real*8,intent(in) :: SLIP_S(54),SLIP_N(54)

      real*8,intent(out) :: DStress(6)

      real*8,intent(in) :: CinS(21)
      
      real*8:: HYDROSTRAIN, DGA(18), DUM1
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
	  
      HYDROSTRAIN=DSTRAN(1)+DSTRAN(2)+DSTRAN(3)
	  
	  
      DO ISLIPS=1,18
	        DGA(ISLIPS)=GammaDot(ISLIPS)*DTIME
      END DO
C XXX	  

      DO ISYS=1,6
        DStress(ISYS)=0.0
        I = VOIGT2FULL(1,ISYS)
        J = VOIGT2FULL(2,ISYS)
        DO K=1,3
        DO L=1,3
          IVAL=FULL2VOIGT(K,L)
          DStress(ISYS)=DStress(ISYS)+
     1	     CinS(MFULL2LIST(I,J,K,L))*
     1	     dSTRAN(IVAL)	 
        END DO
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
          DStress(ISYS)=DStress(ISYS)-
     1	     CinS(MFULL2LIST(I,J,K,L))*
     1       (SLIP_S(K+ICOR)*SLIP_N(L*ICOR))*
     1	     DGA(IVAL)	 
        END DO
        END DO		
        END DO		
      END DO
	  

C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      return
      end subroutine GetDSTRESS