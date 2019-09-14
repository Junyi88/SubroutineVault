      subroutine GetGammaDot(Tau, TauPass, TauCut, V0, RhoM, 
     1 Vs, GammaDot, TauEff, TAUC,	   
     2 CinS)

C Subroutine to calculate forest parallel and mobile dislocations
      
      implicit none
      
      real*8,intent(in) :: Tau(18), TauPass(18), TauCut(18), 
     1	V0(18), RhoM(18), TAUC(18)
      real*8,intent(out) :: Vs(18), GammaDot(18), TauEff(18)

      real*8,intent(in) :: CinS(3)

      integer ISLIPS

C ------------------------------------------------------	
c      DO ISLIPS=1,18
c       Tau(ISLIPS)=0.0 
c       TauPass(ISLIPS)=0.0 
c       TauCut(ISLIPS)=0.0 
c      END DO	  

C ------------------------------------------------------	
      DO ISLIPS=1,18
       TauEff(ISLIPS)=abs(Tau(ISLIPS))-TauPass(ISLIPS)
      IF (TauEff(ISLIPS).LE.0.0) THEN
       TauEff(ISLIPS)=0.0
       Vs(ISLIPS)=0.0
       GAMMADOT(ISLIPS)=0.0	   
      ELSE
       IF (TauEff(ISLIPS).GE.TauC(ISLIPS)) THEN
       Vs(ISLIPS)=V0(ISLIPS)*CINS(1)*
     1	SINH((TauEff(ISLIPS)/TauCut(ISLIPS))**CINS(2))*
     1	SIGN(1.0,TAU(ISLIPS))
	 
       GAMMADOT(ISLIPS)=RhoM(ISLIPS)*CINS(3)*Vs(ISLIPS)
       END IF
      END IF
      END DO
	  
C ------------------------------------------------------	
  
	  
      return
      end subroutine GetGammaDot