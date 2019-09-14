      subroutine GetRhoSSDEvolve(Tau, TauPass, TauCut, V0, RhoM, 
     1 GammaDot, TauEff, SSDDot, RhoSSD, RhoF,   
     2 CinS)

C Subroutine to calculate forest parallel and mobile dislocations
      
      implicit none
      
      real*8,intent(in) :: TauEff(18), Tau(18), TauPass(18) 
      real*8,intent(in) :: TauCut(18), V0(18)
      real*8,intent(in) :: RhoM(18), RhoSSD(18), RhoF(18)
      real*8,intent(in) :: GammaDot(18)
      real*8,intent(out) :: SSDDot(18)

      real*8,intent(in) :: CinS(6)

      integer ISLIPS

C ------------------------------------------------------	
      DO ISLIPS=1,18
      SSDDot(ISLIPS)=0.0
      ENDDO

      DO ISLIPS=1,18
       
      IF (GAMMADOT(ISLIPS).GT.0.0) THEN
	  
      SSDDot(ISLIPS)=SSDDot(ISLIPS)+CinS(1)*sqrt(RhoF(ISLIPS))
      IF (TauEff(ISLIPS).GT.10.0) THEN
      SSDDot(ISLIPS)=SSDDot(ISLIPS)+CinS(2)*RhoM(ISLIPS)/TauEff(ISLIPS)
      ELSE
      SSDDot(ISLIPS)=SSDDot(ISLIPS)+CinS(2)*RhoM(ISLIPS)*0.1
      END IF  
      SSDDot(ISLIPS)=SSDDot(ISLIPS)-CinS(3)*RhoSSD(ISLIPS)
      SSDDot(ISLIPS)=SSDDot(ISLIPS)*abs(GAMMADOT(ISLIPS))
	  
c      SSDDot(ISLIPS)=SSDDot(ISLIPS)-CinS(4)*(RhoSSD(ISLIPS)**2.0)*
c     1 abs(TAU(ISLIPS))*((abs(GAMMADOT(ISLIPS))/CinS(6))**CinS(5))	  
      ELSE
      SSDDot(ISLIPS)=0.0
      END IF  
      END DO

      return
      end subroutine GetRhoSSDEvolve