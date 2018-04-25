      subroutine FlowRule(GammaRate, Tau,
     1 RhoP, RhoF, RhoM, TauPass,TauCut,
     2 tauC, rhoC,
     3 bLambdaNu,Coefficient1)
	 
	 
C Subroutine to calculate forest parallel and mobile dislocations
      
      implicit none
      real*8,intent(out) :: GammaRate(12)
      real*8,intent(inout) :: TauCut(12)
      real*8,intent(in) :: TauPass(12),TauCut(12)
      real*8,intent(in) :: RhoP(12),RhoF(12),RhoM(12)
      real*8,intent(in) :: tauC, rhoC
      real*8,intent(in) :: bLambdaNu,Coefficient1
	  
      integer nI
	  
C ------------------------------------------------------	  
      DO nI=1,12
        IF (RhoM(nI).GE.rhoC) THEN
         IF ((Tau(nI).GE.tauC).AND.(Tau(nI).GE.TauPass(nI))) THEN
         GAMMARATE(nI)=RhoM(nI)*bLambdaNu*SIGN(1.0,TAU(nI))*exp(
     1      -Coefficient1*(1.0-(Tau(nI)-TauPass(nI))/TauCut(nI)))
         ELSE
         GAMMARATE(nI)=RhoM(nI)*bLambdaNu*SIGN(1.0,TAU(nI))*exp(
     1      -Coefficient1)		 
         ENDIF
        ELSE
         GAMMARATE(nI)=0.0	
        END IF
	   
	   
      END DO
	  
      return
      end subroutine FlowRule