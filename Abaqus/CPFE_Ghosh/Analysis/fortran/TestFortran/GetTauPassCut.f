      subroutine GetTauPassCut(TauPass,TauCut,
     1 RhoP, RhoF, RhoM,
     2 CoefficientPass,CoefficientCut)

C Subroutine to calculate forest parallel and mobile dislocations
      
      implicit none
      real*8,intent(out) :: TauPass(12),TauCut(12)
	  
      real*8,intent(in) :: RhoP(12),RhoF(12),RhoM(12)
      real*8,intent(in) :: CoefficientPass,CoefficientCut
	  
      integer nI
	  
C ------------------------------------------------------	  
      DO nI=1,12
       TauPass(nI)=CoefficientPass*sqrt(RhoP(nI)+RhoM(nI))
       TauCut(nI)=CoefficientCut*sqrt(RhoF(nI))
      END DO
	  
      return
      end subroutine GetTauPassCut