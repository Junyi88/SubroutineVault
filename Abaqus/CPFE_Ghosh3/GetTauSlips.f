      subroutine GetTauSlips(RhoP,RhoF,RhoM,
     1 TauPass, TauCut, 
     2 Cin)

C Subroutine to calculate forest parallel and mobile dislocations
      
      IMPLICIT NONE  
      
      REAL*8, INTENT(IN) :: RhoP(18),RhoF(18),RhoM(18)
      REAL*8, INTENT(OUT) :: TauPass(18), TauCut(18)

      REAL*8, INTENT(IN) :: Cin(2)
      INTEGER :: ISLIPS

C ------------------------------------------------------	
      DO ISLIPS=1,18
       TauPass(ISLIPS)= Cin(1)*sqrt(RhoP(ISLIPS)+RhoM(ISLIPS))
       TauCut(ISLIPS)=Cin(2)*sqrt(RhoF(ISLIPS))
      END DO

      return
      end subroutine GetTauSlips