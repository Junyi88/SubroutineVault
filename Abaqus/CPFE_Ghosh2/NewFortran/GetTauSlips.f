      subroutine GetTauSlips(RhoP,RhoF,RhoM,
     1 TauPass, TauCut, V0,	  
     2 Cin)

C Subroutine to calculate forest parallel and mobile dislocations
      
      implicit none
      
      real*8,intent(in) :: RhoP(18),RhoF(18),RhoM(18)
      real*8,intent(out) ::  TauPass(18), TauCut(18), V0(18)

      real*8,intent(in) :: Cin(3)

      integer ISLIPS

C ------------------------------------------------------	
      DO ISLIPS=1,18
       TauPass(ISLIPS)= Cin(1)*sqrt(RhoP(ISLIPS)+RhoM(ISLIPS))
       TauCut(ISLIPS)=Cin(2)*sqrt(RhoF(ISLIPS))
       V0(ISLIPS)=Cin(3)/sqrt(RhoP(ISLIPS)*RhoF(ISLIPS))
      END DO

      return
      end subroutine GetTauSlips