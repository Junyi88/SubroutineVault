      subroutine GetRhoPFM(RhoP,RhoF,RhoM,
     1 RhoSSD,
     2 SLIP_N,SLIP_T,
     5 Coefficient1)

C Subroutine to calculate forest parallel and mobile dislocations
      
      implicit none
      
      real*8,intent(out) :: RhoP(18),RhoF(18),RhoM(18)
      real*8,intent(in) :: RhoSSD(18)
      real*8,intent(in) :: SLIP_N(54),SLIP_T(54)
      real*8,intent(in) :: Coefficient1
      real*8 :: CosProjection, SinProjection
      integer nA, nB, nAX, nStartA, nStartB, NFinA, NFinB

C ------------------------------------------------------	
      DO nA=1,18
       RhoP(nA)=0.0
       RhoF(nA)=0.0
      END DO

C=======================  
      DO nA=1,18	
        nSTARTA=3*(nA-1)+1
        nFinA=3*nA
		
       DO nB=1,18    
        nSTARTB=3*(nB-1)+1
        nFinB=3*nB		
		
        call VectorProjections(CosProjection,SinProjection,
     1 SLIP_N(nSTARTA:nFinA),SLIP_T(NSTARTB:NFINB))
	 
        RhoF(nA)=RhoF(nA)+CosProjection*(RhoSSD(nB))
        RhoP(nA)=RhoP(nA)+SinProjection*(RhoSSD(nB))
       END DO
	   
      END DO
	
C-------------------	
      DO nA=1,18
       RhoM(nA)=Coefficient1*sqrt(RhoF(nA)*RhoP(nA))
      END DO
      return
      end subroutine GetRhoPFM