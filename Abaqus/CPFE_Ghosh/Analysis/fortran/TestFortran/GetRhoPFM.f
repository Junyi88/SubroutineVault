      subroutine GetRhoPFM(RhoP,RhoF,RhoM,
     1 RhoSSD,RHOGNDs,RHOGNDet,RHOGNDen,
     2 SlipsN,SlipsS,SlipsT,
     3 Chi,Coefficient1)

C Subroutine to calculate forest parallel and mobile dislocations
      
      implicit none
      
      real*8,intent(out) :: RhoP(12),RhoF(12),RhoM(12)
      real*8,intent(in) :: RhoSSD(12),RHOGNDs(12)
      real*8,intent(in) :: RHOGNDet(12),RHOGNDen(12)
      real*8,intent(in) :: SlipsN(36),SlipsS(36),SlipsT(36)
      real*8,intent(in) :: Chi, Coefficient1
      real*8 :: CosProjection, SinProjection
      integer nPosA,nPosB, nA, nB
	  
C ------------------------------------------------------	  
      DO nA=1,12
       RhoP(nA)=0.0
       RhoF(nA)=0.0
       RhoM(nA)=0.0
       nPosA=1+3*(nA-1)
		
       DO nB=1,12     
        nPosB=1+3*(nB-1)
		
        call VectorProjections(CosProjection,SinProjection,
     1 SlipsN(nPosA),SlipsT(nPosB))
	 
        RhoF(nA)=RhoF(nA)+Chi*CosProjection*(RhoSSD(nB)+RHOGNDet(nB))
        RhoP(nA)=RhoP(nA)+Chi*SinProjection*(RhoSSD(nB)+RHOGNDet(nB))
		
        call VectorProjections(CosProjection,SinProjection,
     1 SlipsN(nPosA),SlipsS(nPosB))
        RhoF(nA)=RhoF(nA)+Chi*CosProjection*RHOGNDs(nB)
        RhoP(nA)=RhoP(nA)+Chi*SinProjection*RHOGNDs(nB)

        call VectorProjections(CosProjection,SinProjection,
     1 SlipsN(nPosA),SlipsN(nPosB))
        RhoF(nA)=RhoF(nA)+Chi*CosProjection*RHOGNDen(nB)
        RhoP(nA)=RhoP(nA)+Chi*SinProjection*RHOGNDen(nB)
		
       END DO
       RhoM(nA)=Coefficient1*sqrt(RhoF(nA)*RhoP(nA))
      END DO
	  
      return
      end subroutine GetRhoPFM