      subroutine GetRhoPFM(RhoP,RhoF,RhoM,
     1 RhoSSD,
     2 FCC_N,FCC_T,
     4 CUBIC_N,CUBIC_T,	 
     5 Coefficient1)

C Subroutine to calculate forest parallel and mobile dislocations
      
      implicit none
      
      real*8,intent(out) :: RhoP(18),RhoF(18),RhoM(18)
      real*8,intent(in) :: RhoSSD(18)
      real*8,intent(in) :: FCC_N(12,3),FCC_T(12,3)
      real*8,intent(in) :: CUBIC_N(6,3),CUBIC_T(6,3)
      real*8,intent(in) :: Coefficient1
      real*8 :: CosProjection, SinProjection
      integer nA, nB, nAX

C ------------------------------------------------------	
      DO nA=1,18
       RhoP(nA)=0.0
       RhoF(nA)=0.0
      END DO

C=======================  
      DO nA=1,12	
       DO nB=1,12     
		
        call VectorProjections(CosProjection,SinProjection,
     1 FCC_N(nA,1:3),FCC_T(nB,1:3))
	 
        RhoF(nA)=RhoF(nA)+CosProjection*(RhoSSD(nB))
        RhoP(nA)=RhoP(nA)+SinProjection*(RhoSSD(nB))
       END DO
	   
C 3333	   
      DO nB=1,6     
		
        call VectorProjections(CosProjection,SinProjection,
     1 FCC_N(nA,1:3),CUBIC_T(nB,1:3))
	 
        RhoF(nA)=RhoF(nA)+CosProjection*(RhoSSD(nB))
        RhoP(nA)=RhoP(nA)+SinProjection*(RhoSSD(nB))
		
       END DO
      END DO
	

C=======================  
      DO nA=1,6
		nAX=nA+12
       DO nB=1,12     
		
        call VectorProjections(CosProjection,SinProjection,
     1 CUBIC_N(nA,1:3),FCC_T(nB,1:3))
	 
        RhoF(nAX)=RhoF(nAX)+CosProjection*(RhoSSD(nB))
        RhoP(nAX)=RhoP(nAX)+SinProjection*(RhoSSD(nB))
       END DO
	   
C 3333	   
      DO nB=1,6     
		
        call VectorProjections(CosProjection,SinProjection,
     1 CUBIC_N(nA,1:3),CUBIC_T(nB,1:3))
	 
        RhoF(nAX)=RhoF(nAX)+CosProjection*(RhoSSD(nB))
        RhoP(nAX)=RhoP(nAX)+SinProjection*(RhoSSD(nB))
		
       END DO
      END DO
	  
C-------------------	
      DO nA=1,18
       RhoM(nA)=Coefficient1*sqrt(RhoF(nA)*RhoP(nA))
      END DO
      return
      end subroutine GetRhoPFM