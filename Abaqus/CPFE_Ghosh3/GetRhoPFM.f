      subroutine GetRhoPFM(RhoP,RhoF,RhoM,
     1 RhoSSD,
     2 RhoS,RhoET,RhoEN,
     2 FCC_S, FCC_N, FCC_T,  
     2 CUBIC_S, CUBIC_N, CUBIC_T,
     5 Coefficient1)

C Subroutine to calculate forest parallel and mobile dislocations
      
      IMPLICIT NONE
      INTEGER, PARAMETER :: NDOF = 3
      INTEGER, PARAMETER :: NFCC = 12
      INTEGER, PARAMETER :: NCUB = 6
      INTEGER, PARAMETER :: NSLIPS = 18
      
      REAL*8, INTENT(OUT) :: RhoP(18),RhoF(18),RhoM(18)
      
      REAL*8, INTENT(IN) :: RhoSSD(18)
      REAL*8, INTENT(IN) :: RhoS(18),RhoET(18),RhoEN(18)       
      
      REAL*8, INTENT(IN) :: FCC_S(NDOF,NFCC), FCC_N(NDOF,NFCC)
      REAL*8, INTENT(IN) :: FCC_T(NDOF,NFCC)      
      REAL*8, INTENT(IN) :: CUBIC_S(NDOF,NCUB), CUBIC_N(NDOF,NCUB)
      REAL*8, INTENT(IN) :: CUBIC_T(NDOF,NCUB)        
      
      REAL*8, INTENT(IN) :: Coefficient1      
   
      INTEGER :: I,J,K
      REAL*8 :: CosProjection, SinProjection      
      INTEGER :: nA, nB, nAX      
      REAL*8 :: MyFactor=0.05 ! DANGER
      
      REAL*8 :: TSlipNA(3), TSlipNB(3), TSlipSB(3), TSlipTB(3)

C ------------------------------------------------------	
      DO nA=1,18
       RhoP(nA)=0.0
       RhoF(nA)=0.0
      END DO

C=======================  
      DO nA=1,18	
        
        IF (nA.LE.12) THEN
         TSlipNA = FCC_N(:,nA)
        ELSE
         TSlipNA = CUBIC_N(:,nA-12)
        END IF
        
       DO nB=1,18    

        IF (nB.LE.12) THEN
         TSlipNB = FCC_N(:,nB)
         TSlipSB = FCC_S(:,nB)
         TSlipTB = FCC_T(:,nB)
        ELSE
         TSlipNB = CUBIC_N(:,nB-12)
         TSlipSB = CUBIC_S(:,nB-12)
         TSlipTB = CUBIC_T(:,nB-12)         
        END IF   
        
c --- SSD		
        call VectorProjections(CosProjection,SinProjection,
     1 TSlipNA,TSlipTB)

        RhoF(nA)=RhoF(nA)+CosProjection*(RhoSSD(nB))
        RhoP(nA)=RhoP(nA)+SinProjection*(RhoSSD(nB))
c --- RhoET		
        RhoF(nA)=RhoF(nA)+MyFactor*CosProjection*(abs(RhoET(nB)))
        RhoP(nA)=RhoP(nA)+MyFactor*SinProjection*(abs(RhoET(nB)))		
c --- RhoS		
        call VectorProjections(CosProjection,SinProjection,
     1 TSlipNA,TSlipSB)

        RhoF(nA)=RhoF(nA)+MyFactor*CosProjection*(abs(RhoS(nB)))
        RhoP(nA)=RhoP(nA)+MyFactor*SinProjection*(abs(RhoS(nB)))
c --- RhoEN		
        call VectorProjections(CosProjection,SinProjection,
     1 STSlipNA,TSlipNB)

        RhoF(nA)=RhoF(nA)+MyFactor*CosProjection*(abs(RhoEN(nB)))
        RhoP(nA)=RhoP(nA)+MyFactor*SinProjection*(abs(RhoEN(nB)))		
       END DO
	   
      END DO
	
C-------------------	
      DO nA=1,18
       RhoM(nA)=Coefficient1*sqrt(RhoF(nA)*RhoP(nA))
      END DO
      
      
      return
      end 

      include 'VectorProjections.f'