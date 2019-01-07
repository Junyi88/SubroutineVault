      subroutine CalculateTauS(StressTrialMat, 
     1  TAU, TAU_SIGN,
     +  FCC_S, FCC_N)	 

C Subroutine Calculating All Values of ResolveShearStress
      
      implicit none
      
      INTEGER, PARAMETER :: NDOF = 3
      INTEGER, PARAMETER :: NFCC = 12
      INTEGER, PARAMETER :: NCUB = 6
      
      REAL*8, INTENT(IN) :: StressTrialMat(3,3)
      REAL*8, INTENT(IN) :: FCC_S(NDOF,NFCC), FCC_N(NDOF,NFCC)
      REAL*8, INTENT(IN) :: CUBIC_S(NDOF,NCUB), CUBIC_N(NDOF,NCUB)
	  
      REAL*8, INTENT(IN) :: FCC_SPE(NDOF,NFCC), FCC_NPE(NDOF,NFCC)
      REAL*8, INTENT(IN) :: FCC_SSE(NDOF,NFCC), FCC_NSE(NDOF,NFCC)
      REAL*8, INTENT(IN) :: FCC_SCB(NDOF,NFCC), FCC_NCB(NDOF,NFCC)

      REAL*8, INTENT(OUT) :: TAU(12), TAUPE(12), TAUSE(12), TAUCB(12)  
      REAL*8, INTENT(OUT) :: TAU_SIGN(12)

      INTEGER :: ISLIPS, I
      REAL*8 :: tempNorm(3), tempDir(3)
      REAL*8 :: PROD(3)
      
C --------------------------------------------------------------------------
      DO ISLIPS=1,12    
       TAU(ISLIPS)=0.0
      END DO 

      
C --------------------------------------------------------------------------
C Calculate Tau	  
      DO I=1,12
        IF (I.LE.12) THEN            
		 tempNorm = FCC_N(:,I)
         tempDir = FCC_S(:,I)

      
        PROD = matmul(StressTrialMat,tempNorm)
        tau(I)= dot_product(prod,tempDir)
        
        IF(tau(I) .LT. 0.0) THEN
          tau(I) = -1.0*tau(I)
          TAU_SIGN(I) = -1.0
        ELSE 
          TAU_SIGN(I) = 1.0
        END IF
        
      END DO

	  
      return
      end