      subroutine RotateSlipSystems(ROTM,
     1 FCC_S, FCC_N, FCC_T,
     1 FCC_SPE, FCC_NPE, 
     1 FCC_SSE, FCC_NSE,
     1 FCC_SCB, FCC_NCB,
     1 CUBIC_S, CUBIC_N, CUBIC_T)
      
      IMPLICIT NONE     
      INCLUDE 'DeclareParameterSlipsO.f'

      INTEGER :: I,J,K
      INTEGER, PARAMETER :: NDOF = 3
      INTEGER, PARAMETER :: NFCC = 12
      INTEGER, PARAMETER :: NCUB = 6
	  
      REAL*8, INTENT(IN) :: ROTM(NDOF,NDOF)
      REAL*8, INTENT(OUT) :: FCC_S(NDOF,NFCC), FCC_N(NDOF,NFCC)
      REAL*8, INTENT(OUT) :: FCC_T(NDOF,NFCC)
	  
      REAL*8, INTENT(OUT) :: FCC_SPE(NDOF,NFCC), FCC_NPE(NDOF,NFCC)
      REAL*8, INTENT(OUT) :: FCC_SSE(NDOF,NFCC), FCC_NSE(NDOF,NFCC)
      REAL*8, INTENT(OUT) :: FCC_SCB(NDOF,NFCC), FCC_NCB(NDOF,NFCC)
	  
      REAL*8, INTENT(OUT) :: CUBIC_S(NDOF,NCUB), CUBIC_N(NDOF,NCUB)
      REAL*8, INTENT(OUT) :: CUBIC_T(NDOF,NCUB)

      REAL*8 :: TMagS, TVecS(NDOF)	  
      REAL*8 :: TMagN, TVecN(NDOF)
      REAL*8 :: TMagT, TVecT(NDOF)
C ----------------------------------------
      
      DO K = 1, NFCC
        DO I = 1, NDOF
           TVecS(I) = FCC_S0(I,K)
           TVecN(I) = FCC_N0(I,K)
           TVecT(I) = FCC_T0(I,K)
        END DO		
		
        TVecS = matmul(ROTM,TVecS)
        TVecN = matmul(ROTM,TVecN)
        TVecT = matmul(ROTM,TVecT)
		
        TMagS = sqrt(dot_product(TVecS,TVecS))
        TMagN = sqrt(dot_product(TVecN,TVecN))
        TMagT = sqrt(dot_product(TVecT,TVecT))
        
        DO I = 1, NDOF
         FCC_S(I,K) = TVecS(I)/TMagS
         FCC_N(I,K) = TVecN(I)/TMagN
         FCC_T(I,K) = TVecT(I)/TMagT
        END DO		
      END DO		  
	  
c -------------	  
      DO K = 1, NFCC
        DO I = 1, NDOF
           TVecS(I) = FCC_SPE0(I,K)
           TVecN(I) = FCC_NPE0(I,K)
        END DO		
		
        TVecS = matmul(ROTM,TVecS)
        TVecN = matmul(ROTM,TVecN)
		
        TMagS = sqrt(dot_product(TVecS,TVecS))
        TMagN = sqrt(dot_product(TVecN,TVecN))
		
        DO I = 1, NDOF
         FCC_SPE(I,K) = TVecS(I)/TMagS
         FCC_NPE(I,K) = TVecN(I)/TMagN
        END DO		
      END DO	

c -------------	  
      DO K = 1, NFCC
        DO I = 1, NDOF
           TVecS(I) = FCC_SSE0(I,K)
           TVecN(I) = FCC_NSE0(I,K)
        END DO		
		
        TVecS = matmul(ROTM,TVecS)
        TVecN = matmul(ROTM,TVecN)
		
        TMagS = sqrt(dot_product(TVecS,TVecS))
        TMagN = sqrt(dot_product(TVecN,TVecN))
		
        DO I = 1, NDOF
         FCC_SSE(I,K) = TVecS(I)/TMagS
         FCC_NSE(I,K) = TVecN(I)/TMagN
        END DO		
      END DO

c -------------	  
      DO K = 1, NFCC
        DO I = 1, NDOF
           TVecS(I) = FCC_SCB0(I,K)
           TVecN(I) = FCC_NCB0(I,K)
        END DO		
		
        TVecS = matmul(ROTM,TVecS)
        TVecN = matmul(ROTM,TVecN)
		
        TMagS = sqrt(dot_product(TVecS,TVecS))
        TMagN = sqrt(dot_product(TVecN,TVecN))
		
        DO I = 1, NDOF
         FCC_SCB(I,K) = TVecS(I)/TMagS
         FCC_NCB(I,K) = TVecN(I)/TMagN
        END DO		
      END DO
	  
	  
c -------------	 
      DO K = 1, NCUB
        DO I = 1, NDOF
           TVecS(I) = CUBIC_S0(I,K)
           TVecN(I) = CUBIC_N0(I,K)
           TVecT(I) = CUBIC_T0(I,K)
        END DO		
		
        TVecS = matmul(ROTM,TVecS)
        TVecN = matmul(ROTM,TVecN)
        TVecT = matmul(ROTM,TVecT)
		
        TMagS = sqrt(dot_product(TVecS,TVecS))
        TMagN = sqrt(dot_product(TVecN,TVecN))
        TMagT = sqrt(dot_product(TVecT,TVecT))
		
        DO I = 1, NDOF
         CUBIC_S(I,K) = TVecS(I)/TMagS
         CUBIC_N(I,K) = TVecN(I)/TMagN
         CUBIC_T(I,K) = TVecT(I)/TMagT
        END DO		
      END DO	
	  
      return
      end 