      subroutine RotateSlipSystems(ROTM,
     1 FCC_S, FCC_N, FCC_T)
      
      IMPLICIT NONE     
      INCLUDE 'DeclareParameterSlipsO.f'

      INTEGER :: I,J,K
      INTEGER, PARAMETER :: NDOF = 3
      INTEGER, PARAMETER :: NFCC = 12
      INTEGER, PARAMETER :: NCUB = 6
	  
      REAL*8, INTENT(IN) :: ROTM(NDOF,NDOF)
      REAL*8, INTENT(OUT) :: FCC_S(NDOF,NFCC), FCC_N(NDOF,NFCC)
      REAL*8, INTENT(OUT) :: FCC_T(NDOF,NFCC)
	  

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
	  
	  
	  

      return
      end 