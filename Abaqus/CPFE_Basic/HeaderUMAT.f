      INTEGER :: I,J,K, NNN
c  --- For Initialise
      REAL*8 :: ROTM(3,3) 
      INTEGER :: info 

c  --- For ROTATE TENSOR
      REAL*8 :: StiffR(6,6)

c  --- For ROTATE SLIPS
      INTEGER, PARAMETER :: NFCC = 12
      
      REAL*8 :: FCC_S(3,NFCC), FCC_N(3,NFCC)
      REAL*8 :: FCC_T(3,NFCC)
	  
c --- Load FP
      REAL*8 :: FP(3,3)
 
c  --- For Declaration      
      REAL*8 :: kgausscoords, kFp, kcurlFp      
      COMMON/UMPS/kgausscoords(TOTALELEMENTNUM,8,3),
     1 kFp(TOTALELEMENTNUM,8, 9),
     1 kcurlFp(TOTALELEMENTNUM, 8, 9)    
	 
c ---- FOR PrepNewton Iter   
      REAL*8 :: StressV(6), StressVMat(3,3)
      REAL*8 :: StressTrial(6), StressTrialMat(3,3)
      REAL*8 :: Ltot(3,3), DTotStran(6)      

      REAL*8 :: FaiValue = 1.0 + IterAccu
      INTEGER :: ITRNUM= 0

c----  FOR CalculateTauS      
      REAL*8 :: TAU(12)
      REAL*8 :: TAU_SIGN(12)	 
c----  NEWTON CALCULATIONS
      REAL*8 :: PlasticFlag
      REAL*8, PARAMETER :: xIden6(6,6) = reshape([
     1  1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1  0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
     1  0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
     1  0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
     1  0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
     1  0.0, 0.0, 0.0, 0.0, 0.0, 1.0     
     1  ], [6,6]) 
     
      REAL*8 :: plasStrainRate(3,3), plasStrainInc2(6)
      REAL*8 :: xFai(6,6)
      REAL*8 :: xFaiInv(6,6)
      REAL*8 :: dStress(6), Fai(6), FUCK(6)
c --- After NEWTON
      REAL*8 :: tempSys1(3,3),tempSys2(3,3), deter
      REAL*8, PARAMETER :: xI(3,3) = reshape([
     1  1.0, 0.0, 0.0, 
     1  0.0, 1.0, 0.0, 
     1  0.0, 0.0, 1.0 
     1  ], [3,3]) 

	 
      REAL*8 :: dG(12), HH(12)
c      REAL*8 :: TAU_SIGN(12)	 

     
      REAL*8, INTENT(OUT) :: Lp(3,3), GammaDot(12), dGammadTau(6,6)

	  