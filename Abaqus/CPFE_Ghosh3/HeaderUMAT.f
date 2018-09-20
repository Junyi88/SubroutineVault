c  --- For Declaration      
      REAL*8 :: kgausscoords, kFp, kcurlFp,      
      COMMON/UMPS/kgausscoords(TOTALELEMENTNUM,8,3),
     1 kFp(TOTALELEMENTNUM,8, 9),
     1 kcurlFp(TOTALELEMENTNUM, 8, 9)
      
      INTEGER :: I,J,K
c  --- For Initialise
      REAL*8 :: ROTM(3,3) 
      INTEGER :: info 
      REAL*8 :: tempValGamma
c  --- For ROTATE TENSOR
      REAL*8 :: StiffR(6,6)

c  --- For ROTATE SLIPS
      INTEGER, PARAMETER :: NFCC = 12
      INTEGER, PARAMETER :: NCUB = 6
      INTEGER, PARAMETER :: NALLSYS = 18
      
      REAL*8 :: FCC_S(3,NFCC), FCC_N(3,NFCC)
      REAL*8 :: FCC_T(3,NFCC)
	  
      REAL*8 :: FCC_SPE(3,NFCC), FCC_NPE(3,NFCC)
      REAL*8 :: FCC_SSE(3,NFCC), FCC_NSE(3,NFCC)
      REAL*8 :: FCC_SCB(3,NFCC), FCC_NCB(3,NFCC)
	  
      REAL*8 :: CUBIC_S(3,NCUB), CUBIC_N(3,NCUB)
      REAL*8 :: CUBIC_T(3,NCUB)
c --- Load FP
      REAL*8 :: FP(3,3)
      
c ---- FOR GetRhoPFM
      REAL*8 :: RhoP(18),RhoF(18),RhoM(18)
      
c ---- FOR GetTauSlips
      REAL*8 :: TauPass(18), TauCut(18)
      
c ---- FOR PrepNewton Iter   
      REAL*8 :: StressV(6), StressVMat(3,3)
      REAL*8 :: StressTrial(6), StressTrialMat(3,3)
      REAL*8 :: Ltot(3,3), DTotStran(6)      

      REAL*8 :: FaiVal = 1.0 + IterAccu
      INTEGER :: ITERN = 0

c----  FOR CalculateTauS      
      REAL*8 :: TAU(18), TAUPE(12), TAUSE(12), TAUCB(12)  
      REAL*8 :: TAU_SIGN(18)
      
c----  FOR GetCSDHTauC
      REAL*8 :: TAUC(18)
      
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
      REAL*8 :: dStress(6), Fai(6)
c --- After NEWTON
      REAL*8 :: tempSys1(3,3),tempSys2(3,3), deter
      REAL*8, PARAMETER :: xI(3,3) = reshape([
     1  1.0, 0.0, 0.0, 
     1  0.0, 1.0, 0.0, 
     1  0.0, 0.0, 1.0, 
     1  ], [3,3]) 

c----  FOR CalculateSlipRate
      REAL*8 :: Lp(3,3), GammaDot(18), dGammadTau(6,6)
      REAL*8 :: TauEff(18)
      
c ---- For Calculate SSD EVOLVE
      real*8 :: SSDDot(18)
      
c ---- For ROTATION
      real*8 :: FPinv(3,3), Fe(3,3), Le(3,3), LpFeinv(3,3)
      real*8 :: elasspin(3,3), matrix(3,3), update(3,3), ROTMnew(3,3)
      
      