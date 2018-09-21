      subroutine CalculateSlipRate( 
     1  TAU, TAU_SIGN,
     2  TauPass, TauCut, TauC,
     +  FCC_S, FCC_N, CUBIC_S, CUBIC_N,
     2  Lp, GammaDot, dGammadTau, TauEff, 
     3  CinS)
	 
      implicit none
      
      INTEGER, PARAMETER :: NDOF = 3
      INTEGER, PARAMETER :: NFCC = 12
      INTEGER, PARAMETER :: NCUB = 6
      
      REAL*8, INTENT(IN) :: TAU(18), TAU_SIGN(18)
      REAL*8, INTENT(IN) :: TauPass(18), TauCut(18), TAUC(18)
      REAL*8, INTENT(IN) :: FCC_S(NDOF,NFCC), FCC_N(NDOF,NFCC)
      REAL*8, INTENT(IN) :: CUBIC_S(NDOF,NCUB), CUBIC_N(NDOF,NCUB)
      REAL*8, INTENT(IN) :: CinS(2)
      
      REAL*8, INTENT(OUT) :: Lp(3,3), GammaDot(18), dGammadTau(6,6)
      REAL*8, INTENT(OUT) :: TauEff(18)
      
      REAL*8 :: tempDGammaDTau(6,6)
      INTEGER :: I, J, K, ISLIPS
      REAL*8 :: tempNorm(3), tempDir(3)
      REAL*8 :: tempVal, tempStress
      REAL*8 :: xsnt(3,3),xsnv(6),xnsv(6),xsnnst(6,6),xnst(3,3)
c -------------------------------------------
      dGammadTau = 0.0
      Lp = 0.0 
      tempDGammaDTau = 0.0
      

      DO ISLIPS = 1, 18   
        IF (Tau(ISLIPS).GT.TauPass(ISLIPS)) THEN 
            TauEff(ISLIPS) = Tau(ISLIPS) - TauPass(ISLIPS)
        ELSE
            TauEff(ISLIPS) = 0.0
        END IF
        
c ================================     
        IF (TauEff(ISLIPS).GT.TauC(ISLIPS)) THEN
            tempStress = TauEff(ISLIPS) / TauCut(ISLIPS)
            tempVal = tempStress**CinS(2)
            GammaDot(ISLIPS) = CinS(1)*Sinh(tempVal)

c ---- Calculate the differential
        if (ISLIPS.LE.12) THEN
            tempDir = TAU_SIGN(ISLIPS)*FCC_S(:,ISLIPS)
            tempNorm = FCC_N(:,ISLIPS)
        ELSE
            tempDir = TAU_SIGN(ISLIPS)*CUBIC_S(:,ISLIPS-12)
            tempNorm = CUBIC_N(:,ISLIPS-12)
        ENDIF
c ----  
         xsnt = spread(tempDir,2,3)*spread(tempNorm,1,3)
         xnst = spread(tempNorm,2,3)*spread(tempDir,1,3)
         CALL KGMATVEC6(xsnt,xsnv)         
         CALL KGMATVEC6(xnst,xnsv) 
         xsnnst = spread(xsnv,2,6)*spread(xnsv,1,6)          
c ----
         tempDGammaDTau = tempDGammaDTau +
     1       CinS(1)*CinS(2)*(tempStress**(Cins(2)-1.0))*
     1       cosh(tempVal)/TauCut(ISLIPS)
     
         Lp = Lp + gammaDot(ISLIPS)*xsnt           
        ELSE
            GammaDot(ISLIPS) = 0.0
        END IF        
            
      END DO
      
        dGammadTau = 0.5*(tempDGammaDTau+transpose(tempDGammaDTau))   
     
      return
      end 	  