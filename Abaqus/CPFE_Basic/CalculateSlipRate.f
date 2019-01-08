      subroutine CalculateSlipRate( 
     1  TAU, TAU_SIGN,
     +  FCC_S, FCC_N, G, Gamma,
     2  Lp, GammaDot, dGammadTau, 
     3  PROPS,nprops)
	 
      implicit none
      
      INTEGER, PARAMETER :: NDOF = 3
      INTEGER, PARAMETER :: NFCC = 12
      
      REAL*8, INTENT(IN) :: TAU(12), TAU_SIGN(12)
      REAL*8, INTENT(IN) :: FCC_S(NDOF,NFCC), FCC_N(NDOF,NFCC)
      REAL*8, INTENT(IN) :: PROPS(18)
      INTEGER, INTENT(IN) :: nprops

      REAL*8, INTENT(IN) :: G(12), Gamma(12)
      
      REAL*8, INTENT(OUT) :: Lp(3,3), GammaDot(12), dGammadTau(6,6)

      
      REAL*8 :: tempDGammaDTau(6,6), X
      INTEGER :: I, J, K, ISLIPS
      REAL*8 :: tempNorm(3), tempDir(3)
      REAL*8 :: tempVal, tempStress
      REAL*8 :: xsnt(3,3),xsnv(6),xnsv(6),xsnnst(6,6),xnst(3,3)
c -------------------------------------------
      dGammadTau = 0.0
      Lp = 0.0 
      tempDGammaDTau = 0.0
      

      DO ISLIPS = 1, 12   
	    
		X = TAU(ISLIPS)/G(ISLIPS)
c		X = TAU(ISLIPS)/40.0
      !write(201,*) "PROPS"	, PROPS(17),PROPS(18) 
      !write(201,*) "TG"	, ISLIPS,TAU(ISLIPS),G(ISLIPS)  
        GAMMADOT(ISLIPS) = PROPS(17)*(X**(PROPS(18)))
            tempDir = TAU_SIGN(ISLIPS)*FCC_S(:,ISLIPS)
            tempNorm = FCC_N(:,ISLIPS)		
c ----  
         xsnt = spread(tempDir,2,3)*spread(tempNorm,1,3)
         xnst = spread(tempNorm,2,3)*spread(tempDir,1,3)
         CALL KGMATVEC6(xsnt,xsnv)         
         CALL KGMATVEC6(xnst,xnsv) 
         xsnnst = spread(xsnv,2,6)*spread(xnsv,1,6)          
c ----
         tempDGammaDTau = tempDGammaDTau +
     1       xsnnst*(
     1       PROPS(17)*PROPS(18)*(X**(PROPS(18)-1.0)))
     
         Lp = Lp + gammaDot(ISLIPS)*xsnt           

            
      END DO
      
        dGammadTau = 0.5*(tempDGammaDTau+transpose(tempDGammaDTau))   
     
      return
      end 	  