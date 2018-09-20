      SUBROUTINE UMAT(stress,statev,ddsdde,sse,spd,scd,
     1 rpl, ddsddt, drplde, drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
	 
      include 'aba_param.inc'
      CHARACTER*8 CMNAME
      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     2 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     3 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
     
c ---- User Start ------     
      include 'UserParameters.f'     
      include 'HeaderUMAT.f'      
     
     
c ---- INITIALISE
      CALL InitialiseUMAT(
     1 PROPS, STATEV, 
     1 NPROPS, NSTATV,
     1 KSTEP, KINC, NOEL, NPT,
     1 ROTM
     1 )     

C ==== Preliminary Calculations before constitutive =============     
c ---- ROTATE STIFFNESS TENSOR
      CALL RotateStiffnessTensorSimple(ROTM,
     1 PROPS(10:12), StiffR)     
     
c ---- ROTATE Slip Systems
      CALL RotateSlipSystems(ROTM,
     1 FCC_S, FCC_N, FCC_T,
     1 FCC_SPE, FCC_NPE, 
     1 FCC_SSE, FCC_NSE,
     1 FCC_SCB, FCC_NCB,
     1 CUBIC_S, CUBIC_N, CUBIC_T)     
     
c ---- Calculate Rhos    
      CALL GetRhoPFM(RhoP,RhoF,RhoM,
     1 SDV(28:45),
     2 SDV(64:81),SDV(82:99),SDV(100:117),
     2 FCC_S, FCC_N, FCC_T,  
     2 CUBIC_S, CUBIC_N, CUBIC_T,
     5 PROPS(13))
     
c ---- Calculate Tau Slips 
      CALL GetTauSlips(RhoP,RhoF,RhoM,
     1 TauPass, TauCut, PROPS(14:15))

C ==== Initialise Loop For Newton =============       
      CALL PrepNewtonIter(
     1 DFGRD1,DFGRD0, DTIME,
     1 StiffR, STRESS,
     1 Ltot, DTotStran,
     1 StressV, StressVMat,
     1 StressTrial, StressTrialMat
     1 )     
     
      CALL CalculateTauS(StressTrialMat, 
     1  TAU, TAUPE, TAUSE, TAUCB, TAU_SIGN,
     +  FCC_S, FCC_N, CUBIC_S, CUBIC_N,
     +  FCC_SPE, FCC_NPE,
     +  FCC_SSE, FCC_NSE,
     +  FCC_SCB, FCC_NCB)	     
     
      CALL GetCSDHTauC(TAUPEI,TAUSEI,TAUCBI,
     1 SDV(46:63), TAUC, 	   
     2 PROP(16:28))
 
      CALL CalculateSlipRate( 
     1  TAU, TAU_SIGN,
     2  TauPass, TauCut, TauC,
     +  FCC_S, FCC_N, CUBIC_S, CUBIC_N,
     2  Lp, GammaDot, dGammadTau, TauEff, 
     3  PROP(29:30))
      
      PlasticFlag = maxval(TauEff)
      IF (PlasticFlag.GT.UserZero) THEN     
C **** START NEWTON CALCULATIONS ***************
      DO WHILE (FaiVal .GT. IterAccu)        
        ITERN = ITERN + 1
          if(any(dGammadTau /= tmatdGammadTau) .or. 
     1          any(dGammadTau-1 == dGammadTau)) then
              pnewdt = 0.5 ! if sinh( ) has probably blown up then try again with smaller dt
              return
          end if  
c --- Calculate Plastic Strains
      plasStrainRate = (Lp+transpose(Lp))*0.5*dtime
      CALL kmatvec6(plasStrainRate,plasStrainInc2)
      plasStrainInc2(4:6) = 2.0*plasStrainInc2(4:6) 
c --- Calculate Stress Increments
      xFai =  xIden6 + matmul(StiffR,dGammadTau)
      CALL lapinverse(xFai,6,info,xFaiInv)   

      Fai = StressTrial - StressV - matmul(StiffR,plasStrainInc2)
      dStress = matmul(xFaiInv,Fai)
      StressV = StressV + dStress
      CALL kvecmat6(StressV,StressVMat)      
      FaiValue = sqrt(sum(Fai*Fai))
      
c --- Ditch if too big
      IF (ITERN .gt. MaxITER) THEN
          pnewdt = 0.5
      END IF

c --- Now to redo the stuff
      CALL CalculateTauS(StressVMat, 
     1  TAU, TAUPE, TAUSE, TAUCB, TAU_SIGN,
     +  FCC_S, FCC_N, CUBIC_S, CUBIC_N,
     +  FCC_SPE, FCC_NPE,
     +  FCC_SSE, FCC_NSE,
     +  FCC_SCB, FCC_NCB)  
     
      CALL GetCSDHTauC(TAUPEI,TAUSEI,TAUCBI,
     1 SDV(46:63), TAUC, 	   
     2 PROP(16:28))
 
      CALL CalculateSlipRate( 
     1  TAU, TAU_SIGN,
     2  TauPass, TauCut, TauC,
     +  FCC_S, FCC_N, CUBIC_S, CUBIC_N,
     2  Lp, GammaDot, dGammadTau, TauEff, 
     3  PROP(29:30))      
      
      END DO

C **** FINISH NEWTON CALCULATIONS ***************
c -- Calculate Other Stuff



      ELSE
cc --- Elastic      
c      StressVMat = StressVMat
      DDSDDE = StiffR 
      
      END
C ================================================================

c ===

  
      return
      end subroutine UMAT

c ---- Includes ------   
      include 'uexternaldb.f'
      include 'utils.f'
      
      include 'InitialiseUMAT.f'
      include 'RotateStiffnessSimple.f'
      include 'RotateSlipSystems.f'
      
      include 'GetRhoPFM.f'
      include 'GetTauPass.f'
      
      include 'PrepNewtonIter.f'
      include 'CalculateTauS.f'
      include 'GetCSDHTauC.f'
      include 'CalculateSlipRate.f'
      
      