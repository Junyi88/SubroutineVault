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
     
c DBGMSG ---     
C       IF (((NOEL.EQ.1).AND.(NPT.EQ.1)).AND.(KINC.LT.10)) THEN
C       write(6,*), '    '
C       write(6,*), '================================================='
C       write(6,*), KSTEP, KINC, NOEL, NPT 
C       END IF     
c ---- INITIALISE
      CALL InitialiseUMAT(
     1 PROPS, STATEV, 
     1 NPROPS, NSTATV,
     1 KSTEP, KINC, NOEL, NPT,
     1 ROTM, COORDS
     1 )     

C ==== Preliminary Calculations before constitutive ============= 
c ---- LOAD SDV
      DO I = 1,3
      DO J = 1,3
         ROTM(I,J)=STATEV(J+(I-1)*3)
         FP(I,J) = STATEV(9+J+(I-1)*3)
         STATEV(18+J+(I-1)*3) = kcurlFp(noel,J+(I-1)*3,i)
      END DO
      END DO
    
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
     1 STATEV(28:45),
     2 STATEV(64:81),STATEV(82:99),STATEV(100:117),
     2 FCC_S, FCC_N, FCC_T,  
     2 CUBIC_S, CUBIC_N, CUBIC_T,
     5 PROPS(13))
     
c DBGMSG ---     
C       IF (((NOEL.EQ.1).AND.(NPT.EQ.1)).AND.(KINC.LT.10)) THEN
C       write(6,*), '    '
C       write(6,*), '-  CP1  -----------------------'
C       write(6,*), 'PlasticFlag = ', PlasticFlag
C       write(6,*), 'RhoP = ', RhoP
C       write(6,*), '    '
C       write(6,*), 'RhoF = ', RhoF
C       write(6,*), '    '
C       write(6,*), 'RhoM = ', RhoM
C       write(6,*), '    '
C       write(6,*), 'FCC_S = ', FCC_S
C       write(6,*), '    '
C       write(6,*), 'FCC_N = ', FCC_N
C       write(6,*), '    '
C       write(6,*), 'FCC_T = ', FCC_T
C       write(6,*), '    '
C       END IF
      
c ---- Calculate Tau Slips 
      CALL GetTauSlips(RhoP,RhoF,RhoM,
     1 TauPass, TauCut, PROPS(14:15))

C ==== Initialise Loop For Newton =============       
      CALL PrepNewtonIter(
     1 DFGRD1,DFGRD0, DTIME,
     1 StiffR, STRESS,
     1 Ltot, DTotStran,
     1 StressV, StressVMat,
     1 StressTrial, StressTrialMat, F
     1 )     
     
      CALL CalculateTauS(StressTrialMat, 
     1  TAU, TAUPE, TAUSE, TAUCB, TAU_SIGN,
     +  FCC_S, FCC_N, CUBIC_S, CUBIC_N,
     +  FCC_SPE, FCC_NPE,
     +  FCC_SSE, FCC_NSE,
     +  FCC_SCB, FCC_NCB)	     
     
      CALL GetCSDHTauC(TAUPE,TAUSE,TAUCB,
     1 STATEV(46:63), TAUC, 	   
     2 PROPS(16:28),PROPS(39))
c   --- DBG
C        TAUC=0.0

      CALL CalculateSlipRate( 
     1  TAU, TAU_SIGN,
     2  TauPass, TauCut, TauC,
     +  FCC_S, FCC_N, CUBIC_S, CUBIC_N,
     2  Lp, GammaDot, dGammadTau, TauEff, 
     3  PROPS(29:30))
      
      PlasticFlag = maxval(TauEff)
      
c DBGMSG ---     
C       IF (((NOEL.EQ.1).AND.(NPT.EQ.1)).AND.(KINC.LT.10)) THEN
C       write(6,*), '    '
C       write(6,*), '-  CP1  -----------------------'
C       write(6,*), 'PlasticFlag = ', PlasticFlag
C       write(6,*), 'Tau = ', Tau
C       write(6,*), '    '
C       write(6,*), 'TauPass = ', TauPass
C       write(6,*), '    '
C       write(6,*), 'TauCut = ', TauCut
C       write(6,*), '    '
C       write(6,*), 'TauEff = ', TauEff
C       write(6,*), '    '
C       write(6,*), 'StressTrialMat = ', StressTrialMat
C       write(6,*), 'PNEWDT = ', PNEWDT
C       END IF
      ITERN = 0
      FaiValue = 1.0

      IF (PlasticFlag.GT.UserZero) THEN     
C **** START NEWTON CALCULATIONS ***************
C       IF ((NPT.EQ.1).AND.(STATEV(200).LT.0.5)) THEN
C       write(6,*), '   '
C             write(6,*), '-------------------------'
C             write(6,*), KINC
C             write(6,*), 'stress trial'
C             write(6,*), StressTrial
C       ENDIF 
      IF (((NOEL.EQ.10).AND.(NPT.EQ.1))) THEN
      write(6,*),''
      write(6,*), KINC, '-----------------------------'
      END IF   
      
      DO WHILE (FaiValue .GT. IterAccu)        
        ITERN = ITERN + 1
c DBGMSG ---     
c     DBG
      IF (((NOEL.EQ.10).AND.(NPT.EQ.2))) THEN
      write(6,*), 'CP1', ITERN, Faivalue, NOEL, NPT, KINC
      END IF  
      IF (((NOEL.EQ.10).AND.(NPT.EQ.2))) THEN
      write(6,*), 'CP1a', ITERN, Faivalue, NOEL, NPT, KINC
      END IF  
      IF (((NOEL.EQ.10).AND.(NPT.EQ.2))) THEN
      write(6,*), 'CP1b', ITERN, Faivalue, NOEL, NPT, KINC
      END IF  
      IF (((NOEL.EQ.10).AND.(NPT.EQ.2))) THEN
      write(6,*), ITERN, Faivalue, NOEL, NPT, KINC
      write(6,*),'Fai'
      write(6,*), Fai      
      write(6,*), ''
      write(6,*),'StressV'
      write(6,*), StressV     
      write(6,*), ''
      write(6,*),'dStress'
      write(6,*), dStress     
      write(6,*), ''
      END IF       
c --- Now to redo the stuff
      IF (ITERN.EQ.1) THEN
      ELSE
      CALL CalculateTauS(StressVMat, 
     1  TAU, TAUPE, TAUSE, TAUCB, TAU_SIGN,
     +  FCC_S, FCC_N, CUBIC_S, CUBIC_N,
     +  FCC_SPE, FCC_NPE,
     +  FCC_SSE, FCC_NSE,
     +  FCC_SCB, FCC_NCB)  
      IF (((NOEL.EQ.10).AND.(NPT.EQ.2))) THEN
      write(6,*), 'CP2', ITERN, Faivalue, NOEL, NPT, KINC
      END IF  
     
      CALL GetCSDHTauC(TAUPE,TAUSE,TAUCB,
     1 STATEV(46:63), TAUC, 	   
     2 PROPS(16:28),PROPS(39))
c   --- DBG
C        TAUC=0.0 
      CALL CalculateSlipRate( 
     1  TAU, TAU_SIGN,
     2  TauPass, TauCut, TauC,
     +  FCC_S, FCC_N, CUBIC_S, CUBIC_N,
     2  Lp, GammaDot, dGammadTau, TauEff, 
     3  PROPS(29:30))      
       END IF       
       IF (((NOEL.EQ.10).AND.(NPT.EQ.2))) THEN
      write(6,*), 'CP3', ITERN, Faivalue, NOEL, NPT, KINC
      END IF         
          if(any(dGammadTau /= dGammadTau) .or. 
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
       IF (((NOEL.EQ.10).AND.(NPT.EQ.2))) THEN
      write(6,*), 'CP4', ITERN, Faivalue, NOEL, NPT, KINC
      END IF      
C     DBG
C       IF ((NPT.EQ.1).AND.(STATEV(200).LT.0.5)) THEN
C           write(6,*), ITERN, FAIVALUE
C           write(6,*), 'STRESSV'
C           write(6,*), STRESSV
C           write(6,*), 'dSTRESS'
C           write(6,*), dSTRESS
C           write(6,*), 'Fai'
C           write(6,*), Fai
C           write(6,*), '  '
C       ENDIF      
      
      CALL kvecmat6(StressV,StressVMat)      
      FaiValue = sqrt(sum(Fai*Fai))
      
      
    
      
      
c --- Ditch if too big
      IF (ITERN .gt. MaxITER) THEN
      
        IF (FAIVALUE.GT.FuzzyAccept) THEN
           pnewdt = 0.1       
           return           
         ELSE
            FaiValue = 0.0
         ENDIF
C           write(6,*) 'HIGH ITER', FaiValue, ITERN, KINC, NPT
C            return
C           write(6,*) 'FAI  = ', Fai
C           write(6,*) 'StressV  = ', StressV
C           write(6,*) 'StressTrial  = ', StressTrial
c DBG  
c       StressV = StressV - matmul(StiffR,plasStrainInc2)
c      StressV = StressTrial - matmul(StiffR,plasStrainInc2)
C       FaiValue = 0.0
C             STATEV(200) = 1.0
      END IF

c DBG      
C       DO I = 1,18
        STATEV(150) = FaiValue
        STATEV(149) = ITERN        
C       END DO
      

      
      END DO

C **** FINISH NEWTON CALCULATIONS ***************
c -- Calculate Other Stuff
      DDSDDE =  matmul(xFaiInv,StiffR)  
      plasStrainrate=(Lp+transpose(Lp))*0.5  
      
C     *** UPDATE PLASTIC DEFORMATION GRADIENT    
      tempSys1 = 0.; tempSys2 = 0.
      tempSys1 = xI - Lp*dtime      
      CALL kdeter(tempSys1,deter)      
      IF (deter /= 0.0) THEN
         CALL lapinverse(tempSys1,3,info,tempSys2)
         fp = matmul(tempSys2,fp)
      ELSE
         fp = fp
      END IF   
      
      ELSE
cc --- Elastic      
      StressVMat = StressTrialMat
      DDSDDE = StiffR 
      plasStrainrate= 0.0 
      ENDIF
C ================================================================

c --   Calculate Other Stuff
      CALL GetSSDEvolve(RhoM, RhoF,STATEV(28:45),
     1 GammaDot, TauEff, SSDDot, TAU,   
     2 PROPS(31:37))
     
      CALL kmatvec6(StressVMat,STRESS) !output stress
c ===
c DBGMSG ---     
C       IF (((NOEL.EQ.1).AND.(NPT.EQ.1)).AND.(KINC.LT.10)) THEN
C       write(6,*), '-  CP2 '
C       write(6,*), 'Stress = ', Stress
C       write(6,*), 'PNEWDT = ', PNEWDT
C       END IF
      
C ====  NOW TO CALCULATE GNDs ===================================================== 
      IF (NPT.EQ.8) THEN 
         CALL PerformCurl(NOEL,NPT)
      END IF
      
C ====  NOW TO UPDATE EVERYTHING ===================================================== 
C     STRESS = DONE
C     DDSDDE = DONE

      DO I = 1,3
      DO J = 1,3
         STATEV(9+J+(I-1)*3) = FP(I,J) 
      END DO
      END DO
 
      call MutexLock( 5 )      
      DO i=1, 9
          kFp(noel,NPT,i) = STATEV(9+J+(I-1)*3)
      END DO
      call MutexUnlock( 5 )   
      
C     RHO SSD
C       DO I = 1,18
C         STATEV(27+I) = STATEV(27+I)+SSDDot(I)*DTIME
C       END DO

C     RHO CSD DONE
C     GNDS

C DBG
C       Call kcalcGND(FCC_S, FCC_N, FCC_T, 
C      + CUBIC_S, CUBIC_N, CUBIC_T,
C      + STATEV(64:81),STATEV(82:99),STATEV(100:117), STATEV(19:27),
C      + PROPS(38))

c     CUMULATIVE GAMMA
      STATEV(140) = 0.0
      STATEV(141) = 0.0
      STATEV(142) = 0.0      
      
      DO I=1,18
        tempValGamma = GAMMADOT(I)*DTIME
        STATEV(117+I) = STATEV(117+I)+tempValGamma
        STATEV(139) = STATEV(139)+tempValGamma

      STATEV(140) = STATEV(140) +  STATEV(27+I)  
      STATEV(141) = STATEV(141) +  STATEV(45+I) 
      
      STATEV(142) = STATEV(142) +  STATEV(63+I)*STATEV(63+I)
      STATEV(142) = STATEV(142) +  STATEV(81+I)*STATEV(81+I)
      STATEV(142) = STATEV(142) +  STATEV(99+I)*STATEV(99+I)      
      END DO      
      STATEV(142) = sqrt(STATEV(142))
      
c DBGMSG ---     
C       IF (((NOEL.EQ.1).AND.(NPT.EQ.1)).AND.(KINC.LT.10)) THEN
C       write(6,*), '-  CP3 '
C       write(6,*), 'FP = ', FP
C       write(6,*), 'F = ', F
C       END IF      
      
C ===== ROTATE ========================
      CALL kdeter(Fp,deter)            
      IF (deter /= 0.) THEN
         Fpinv = 0.
         CALL lapinverse(Fp,3,info,Fpinv)
!         IF(info5 /= 0) write(6,*) "inverse failure: print3 in kmat"
         Fe = matmul(F,Fpinv)  ! SHOULD BE THE NEW F RIGHT?
      ELSE
         write(6,*) "Error in orientation update: finding inv(Fp)",noel,
     +    npt, kinc
         call XIT       
      END IF 
      
      CALL kdeter(Fe,deter)            
      IF (deter /= 0.) THEN
         Feinv = 0.
         CALL lapinverse(Fe,3,info5,Feinv)
!         IF(info5 /= 0) write(6,*) "inverse failure: print3 in kmat"         
      ELSE
          write(6,*) "Error in orientation update: finding inv(Fe)",noel
     +     ,npt, kinc
         call XIT      
      END IF    
      
      LpFeinv = 0.; 
      LpFeinv = matmul(Lp, Feinv)
      Le = L - matmul(Fe,LpFeinv)        
      elasspin=(Le-transpose(Le))*0.5
      matrix = xI - elasspin*dtime      
      CALL kdeter(matrix,deter)         
      IF (deter /= 0.) THEN
         update = 0.
         CALL lapinverse(matrix,3,info,update)
         IF(info /= 0) write(6,*) "inverse failure: print3 in kmat"
         !print3 = 0.
         !print3 = gmatinv + dtime*matmul(elasspin,gmatinv)
         ROTMnew = matmul(update,ROTM)                        
         !write(6,*) "gmatinv, print3", gmatinv, print3      
      ELSE         
         ROTMnew = ROTM
      write(6,*) "WARNING gmatinv not updated at noel,npt, kinc:", noel,
     + npt, kinc
      END IF    

c DBGMSG ---     
C       IF (((NOEL.EQ.1).AND.(NPT.EQ.1)).AND.(KINC.LT.10)) THEN
C       write(6,*), 'ROTM = ', ROTM
C       write(6,*), 'ROTMNew = ', ROTMnew
C       write(6,*), 'F = ', F
C       write(6,*), 'Fp = ', Fp
C       write(6,*), 'Fe = ', Fe
C       END IF   

      
      DO I = 1,3
        DO J = 1,3 
           STATEV(J+(I-1)*3) = ROTMnew(I,J)
        END DO
      END DO        
       if (maxval(ROTMnew) > 1) then
          write(6,*) "something very wrong with gmatinv"
          call XIT
       end if 
C   DEBUG VALUES -------------
C       DO I = 1,18
C         STATEV(200+I) = TauPass(I)
C         STATEV(220+I) = TauCut(I)        
C       END DO

c =========== 
C       IF ((NOEL.EQ.1).AND.(NPT.EQ.1)) THEN
C       write(6,*), 'PNEWDT = ', PNEWDT
C       END IF
      
c --------------------------------------
C       DO ISLIPS=1,6
C        IF ((ABS(DStress(ISLIPS)).LT.5.0e1)) THEN
C        ELSE
C          PNEWDT=0.5
C          write(6,*) "LargeStresses", noel,
C      + npt, kinc
C        END IF	   
C       END DO		
	  
c ------------------------------------------------	
      return
      end subroutine UMAT

c ---- Includes ------   
      include 'uexternaldb.f'
      include 'utils.f'
      include 'ksvd2.f'
      
      include 'InitialiseUMATDBG.f'
      include 'RotateStiffnessTensorSimple.f'
      include 'RotateSlipSystems.f'
      
      include 'GetRhoPFM.f'
      include 'GetTauSlips.f'
      
      include 'PrepNewtonIter.f'
      include 'CalculateTauS.f'
      include 'GetCSDHTauC.f'
       include 'CalculateSlipRate.f'
C       include 'CalculateSlipRateLinear.f' 
      
      include 'GetSSDEvolve.f'
      include 'PerformCurl.f'
      include 'kCalcGND.f'
      
      