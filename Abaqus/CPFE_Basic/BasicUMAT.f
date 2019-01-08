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
      open (201, file = '/home/jl1908/projects/scrap/X/LogLog.txt', 
     1 ACCESS='APPEND', ACTION='WRITE') 	  
      !write(201,*) "C0"	 
c ------------------------------------------------	
c ---- INITIALISE
      CALL InitialiseUMAT(
     1 PROPS, STATEV, 
     1 NPROPS, NSTATV,
     1 KSTEP, KINC, NOEL, NPT,
     1 ROTM, COORDS
     1 )     
	 
        DO I = 31, 42
            STATEV(I) = PROPS(13)
        END DO	 
      !write(201,*) "C1"	 	  
C ==== Preliminary Calculations before constitutive ============= 
c ---- LOAD SDV
      DO I = 1,3
      DO J = 1,3
         ROTM(I,J)=STATEV(J+(I-1)*3)                                   
         FP(I,J) = STATEV(9+J+(I-1)*3)
c         STATEV(18+J+(I-1)*3) = kcurlFp(NOEL, NPT,J+(I-1)*3) !?
      END DO
      END DO
      !write(201,*) "C2"	 	  
c ---- ROTATE STIFFNESS TENSOR
      CALL RotateStiffnessTensorSimple(ROTM,
     1 PROPS(10:12), StiffR)  

      !write(201,*) "C3"	 	 
      CALL RotateSlipSystems(ROTM,
     1 FCC_S, FCC_N, FCC_T)   
      !write(201,*) "C4"	 	 
c -------------------------------	  
C ==== Initialise Loop For Newton =============       
      CALL PrepNewtonIter(
     1 DFGRD1,DFGRD0, DTIME,
     1 StiffR, STRESS,
     1 Ltot, DTotStran,
     1 StressV, StressVMat,
     1 StressTrial, StressTrialMat, F
     1 )     
      !write(201,*) "C5"	 
      call CalculateTauS(StressTrialMat, 
     1  TAU, TAU_SIGN,
     +  FCC_S, FCC_N)	 	 
      !write(201,*) "C6"	 	 
      call CalculateSlipRate( 
     1  TAU, TAU_SIGN,
     +  FCC_S, FCC_N, STATEV(31:42), STATEV(43:54),
     2  Lp, GammaDot, dGammadTau,  
     3  PROPS,nprops)
      !write(201,*) "C7"	 
c     SDV(31->42) :: G
c     SDV(43->54) :: gamma	 
       IF (KINC>1) THEN
c ------------------------------------------------	
      ITRNUM= 0
      FaiValue = 1.0
      DO WHILE (FaiValue  .GT. IterAccu)        
        ITRNUM= ITRNUM+ 1
      !write(201,*) "ITRNUM------"	,ITRNUM , NOEL, NPT, KINC
      !write(201,*) dGammadTau 	
	        !write(201,*) "G----"
      !write(201,*) STATEV(31:42) 	
      !write(201,*) "TAU" , TAU	  
          if(any(dGammadTau /= dGammadTau) .or. 
     1          any(dGammadTau-1 == dGammadTau)) then

              pnewdt = 0.5 ! if sinh( ) has probably blown up then try again with smaller dt
              return
          end if  

       IF (ITERN.GT.1) THEN
c --- Now to redo the stuff
      call CalculateTauS(StressTrialMat, 
     1  TAU, TAU_SIGN,
     +  FCC_S, FCC_N)	 	 
	 
      call CalculateSlipRate( 
     1  TAU, TAU_SIGN,
     +  FCC_S, FCC_N, STATEV(31:42), STATEV(43:54),
     2  Lp, GammaDot, dGammadTau, 
     3  PROPS,nprops)
      END IF		  
c --- Calculate Plastic Strains
      plasStrainRate = (Lp+transpose(Lp))*0.5*dtime
      CALL kmatvec6(plasStrainRate,plasStrainInc2)
      plasStrainInc2(4:6) = 2.0*plasStrainInc2(4:6) 
c --- Calculate Stress Increments
      dGammadTau=dGammadTau*DTIME
      xFai =  xIden6 + matmul(StiffR,dGammadTau)
      CALL lapinverse(xFai,6,info,xFaiInv)   
        !write(201,*) "GammaDot" , GammaDot
		!write(201,*) "dGammadTau" , dGammadTau
	        !write(201,*) "LP" , LP
	        !write(201,*) "xfai" , xfai
			!write(201,*) "xFaiInv" , xFaiInv
      Fai = StressTrial - StressV - matmul(StiffR,plasStrainInc2)
	   !write(201,*) "FAI" , FAI
      dStress = matmul(xFaiInv,Fai)
	  !write(201,*) "dStress" , dStress
      StressV = StressV + dStress
      CALL kvecmat6(StressV,StressVMat)      
      FaiValue = sqrt(sum(Fai*Fai))
      
c --- Ditch if too big
      IF (ITRNUM.gt. MaxITER) THEN
        IF (FAIVALUE.GT.FuzzyAccept) THEN
           pnewdt = 0.1       
           return           
         ELSE
            FaiValue = 0.0
         ENDIF
      END IF

        STATEV(55) = FaiValue
        STATEV(56) = ITRNUM                    

      
      END DO
      !write(201,*) "C8"	 
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
      CALL kmatvec6(StressVMat,STRESS) !output stress	  
C **** FINISH NEWTON CALCULATIONS ***************  
c PROP(13): g0
c PROP(14): tau0
c PROP(15): taus
c PROP(16): h0
c PROP(17): adot
c prOP(18): N
      !write(201,*) "C9"	 
c     SDV(31->42) :: G
c     SDV(43->54) :: gamma
      do I =1 ,12
       STATEV(42+I) = STATEV(42+I) + GAMMADOT(I)*DTIME
       HH(I) = 1.0/cosh(abs(PROPS(16)*STATEV(42+I))/(
     1 PROPS(15)-PROPS(14)))
	   dG(I) = (HH(I)**2.0)*DTIME*PROPS(16)*GAMMADOT(I)
	   STATEV(30+I)= STATEV(30+I) +dG(I)
      END DO

C **** FINISH NEWTON CALCULATIONS ***************  
	        ELSE
cc --- Elastic      
      StressVMat = StressTrialMat
      DDSDDE = StiffR 
      plasStrainrate= 0.0 

      IF ((KSTEP.GT.1).AND.(KINC.LE.1)) THEN
      DO I = 1,3
      DO J = 1,3
         FP(I,J) = STATEV(9+J+(I-1)*3)  
      END DO
      END DO
      ENDIF
      CALL kmatvec6(StressVMat,STRESS) !output stress
      ENDIF
C ====  NOW TO UPDATE EVERYTHING ===================================================== 
C     STRESS = DONE
C     DDSDDE = DONE
      !write(201,*) "C10"	 
      DO I = 1,3
      DO J = 1,3
         STATEV(9+J+(I-1)*3) = FP(I,J) 
      END DO
      END DO
 
      call MutexLock( 5 )      
      DO i=1, 9
          kFp(noel,NPT,i) = STATEV(I)
      END DO
      call MutexUnlock( 5 )  
      !write(201,*) "C11"	 
C ===== ROTATE ========================
      CALL kdeter(Fp,deter)            
      IF (deter /= 0.) THEN
         Fpinv = 0.
         CALL lapinverse(Fp,3,info,Fpinv)
!         IF(info5 /= 0) !write(201,*) "inverse failure: print3 in kmat"
         Fe = matmul(F,Fpinv) ! SHOULD BE THE NEW F RIGHT?         
      ELSE
         !write(201,*) "Error in orientation update: finding inv(Fp)",noel,
     +    npt, kinc, plasStrainrate
         !write(201,*) Fp
         call XIT       
      END IF 
      
      CALL kdeter(Fe,deter)            
      IF (deter /= 0.) THEN
         Feinv = 0.
         CALL lapinverse(Fe,3,info5,Feinv)
!         IF(info5 /= 0) !write(201,*) "inverse failure: print3 in kmat"         
      ELSE
          !write(201,*) "Error in orientation update: finding inv(Fe)",noel
     +     ,npt, kinc, deter
          !write(201,*) Fe
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
         IF(info /= 0) !write(201,*) "inverse failure: print3 in kmat"
         !print3 = 0.
         !print3 = gmatinv + dtime*matmul(elasspin,gmatinv)
         ROTMnew = matmul(update,ROTM)                        
         !!write(201,*) "gmatinv, print3", gmatinv, print3      
      ELSE         
         ROTMnew = ROTM
      !write(201,*) "WARNING gmatinv not updated at noel,npt, kinc:", noel,
     + npt, kinc
      END IF    
      
      DO I = 1,3
        DO J = 1,3 
        STATEV(J+(I-1)*3) = ROTMnew(I,J)
        END DO
      END DO        
       if (maxval(ROTMnew) > 1.0) then
          !write(201,*) "something very wrong with gmatinv"
          call XIT
       end if 
      !write(201,*) "CE"	 
	  
      close(201)
      return
      end subroutine UMAT
	  
	  
c ---- Includes ------   
      include 'uexternaldb.f'
      include 'utils.f'
      include 'ksvd2.f'
      
      include 'InitialiseUMAT.f'
      include 'RotateStiffnessTensorSimple.f'
      include 'RotateSlipSystems.f'

	  
      include 'PrepNewtonIter.f'
      include 'CalculateTauS.f'
      include 'CalculateSlipRate.f'