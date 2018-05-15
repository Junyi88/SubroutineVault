      SUBROUTINE UMAT(stress,statev,ddsdde,sse,spd,scd,
     1 rpl, ddsddt, drplde, drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
	 
      include 'aba_param.inc'
c
      CHARACTER*8 CMNAME
      EXTERNAL F

      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     2 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     3 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
	 
      include 'DeclareParameterSlipsO.f'
      
      INTEGER:: ISLIPS, I, J, NDUM1, NA, NB
      real*8 :: TAU(18), TAUPE(12), TAUSE(12), TAUCB(12), SLIP_T(54)	  
      real*8 :: RhoP(18),RhoF(18),RhoM(18),RhoSSD(18)
      real*8 :: TauPass(18), TauCut(18), V0(18)
      real*8 :: H(12), RhoCSD(12), TAUC(18) 
      real*8 :: Vs(18) , GammaDot(18) , TauEff(18), SSDDot(18)
      real*8 :: DStress(6) 
	  
      real*8 :: ORI_ROT(3,3), SPIN_TENSOR(3,3)
      real*8:: dFP(9), dRhoS(18),dRhoET(18),dRhoEN(18)
c ------------------------------------------------	  
C
C     CALCULATE VELOCITY GRADIENT FROM DEFORMATION GRADIENT.
C     REFERENCE: Li & al. Acta Mater. 52 (2004) 4859-4875
C     
      Real*8:: FTINV(3,3),STRATE(3,3),VELGRD(3,3),AUX1(3,3),ONEMAT(3,3)
      PARAMETER (ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,SIX=6.0D0)
      DATA NEWTON,TOLER/10,1.D-6/
      Real*8:: gausscoords(3,8)
	  real*8 :: kgausscoords, kFp, kcurlFp
c XDANGER
      COMMON/UMPS/kgausscoords(1,8,3),kFp(1,8, 3),
     1 kcurlFp(1, 8, 3)
	  
c      print *, '*****************************************'
	
      CALL ONEM(ONEMAT)     
      CALL ZEROM(FTINV)
      CALL ZEROM(AUX1)
      CALL ZEROM(VELGRD)
      CALL M3INV(DFGRD0,FTINV)
      CALL MPROD(DFGRD1,FTINV,AUX1)
      DO 231 I=1,3
        DO 231 J=1,3
          VELGRD(I,J) = (AUX1(I,J)-ONEMAT(I,J))
231   CONTINUE	
 
      DO I=1,3
       DO J=1,3
	      SPIN_TENSOR(I,J)=0.5*(VELGRD(I,J)-VELGRD(J,I))
       END DO
      END DO	
c ------------------------------------------------		
c Perform Initialisation

      IF (KINC.LE.1) THEN
	  
       DO ISLIPS=1,nstatv
          STATEV(ISLIPS)=0.0
       END DO

       NDUM1=0
       DO I=1,3
       DO J=1,3
          NDUM1=NDUM1+1
          ORI_ROT(I,J)=PROPS(NDUM1)
       END DO
       END DO
	   
       DO ISLIPS=1,18
          STATEV(ISLIPS+108)=PROPS(ISLIPS+9)
       END DO

c ---- S_SCHMID
       DO ISLIPS=1,12
        NDUM1=(ISLIPS-1)*3
        NA=NDUM1+1
        NB=NDUM1+3		        
        call ROTATE_Vec(ORI_ROT,FCC_S0(1:3,ISLIPS),STATEV(NA:NB))
       END DO
       DO ISLIPS=1,6
        NDUM1=(ISLIPS+11)*3
        NA=NDUM1+1
        NB=NDUM1+3		        
        call ROTATE_Vec(ORI_ROT,CUBIC_S0(1:3,ISLIPS),STATEV(NA:NB))
       END DO

c ---- N_SCHMID
       DO ISLIPS=1,12
        NDUM1=(ISLIPS-1)*3+54
        NA=NDUM1+1
        NB=NDUM1+3		        
        call ROTATE_Vec(ORI_ROT,FCC_N0(1:3,ISLIPS),STATEV(NA:NB))
       END DO
       DO ISLIPS=1,6
        NDUM1=(ISLIPS+11)*3+54
        NA=NDUM1+1
        NB=NDUM1+3		        
        call ROTATE_Vec(ORI_ROT,CUBIC_N0(1:3,ISLIPS),STATEV(NA:NB))
       END DO	
	   
c ---- S_PE
       DO ISLIPS=1,12
        NDUM1=(ISLIPS-1)*3+184
        NA=NDUM1+1
        NB=NDUM1+3		        
        call ROTATE_Vec(ORI_ROT,FCC_SPE0(1:3,ISLIPS),STATEV(NA:NB))
       END DO
c ---- N_PE
       DO ISLIPS=1,12
        NDUM1=(ISLIPS-1)*3+196
        NA=NDUM1+1
        NB=NDUM1+3		        
        call ROTATE_Vec(ORI_ROT,FCC_NPE0(1:3,ISLIPS),STATEV(NA:NB))
       END DO

c ---- S_SE
       DO ISLIPS=1,12
        NDUM1=(ISLIPS-1)*3+208
        NA=NDUM1+1
        NB=NDUM1+3		        
        call ROTATE_Vec(ORI_ROT,FCC_SSE0(1:3,ISLIPS),STATEV(NA:NB))
       END DO
c ---- N_SE
       DO ISLIPS=1,12
        NDUM1=(ISLIPS-1)*3+220
        NA=NDUM1+1
        NB=NDUM1+3		        
        call ROTATE_Vec(ORI_ROT,FCC_NSE0(1:3,ISLIPS),STATEV(NA:NB))
       END DO	   
	   
c ---- S_CB
       DO ISLIPS=1,12
        NDUM1=(ISLIPS-1)*3+232
        NA=NDUM1+1
        NB=NDUM1+3		        
        call ROTATE_Vec(ORI_ROT,FCC_SCB0(1:3,ISLIPS),STATEV(NA:NB))
       END DO
c ---- N_CB
       DO ISLIPS=1,12
        NDUM1=(ISLIPS-1)*3+244
        NA=NDUM1+1
        NB=NDUM1+3		        
        call ROTATE_Vec(ORI_ROT,FCC_NCB0(1:3,ISLIPS),STATEV(NA:NB))
       END DO	      
c ---- Rotate STIFFNESS TENSOR
        call ROTATE_COMTEN(ORI_ROT,PROPS(28:48),STATEV(164:184))

c--- Do Stuff
       STATEV(257)=1.0
       STATEV(261)=1.0
       STATEV(265)=1.0

c XDANGER
       DO ISLIPS=1,12
        STATEV(265+ISLIPS)=1.0e5
        STATEV(283+ISLIPS)=1.0e5
        STATEV(301+ISLIPS)=1.0e5    
       END DO	 
       DO ISLIPS=1,6
        STATEV(271+ISLIPS)=0.0
        STATEV(289+ISLIPS)=0.0
        STATEV(307+ISLIPS)=0.0   
       END DO	
	   
      call MutexLock( 2 )      ! lock Mutex #2      
      ! use original co-ordinates X     
      do i =1,3
          kgausscoords(noel,npt,i) = coords(i)
      end do
      call MutexUnlock( 2 )   ! unlock Mutex #2
		
      ENDIF
c --------------------------------
C Calculate som common values
        call Get_TfromSN(STATEV(1:54),STATEV(55:108),SLIP_T)
	  
c ------------------------------------------------	
c NOW FOR THE MAIN STUFF
        call CalculateTauS(STRESS, TAU, TAUPE, TAUSE, TAUCB,
     +  STATEV(1:54), STATEV(55:108),
     +  STATEV(185:196), STATEV(197:208),
     +  STATEV(209:220), STATEV(221:232),
     +  STATEV(233:244), STATEV(245:256))	

       call GetRhoPFMGND(RhoP,RhoF,RhoM,
     1 STATEV(109:126),
     2 STATEV(1:54),STATEV(55:108),SLIP_T,
     2 RhoS,RhoET,RhoEN,
     5 PROPS(50))
	 
        call GetTauSlips(RhoP,RhoF,RhoM,
     1 TauPass, TauCut, V0,	  
     2 PROPS(51:53))
	 
      call GetCSDHTauC(TAUPE,TAUSE,TAUCB,
     1 H, RhoCSD, TAUC, 	   
     2 PROPS(54:66))	 

      call GetGammaDot(Tau, TauPass, TauCut, V0, RhoM, 
     1 Vs, GammaDot, TauEff, TAUC,	    	   
     2 PROPS(67:69))	 	 

      call GetRhoSSDEvolve(Tau, TauPass, TauCut, V0, RhoM, 
     1 GammaDot, TauEff, SSDDot, STATEV(109:126), RhoF,      	   
     2 PROPS(70:75))

      call GetDSTRESSFP(DStress,GammaDot,dstran,Stress,dTIME, 
     1 STATEV(1:54),STATEV(55:108), dFP,
     2 PROPS(28:48))

      call GetDDSDDE(DDSDDE,Stress,	   
     2 PROPS(28:48))
	 
c ------------------------------------------------	
c UPDATE ALL
      DO ISLIPS=1,18
       STATEV(ISLIPS+108)=STATEV(ISLIPS+108)+DTIME*SSDDot(ISLIPS)	
       STATEV(ISLIPS+144)=STATEV(ISLIPS+144)+DTIME*GammaDot(ISLIPS)
       STATEV(163)=STATEV(163)+DTIME*GammaDot(ISLIPS)
      END DO
      DO ISLIPS=1,6
       Stress(ISLIPS)=Stress(ISLIPS)+DStress(ISLIPS)	
      END DO
      DO ISLIPS=1,12
       STATEV(ISLIPS+126)=RhoCSD(ISLIPS)	
      END DO
c ------------------------------------------------	
c Rotate The Slip Systems

c ---- S
      call RotateSlipSystems(GammaDot,dTIME,DSTRAN, SPIN_TENSOR,
     1 STATEV(1:54),STATEV(55:108),
     2 STATEV(185:196),STATEV(197:208),
     2 STATEV(209:220),STATEV(221:232),
     2 STATEV(233:244),STATEV(245:256))



c ------------------------------------------------	
         DO kint =1,8 
             DO i=1,3         
                 gausscoords(i,kint) = kgausscoords(noel,kint,i)                          
             END DO 
         END DO

      call CalculateDRho(STATEV(257:265),dtime,
     1 gammadot,STATEV(1:54),STATEV(55:108),STATEV(109:126),
     1 dRhoS,dRhoET,dRhoEN,PROPS(69),gausscoords,noel,npt)
		 
         DO ISLIPS=1,9
		     STATEV(256+ISLIPS)=STATEV(256+ISLIPS)+dFP(ISLIPS)
         END DO		 

         DO ISLIPS=1,18
		   STATEV(265+ISLIPS)=STATEV(265+ISLIPS)+dRhoS(ISLIPS)
		   STATEV(283+ISLIPS)=STATEV(283+ISLIPS)+dRhoET(ISLIPS)
		   STATEV(301+ISLIPS)=STATEV(301+ISLIPS)+dRhoEN(ISLIPS)		   
         END DO		
		 
c ------------------------------------------------	 
      return
      end subroutine UMAT

      include 'UTILS1.f'
      include 'StiffnessTensorTools.f'
      include 'CalculateTauS.f'
      include 'GetRhoPFMGND.f'
      include 'GetTauSlips.f'
      include 'GetCSDHTauC.f' 
      include 'GetGammaDot.f'
      include 'GetRhoSSDEvolve.f'
      include 'GetDSTRESS2FP.f'	
      include 'GetDDSDDEN.f'
      include 'VectorProjections.f'
      include 'RotateSlipSystems.f'
	  
      include 'VectorCurl.f'	  	  
      include 'CalculateDRho.f'	  
      include 'kshapes.f'
      include 'utils.f'