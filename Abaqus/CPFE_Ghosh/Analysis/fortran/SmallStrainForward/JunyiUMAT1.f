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
	 
      include 'DeclareParametersForSlips.f'
      
      INTEGER:: ISLIPS
      real*8 :: TAU(18), TAUPE(12), TAUSE(12), TAUCB(12)	  
      real*8 :: RhoP(18),RhoF(18),RhoM(18),RhoSSD(18)
      real*8 :: TauPass(18), TauCut(18), V0(18)
      real*8 :: H(12), RhoCSD(12), TAUC(18) 
      real*8 :: Vs(18) , GammaDot(18) , TauEff(18), SSDDot(18)
      real*8 :: DStress(6)
	  
	  
      IF (KINC.LE.1) THEN
      DO ISLIPS=1,18
       RhoSSD(ISLIPS)=PROPS(NPROPS)
       STATEV(ISLIPS)=RhoSSD(ISLIPS)
      END DO
      ENDIF

c ------------------------------------------------	 
      DO ISLIPS=1,18
       RhoSSD(ISLIPS)=STATEV(ISLIPS)
      END DO
	  
c ------------------------------------------------	 	  
      call CalculateTauS(STRESS, TAU, TAUPE, TAUSE, TAUCB,
     +  FCC_N,FCC_S,
     +  FCC_NPE,FCC_SPE,
     +  FCC_NSE,FCC_SSE,
     +  FCC_NCB,FCC_SCB,
     +  CUBIC_N,CUBIC_S)
	 
	 
      DO ISLIPS=1,18
       STATEV(ISLIPS+18)=TAU(ISLIPS)
       STATEV(ISLIPS+36)=TAUPE(ISLIPS)		
       STATEV(ISLIPS+54)=TAUSE(ISLIPS)	
       STATEV(ISLIPS+72)=TAUCB(ISLIPS)	   
      END DO
	  
	  
      call GetRhoPFM(RhoP,RhoF,RhoM,
     1 RhoSSD,
     2 FCC_N,FCC_T,
     4 CUBIC_N,CUBIC_T,	 
     5 PROPS(1))
	 
      call GetTauSlips(RhoP,RhoF,RhoM,
     1 TauPass, TauCut, V0,	  
     2 PROPS(2),PROPS(3),PROPS(4))

      DO ISLIPS=1,18
       STATEV(ISLIPS+90)=RhoP(ISLIPS)
       STATEV(ISLIPS+108)=RhoF(ISLIPS)		
       STATEV(ISLIPS+126)=RhoM(ISLIPS)	
       STATEV(ISLIPS+144)=TauPass(ISLIPS)	  
       STATEV(ISLIPS+162)=TauCut(ISLIPS)
       STATEV(ISLIPS+180)=V0(ISLIPS)	   
      END DO

      call GetCSDHTauC(TAUPE,TAUSE,TAUCB,
     1 H, RhoCSD, TAUC, 	   
     2 PROPS(5:17))
	 
      DO ISLIPS=1,18
       STATEV(ISLIPS+198)=H(ISLIPS)
       STATEV(ISLIPS+216)=RhoCSD(ISLIPS)		
       STATEV(ISLIPS+234)=TAUC(ISLIPS)	   
      END DO
	 
      call GetGammaDot(Tau, TauPass, TauCut, V0, RhoM, 
     1 Vs, GammaDot, TauEff, TAUC,	    	   
     2 PROPS(18:20))	 
	 
      DO ISLIPS=1,18
       STATEV(ISLIPS+252)=GammaDot(ISLIPS)
       STATEV(ISLIPS+270)=Vs(ISLIPS)		
       STATEV(ISLIPS+288)=TauEff(ISLIPS)	   
      END DO
	 
      call GetRhoSSDEvolve(Tau, TauPass, TauCut, V0, RhoM, 
     1 GammaDot, TauEff, SSDDot, RhoSSD, RhoF,      	   
     2 PROPS(21:26))
	 
      call GetDSTRESS(DStress,GammaDot,dstran,Stress,dTIME, 
     1 FCC_Mu,FCC_Ohm,Cubic_Mu,Cubic_Ohm,	   
     2 PROPS(27:29))
 	  
      call GetDDSDDE(DDSDDE,Stress,	   
     2 PROPS(27:29))
c ------------------------------------------------	
c UPDATE ALL
      DO ISLIPS=1,18
       STATEV(ISLIPS)=STATEV(ISLIPS)+DTIME*SSDDot(ISLIPS)	
      END DO
      DO ISLIPS=1,6
       Stress(ISLIPS)=Stress(ISLIPS)+DStress(ISLIPS)	
      END DO
	  
c ------------------------------------------------		  
      return
      end subroutine UMAT
	  
c==========================================================	 
      include 'CalculateTauS.f'
      include 'VectorProjections.f'
      include 'GetRhoPFM.f'	  
      include 'GetTauSlips.f'
      include 'GetCSDHTauC.f' 
      include 'GetGammaDot.f'
      include 'GetRhoSSDEvolve.f'
      include 'GetDSTRESSN.f'	  
      include 'GetDDSDDEN.f'
	  




