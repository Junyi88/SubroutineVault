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
      END DO
      ENDIF

c ------------------------------------------------	 
      DO ISLIPS=1,18
       RhoSSD(ISLIPS)=STATEV(ISLIPS)
      END DO
	  
c ------------------------------------------------	 	  
      call CalculateTauS(STRESS, TAU, TAUPE, TAUSE, TAUCB,
     +  FCC_N,FCC_S
     +  FCC_NPE,FCC_SPE
     +  FCC_NSE,FCC_SSE
     +  FCC_NCB,FCC_SCB
     +  CUBIC_N,CUBIC_S)
	 
      call GetRhoPFM(RhoP,RhoF,RhoM,
     1 RhoSSD,
     2 FCC_N,FCC_T,
     4 CUBIC_N,CUBIC_T,	 
     5 PROPS(1))
	 
      call GetTauSlips(RhoP,RhoF,RhoM,
     1 TauPass, TauCut, V0,	  
     2 PROPS(2),PROPS(3),PROPS(4))
	 
      call GetCSDHTauC(TAUPE,TAUSE,TAUCB,
     1 H, RhoCSD, TAUC, 	   
     2 PROPS(5:17))
	 
      call GetGammaDot(Tau, TauPass, TauCut, V0, RhoM, 
     1 Vs, GammaDot, TauEff, 	   
     2 PROPS(18:20))	 
	 
      call GetRhoSSDEvolve(Tau, TauPass, TauCut, V0, RhoM, 
     1 Vs, GammaDot, TauEff, SSDDot	    	   
     2 PROPS(21:26))
	 
      call GetDSTRESS(DStress,GammaDot,DStrains,Stress,dTIME, 
     1 FCC_Mu,FCC_Ohm,Cubic_Mu,Cubic_Ohm 	   
     2 PROPS(27:29))
 	  
      call GetDDSDDE(DDSDDE,Stress,	   
     2 PROPS(27:29))
c ------------------------------------------------	
c UPDATE ALL
      DO ISLIPS=1,18
       STATEV(ISLIPS)=STATEV(ISLIPS)+DTIME*SSDDot	
      END DO
      DO ISLIPS=1,6
       Stress(ISLIPS)=Stress(ISLIPS)+DStress	
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
      include 'GetDSTRESS.f'	  
      include 'GetDDSDDE.f'
	  




