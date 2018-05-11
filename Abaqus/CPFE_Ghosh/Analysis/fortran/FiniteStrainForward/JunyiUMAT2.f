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
	 
      include 'DeclareParametersSlipsO.f'
      
      INTEGER:: ISLIPS, I, J, NDUM1
      real*8 :: TAU(18), TAUPE(12), TAUSE(12), TAUCB(12)	  
      real*8 :: RhoP(18),RhoF(18),RhoM(18),RhoSSD(18)
      real*8 :: TauPass(18), TauCut(18), V0(18)
      real*8 :: H(12), RhoCSD(12), TAUC(18) 
      real*8 :: Vs(18) , GammaDot(18) , TauEff(18), SSDDot(18)
      real*8 :: DStress(6) 
	  
      real*8 :: ORI_ROT(3,3)
	  
c ------------------------------------------------	  
C
C     CALCULATE VELOCITY GRADIENT FROM DEFORMATION GRADIENT.
C     REFERENCE: Li & al. Acta Mater. 52 (2004) 4859-4875
C     
      Real*8:: FTINV(3,3),STRATE(3,3),VELGRD(3,3),AUX1(3,3),ONEMAT(3,3)
      PARAMETER (ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,SIX=6.0D0)
      DATA NEWTON,TOLER/10,1.D-6/
	  
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
       DO ISLIPS=1,18
          STATEV(ISLIPS+108)=PROPS(NPROPS+9)
       END DO

       DO ISLIPS=1,18
          STATEV(ISLIPS+108)=PROPS(NPROPS+9)
       END DO
	   
      ENDIF





c ------------------------------------------------	 
      return
      end subroutine UMAT

      include 'UTILS1.f'