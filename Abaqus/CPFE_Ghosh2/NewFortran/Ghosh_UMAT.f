      SUBROUTINE UMAT(stress,statev,ddsdde,sse,spd,scd,
     1 rpl, ddsddt, drplde, drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
	 
      include 'aba_param.inc'
c
#include <SMAAspUserSubroutines.hdr>
      CHARACTER*8 CMNAME
c      EXTERNAL F

      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     2 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     3 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
	 
      include 'DeclareParameterSlipsO.f'
      
      INTEGER:: ISLIPS, I, J, NDUM1, NA, NB, ICOR, ISL
      real*8 :: TAU(18), TAUPE(12), TAUSE(12), TAUCB(12)
      real*8 :: SLIP_T(54), IBURG, CFP(3,3)	  
      real*8 :: RhoP(18),RhoF(18),RhoM(18),RhoSSD(18)
      real*8 :: TauPass(18), TauCut(18), V0(18)
      real*8 :: H(12), RhoCSD(12), TAUC(18) 
      real*8 :: Vs(18) , GammaDot(18) , TauEff(18), SSDDot(18)
      real*8 :: DStress(6) , KCURLLOCAL(6)
      real*8 :: MXSLIP=1.0e-3
      real*8 :: ORI_ROT(3,3), SPIN_TENSOR(3,3)
      real*8:: dFP(9), dRhoS(18),dRhoET(18),dRhoEN(18)
c ------------------------------------------------	  
C
C     CALCULATE VELOCITY GRADIENT FROM DEFORMATION GRADIENT.
C     REFERENCE: Li & al. Acta Mater. 52 (2004) 4859-4875
C     
      real*8,parameter  :: zero=1.0e-16,xgauss = 0.577350269189626
      real*8,parameter  :: xweight = 1.0
      integer, parameter :: TOTALELEMENTNUM=159520
c  1728 853200
      Real*8:: FTINV(3,3),STRATE(3,3),VELGRD(3,3),AUX1(3,3),ONEMAT(3,3)
      PARAMETER (ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,SIX=6.0D0)
      DATA NEWTON,TOLER/10,1.D-6/
      Real*8:: gausscoords(3,8)
      real*8 :: kgausscoords, kFp, kcurlFp, kDGA
      real*8:: xnat(20,3),xnat8(8,3),gauss(8,3), DGA(18)
      real*8:: svars(144)
c XDANGER
      COMMON/UMPS/kgausscoords(TOTALELEMENTNUM,8,3),
     1 kFp(TOTALELEMENTNUM,8, 9),
     1 kcurlFp(TOTALELEMENTNUM, 8, 9)
	 
c -------------------------------------------------
c Initialisation
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
        NDUM1=(ISLIPS-1)*3+220
        NA=NDUM1+1
        NB=NDUM1+3		        
        call ROTATE_Vec(ORI_ROT,FCC_NPE0(1:3,ISLIPS),STATEV(NA:NB))
       END DO

c ---- S_SE
       DO ISLIPS=1,12
        NDUM1=(ISLIPS-1)*3+256
        NA=NDUM1+1
        NB=NDUM1+3		        
        call ROTATE_Vec(ORI_ROT,FCC_SSE0(1:3,ISLIPS),STATEV(NA:NB))
       END DO
c ---- N_SE
       DO ISLIPS=1,12
        NDUM1=(ISLIPS-1)*3+292
        NA=NDUM1+1
        NB=NDUM1+3		        
        call ROTATE_Vec(ORI_ROT,FCC_NSE0(1:3,ISLIPS),STATEV(NA:NB))
       END DO	   
	   
c ---- S_CB
       DO ISLIPS=1,12
        NDUM1=(ISLIPS-1)*3+328
        NA=NDUM1+1
        NB=NDUM1+3		        
        call ROTATE_Vec(ORI_ROT,FCC_SCB0(1:3,ISLIPS),STATEV(NA:NB))
       END DO
c ---- N_CB
       DO ISLIPS=1,12
        NDUM1=(ISLIPS-1)*3+364
        NA=NDUM1+1
        NB=NDUM1+3		        
        call ROTATE_Vec(ORI_ROT,FCC_NCB0(1:3,ISLIPS),STATEV(NA:NB))
       END DO	      
c ---- Rotate STIFFNESS TENSOR
        call ROTATE_COMTEN(ORI_ROT,PROPS(28:48),STATEV(164:184))

c--- Do Stuff
       STATEV(401)=1.0
       STATEV(405)=1.0
       STATEV(409)=1.0

c XDANGER
c       DO ISLIPS=1,18
c        STATEV(409+ISLIPS)=(1.0e9)*(1.0e-12)
c        STATEV(429+ISLIPS)=(1.0e9)*(1.0e-12)
c        STATEV(447+ISLIPS)=(1.0e9)*(1.0e-12) 
c       END DO	 
       DO ISLIPS=1,18
        STATEV(409+ISLIPS)= 0.0
        STATEV(429+ISLIPS)= 0.0
        STATEV(447+ISLIPS)= 0.0
       END DO	 
c      DO ISLIPS=1,6
c        STATEV(271+ISLIPS)=0.0
c        STATEV(289+ISLIPS)=0.0
c        STATEV(307+ISLIPS)=0.0   
c       END DO	
	   


      IF (KSTEP.LE.1) THEN
        call MutexLock( 3 )      ! lock Mutex #2      
      ! use original co-ordinates X     
        do i =1,3
          kgausscoords(noel,npt,i) = coords(i)
          statev(480+I) = coords(i)
        end do
        DO kint =1,8 
          DO i=1,3         
           gausscoords(i,kint) = kgausscoords(noel,kint,i)                          
          END DO 
         END DO	  
	  
        kfP(noel,npt,1)=1.0
        kfP(noel,npt,5)=1.0
        kfP(noel,npt,9)=1.0
        call MutexUnlock( 3 )   ! unlock Mutex #2
      ELSE
        call MutexLock( 3 )      ! lock Mutex #2      
      ! use original co-ordinates X     
        do i =1,3
          kgausscoords(noel,npt,i) = statev(480+I)
        end do
        DO kint =1,8 
          DO i=1,3         
           gausscoords(i,kint) = kgausscoords(noel,kint,i)                          
          END DO 
         END DO	  
	  
        kfP(noel,npt,1)=1.0
        kfP(noel,npt,5)=1.0
        kfP(noel,npt,9)=1.0
        call MutexUnlock( 3 )   ! unlock Mutex #2
      ENDIF	  
	  
      ENDIF
c ---


c -------------------------------------------------
      CALL UMATTEMPLATE(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)	 
	 
      return
      end subroutine UMAT
	  
c --------------------
      include 'uexternaldb.f'
      include 'UMAT_TemplateTest.f'
      include 'StiffnessTensorTools.f'
      include 'RotateSlipSystems.f'
