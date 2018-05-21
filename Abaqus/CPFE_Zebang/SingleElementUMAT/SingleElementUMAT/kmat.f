**************************************************
*           USER SUBROUTINE KMAT                 *
**************************************************  
C
C ** HCP, BCC, FCC ALL IN ONE: LARGE DEFORMATION **
C ** SLIP DIRECTIONS AND NORMALS NOT UPDATED **
C
C
      SUBROUTINE kmat(dtime,nsvars,usvars,xI,jelem,kint,knsdv,time,F,
     + L,iphase,C,stressV,dstressinc,totstran,dtotstran,
     + TEMP,DTEMP,vms,pdot, pnewdt)
C
      INCLUDE 'ABA_PARAM.INC'
      
      INTEGER, parameter:: M=3,N=3,L0=24,L1=24,L2=12,KM=6,KN=6
      REAL*8,parameter  :: zero=1.0e-8
C
      !scalars
      INTEGER,intent(in) :: nsvars,jelem,kint,knsdv,iphase
      REAL*8,intent(in) :: dtime
      
      !arrays
      REAL*8,intent(inout):: usvars(nsvars)
      REAL*8,intent(in) :: time(2),F(3,3),L(3,3),xI(3,3)
      REAL*8,intent(out) :: C(6,6)     
      



C     *** USER DEFINED ARRAYS ***
      REAL*8 :: stressM(3,3),plasStrainRate(3,3), totalStrainInc(6),
     +  prod(M),tempNorm(M), tempDir(M), 
     + xRot(M,N),VRES1(6),VRES2(6),plasStrainInc2(6),Lp(3,3),
     + totplasstran(6),Le(3,3),
     + tSigma(6,6),tStran(6,6),tSigmainv(6,6),prod6(6,6),
     + xsnt(M,N),xsnv(KM),xnsv(KM),xsnnst(KM,KN),xIden6(KM,KN),
     + xnst(M,N),tmat(KM,KN),
     + trialstress(6),xjfai(KM,KN),
     + xjfaiinv(KM,KN),stressV(6),fai(6),trialstressM(3,3),
     + dstressinc(6),xu(3,3),xuinver(3,3),
     + xstressmdef(M,N),xstressdef(6),
     + xStiff(6,6),xStiffdef(6,6),elasspin(3,3),plasspin(3,3),
     + tSig(6,6),tStr(6,6),tSiginv(6,6), T(6,6),Tinv(6,6),
     + tStraninv(6,6),gmatinv(3,3),devstress(3,3),
     + compliance(6,6), print1(3,3),print2(3,3),tmpvec(6),print3(3,3),
     + print4(3,3),xrote(3,3),dtotstran(6),totstran(6),
     + totStrainRate(3,3),spin(3,3), tempstrain(3,3),curlfp(3,3),
     + curlfe(3,3),fp(3,3),fe(3,3),feinv(3,3),fpinv(3,3),xfinv(3,3)
     + ,tempNormGND(3), rhogndsys(30), 
     + temp3X3(3,3),temp3X3a(3,3), temp6x1(6), temp6x1a(6), temp3x1(3),
     + temp6x6(6,6), gv(9), rhogndnew(12),
     + tauR(12), sslip(3),nnorm(3),tnorm(3),burgX(M),
     + screw(3,3),edgeN(3,3),edgeT(3,3),scol(9),encol(9),etcol(9)
       
      REAL*8 ::  vresth(6), dstranth(6),thermat(3,3),expanse33(3,3),
     + sigthv(6),sigthm(3,3), gmatinvnew(3,3) 
        
      REAL*8,dimension(:,:),allocatable :: xNorm,xDir
      REAL*8,dimension(:),allocatable :: tau, gammadot, gndall, gndold,
     +  burgerv, tauc, tau2 ,S,gndtotal, rhoOutput
      INTEGER,dimension(:),allocatable :: ids
      INTEGER MZ, NZ, slipsys
      DOUBLE PRECISION,dimension(:,:),allocatable :: V,U,Sinv,A,AM,
     + Ainv,diagm,tempA, UT      
      
      REAL*8 rho, TAUCSSD 
      REAL*8 LpFeinv(3,3), matrix(3,3), update(3,3)
!      REAL(selected_real_kind(8,70)) :: rhoc

      REAL*8 gammast,rhoGND,rhoBasal,rhoPrism,rhoAPyrm,rhoCAP,rhoCAS,
     + gndedge, gndscrew
      
      REAL*8 gama(L0)
      
      character(len=*),parameter :: fmt2 = "(24(' ',(I2,1X)))",
     + fmt3="(3(' ',(ES11.3,1X)))",fmt24 = "(24(' ',(ES11.3,1X)))"
      
     
      
      
      !common/therm/ytemp,ytemprate


!----------------Note for use---------------
C     Tempature must larger than 273K since the material properties are interploted 
C     for only those great than room temperature, i.e. 20 degree C
!--------------------------------------------------------------------------------
!   
C removed (kint-1)*knsdv+ from all the SDV allocation, SDV reduced from 712 to 89, with SDV(1) = sum over k=1,8 of SDV(1+89*(k-1))
! individual integration point values of Fp stored in common block   ET updated 10/04/2015

C     *** MATERIAL CONSTANTS, ETC ***       
      SELECT CASE(iphase)
      CASE(0) !hcp
      nSys = L0
      
      TAUCSSD = 280  
      crssratio = 3.01
      
      
C      write(6,*)'TAUCSSD,crssratio',TAUCSSD,crssratio
!----------------------------------------------------------------------------------------
      
      burger1 = 2.95E-4 ! Metallic atomic radius of 160pm - 
                        ! converted to microns and DOubled.
      caratio = 3.01 !Updated for Zr
      
      gammast = 0.0 !Coefficient for RhoSSD evolution   !Changed 24/03/14
	  
	  
	  E1 = 84745; E3 = 119789
      v12 = 0.46
      v13 = 0.32     
      G12 = E1/(2.0*(1.0+v12)); G13 = 40000   

      rhossd = usvars(54)
      rhoGND= usvars(26)
      rho=rhoGND+rhossd
      XTAUC1 = TAUCSSD + (1.0 * G12 * burger1 * sqrt(rho))
      !XTAUC1 = 240 
	  
      XTAUC0 = 1.125*XTAUC1
      XTAUC2 = crssratio*XTAUC1 
      burger2 = caratio*burger1
        
      alpha1 = 9.5D-6; alpha2 = alpha1; alpha3 = 0.5895*alpha1


      allocate(xNorm(L0,M),xDir(L0,M),tau(L0),gammadot(L0),gndall(L0),
     + gndold(L0),burgerv(L0),tauc(L0),ids(L0),tau2(L0),STAT=ialloc)
     
      burgerv(1:12) = burger1; burgerv(13:L0) = burger2
      tauc(1:3) = XTAUC0; tauc(4:12) = XTAUC1; tauc(13:L0) = XTAUC2        
            
      case(1) !bcc
      nSys = L1
      XTAUC = 280.0
      burger = 2.86E-4
	  
      !bcc properties Kim and Rokhlin (2009) J.Acoust.Soc.Am.
      E1 = 32024; E3 = E1 
      v12 = 0.4556; v13 = v12
      G12 = 54900; G13 = G12
      
      alpha1 = 9.5e-6; alpha2 = alpha1; alpha3 = 0.5895*alpha1
      
      allocate(xNorm(L1,M),xDir(L1,M),tau(L1),gammaDot(L1),gndall(L1),
     + gndold(L1),burgerv(L1),tauc(L1),ids(L1),tau2(L1),STAT=ialloc)
     
      burgerv = burger
      tauc = xtauc
     
      case(2) !fcc
      nSys = L2
      XTAUC = 230.0 ! (563.4/sqrt(6) = 230.0)
      burger = 3.5072e-4
      
   !   E1 = 207.0E+3; E3 = E1
   !   v12 = 0.28; v13 = v12
   !   G12 = E1/(2.0*(1.0+v12)); G13 = G12
      
      E1 = 207.0E+3; E3 = E1
      v12 = 0.28; v13 = v12
      G12 = E1/(2.0*(1.0+v12)); G13 = G12
      
      alpha1 = 13.0e-6; alpha2 = alpha1; alpha3 = alpha1
      
      gammast = 0.05
     
      
      allocate(xNorm(L2,M),xDir(L2,M),tau(L2),gammaDot(L2),gndall(L2),
     + gndold(L2),burgerv(L2),tauc(L2),ids(L2),tau2(L2),STAT=ialloc)
     
      burgerv = burger
      tauc = xtauc
      case(3) !carbide
      nSys = L2
      XTAUC = 230.0D10
      burger = 3.5072e-4
      
      E1 = 207.0E+4; E3 = E1
      v12 = 0.28; v13 = v12
      G12 = E1/(2.0*(1.0+v12)); G13 = G12
      
      alpha1 = 4.5e-6; alpha2 = alpha1; alpha3 = alpha1
      
      allocate(xNorm(L2,M),xDir(L2,M),tau(L2),gammaDot(L2),gndall(L2),
     + gndold(L2),burgerv(L2),tauc(L2),ids(L2),tau2(L2),STAT=ialloc)
     
      burgerv = burger
      tauc = xtauc
      
      case default
      WRITE(6,*)
      WRITE(6,*)"Not sure what crystal type. Material constants."
      END SELECT

      
      h = 0.0      
      kslip = 6
      !0 = original slip rule with no GND coupling.
      !5 = original slip rule with GND coupling, with a new parameter to make rhoc more physically reasonable, slip included.
      !6 = Slip rule with constant Activation Volume
      
      gndon = 0 !1=on, 0=off



C     *** ZERO ARRAYS ***

      result = 0.0; totstran = 0.0; totplasstran = 0.0;devstress=0.
      plasStrainInc2=0.
      xStiff=0.0; xStiffdef=0.0; C=0.0; trialstressM=0.0
      tmat=0.0;xjfai=0.
      xjfaiinv=0.;plasStrainRate=0.;trialstress=0.; stressV=0.
      fai=0.;dstressinc=0.;xIden6=0.; xu=0.
      xuinver=0.;xRot=0.;prod=0.
      T=0; Tinv=0 !NEW ADDITIONS
      tSig=0.;tStr=0.;tSiginv=0.;tStraninv=0.
      tSigmainv=0.;gmatinv=0.;Lp = 0.;compliance=0.;Le = 0.
      print1=0.;print2=0.; tmpvec = 0.;print3 = 0.;print4 = 0.
      thermat =0.
      vresth=0.; dstranth=0.;expanse33=0.;sigthv=0.;sigthm=0.
      xrote=0.;tempstrain=0.;spin=0.;totStrainRate=0.
      curlfp=0.;curlfe=0.;fp=0.;fe=0.;feinv=0.;fpinv=0.;xfinv=0.
      gmatinvnew = 0.
      
      
      DO I=1,KM; xIden6(I,I)=1.; END DO      
      DO I=1,M; xRot(I,I) = 1.; END DO
      
      
      !Zero allocate arrays from above!
      xNorm=0.; xDir=0.; tau=0.; gammadot=0.
      gndall=0.; gndold=0.; ids=0; tau2=0.
      gndcas=0.;gndcap=0.;gndapy=0.;gndapr=0.;gndab=0.;rhognd=0.
      
      !Ben's experimental data storage groups
      capyramedge=0.;capyramscrew=0.;apyramedge=0.;aprismedge=0.
      abasedge=0.;ascrew=0.  
      
      thermat(1,1) = alpha1; thermat(2,2)=alpha2; thermat(3,3) = alpha3
      
C     *** SET UP ELASTIC STIFFNESS MATRIX IN LATTICE SYSTEM ***   
      compliance(1,1:3) = (/1./E1,-v12/E1,-v13/E1/)
      compliance(2,2:3) =         (/1./E1,-v13/E1/)
      compliance(3,3:3) =                 (/1./E3/)
      compliance(4,4:4) =                       (/1./G12/)
      compliance(5,5:5) =                       (/1./G13/)
      compliance(6,6:6) =                       (/1./G13/)
C
      DO i=2,6
         DO j=1,i-1
            compliance(i,j)=compliance(j,i)
         END DO
      END DO
      
      CALL lapinverse(compliance,6,info,xStiff)
     
    
C     *** INITIALIZE USER ARRAYS ***

      DO i=1,3
        DO j=1,3
         gmatinv(i,j) = usvars(j+(i-1)*3)
        END DO
      END DO
C
      p = usvars(10)
C
      DO i=1,6
        totplasstran(i) = usvars(10+i)
      END DO
C
      DO i=1,6
        totstran(i) = usvars(16+i)
      END DO
C        
      DO i=1,6
        xstressdef(i) = usvars(47+i)
      END DO
C
      rhognd = usvars(37)  
C
      !r= usvars((kint-1)*knsdv+56)       
      rhossd = usvars(54)
      
      r = usvars(56)
                
      DO i=1,nSys
        gndold(i) = usvars(56+i)
      END DO
C
      DO i=1,3
        DO j=1,3
          fp(i,j) = usvars(86+j+((i-1)*3))
        END DO
      END DO              
C
      DO i=1,3
        DO j=1,3 
         curlfp(i,j) = usvars(37+j+(i-1)*3)
        END DO
      END DO
       
C    -------------------------------------     
      DO i=1,L0
        gama(i) = usvars(95+i)
      END DO         
C    -------------------------------------
      
      
C     *** DIRECTIONS FROM LATTICE TO DEFORMED SYSTEM ***
        
      CALL kdirns(gmatinv,iphase,nSys,xDir,xNorm)
        

C     *** STIFFNESS FROM LATTICE TO DEFORMED SYSTEM ***

      CALL rotord4sig(gmatinv,tSig)
      CALL rotord4str(gmatinv,tStr)
      CALL lapinverse(tSig,6,info2,tSiginv)

      prod6 = matmul(tSiginv,xStiff)      
      xStiffdef = matmul(prod6,tStr)

      expanse33 = matmul(matmul(gmatinv,thermat),transpose(gmatinv))
C      expanse33 = expanse33*ytemprate*dtime !dstrain = alpha*dT
      expanse33 = expanse33*DTEMP !dstrain = alpha*dT
      
      CALL kmatvec6(expanse33,dstranth)
      dstranth(4:6) = 2.0*dstranth(4:6)
            

    
C     *** DETERMINE INCREMENT IN TOTAL STRAIN (6X1 ***     

      tempstrain=(L+transpose(L))*0.5*dtime
      spin=(L-transpose(L))*0.5 

      CALL kmatvec6(tempstrain,dtotstran)
      dtotstran(4:6) = 2.0*dtotstran(4:6)


C     *** COMPUTE TRIAL STRESS ***

      stressV = xstressdef ! old stress
      trialstress = stressV+matmul(xStiffdef,dtotstran)-
     + matmul(xStiffdef,dstranth)            
      CALL kvecmat6(trialstress,trialstressM) 
            
      CALL kvecmat6(stressV,stressM) 
      trialstressM = trialstressM + (matmul(spin,stressM) - 
     + matmul(stressM,spin))*dtime 


 
C     *** CALCULATE RESOLVED SHEAR STRESS ON A SLIP SYSTEMS ***


      DO I=1,nSys
        tempNorm = xNorm(I,:); tempDir = xDir(I,:)
        prod = matmul(trialstressM,tempNorm)
        tau(I)= dot_product(prod,tempDir)
        IF(tau(I) .LT. 0.0) THEN
          tau(I) = -1.E0*tau(I)
          DO K=1,3
            xDir(I,K)=-1.E0*xDir(I,K)
          END DO
        END IF
      END DO
        
          
      xtau = maxval(tau/tauc)  


C     *** PLASTIC DEFORMATION ***

      IF (xtau >= 1.0 ) THEN
      
C      IF (kint == 1) THEN
C        IF (xtau .lt. 1.001) THEN
C            PRINT *, 'Time = ' , time(1)
C        END IF
C      END IF 
      debug = 0
      do while (debug == 1)
          !paused for debuggin
          end do
      
      faivalue=1.
      xacc=1.e-8
      iter=0


C     *** USE NEWTON METHOD TO DETERMINE STRESS INCREMENT ***

      DO WHILE (faivalue .gt. xacc)      
      
      iter=iter+1


      !============================================================================   
      !  Slip rule:
      !  Returns Lp and tmat required to define the material jacobian.  
      !============================================================================  

      IF (kslip == 0) THEN
      !Original slip rule with no GND coupling i.e. using alpha and beta
      CALL kslip0(xNorm,xDir,tau,tauc,caratio,dtime,nSys,r,iphase,Lp,
     + tmat)
     
      ELSE IF (kslip == 5) THEN
      !Original slip rule with GND coupling            
      CALL kslip5(xNorm,xDir,tau,tauc,burgerv,rhossd,gndold,caratio,
     + dtime,nSys,r,iphase,Lp,tmat)
	 
      ELSE
      !Original slip rule with no GNDs eg no hardending
      CALL kslip6(xNorm,xDir,tau,tauc,burgerv,caratio,      
     + dtime,nSys,r,iphase,Lp,tmat,TEMP,gammaDot)
           
      END IF
      
      if(any(tmat /= tmat) .or. any(tmat-1 == tmat)) then
          pnewdt = 0.5 ! if sinh( ) has probably blown up then try again with smaller dt
          write(*,*) "W! tmat=NaN or inf: jelem, kint, time: ",
     +     jelem, kint, time
          return
      end if
  
      
      
      !============================================================================  


C     *** DETERMINE PLASTIC STRAIN INCREMENETS FOR UEL

      plasStrainRate = (Lp+transpose(Lp))*0.5*dtime
      CALL kmatvec6(plasStrainRate,plasStrainInc2)
      plasStrainInc2(4:6) = 2.0*plasStrainInc2(4:6)            



C     *** CALCULATE THE STRESS INCREMENT ***

      xjfai =  xIden6 + matmul(xStiffdef,tmat)
      CALL lapinverse(xjfai,6,info3,xjfaiinv)
!      IF(info3 /= 0) write(6,*) "inverse failure: xjfai in kmat"
      fai = trialstress - stressV - matmul(xStiffdef,plasStrainInc2)
      dstressinc = matmul(xjfaiinv,fai)
      stressV = stressV + dstressinc
      CALL kvecmat6(stressV,stressM)      
      faivalue = sqrt(sum(fai*fai))        


C     *** UPDATE RESOLVED SHEAR STRESS ACCORDING TO NEW STRESS ***
       
      DO I=1,nSys
      
          tempNorm = xNorm(I,:); tempDir = xDir(I,:)    
          prod = matmul(stressM,tempNorm)
          tau(I)= dot_product(prod,tempDir)

          IF(tau(I) < 0.0) THEN
            tau(I) = -1.E0*tau(I)
            DO K=1,3
              xDir(I,K)=-1.E0*xDir(I,K)
            END DO
          END IF
          !=============================== ids
          IF(tau(i)/tauc(i) >= 1.0) THEN
            ids(i) = i; tau2(i) = tau(i)/tauc(i)
          END IF
          !===============================
          
      END DO
        
      xtau = maxval(tau/tauc) 
          
      
      IF (iter .gt. 50) THEN
          WRITE(*,*) "WARNING NEWTON LOOP NOT CONVERGED: jelem, kint, 
     +     time:", jelem, kint, time
          pnewdt = 0.5
          return
          !CALL XIT 
      END IF
                 
      !*** THE END OF NEWTON ITERATION ***
      END DO



C     *** NOW CALCULATE THE JACOBIAN***

      C = matmul(xjfaiinv,xStiffdef)
       
       
C     *** ROTATE STRESS BACK TO GLOBAL SYSTEM *** 

      xstressmdef = stressM


C     *** UPDATE OUTPUT VARIABLES ***

      plasStrainrate=(Lp+transpose(Lp))*0.5       
      pdot=sqrt(2./3.*sum(plasStrainrate*plasStrainrate))
      
      p = p + pdot*dtime
      r = r + h*pdot*dtime  
      
C      IF (kint == 1) THEN
C        PRINT *, pdot
C      END IF 
      
C     *** UPDATE PLASTIC DEFORMATION GRADIENT
    
      print2 = 0.; print3 = 0.
      print2 = xI - Lp*dtime      
      CALL kdeter(print2,deter)      
      IF (deter /= 0.0) THEN
         CALL lapinverse(print2,3,info4,print3)
         fp = matmul(print3,fp)
      ELSE
         fp = fp
      END IF  
      
      

      
    !=========================================================================
    ! SSD Evolution 
     
      rhossd = rhossd + (gammast*pdot*dtime)
      
    !=========================================================================
               
C     *** ELASTIC DEFORMATION ***     
      ELSE
      xstressmdef = trialstressM
      C = xStiffdef      
      END IF
      
          
      
      
!     CALL kvecmat6(xstressdef,stressM) 
!     xstressmdef = xstressmdef + (matmul(spin,stressM) - matmul(stressM,spin))*dtime
      
      CALL kmatvec6(xstressmdef,xstressdef) !output stress
      devstress = xstressmdef - 1./3.*trace(xstressmdef)*xI
      vms = sqrt(3./2.*(sum(devstress*devstress))) !von mises stress 

    !   call MutexLock( 1 )      ! lock Mutex #1   
      !write(*,*) gmatinv
      !write(*,*) kint, jelem 
    !  call MutexUnLock( 1 )      ! lock Mutex #1   





    !=========================================================================
    ! *** DETERMINE DENSITY OF GNDs
    ! Update P.Ashton November 2015
    !=========================================================================

      
      
      IF(MAXVAL(ABS(curlfp)) <= 1.0e-8) THEN 
        rhoGND = 0.;gndab  = 0.;gndapr = 0.;gndapy = 0.
      
      ELSE

      IF (gndon == 0) THEN !Switching GND evolution on and off
        rhoGND = 0.;gndab  = 0.;gndapr = 0.;gndapy = 0.
      
      ELSE   
      
      MZ=9;NZ=3*nSys
      
      ALLOCATE(A(mz,nz),AM(mz,nz),V(nz,nz),S(mz),U(mz,mz),Sinv(nz,mz),
     + Ainv(nz,mz),tempA(nz,mz),UT(mz,mz),rhoOutput(nz),STAT=ialloc)
     
      A=0.;AM=0.;V=0.;S=0.;U=0.;Sinv=0.;Ainv=0.;tempA=0.;UT=0.
      
      gv = reshape(curlfp,(/9/))   
      

           
      tauR = tau2(1:nSys)
        
      DO i=1,12
        IF(tauR(i) > 1) THEN
           tauR(i) = tauR(i) - 1
        END IF 
      END DO 
      
      tausum = sum(tauR)
      tauBas = sum(tauR(1:3))
      tauPris= sum(tauR(4:6))
      tauPyr = sum(tauR(7:12))
       
  !    DO slipSys=1,nsys ! Cycle through each system
      
  !    IF (ids(slipSys) > 0) THEN
   
        DO I=1,nSys ! 1-12 for <a> type
        
            burgX = xDir(I,:)
            burgX = burgX*burger1 
         
            sslip = xDir(I,:); nnorm = xNorm(I,:) 
            CALL CrossProd(sslip,nnorm,tnorm) 
         
            CALL DyadicProd(sslip,burgX,screw)
            CALL DyadicProd(nnorm,burgX,edgeN)
            CALL DyadicProd(tnorm,burgX,edgeT)
         
            CALL Convert2Col(screw,scol)
            CALL Convert2Col(edgeN,encol)
            CALL Convert2Col(edgeT,etcol)
         
            DO J=1,9
                A(J,I)    = scol(J)
                A(J,I+nSys) = encol(J)
                A(J,I+(2*nSys)) = etcol(J)
            END DO
                
        END DO   ! Generate current A-matrix loop
      


  ! ***********************************************************************    
      ! Matrix inversion by singular value decomposition: A+ = [V][S+][U*]
  ! ***********************************************************************   
        CALL SVD(A,U,S,V,MZ,NZ) 
      
        DO i = 1, ubound(S,1)
            IF (S(i) > 1e-6) THEN
                Sinv(i,i)= 1.0 / S(i)
            END IF
        END DO 
        
            
        UT=transpose(U)
      
        NCOL=NZ;MX=NZ;NX=MZ 
        CALL KMLTM (V,Sinv,tempA,NCOL,MX,NX) 

      
        NCOL=MZ;MX=NZ;NX=MZ 
        CALL KMLTM (tempA,UT,Ainv,NCOL,MX,NX) 

        DO i=1,nz
        x = 0.
          DO j=1,9
           x = x + Ainv(i,j)*gv(j)
          END DO
        
        rhoOutput(i) = x
        END DO
 
        rhos = sum(rhoOutput(1:nSys))
        rhoen = sum(rhoOutput((nSys+1):(2*nSys)))
        rhoet = sum(rhoOutput((2*nSys+1):(3*nSys)))
        rhofinal = sqrt((rhos*rhos) +(rhoen*rhoen) +(rhoet*rhoet))

     
      rhoGND = rhofinal      



      

      END IF ! End of GND on/off switch   
      
      END IF ! if curl(Fp) < 1e-8


C     *** ORIENTATION UPDATE ***
      !Assuming that all rigid body rotatoin is lumped into Fe and that the elastic strians are small 
!     then the elastic spin is We = d(Fe)/dt inv(Fe)
      !L = We + Fe Lp inv(Fe) therefore 
      !We = L - Fe Lp inv(Fe)
      ! G(t+dt) = G(t) + We G(t)dt dt or an implicit update is G(t+dt)  = G(t)exp[We(t+dt)dt]  ~ inv[I - We(t+dt) dt] G(t) 
      
      ! We need Fe and inv(Fe) using F = Fe Fp gives Fe = F.inv(Fp)
      CALL kdeter(Fp,deter)      
      
      IF (deter /= 0.) THEN
         Fpinv = 0.
         CALL lapinverse(Fp,3,info5,Fpinv)
!         IF(info5 /= 0) write(6,*) "inverse failure: print3 in kmat"
         Fe = matmul(F,Fpinv)          
      ELSE
         write(*,*) "Error in orientation update: finding inv(Fp)",noel,
     +    npt, kinc
         call XIT 
      
      END IF  
      
      
      CALL kdeter(Fe,deter)      
      
      IF (deter /= 0.) THEN
         Feinv = 0.
         CALL lapinverse(Fe,3,info5,Feinv)
!         IF(info5 /= 0) write(6,*) "inverse failure: print3 in kmat"         
      ELSE
          write(*,*) "Error in orientation update: finding inv(Fe)",
     +     jelem,kint, time
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
         CALL lapinverse(matrix,3,info5,update)
         IF(info5 /= 0) write(*,*) "inverse failure: print3 in kmat"
         !print3 = 0.
         !print3 = gmatinv + dtime*matmul(elasspin,gmatinv)
         gmatinvnew = matmul(update,gmatinv)                  
       
         !write(*,*) "gmatinv, print3", gmatinv, print3
      
      ELSE         
         gmatinvnew = gmatinv
      write(*,*) "WARNING gmatinv not updated at noel,npt, kinc:",jelem,
     + kint, time
      END IF      

      gmatinv =  gmatinvnew            
      
       if (maxval(gmatinv) > 1) then
          write(*,*) "something very wrong with gmatinv"
          call XIT
       end if

C     ******************************************
C     Heat                     !Changed 19/06/2014
C      rho=7.8e-3
C      specHeat=460
C      beta=0.95
C      ytemprate=vms*pdot*beta/(rho*specHeat)
C     ******************************************     


C************************************************
C     Stored Energy Criterion David W.
C************************************************
C     Update Stored Energy
      Gdot=0
      DO i = 1, 6
          Gdot = Gdot + ABS(0.05 * xstressdef(i) * plasStrainInc2(i))
      END DO

C Combine SSDs and GNDs      
      DislDens = rhognd + rhossd
      DislDens = ABS(SQRT(DislDens))
      
      if (DislDens == 0) then
          write(*,*) "Dislocation Density Wrong!"
          DislDens = 0.01
       end if
      
      Gdot = Gdot / DislDens
      
      Gstored = usvars(126) + Gdot
      
C************************************************

C     *** UPDATE STATE VARIABLES *** ! Free: None

      DO i=1,3
        DO j=1,3
        usvars(j+(i-1)*3) = gmatinv(i,j)
        END DO
      END DO

      usvars(10) = p
 
      DO i=1,6
        usvars(10+i) = totplasstran(i) + 
     +                                  plasStrainInc2(i)
      END DO
C
      DO i=1,6
        usvars(16+i) = totstran(i) + dtotstran(i)
      END DO
      
        
      !accumulated slip
      DO i=1,L0
         gama(i)= gama(i)+gammaDot(i)*dtime
      ENDDO    
      
      
      !23-25 are rotations stored in UEL
    
C     usvars((kint-1)*knsdv+34) = xtau
     
      usvars(26) = rhognd !all
      usvars(27) = gndab  !a basal
      usvars(28) = gndapr !a prismatic  
      usvars(29) = gndapy !a pyramidal
      usvars(30) = gndcap !c+a primary    
      usvars(31) = gndcas !c+a secondary
     
      !Ben's experimental data storage groups
      !usvars((kint-1)*knsdv+32) = capyramedge !c+a pyramidal edge group
      !usvars((kint-1)*knsdv+33) = capyramscrew !c+a pyramidal screw group
      !usvars((kint-1)*knsdv+34) = apyramedge !a pyramidal edge
      usvars(35) = aprismedge !a prismatic edge
      usvars(36) = abasedge  !a basal edge
      usvars(37) = ascrew   !a all screw    

      !38-46 are curlfp terms calulated in gradient routine (kcurl)
      
      DO i=1,6
       usvars(47+i) = xstressdef(i)
      END DO

      usvars(32) = maxval(plasStrainrate)
      usvars(33) = pdot
      usvars(34) = xtau !XTAUC1 
      usvars(54) = rhossd
      usvars(55) = vms
      usvars(56) = r
     
      !GNDs on indiviual systems
      !max(nSys) is currently limited to 24. IF all 48 of bcc is needed, storage should be raised to match that!
      DO i=1,nSys
       usvars(56+i) = gndold(i)
      END DO
!
      
      DO i=1,3
       DO j=1,3 
        usvars(86+j+(i-1)*3) = fp(i,j)
       END DO
      END DO   
      
C     ----------------------------------------------      
      DO i=1,L0
        usvars(95+i) = gama(i)
      END DO         
C     ----------------------------------------------      

      usvars(126) = GStored

!C     *** WRITE RESULTS TO DUMMY UMAT!!
!
!      DO i=1,knsdv
!        ksdv(kint,i) = usvars(i)
!      END DO

           
      DEALLOCATE(xNorm,xDir,tau,gammadot,gndall,gndold,burgerv,
     + tauc,ids,tau2)

      RETURN

      END

