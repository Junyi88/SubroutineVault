*****************************************************************************
**  UMAT FOR ABAQUS/STANDARD INCORPORATING ELASTIC BEHAVIOUR  FOR PLANE    **
**  STRAIN AND AXI-SYMMETRIC ELEMENTS.                                     **
*****************************************************************************
*****************************************************************************
**
**
**
*USER SUBROUTINE
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C
C
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
C
C      PARAMETER (M=3,N=3,ID=3,ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,
C     +          SIX=6.D0, NINE=9.D0, TOLER=0.D-6)

      PARAMETER (M=3,N=3,ID=3,ZERO=0.,ONE=1.,TWO=2.,THREE=3.,
     +          SIX=6., NINE=9., TOLER=0.E-6)  
C
      DIMENSION DSTRESS(6), DDS(6,6),XDSTRAN(6)
C
C      DOUBLE PRECISION TimeCount

      DIMENSION xb(6,60), xj(3,3), xjinv(3,3), dndx(20,3),xprod8(60,6),
     + dndloc(20,3),xbtrans(60,6), d(60,60), duv(60),
     + dfplas(60),gauss(8,3),dfplastmp(60),weight(8),
     + stressx(6),dtempx(60,60),
     + currentcoords(3,20),gausscoords(3,8),
     + xjcurrent(3,3),xjinvcurrent(3,3),
     + xstressvc(9,9),xstressvg(6,6),xstresssr(6,6),stressmat(3,3),
     + dndxcurrent(20,3),xg(9,60),xgtrans(60,9),xprod9(60,9),
     + dtemp2(60,60),temprot(3,3),
     + xn(20),xh(6,9),sigltp(6,6),sigutp(6,6),sigtp(6,6),xI(3,3),
     + xnmat(3,60),xnmatI(60,3),xrotgrad(3,3),xlatrot(3),
     + curv(3,3),xnat(20,3),curlfp(3,3),tmpv(3),xnat8(8,3),
     + xtime(10000)
     
C *** USER DEFINED ARRAYS ***
      real*8 :: spin(3,3), tempstrain(3,3) ,stressV(6),dstressinc(6),
     + dtotstran(6),totstran(6),dstressth(6),dstranth(6), SVARS(NSTATV*8),
     + Fdot(3,3), invF(3,3), F(3,3), L(3,3), LETEST1(3,3),LETEST2(3,3),LETESTV(6)
 
      integer debug
      character(len=*),parameter :: fmt20 = "(' ',20(F6.4,1X))",
     + fmt3="(3(' ',(ES11.3,1X)))",fmt6 = "(' ',6(F12.3,1X))",
     + fmt5 = "(' ',5(F12.8,1X))",fmt7 = "(' ',7(F12.8,1X))"

       PARAMETER ( xgauss = 0.577350269189626,  xweight = 1.0) 
     
       include 'mycommon.f'
      

       
c #include <SMAAspUserSubroutines.hdr>
 
    
!            
      debug = 0
      do while (debug == 1 .and. npt == 1 .and. noel == 1 .and. kinc ==1)
          ! wait here to attach debugger      
          
      end do
      !
      !! create a local array with ID=1 and SIZE=100

 
          

          
      
C*** ZERO ARRAYS ***
      tempstrain=0.;spin=0.
      

      
    
C     SPECIFY NUMBER OF USER STATE VARIABLES (nsvars = 8 x knsdv ---- 8 integration points) 
       
      knsdv = nsdv     
C      knsdv = 89     
      iphase = int(props(1))

C     WRITE PROPS INTO SVARS TO INITIALIZE (ONCE ONLY), 

      if (kinc == 1 .and. kstep==1) then
                

  
        do i=1,NSTATV
         STATEV(i) = 0. 
        end do 
      
   
         do i=1,3
         do j=1,3
          STATEV(j+(i-1)*3) = props(j+1+((i-1)*3)) !gmat
         end do
         end do
         

      !fp     
      STATEV(81) = 1.0
      STATEV(85) = 1.0
      STATEV(89) = 1.0

      STATEV(54) = 0.01 ! initial sessile SSD density
  
      svars = 0
      
  
       
      call MutexLock( 1 )      ! lock Mutex #1      
      do i=1,9      
          kFp(noel,npt,i) = 0.0
      end do
      call MutexUnlock( 1 )   ! unlock Mutex #1
              
      call MutexLock( 2 )      ! lock Mutex #2      
      ! use original co-ordinates X     
      do i =1,3
          kgausscoords(noel,npt,i) = coords(i)
      end do
      call MutexUnlock( 2 )   ! unlock Mutex #2
      
       call MutexLock( 3 )      ! lock Mutex #3     
      do i=1,9      
          kcurlFp(noel,npt,i) = 0.0
      end do
      call MutexUnlock( 3 )   ! unlock Mutex #3
      
      end if    

C    ZERO ARRAYS
      xI=0.            
      DO I=1,3
         xI(I,I)=1.
      END DO

100   FORMAT (4(2X, E20.5)) 
150   FORMAT (4(2X, E20.5))     

CC    LOOP OVER EIGHT INTEGRATION POINTS
       stressx = 0.; 
       dtempx = 0.0;
       F=0.  
       invF = 0. 
       Fdot = 0.
       L=0.
C   DETERMINE DEFORMATION AND VELOCITY GRADIENTS - CORRECTED ET 20/05/15     
       F = DFGRD0   
       Fdot = (DFGRD1-DFGRD0)/DTIME 
       CALL lapinverse(F,3,info,invF) 
   
C   CALL KMAT FOR MATERIAL BEHAVIOUR - Global stiffness matrix C
C      and stress 
       L = matmul(Fdot,invF)

   !   call KMLT(Fdot,invF,L)
      
 
      DO i=1,9
          statev(37+i) = kcurlfp(noel,npt,i)
      END DO
      
      write(*,*) 'QZ1'
      call kmat(dtime,NSTATV,STATEV,xI,NOEL,NPT,knsdv,time,F,L,iphase,
     +    DDSDDE,stressV,dstressinc,totstran,dtotstran,
     +    TEMP,DTEMP,vms,pdot,pnewdt)
      write(*,*) 'QZ6' 
   
C    RECOVER stress for calculating residual force
      DO K=1,6
           stress(K)=STATEV(47+K)
      END DO
      
      ! an attempt to calculate the same LE that abaqus calculates
      ! only seems to agrees when no rotation 
      !LETEST1 = 0.; LETEST2 = 0.;   
      !LETEST1 = 0.5*Log(matmul(F,transpose(F)))
      !LETEST2 = -0.5*Log(matmul(transpose(invF),invF))
      !
      !LETESTV=0.;   
      !CALL kmatvec6(LETEST1,LETESTV)
      !do k =1, 6
      !    
      !    statev(56+k) = LETESTV(k)    
      !end do
      !LETESTV=0.;
      ! CALL kmatvec6(LETEST2,LETESTV)
      !do k =1, 6
      !    
      !    statev(62+k) = LETESTV(k)    
      !end do
      
      call MutexLock( 1 )      ! lock Mutex #1 
      DO i=1,9                                                      
          kFp(noel,npt,i)= statev(80+i)
      END DO
      call MutexUnlock( 1 )      ! lock Mutex #1 
      
      
      write(*,*) 'QZ7'       
      
      IF (npt == 8 ) THEN ! update curl Fp
      
             
C=======================================================================   
C   SPECIFY GAUSS POINT LOCAL COORDS, AND WEIGHTING CONSTANTS 
      
      INCLUDE 'kgauss.f'     
      xnat8 = xnat(1:8,:) 


C=======================================================================          
         nsvars = nintpts * nsdv ! 8 x 89 = 712       
         
          
         DO kint =1,8
             
             DO i=1,3         
                 gausscoords(i,kint) = kgausscoords(noel,kint,i)                          
             END DO
         
             DO i=1,9          
                 svars(80 + i + 89*(kint-1)) = kFp(noel,kint,i)         
             END DO
         END DO
   
         
C=======================================================================   
C   A FULL INTEGRATION GRADIENT SCHEME
      write(*,*) 'QZ8' 
      CALL kcurl(nsvars,svars,knsdv,xnat8,gauss,gausscoords)
      write(*,*) 'QZ9'       
      call MutexLock( 3 )      ! lock Mutex #1 
      DO kint =1, 8
          DO i=1, 9
              kcurlFp(noel,kint,i) = svars(37+i + 89*(kint-1))
          END DO
      END DO
      call MutexUnlock( 3 )      ! lock Mutex #1 
    
            
      END IF
C======================================================================= 
          
      RETURN
      END
**

       
 
      include 'kmat.f'
      include 'uexternaldb.f'
C      include 'umatetc.f' 
C      include 'kumatht.f'     
      !include 'khetval.f'
      include 'kdirns.f'
      include 'kslip0.f'
      include 'kslip5.f'
      include 'kslip6.f'
      include 'kgndl2.f'
      include 'ksvd2.f'
      include 'kcurl.f'
      include 'kshapes.f'
      include 'utils.f'

            
            