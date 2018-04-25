INCLUDE './minpack/r1mpyq.f'
INCLUDE './minpack/r1updt.f'
INCLUDE './minpack/dogleg.f'
INCLUDE './minpack/qform.f'
INCLUDE './minpack/qrfac.f'
INCLUDE './minpack/enorm.f'
INCLUDE './minpack/dpmpar.f'
INCLUDE './minpack/hybrj.f'
INCLUDE './minpack/hybrj1.f'
INCLUDE './minpack/hybrd.f'
INCLUDE './minpack/hybrd1.f'
INCLUDE './minpack/fdjac1.f'
INCLUDE './Module3.f'
INCLUDE './uexternaldb.f'
!
!=====================================================================72
!
!  UMAT subroutine for crystal plasticity - Originaly form SNL & MSU, Tayloir Hytpothesis,
!  Modified by Xiang Zhang, Robby and Tung Phan from Vanderbilt University to 
!   1) incorporate new flow rule and evolution functions, 
!   2) can be run in abaqus in parallel.
!   3) incorporate a flow rule of dislocation climb
!   4) combine uel subroutine for cohesive zone model
!
!  Last edited by Xiang Zhang (03/01/2016) 
!	and Van-Tung Phan (02/01/2017)
!=====================================================================72
!
      subroutine umat ( stress,  statev,  ddsdde,  sse,     spd, &
                        scd,     rpl,     ddsddt,  drplde,  drpldt, &
                        strain,  dstrain, time,    dtime,   temp, &
                        dtemp,   predef,  dpred,   cmname,  ndi, &
                        nshr,    ntens,   nstatv,  props,   nprops, &
                        coords,  drot,    pnewdt,  celent,  dfgrd0, &
                        dfgrd1,  noel,    npt,     layer,   kspt, &
                        kstep,   kinc )
!
      use     dummyvar
      use     datatype
      use     dummyvar
      use     InitialState
!      use     InitialLocalRotMatrix
      use     SlipData
      use     DataType
      !include 'ABA_PARAM.INC'
      implicit none

      character*30 cmname
!
      integer ntens, nprops, nstatv, ndi, nshr, noel, npt, layer, kspt, kstep, kinc
      real*8  sse, spd, scd, rpl
      real*8  drpldt(ntens), temp, dtemp, dtime
      real*8  celent, pnewdt
!---------------------------------------------------------------------72
!---- Dimension arrays passed into the UMAT sub
!---------------------------------------------------------------------72
!
        real*8         &
!       dimension  &
        coords(3),     & ! Coordinates of Gauss pt. being evaluated
        ddsdde(ntens,ntens), & ! Tangent Stiffness Matrix
        ddsddt(ntens), & ! Change in stress per change in temperature
        dfgrd0(3,3),  & ! Deformation gradient at t_n
        dfgrd1(3,3),  & ! Deformation gradient at t_(n+1)
        dpred(1),     & ! Change in predefined state variables
        drplde(ntens), &! Change in heat generation per change in strain
        drot(3,3),     &! Incremental rotation matrix
        dstrain(ntens),&! Strain increment tensor (vector form)
        predef(1),     &! Predefined state vars dependent on field variables
        props(nprops), &! Material properties
        statev(nstatv),&! State variables
        strain(ntens), &! Strain tensor (vector form)
        stress(ntens), &! Cauchy stress (vector form)
        time(2)         ! Time Step and Total Time
!
!---------------------------------------------------------------------72
!---- Contents of the props array (nprops = 21)
!---------------------------------------------------------------------72
!
!------- crystal type and # of grain per ip
!          props(1)  crystalID : 1-FCC, 2-BCC, 3HCP
!          props(2)  numgrn
!
!------- crystal elasticity
!          props(3)  elastID : 1-ISO, 2-ANISO
!                    ISO      ANISO
!                          CUBIC  HCP
!          props(4)  eMod   C11   C11
!          props(5)  eMu    C12   C12
!          props(6)         C44   C13
!          props(7)               C33
!          props(8)               C44
!
!------- slip system kinetic equation: power law
!          props( 9)  xm     : strain-rate sensitivity (m)
!          props(10)  gam0   : reference shear rate (d_0)
!
!------- slip system hardening law : Voce's type
!          props(11)  h0     : initial work hardening rate (Theta_0)
!          props(12)  tausi  : slip system strength (tau_0)
!          props(13)  taus0  : saturation threshold strength (tau_s0)
!          props(14)  xms    : strain-rate sensitivity - (m')
!          props(15)  gamss0 : ref deformation rate - g_s^star
!          props(16)  kappa0 : initial slip system strength / slip
!                                  taus0 .ge. kappa0
!------- Lattice orientation codes
!          props(17)  kODF
!          props(18)  kODFOut
!
!------- iterations/tolerance data for state and Newton's method
!          props(19)  maxIterState
!          props(20)  tolerState
!          props(21)  maxIterNewt
!          props(22)  tolerNewt
!
      integer numel_aba, numqpt_aba, grainid
      real*8   gcrot0(3,3,1,1,1) !Rotation tensor C0
      type (xtalVars) crystalvar
      type (xtalVars_n) crystalvar_n
!
!---------------------------------------------------------------------72
!---- Initialize crystal parameters
!---------------------------------------------------------------------72
!
      if (mod(kinc,10)==1 .and. noel==1 .and. npt==1) write(*,'(A5, 2x, i8, 5x, A5, 2x, i4)') 'kinc=', kinc, 'kstep=', kstep

!
!-----Obtain grain ID
!
!     if ((index(cmname, 'ELSET') + index(cmname, 'elset')) <1  .or. (index(cmname, '_DAMAGE')+index(cmname, '_damage') ) <1 )  &
!         write (*,*) 'Make sure your material name follows the Elset_grainid_DAMAGE'
!     %get the grain id from the material name
!     read(cmname(6:(index(cmname, '_')-1)), *) grainid
     
     if ((index(cmname, 'Material-') + index(cmname, 'MATERIAL-')) <1   )  &
         write (*,*) 'Make sure your material name follows the Material-grainid'
!     %get the grain id from the material name
	     read(cmname((index(cmname, '-')+1):len(cmname)), *) grainid

      numel_aba  = nint (props(1)) ! number of elements in Abaqus model
      numqpt_aba = nint (props(2)) ! number of gauss points
!
!     Setupcrystalprop is called within uexternaldb.f

      !write(*,*) 'Check 2'
      if ( (time(2) .eq. 0.0) .and. (kinc .gt. 0) .and. (kstep==1)) then
!
!------- initialize arrays
!
	  	!write(*,*) "Check in CIA"
	  	!write(*,*) "no. of elements:", numel_aba
	  	!write(*,*) "no. of GPs:", numqpt_aba
	  	  
      	call CrystalInitializeArrays(numqpt, numel,  &
                                   numgrn, numslip,  &
                                   slipr0,  &
                                   backst0,  &
                                   backst0c,  &
                                   gcrot0,  &
                                   crystalvar%gstress,  &
                                   crystalvar_n%gstress_n,  &
                                   crystalvar%gestran,  &
                                   crystalvar_n%gestran_n,  &
                                   crystalvar%gslipr,  &
                                   crystalvar_n%gslipr_n,  &
                                   crystalvar%gbackst,  &
                                   crystalvar_n%gbackst_n,  &
                                   crystalvar%gbackstc,  &
                                   crystalvar_n%gbackstc_n,  &
                                   crystalvar%gstatev,  &
                                   crystalvar_n%gstatev_n,  &
                                   crystalvar%geqvalues,  &
                                   crystalvar%ggamdot,  &
                                   crystalvar%ggamdot2,  &
                                   crystalvar%gcrot,  &
                                   crystalvar_n%gcrot_n,  &
                                   crystalvar%grrot,  &
                                   crystalvar_n%grrot_n)
        !write(*,*) "Check ***"
        !write(*,*) "grainid", grainid

      	call SetUpInitialState(grainid, gcrot0, numqpt, numel, noel, npt,  &
                                nstatv, statev)

         !write(*,*) 'SetUpInitialState is called!'
      endif
!

!---------------------------------------------------------------------72
!---- Evolve material state - compute stresses and algorithmic moduli
!---------------------------------------------------------------------72
!
      !write(*,*) 'Check 3: nstatv=', nstatv
      if ( (time(2) .ge. 0d0) .and. (kinc .ge. 1) ) then
!
            call EvolveMatlState(stress, statev, ddsdde, strain,  &
               dstrain, time, dtime, temp, dtemp, & ! ndi, nshr, ntens, 
               nstatv, drot, dfgrd0, dfgrd1, noel, npt,  &
               kinc, kstep, numel, numqpt, iprint, pnewdt, numel_aba, numqpt_aba, gcrot0, grainid)
!
      endif
      !write(*,*) 'Check 4'

!

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SetUpCrystalProps(  )
!
      use     SlipOverStress
      use     WriteControl
!      use     InitialLocalRotMatrix
      use     IterData
      use     InitialState
      use     SlipData
      use     InitialSlipSys
      use     OriData
      use     EProp
      use     XtalPar
      use     WorkDir
      use     DummyVar
      implicit none
      include 'params_xtal.inc'

      integer mprops
      !real*8  props(mprops)
      real*8  props(2)
!
!
      real*8  gstress    (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstress_n  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran    (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gslipr     (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gslipr_n   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackst    (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackst_n  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackstc   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackstc_n (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev    (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n  (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues  (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot2   (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot      (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot      (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
!
!---------------------------------------------------------------------72
!
!-----------This is a dummy variable that is not necessary
!-----------leaving it here just becauase i do not want to change the
!           afterward routines to eliminated
           mprops=2
           props=1.d0
!--------------------------------------------------------------------
!---------------------------------------------------------------------72
!
!------- open input/output files
!

      call CrystalOpenIOFiles(filePath, fileRoot,  &
                              XTAL_I, XTAL_E, XTAL_O)
      call CrystalOpenIOFiles_2(filePath, fileRoot,  &
                                XTAL_STRESS_O, XTAL_STRAIN_O,  &
                                XTAL_EFFSS_O, XTAL_TRUESS_O,  &
                                XTAL_ITER_O, AGG_EFFSS_O)
!
!------- read/echo material data
!		 modified for dislocation climb: add kBar0
!
      call CrystalModelData(props, mprops, numqpt, numel,  &
                            numgrn, numslip, numvtx, kODFout,  &
                            !kappa0, fCeDevVol, fCeVol,  &
                            slipr0, backst0, backst0c, fCeDevVol, fCeVol,  &
                            matProp, tauSlip,  &
                            fCeDev, fCeiDev,  &
                            zBar0,  &
                            pBar0, qBar0,  &
                            pBar0Vec, qBar0Vec,  &
                            ppTBar0, KBar0,  &
                            KpBar0, KqBar0,  &
                            KpBar0Vec, KpBar0devVec, KqBar0Vec,  &
                            KpKpTBar0,  &
!                            gcrot0,  &
                            sigfs,  &
                            overstress,  &
                            hardmtx,  &
                            maxIterState, maxIterNewt,  &
                            tolerState, tolerNewt,  &
                            XTAL_I, XTAL_E, XTAL_O, filePath,  &
                            kODF, numor, seed, angles,  &
                            XTAL_TXT_OUT)
!
!------- initialize useful matrices to compute algorithmic moduli Cep
!
      call CrystalInitializeMatrxCep( )
!
!------- parameters for global iterations
!
      call GlobalIterationParams( )

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SetUpInitialState( grainid, gcrot0, &
         numqpt, numel, noel, npt, nstatv, statev &
         )
!
!      use     InitialLocalRotMatrix
      use     InitialState
      use     SlipData
      use     OriData
      implicit none
      include  'params_xtal.inc'
!
!      character*30 cmname
      integer numel, numqpt, noel, npt, grainid
      integer nstatv
      real*8  statev(nstatv)
      real*8  gcrot0(DIMS, DIMS, 1, NUMQPT_T, NUMEL_T)

!

!
!---------------------------------------------------------------------72
!
      call AssignCrystalODF(grainid, numel, numqpt, noel, npt, &
                            numgrn, numslip, kODF, &
                            numor, seed, &
                            gcrot0, &
                            angles, &
!                            euler,  &
                            XTAL_O, XTAL_TXT_OUT)

      call SetUpStateVars(nstatv, statev,  &
                          numgrn, numslip,  &
                          slipr0, backst0, backst0c, &
                          !kappa0,  &
                          gcrot0,  &
                          XTAL_O)
      !write(*,*) "Called!"                          

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SetUpStateVars(  &
         nstatv, statev,  &
         numgrn, numslip,  &
         slipr0, backst0, backst0c, &
         !kappa0,  &
         gcrot0,  &
         FILE_O  &
         ) 

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer nstatv
      real*8  statev(nstatv)

      integer numgrn, numslip, FILE_O
      real*8  slipr0(NKAPP), backst0(NKAPP), backst0c(NKAPP)
      !real*8  kappa0(NKAPP)
      real*8  gcrot0 (DIMS, DIMS, 1, NUMQPT_T, NUMEL_T)

      integer dex, ig, varsPerGrn
!
!---------------------------------------------------------------------72
!
!---- initialize state variable vector per ip
!
      varsPerGrn = 2*NVECS         & ! stress, estran 2*5=10
!     &             + NKAPP          ! kappa 
                   !+ numslip         & ! kappa
                   + 3*numslip       & !slipr, backst, backstc 3*12=36
                   + NSTAV           & ! ee_v, p, detVe, WpNorm, ssact 5
                   + NEQVA/2         & ! eqps, eqstr, gam_star, gamtot 8/2=4
                   + 3*DIMS*DIMS     & ! C_0, C, R 3*3*3=27
                   + 2*numslip        ! gamdot 2*12=24
	  !write(*,*) 'setup: number of vars per grain', varsPerGrn	 
	  
      
      if (nstatv .ne. numgrn*varsPerGrn)  &
         call RunTimeError(FILE_O, 'nstatv .ne. numgrn*varsPerGrn')
      
      !write(*,*) "nstatv, numgrn, varspergrn=", nstatv, numgrn, varsPerGrn
      do ig = 1, numgrn

		 dex = (ig-1)*varsPerGrn + 1           ! stress (xtal)
         call SetTensor(statev(dex), pzero, NVECS)
         dex = dex + NVECS                      ! estrain

         call SetTensor(statev(dex), pzero, NVECS)
 		 dex = dex + NVECS                      ! kappa
 		 
!         call SetTensor(statev(dex), kappa0, NKAPP)
         !call EqualTensors(kappa0, statev(dex), numslip)
         call EqualTensors(slipr0, statev(dex), numslip)
         dex = dex + numslip
         
         call EqualTensors(backst0, statev(dex), numslip)
!        dex = dex + NKAPP              ! ee_v, p, detVe, WpNorm, ssact
         dex = dex + numslip

         call EqualTensors(backst0c, statev(dex), numslip)
         dex = dex + numslip         
         
         call SetTensor(statev(dex), pzero, NSTAV)
         dex = dex + NSTAV               ! eqps, eqstr, gam_star, gamtot

         call SetTensor(statev(dex), pzero, NEQVA/2)
         dex = dex + NEQVA/2                    ! C_0

         call EqualTensors(gcrot0(1,1,ig,1,1), statev(dex), DIMS*DIMS)
         dex = dex + DIMS*DIMS                  ! C

         call EqualTensors(gcrot0(1,1,ig,1,1), statev(dex), DIMS*DIMS)
         dex = dex + DIMS*DIMS                  ! R

         call EqualTensors(Ident2nd, statev(dex), DIMS*DIMS)
         dex = dex + DIMS*DIMS                  ! gamdot

         call SetTensor(statev(dex), pzero, numslip)
         dex = dex + numslip                    ! gamdotC

         call SetTensor(statev(dex), pzero, numslip)
         
         !write(*,*) "dex and nstat", dex, nstatv

      enddo
	  
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE EvolveMatlState(  &
         stress, statev, ddsdde, strain, dstrain, time,  &
         dtime, temp, dtemp,  & ! ndi, nshr, ntens, 
         nstatv, drot, dfgrd0, dfgrd1, noel, npt,  &
         incr, kstep, numel, numqpt, iprint, pnewdt, numel_aba, numqpt_aba, gcrot0, grainid  &
         )
!
!      use     InitialLocalRotMatrix
      use     SlipData
      use     DataType
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer nstatv , numel_aba, numqpt_aba, grainid
      integer noel, npt, incr, numel, numqpt, iprint, kstep
      real*8  dtime, temp, dtemp, pnewdt
      real*8  stress(ntens), statev(nstatv), ddsdde(ntens, ntens)
      real*8  strain(ntens), dstrain(ntens), time(2)
      real*8  dfgrd0(3,3), dfgrd1(3,3), drot(3,3)

      integer statusFlag, numIncrs
      real*8  gstress   (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran   (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      !real*8  gkappa    (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gslipr    (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  gbackst   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackstc  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  gstatev   (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot     (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
!
      real*8  gcrot0(DIMS, DIMS, 1, NUMQPT_T, NUMEL_T)
!
      real*8  grrot     (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot   (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot2  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  gstress_n (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      !real*8  gkappa_n  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gslipr_n  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  gbackst_n  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackstc_n  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      
      real*8  gstatev_n (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

      real*8  vgrad(DIMS, DIMS)
      real*8  savg_ij(DIMS, DIMS, NUMQPT_T) 
      real*8  cavg_ijkl(DIMS2, DIMS2, NUMQPT_T)
      real*8  epsxx(NUMQPT_T), epsyy(NUMQPT_T), epsth(NUMQPT_T)
      real*8  epsxy(NUMQPT_T), epsum(NUMQPT_T), spinn(NUMQPT_T)
      real*8  qptpo(NUMQPT_T), qptp(NUMQPT_T)
      real*8  sigxx(NUMQPT_T), sigyy(NUMQPT_T)
      real*8  sigth(NUMQPT_T), sigxy(NUMQPT_T)

      real*8  epsxz(NUMQPT_T),  epsyz(NUMQPT_T)
      real*8  spinxz(NUMQPT_T), spinyz(NUMQPT_T)
      real*8  sigxz(NUMQPT_T),  sigyz(NUMQPT_T)
      integer ielem
      type (xtalVars) crystalvar
      type (xtalVars_n) crystalvar_n
!
!---------------------------------------------------------------------72
!
!-------- kinematic quantities - deformation path
!
     !write(*,*) '1noel:', noel, '1grainid', grainid
      !write(*,*) 'check a'
      call VelocityGrad(dstrain, drot, dfgrd0, dfgrd1, vgrad, dtime,  &
                        noel, npt, incr)
     !write(*,*) '2noel:', noel, '2grainid', grainid
      
      call DeformPath(epsxx, epsyy, epsth, epsxy, epsum, spinn,  epsxz, epsyz, spinxz, spinyz, sigxz, sigyz, &
            qptpo, qptp, vgrad, temp, (temp+dtemp), numel, numqpt)
     !write(*,*) '3noel:', noel, '3grainid', grainid
!
!-------- fetch state variables from abaqus vector array
!
      !write(*,*) 'check b! need to solve RecoverStateVars'
      !write(*,*) "crystalvar_n%gstress_n", crystalvar_n%gstress_n
      !write(*,*) "crystalvar_n%gestran_n", crystalvar_n%gestran_n
      !write(*,*) "crystalvar_n%gslipr_n", crystalvar_n%gslipr_n
      !write(*,*) "crystalvar_n%gbackst_n", crystalvar_n%gbackst_n
      !write(*,*) "crystalvar_n%gstatev_n", crystalvar_n%gstatev_n
      !write(*,*) "crystalvar_n%gcrot_n", crystalvar_n%gcrot_n
      !write(*,*) "crystalvar_n%grrot_n", crystalvar_n%grrot_n
      !write(*,*) "gcrot0", gcrot0
      !write(*,*) "crystalvar%ggamdot", crystalvar%ggamdot
      !write(*,*) "crystalvar%ggamdot2", crystalvar%ggamdot2
      !write(*,*) "crystalvar%geqvalues", crystalvar%geqvalues
      !write(*,*) "Check in RecoverStateVars:"
      !write(*,*) "no. of statev:", nstatv
      !write(*,*) "no. of elements:", numel_aba
	  !write(*,*) "no. of grn:", numgrn
	  !write(*,*) "no. of GPs:", numqpt_aba
      call RecoverStateVars(statev, nstatv,  &
                            crystalvar_n%gstress_n, crystalvar_n%gestran_n,  &
                            crystalvar_n%gslipr_n, crystalvar_n%gbackst_n, crystalvar_n%gbackstc_n,  &
                            crystalvar_n%gstatev_n, crystalvar_n%gcrot_n, crystalvar_n%grrot_n, gcrot0,  & 
                            crystalvar%ggamdot,  &  ! glide rate
                            crystalvar%ggamdot2,  & ! climb rate
                            crystalvar%geqvalues,  &
                            numgrn, numslip)
!
     !write(*,*) '4noel:', noel, '4grainid', grainid
!
!-------- compute state
!
      !write(*,*) 'check c'
      numIncrs = 1000
      ielem    = 1
      statusFlag = XTAL_CONVERGED
      call StressCrystal(sigxx, sigyy, sigth, sigxy,  &
                         epsxx, epsyy, epsth, epsxy, epsum,  &
                         spinn, epsxz, epsyz, spinxz, spinyz, sigxz, sigyz,  &
                         qptpo, qptp, dtime, time,  &
                         ielem, incr, kstep, numel, numqpt, iprint, savg_ij, cavg_ijkl,  &
                         statusFlag, numIncrs, noel, npt, numel_aba, numqpt_aba,  &
                         crystalvar, crystalvar_n, gcrot0, grainid)
!                         
      if (statusFlag .ne. XTAL_CONVERGED) then
         write(XTAL_O, *) ' ** Umat did not converged       **'
         write(XTAL_O, *) ' ** re-scaling time step by 0.75 **'
         pnewdt = 0.75
         return
      endif
!
     !write(*,*) '5noel:', noel, '5grainid', grainid
!-------- save state variables in abaqus vector array
!
      !write(*,*) 'check d'
      !call SaveStateVars(statev, nstatv, gstress, gestran, gkappa,  &
      !write(*,*) "nstat=", nstatv      
      call SaveStateVars(statev, nstatv, crystalvar%gstress,   &
            crystalvar%gestran, crystalvar%gslipr, crystalvar%gbackst, crystalvar%gbackstc,  &
            crystalvar%gstatev, crystalvar%gcrot, crystalvar%grrot, gcrot0, crystalvar%ggamdot,  &
            crystalvar%ggamdot2,  &
            crystalvar%geqvalues, numgrn,  &
            numslip)
      !write(*,*) 'check e'
     !write(*,*) '60noel:', noel, '60grainid', grainid
!
!-------- stresses and algorithmic moduli in abaqus format
!
      call SaveStressModuli(stress, ddsdde, savg_ij(1, 1, 1),  &
                            cavg_ijkl(1, 1, 1))
     !write(*,*) '70noel:', noel, '70grainid', grainid

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE VelocityGrad(  &
         dstran, drot_aba, dfgrd0, dfgrd1, vgrad, dtime, ielem, iqpt,  &
         incr  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

!      integer ntens, nshr
      integer ielem, iqpt, incr
      real*8  dtime
      real*8  dstran(NTENS)
      real*8  drot_aba(3,3), dfgrd0(3,3), dfgrd1(3,3), vgrad(3,3)

      integer i, j
      real*8  det
      real*8  d_aba(3,3)
      real*8  d_umat(3,3), w_umat(3,3), drot_umat(3,3), spin(3)
      real*8  matx_1(3,3), matx_2(3,3), matx_3(3,3), rdfgrd(3,3)
!
!---------------------------------------------------------------------72
!
!-------- relative deformation gradient: f=F_(n+1)*F_n^{-1}
!
      call InverseOfMat3x3(dfgrd0, matx_1, det)
      call MultAxB(dfgrd1, matx_1, rdfgrd, 3)
!
!-------- approximation to velocity gradient: L=2/dt*(f-I)*(f+I)^{-1}
!
      call AddTensors(pone, rdfgrd, -pone, Ident2nd, matx_2, 3*3)
      call AddTensors(pone, rdfgrd, pone, Ident2nd, matx_1, 3*3)
      call InverseOfMat3x3(matx_1, matx_3, det)

      call MultAxB(matx_2, matx_3, vgrad, 3)
      call SetToScaledtensor(2.0/dtime, vgrad, vgrad, 3*3)
!
!-------- rate of deformation and spin: d_umat=sym(L), w_umat=skw(L)
!
      call SetTensor(d_umat, pzero, DIMS*DIMS)
      call SetTensor(w_umat, pzero, DIMS*DIMS)

      call SymmetrizeTensor(vgrad, d_umat, 3)
      call SkewSymmetrizeTensor(vgrad, w_umat, 3)
!
!-------- incremental rotation tensor: drot_umat = exp(dt*w_umat)
!
      call Mat3x3ToVec3x1Skew(w_umat, spin, 3)
      call IncrementalRotTensor(drot_umat, spin, dtime)
!
!-------- rate of deformation from ABAQUS
!
      call SetTensor(d_aba, pzero, DIMS*DIMS)
      d_aba(1,1) = dstran(1)/dtime
      d_aba(2,2) = dstran(2)/dtime
      d_aba(3,3) = dstran(3)/dtime
      d_aba(1,2) = dstran(4)/ptwo/dtime
      d_aba(2,1) = d_aba(1,2)

      if (NSHR .gt. 1) then
         d_aba(1,3) = dstran(5)/ptwo/dtime
         d_aba(2,3) = dstran(6)/ptwo/dtime
         d_aba(3,1) = d_aba(1,3)
         d_aba(3,2) = d_aba(2,3)
      endif
!
!-------- print values for comparison
!
!      if (ielem .eq. kPRINT_ELEM .and. iqpt .eq. kPRINT_QPT) then
      if (ielem .eq. kPRINT_ELEM .and. iqpt .eq. 100) then
         if (incr .eq. 1) write(XTAL_E, 1000) ielem, iqpt
         write(XTAL_E, '(/i5)') incr
         do i = 1, 3
            write(XTAL_E, 5000) (vgrad(i,j), j=1,3)
         enddo
         write(XTAL_E, *)
         do i = 1, 3
            write(XTAL_E, 5000) (d_aba(i,j), j=1,3),  &
                                (d_umat(i,j),j=1,3),  &
                                (w_umat(i,j),j=1,3)
         enddo
         write(XTAL_E, *)
         do i = 1, 3
            write(XTAL_E, 5000) (drot_aba(i,j), j=1,3),  &
                                (drot_umat(i,j),j=1,3)
         enddo
      endif

1000  format(' INCR',/,18x,' L_ij',/,  &
             18x,' D_aba',32x,'D_umat',30x,'W_umat',/,  &
             18x,' Drot_aba',32x,'Drot_umat',/,  &
             38x,' (elem # ', i5, ',  qpt # ', i2, ')')
5000  format((3(2x,3(1x,e11.4))))
!
!-------- recompute velocity gradient: L = d_aba + w_umat
!
      call AddTensors(pone, d_aba, pone, w_umat, vgrad, DIMS*DIMS)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE DeformPath(  &
         epsxx, epsyy, epsth, epsxy, epsum, spinn,  &
         epsxz, epsyz, spinxz, spinyz,  &
         sigxz, sigyz, qptpo, qptp, vgrad,  &
         theta_n, theta, numel, numqpt  &
         )

      implicit none
      include 'params_xtal.inc'

      integer numel, numqpt
      real*8  theta_n, theta
      real*8  vgrad(DIMS, DIMS)
      real*8  epsxx(numqpt), epsyy(numqpt), epsth(numqpt)
      real*8  epsxy(numqpt), epsum(numqpt), spinn(numqpt)
      real*8  qptpo(numqpt), qptp(numqpt)

      real*8  epsxz(NUMQPT_T),  epsyz(NUMQPT_T)
      real*8  spinxz(NUMQPT_T), spinyz(NUMQPT_T)
      real*8  sigxz(NUMQPT_T),  sigyz(NUMQPT_T)

!
!---------------------------------------------------------------------72
!
!------- deformation rate: deviatoric and volumetric 
!
      epsxx(1) = vgrad(1, 1)
      epsyy(1) = vgrad(2, 2)
      epsth(1) = vgrad(3, 3)

      epsum(1) = epsxx(1) + epsyy(1) + epsth(1)

      epsxx(1) = epsxx(1) - epsum(1) / 3.0d0
      epsyy(1) = epsyy(1) - epsum(1) / 3.0d0
      epsth(1) = epsth(1) - epsum(1) / 3.0d0

      epsxy(1) = (vgrad(2, 1) + vgrad(1, 2))/2.d0
      epsxz(1) = (vgrad(3, 1) + vgrad(1, 3))/2.d0
      epsyz(1) = (vgrad(3, 2) + vgrad(2, 3))/2.d0
!
!------- axial vector of spin
!
      spinn(1)  = (vgrad(1, 2) - vgrad(2, 1))/2.d0
      spinxz(1) = (vgrad(1, 3) - vgrad(3, 1))/2.d0
      spinyz(1) = (vgrad(2, 3) - vgrad(3, 2))/2.d0
!
!------- temperature (K)
!
      qptpo(1) = theta_n
      qptp(1)  = theta

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE RecoverStateVars(  &
         !statev, nstatv, gstress_n, gestran_n, gkappa_n, gstatev_n,  &
         statev, nstatv,  &
         gstress_n, gestran_n, gslipr_n,  &
         gbackst_n, gbackstc_n, gstatev_n,  &
         gcrot_n, grrot_n, gcrot0,  &
         ggamdot,  &
         ggamdot2,  &
         geqvalues,  &
         numgrn, numslip)

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
      
      integer nstatv, numgrn, numslip
      real*8  statev(nstatv)

      real*8  gstress_n (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      !real*8  gkappa_n  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gslipr_n  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackst_n (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackstc_n(NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)      
      real*8  gstatev_n (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot0    (DIMS, DIMS, 1, NUMQPT_T, NUMEL_T)

      real*8  geqvalues (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot   (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot2  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)

      integer varsPerGrn, ig, dex
!
!---------------------------------------------------------------------72
!
!---- number of state variables per ip, inputfile: ... *Depvar 106
!
      varsPerGrn = 2*NVECS        &  ! stress, estran 2*5=10
!     &             + NKAPP          ! kappa
                   !+ numslip         & ! kappa
                   + 3*numslip        & ! slipr, backst and backstc 3*12=36
                   + NSTAV           & ! ee_v, p, detVe, WpNorm, ssact 5
                   + NEQVA/2         & ! eqps, eqstr, gam_star, gamtot 8/2=4
                   + 3*DIMS*DIMS     & ! C_0, C, R 3*3*3=27
                   + 2*numslip         ! gam_dot and gam_dotC 2*12
                   
	  !write(*,*) "RecoverStatev: number of vars per grain=", varsPerGrn
!
!---- recover state variables from abaqus vector array
!
      do ig = 1, numgrn

         dex = (ig-1)*varsPerGrn + 1           ! stress

         call EqualTensors(statev(dex), gstress_n(1,ig,1,1), NVECS)
         dex = dex + NVECS                      ! estrain
         
         call EqualTensors(statev(dex), gestran_n(1,ig,1,1), NVECS)
         dex = dex + NVECS                      ! kappa
         
!         call EqualTensors(statev(dex), gkappa_n(1,ig,1,1), NKAPP)
         !call EqualTensors(statev(dex), gkappa_n(1,ig,1,1), numslip)
         
         call EqualTensors(statev(dex), gslipr_n(1,ig,1,1), numslip)
         dex = dex + numslip             ! ee_v, p, detVe, WpNorm, ssact
         
         call EqualTensors(statev(dex), gbackst_n(1,ig,1,1), numslip)
!         dex = dex + NKAPP              ! ee_v, p, detVe, WpNorm, ssact
         dex = dex + numslip             

         call EqualTensors(statev(dex), gbackstc_n(1,ig,1,1), numslip)
         dex = dex + numslip             

         call EqualTensors(statev(dex), gstatev_n(1,ig,1,1), NSTAV)
         dex = dex + NSTAV               ! eqps, eqstr, gam_star, gamtot

         call EqualTensors(statev(dex), geqvalues(1,ig,1,1), NEQVA/2)
         dex = dex + NEQVA/2                    ! C_0

         call EqualTensors(statev(dex), gcrot0(1,1,ig,1,1), DIMS*DIMS)
         dex = dex + DIMS*DIMS                  ! C

         call EqualTensors(statev(dex), gcrot_n(1,1,ig,1,1), DIMS*DIMS)
         dex = dex + DIMS*DIMS                  ! R

         call EqualTensors(statev(dex), grrot_n(1,1,ig,1,1), DIMS*DIMS)
         dex = dex + DIMS*DIMS                  ! gamdot

         call EqualTensors(statev(dex), ggamdot(1,ig,1,1), numslip)
         dex = dex + numslip                  ! gamdotC

         call EqualTensors(statev(dex), ggamdot2(1,ig,1,1), numslip)
         
      enddo

      
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SaveStateVars(  &
         !statev, nstatv, gstress, gestran, gkappa, gstatev, gcrot,  &
         statev, nstatv, gstress, gestran, gslipr, gbackst, gbackstc,  &
         gstatev, gcrot,  &
         grrot, gcrot0, ggamdot,  &
         ggamdot2,  &
         geqvalues, numgrn, numslip  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
      
      integer nstatv, numgrn, numslip
      real*8  statev(nstatv)

      real*8  gstress   (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran   (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gslipr    (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackst   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackstc  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)      
      !real*8  gkappa    (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev   (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot     (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot     (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot0    (DIMS, DIMS, 1, NUMQPT_T, NUMEL_T)

      real*8  geqvalues (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot   (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot2  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)

      integer varsPerGrn, ig, dex
!
!---------------------------------------------------------------------72
!
!---- number of state variables per ip
!
      varsPerGrn = 2*NVECS        &  ! stress, estran
!     &             + NKAPP          ! kappa
                   !+ numslip         & ! kappa
                   + 3*numslip        & ! slipr, backst and backstc
                   + NSTAV           & ! ee_v, p, detVe, WpNorm, ssact
                   + NEQVA/2         & ! eqps, eqstr, gam_star, gamtot 
                   + 3*DIMS*DIMS     & ! C_0, C, R
                   + 2*numslip        ! gam_dot + gam_dotC
                   
	  !write(*,*) "Save: no. of vars=", varsPerGrn
!
!---- save state variables in abaqus vector array
!
      do ig = 1, numgrn

         dex = (ig-1)*varsPerGrn + 1           ! stress

         call EqualTensors(gstress(1,ig,1,1), statev(dex), NVECS)
         dex = dex + NVECS                      ! estrain

         call EqualTensors(gestran(1,ig,1,1), statev(dex), NVECS)
         dex = dex + NVECS                      ! kappa

!        call EqualTensors(gkappa(1,ig,1,1), statev(dex), NKAPP)
!		 call EqualTensors(gkappa(1,ig,1,1), statev(dex), numslip)
         
         call EqualTensors(gslipr(1,ig,1,1), statev(dex), numslip)
         dex = dex + numslip
         
         call EqualTensors(gbackst(1,ig,1,1), statev(dex), numslip)
!         dex = dex + NKAPP              ! ee_v, p, detVe, WpNorm, ssact
         dex = dex + numslip

         call EqualTensors(gbackstc(1,ig,1,1), statev(dex), numslip)
         dex = dex + numslip

         call EqualTensors(gstatev(1,ig,1,1), statev(dex), NSTAV)
         dex = dex + NSTAV               ! eqps, eqstr, gam_star, gamtot
         
         call EqualTensors(geqvalues((NEQVA/2+1),ig,1,1), statev(dex),  &
                                                              NEQVA/2)
         dex = dex + NEQVA/2                    ! C_0
         !call EqualTensors(gcrot0(1,1,ig,1,1), statev(dex), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! C
         call EqualTensors(gcrot(1,1,ig,1,1), statev(dex), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! R
         call EqualTensors(grrot(1,1,ig,1,1), statev(dex), DIMS*DIMS)

         dex = dex + DIMS*DIMS                  ! gamdot
         call EqualTensors(ggamdot(1,ig,1,1), statev(dex), numslip)
         
         dex = dex + numslip                  ! gamdotC
         call EqualTensors(ggamdot2(1,ig,1,1), statev(dex), numslip)

      enddo
	  !write(*,*) "dex and nstat", dex, nstatv
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE DeformationRate_n(  &
         dvec_n, wvec_n, d_kk_n, iqpt  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer iqpt
      real*8  d_kk_n
      real*8  dvec_n(NVECS), wvec_n(3)
!
!---------------------------------------------------------------------72
!
!------ zeros arrays
!
      d_kk_n = 0.0
      call SetTensor(dvec_n, pzero, NVECS)
      call SetTensor(wvec_n, pzero, 3)
!
!------ stop program if called during MPS runs ??
!
!      print *, 'DeformationRate_n called ... Exiting'
!      call RunTimeError(XTAL_O, 'DeformationRate_n called')

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SaveStressModuli(  &
         stress, ddsdde, savg_ij, cavg_ijkl  &
         )

      implicit none
      include 'params_xtal.inc'

      real*8  stress(NTENS), ddsdde(NTENS, NTENS)
      real*8  savg_ij(DIMS, DIMS), cavg_ijkl(DIMS2, DIMS2)

      integer i, j
!
!---------------------------------------------------------------------72
!
!---- stresses
!
      stress(1) = savg_ij(1,1)
      stress(2) = savg_ij(2,2)
      stress(3) = savg_ij(3,3)
      stress(4) = savg_ij(1,2)

      if (NSHR .gt. 1) then
         stress(5) = savg_ij(1,3)
         stress(6) = savg_ij(2,3)
      endif
!
!---- algorithmic moduli
!
      do i = 1, NTENS
         do j = 1, NTENS
            ddsdde(i,j) = cavg_ijkl(i,j)
         enddo
      enddo

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE CrystalOpenIOFiles(  &
         filePath, fileRoot,  &
         FILE_I, FILE_E, FILE_O  &
         )

      implicit none
    
      character*80 filePath, fileRoot
      integer      FILE_I, FILE_E, FILE_O

      integer      length1, length2
      character*80 dataFile, filename
!
!---------------------------------------------------------------------72
!
!------- root name of input/output files
!
!      write(*,*) ' Root name of crystal-model input/output files ?'
!      read 1000, dataFile
!
!------- open files
!
      length1 = index(filePath,' ') - 1
      length2 = index(fileRoot,' ') - 1
	  
!filename = filePath(1:length1)//fileRoot(1:length2)//'.xtali'
      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtali'
      open(unit=FILE_I, file=filename, status='unknown',  &
                                                 access='sequential')
      rewind(FILE_I)
	  
!      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtale'
      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtale'
      open(unit=FILE_E, file=filename, status='unknown')
      rewind(FILE_E)

!       filename = filePath(1:length1)//fileRoot(1:length2)//'.xtalo'	  
      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtalo'
      open(unit=FILE_O, file=filename, status='unknown')

!
!------- formats
!
1000  format(a80)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE CrystalOpenIOFiles_2(  &
         filePath, fileRoot,  &
         FILE_STRESS_O, FILE_STRAIN_O,  &
         FILE_EFFSS_O, FILE_TRUESS_O,  &
         FILE_ITER_O, FILE_AGG_EFFSS_O  &
         )

      implicit none

      character*80 filePath, fileRoot
      integer FILE_STRESS_O, FILE_STRAIN_O,  &
              FILE_EFFSS_O, FILE_TRUESS_O,  &
              FILE_ITER_O, FILE_AGG_EFFSS_O
    
      integer      length1, length2
      character*80 dataFile, filename
!
!---------------------------------------------------------------------72
!
!------- root name of i/o files
!
!      write(*,*) ' Root name of additional output files ?'
!      read 1000, dataFile
!
!------- open files for output at specified elem/iqpt/grain
!
      length1 = index(filePath,' ') - 1
      length2 = index(fileRoot,' ') - 1

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtal.strs'
      open(unit=FILE_STRESS_O, file=filename, status='unknown')

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtal.strn'
      open(unit=FILE_STRAIN_O, file=filename, status='unknown')

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtal.efss'
      open(unit=FILE_EFFSS_O, file=filename, status='unknown')

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtal.trss'
      open(unit=FILE_TRUESS_O, file=filename, status='unknown')

      filename = filePath(1:length1)//fileRoot(1:length2)//'.xtal.iter'
      open(unit=FILE_ITER_O, file=filename, status='unknown')
!
!------- open files for aggregate output at specified elem/iqpt
!
      filename = filePath(1:length1)//fileRoot(1:length2)//'.agg.efss'
      open(unit=FILE_AGG_EFFSS_O, file=filename, status='unknown')
!
!------- formats
!
1000  format(a80)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE CrystalModelData(  &
         props, mprops, numqpt, numel,  &
         numgrn, numslip, numvtx, kODFout,  &
         !kappa0, fCeDevVol, fCeVol,  &
         slipr0, backst0, backst0c, fCeDevVol, fCeVol,  &
         matProp, tauSlip,  &
         fCeDev, fCeiDev,  &
         zBar0, &
         pBar0, qBar0,  &
         pBar0Vec, qBar0Vec,  &
         ppTBar0, KBar0,  &
         KpBar0, KqBar0,  &
         KpBar0Vec, KpBar0devVec, KqBar0Vec,  &
         KpKpTBar0,  &
!         gcrot0,  &
         sigfs,  &
         overstress,  &
         hardmtx,  &
         maxIterState, maxIterNewt,  &
         tolerState, tolerNewt,  &
         FILE_I, FILE_E, FILE_O, filePath,  &
         kODF, numor, seed, angles,   &
         FILE_TXT_OUT  &
         )

      implicit none
      include 'params_xtal.inc'

      character*80 filePath
      integer mprops, numqpt, numel
      real*8  props(mprops)

      integer numgrn, numslip, numvtx, kODFout
      !real*8  kappa0(NKAPP), fCeDevVol(NVECS), fCeVol
      real*8  slipr0(NKAPP), backst0(NKAPP), backst0c(NKAPP) , fCeDevVol(NVECS), fCeVol
      real*8  matProp(NPROPS, MAX_SLIP), tauSlip(MAX_SLIP)
      real*8  fCeDev(NVECS, NVECS), fCeiDev(NVECS, NVECS)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), qBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0Vec(NVECS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  KBar0(DIMS, DIMS, MAX_SLIP)      
      real*8  KpBar0(DIMS, DIMS, MAX_SLIP), KqBar0(DIMS, DIMS, MAX_SLIP)
      real*8  KpBar0Vec(NVECS, MAX_SLIP), KpBar0devVec(NVECS, MAX_SLIP), KqBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      real*8  KpKpTBar0(NVECS, NVECS, MAX_SLIP)
!      real*8  gcrot0(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  sigfs(5, MAX_VTX)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      real*8  overstress(MAX_SLIP)

      integer maxIterState, maxIterNewt
      real*8  tolerState, tolerNewt

      integer FILE_I, FILE_E, FILE_O
      integer FILE_TXT_OUT

      integer kODF, numor, seed
      real*8  angles(DIMS, MAX_ORIEN)
!      real*8  euler(DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
!
!---------------------------------------------------------------------72
!
!------- input material data for single crystal
!
      call CrystalMaterialData(props, mprops, numqpt, numel,  &
                            numgrn, numslip, numvtx, kODFout,  &
                            !kappa0, fCeDevVol, fCeVol,  &
                            slipr0, backst0, backst0c, fCeDevVol, fCeVol,  &
                            matProp, tauSlip,  &
                            fCeDev, fCeiDev,  &
                            zBar0,  &
                            pBar0, qBar0,  &
                            pBar0Vec, qBar0Vec,  &
                            ppTBar0, KBar0,  &
                            KpBar0, KqBar0,  &
                            KpBar0Vec, KpBar0devVec, KqBar0Vec,  &
                            KpKpTBar0,  &
                            !gcrot0,  &
                            sigfs,  &
                            overstress,  &
                            hardmtx,  &
                            FILE_I, FILE_E, FILE_O, filePath,  &
                            kODF, numor, seed, angles,  &
                            FILE_TXT_OUT)
!
!------- input convergence data (state iterations & constitutive solver)
!
      call CrystalSolverData(props, mprops,  &
                            maxIterState, maxIterNewt,  &
                            tolerState, tolerNewt,  &
                            FILE_I, FILE_E)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE CrystalMaterialData(  &
         props, mprops, numqpt, numel,  &
         numgrn, numslip, numvtx, kODFout,  &
         !kappa0, fCeDevVol, fCeVol,  &
         slipr0, backst0, backst0c, fCeDevVol, fCeVol,  &
         matProp, tauSlip,  &
         fCeDev, fCeiDev,  &
         zBar0,  &
         pBar0, qBar0,  &
         pBar0Vec, qBar0Vec,  &
         ppTBar0, KBar0,  &
         KpBar0, KqBar0,  &
         KpBar0Vec, KpBar0devVec, KqBar0Vec,  &
         KpKpTBar0,  &
         !gcrot0,  &
         sigfs,  &
         overstress,  &
         hardmtx,  &
         FILE_I, FILE_E, FILE_O, filePath,  &
         kODF, numor, seed, angles,   &
         FILE_TXT_OUT  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      character*80 filePath
      integer mprops, numqpt, numel
      real*8  props(mprops)

      integer numgrn, numslip, numvtx, kODFout
      !real*8  kappa0(NKAPP), fCeDevVol(NVECS), fCeVol
      real*8  slipr0(NKAPP), backst0(NKAPP), backst0c(NKAPP), fCeDevVol(NVECS), fCeVol
      real*8  matProp(NPROPS, MAX_SLIP), tauSlip(MAX_SLIP)
      real*8  fCeDev(NVECS, NVECS), fCeiDev(NVECS, NVECS)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), qBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0Vec(NVECS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  KBar0(DIMS, DIMS, MAX_SLIP)      
      real*8  KpBar0(DIMS, DIMS, MAX_SLIP), KqBar0(DIMS, DIMS, MAX_SLIP)
      real*8  KpBar0Vec(NVECS, MAX_SLIP), KpBar0devVec(NVECS, MAX_SLIP), KqBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      real*8  KpKpTBar0(NVECS, NVECS, MAX_SLIP)      
!      real*8  gcrot0(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  sigfs(5, MAX_VTX)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      real*8  overstress(MAX_SLIP)

      integer FILE_I, FILE_E, FILE_O
      integer FILE_TXT_OUT

      integer kODF, numor, seed
      real*8  angles(DIMS, MAX_ORIEN)
!      real*8  euler(DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  xtalProp(NPROPS)

      integer is, i
      character*16 scysFile
      integer crystalID

      integer SXFILE_I
      parameter(SXFILE_I=91)

      integer    length1, length2
      character  filename*80, SingleXtalFile*80, prosa*80
!
!---------------------------------------------------------------------72
!
!---- initialize some arrays that store material data / slip system
!
      !call SetTensor(kappa0,  pzero, NKAPP)
      call SetTensor(slipr0,  pzero, NKAPP)
      call SetTensor(backst0,  pzero, NKAPP)
      call SetTensor(backst0c, pzero, NKAPP)
      call SetTensor(matProp, pzero, NPROPS*MAX_SLIP)
      call SetTensor(hardmtx, pzero, MAX_SLIP*MAX_SLIP)
  
      call SetTensor(xtalProp, pzero, NPROPS)
!
!---- crystal type and number of grains per IP
!
      call SetCrystalType(crystalID, numgrn, props, mprops,  &
                          FILE_I, FILE_E, FILE_O)
!
!---- open single crystal filename

      read(FILE_I,'(a18)') SingleXtalFile

      length1 = index(filePath,' ') - 1
      length2 = index(SingleXtalFile,' ') - 1

      filename = filePath(1:length1)//SingleXtalFile(1:length2)
      open(unit=SXFILE_I, file=filename, status='unknown',  &
                                         access='sequential')
      rewind(SXFILE_I)
!
!---- elastic stiffness matrix of single crystal
!
      call SetCrystalElasticity(crystalID, fCeDev, fCeiDev, xtalProp,  &
                                fCeDevVol, fCeVol, props, mprops,  &
                                SXFILE_I, FILE_E, FILE_O)
!
!---- thermal expansion coefficients
!
!      call SetCrystalThermal(xtalProp, SXFILE_I, FILE_E)
       read(SXFILE_I,'(a)') prosa
!
!---- set crystal geometry: slip / twinning data
!
      call SetCrystalGeometry(crystalID, numslip, numvtx, scysFile,  &
         zBar0, pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0, &
         KBar0, KpBar0, KqBar0, KpBar0Vec, KpBar0devVec, KqBar0Vec,  &
         KpKpTBar0, xtalProp,  &
         !matProp, hardmtx, kappa0, tauSlip, props, mprops, SXFILE_I,  &
         matProp, hardmtx, slipr0, backst0, backst0c, tauSlip, props, mprops, SXFILE_I,  &
         FILE_E)
!
!---- close single crystal filename
!
      close(SXFILE_I)
!
!---- crystal orientations
!
      call SetCrystalLatticeOrient(numgrn, numqpt, numel, kODFout,  &
                                   props, mprops,  &
                                   kODF, numor,seed, angles,   &
                                   FILE_I, FILE_E, FILE_O, filePath,  &
                                   FILE_TXT_OUT)
!
!---- vertices of rate independent yield surface (single crystal)
!
      if (crystalID .eq. kFCC .or. crystalID .eq. kBCC) then
         !call VerticesSCYS_Cubic(numvtx, kappa0(1), sigfs, scysFile,  &
         call VerticesSCYS_Cubic(numvtx, slipr0(1), backst0(1), backst0c(1), sigfs, scysFile,  &
                                 FILE_E)
      else  ! crystalID=kHCP
         !call VerticesSCYS_HCP(numvtx, kappa0(1), sigfs, scysFile,  &
         call VerticesSCYS_HCP(numvtx, slipr0(1), backst0(1), backst0c(1), sigfs, scysFile,  &
                               FILE_E, filePath)
      endif
!
!---- For HCP crystals, the overstress will keep the same initial 
!---- difference in slip system hardness between basal/prismatic
!---- and pyramidal slip systems. 
!---- For FCC/BCC crystals, it won't have any effect.
!
      do is = 1, numslip
         overstress(is) = (tauSlip(is) - 1.0) * slipr0(1)
      enddo

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE CrystalSolverData(  &
         props, mprops,  &
         maxIterState, maxIterNewt,  &
         tolerState, tolerNewt,  &
         FILE_I, FILE_E  &
         )

      implicit none

      integer mprops
      real*8  props(mprops)
     
      integer maxIterState, maxIterNewt
      real*8  tolerState, tolerNewt

      integer FILE_I, FILE_E
!
!---------------------------------------------------------------------72
!
!------- number iterations and tolerance for state iterations
!
      read(FILE_I, *) maxIterState, tolerState
!       maxIterState = nint (props(19))
!       tolerState   = props(20)
!
!------- number iterations and tolerance for newton method
!
      read(FILE_I, *) maxIterNewt, tolerNewt
!       maxIterNewt = nint (props(21))
!       tolerNewt   = props(22)
!
!------- echo input data
!
      write(FILE_E, 1000) maxIterState, maxIterNewt,  &
                          tolerState, tolerNewt      
!
!------- format
!
1000  format(/'*-----   Local Convergence Control Parameters  -----*'/,  &
              7x,'  max iters State  = ',i5 /,  &
              7x,'  max iters Newton = ',i5 /,  &
              7x,'  tolerance State  = ',e12.5 /,  &
              7x,'  tolerance Newton = ',e12.5)
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE CrystalInitializeArrays(  &
         numqpt, numel,  &
         numgrn, numslip,  &
         !kappa0,  &
         slipr0,  &
         backst0,  &
         backst0c,  &
         gcrot0,  &
         gstress,  &
         gstress_n,  &
         gestran,  &
         gestran_n,  &
         !gkappa,  &
         !gkappa_n,  &
         gslipr,  &
         gslipr_n,  &
         gbackst,  &
         gbackst_n,  &
         gbackstc,  &
         gbackstc_n,  &
         gstatev,  &
         gstatev_n,  &
         geqvalues,  &
         ggamdot,  &
         ggamdot2,  &
         gcrot,  &
         gcrot_n,  &
         grrot,  &
         grrot_n  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numqpt, numel
      integer numgrn, numslip
      !real*8  kappa0(NKAPP)
      real*8  slipr0(NKAPP)
      real*8  backst0(NKAPP), backst0c(NKAPP)
      real*8  gcrot0 (DIMS, DIMS, 1, NUMQPT_T, NUMEL_T)

      real*8  gstress    (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstress_n  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran    (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      !real*8  gkappa     (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      !real*8  gkappa_n   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gslipr     (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gslipr_n   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      
      real*8  gbackst     (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackst_n   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackstc    (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackstc_n  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      
      real*8  gstatev    (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n  (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues  (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot    (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot2   (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot      (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot      (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

      integer ie, ip, ig
!
!---------------------------------------------------------------------72
!
!------- initialize arrays used in XTAL constitutive integration method
!
      do ie = 1, numel
        do ip = 1, numqpt
          do ig = 1, numgrn

            call SetTensor(gstress  (1,ig,ip,ie), pzero, NVECS)
            call SetTensor(gestran  (1,ig,ip,ie), pzero, NVECS)
            !call SetTensor(gkappa   (1,ig,ip,ie), pzero, NKAPP)
            call SetTensor(gslipr   (1,ig,ip,ie), pzero, NKAPP)
            call SetTensor(gbackst  (1,ig,ip,ie), pzero, NKAPP)
            call SetTensor(gbackstc (1,ig,ip,ie), pzero, NKAPP)
            call SetTensor(gstatev  (1,ig,ip,ie), pzero, NSTAV)
            call SetTensor(gstress_n(1,ig,ip,ie), pzero, NVECS)
            call SetTensor(gestran_n(1,ig,ip,ie), pzero, NVECS)
            !call SetTensor(gkappa_n (1,ig,ip,ie), pzero, NKAPP)
            call SetTensor(gslipr_n (1,ig,ip,ie), pzero, NKAPP)
            call SetTensor(gbackst_n(1,ig,ip,ie), pzero, NKAPP)
            call SetTensor(gbackstc_n(1,ig,ip,ie), pzero, NKAPP)
            call SetTensor(gstatev_n(1,ig,ip,ie), pzero, NSTAV)
            call SetTensor(geqvalues(1,ig,ip,ie), pzero, NEQVA)
            call SetTensor(ggamdot  (1,ig,ip,ie), pzero, MAX_SLIP)
            call SetTensor(ggamdot2 (1,ig,ip,ie), pzero, MAX_SLIP)

            call SetTensor(gcrot  (1,1,ig,ip,ie), pzero, DIMS*DIMS)
            call SetTensor(gcrot_n(1,1,ig,ip,ie), pzero, DIMS*DIMS)
            call SetTensor(grrot  (1,1,ig,ip,ie), pzero, DIMS*DIMS)
            call SetTensor(grrot_n(1,1,ig,ip,ie), pzero, DIMS*DIMS)

          enddo
        enddo
      enddo
!
!------- initialize state variables and rotation tensors
!
      do ie = 1, numel
        do ip = 1, numqpt
          do ig = 1, numgrn

            !call EqualTensors(kappa0, gkappa_n(1,ig,ip,ie), NKAPP)
            call EqualTensors(slipr0, gslipr_n(1,ig,ip,ie), NKAPP)
            call EqualTensors(backst0, gbackst_n(1,ig,ip,ie), NKAPP)
            call EqualTensors(backst0c, gbackstc_n(1,ig,ip,ie), NKAPP)
            call EqualTensors(gcrot0(1,1,ig,ip,ie),  &
                               gcrot_n(1,1,ig,ip,ie), DIMS*DIMS)
            call EqualTensors(Ident2nd,  &
                               grrot_n(1,1,ig,ip,ie), DIMS*DIMS)

          enddo
        enddo
      enddo
	  !write(*,*) "CIArrays: done"
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE CrystalCloseIOFiles( )

      implicit none
      include 'params_xtal.inc'
!
!---------------------------------------------------------------------72
!
!------- close files
!
      close(XTAL_I)
      close(XTAL_O)
      close(XTAL_E)
!
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE CrystalCloseIOFiles_2( )

      implicit none
      include 'params_xtal.inc'
!
!---------------------------------------------------------------------72
!
!------- close files
!
      close(XTAL_STRESS_O)
      close(XTAL_STRAIN_O)
      close(XTAL_EFFSS_O)
      close(XTAL_TRUESS_O)
      close(XTAL_ITER_O)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE CrystalInitializeMatrxCep(  &
         )
!
      use     TransData
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer i, j
      real*8  unit4(DIMS2, DIMS2), unitDev4(DIMS2, DIMS2)

!
!---------------------------------------------------------------------72
!
!----- zero 4th order unit tensors and transformation matrices

      call SetTensor(unit4, pzero, DIMS2*DIMS2)
      call SetTensor(unitDev4, pzero, DIMS2*DIMS2)

      call SetTensor(fDevMat5x6, pzero, NVECS*DIMS2)
      call SetTensor(fDevMat6x5, pzero, DIMS2*NVECS)
      call SetTensor(fMatTId5x6, pzero, NVECS*DIMS2)
!
!----- set 4th order unit tensor (I)
!
      do i = 1, DIMS2
         unit4(i, i) = pone
      enddo
!
!----- set 4th order deviatoric unit tensor (I_dev) 
!
      do i = 1, DIMS2/2
         unitDev4(i, i) = pone
         do j = 1, DIMS2/2
            unitDev4(i, j) = unitDev4(i, j) - pthird
         enddo
      enddo
 
      do i = DIMS2/2+1, DIMS2
         unitDev4(i, i) = pone
      enddo
!
!----- set tranformation matrix [T] = [tDevMat]_5x6:
!-----    {}'_5x1 = [tDevMat]_5x6 {}'_6x1
!
      fDevMat5x6(1,1) =  pone/sqr2
      fDevMat5x6(1,2) = -pone/sqr2
      fDevMat5x6(2,3) = sqr32
      fDevMat5x6(3,4) = sqr2
      fDevMat5x6(4,5) = sqr2
      fDevMat5x6(5,6) = sqr2
!
!----- set tranformation matrix [J] = [tDevMat]_6x5:
!-----    {}'_6x1 = [tDevMat]_6x5 {}'_5x1
!
      fDevMat6x5(1,1) =  pone/sqr2
      fDevMat6x5(1,2) = -pone/(sqr2*sqr3)
      fDevMat6x5(2,1) = -pone/sqr2
      fDevMat6x5(2,2) = -pone/(sqr2*sqr3)
      fDevMat6x5(3,2) = sqr23
      fDevMat6x5(4,3) = pone/sqr2
      fDevMat6x5(5,4) = pone/sqr2
      fDevMat6x5(6,5) = pone/sqr2
!
!----- needed product: [TId]_5x6 = [T]_5x6 [I_dev]_6x6
!
      call MultAxB_G(fDevMat5x6, unitDev4, fMatTId5x6, NVECS, DIMS2,  &
                     DIMS2) 

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE GlobalIterationParams(  &
         )

      implicit none
      include 'params_xtal.inc'

      integer maxIterN
      real*8  tolerN, divTolerN, zeroTolerN


      logical alwaysLineSearch
      integer searchIters
      real*8  maxStepSize, orthoToler


      integer maxIncrCuts, quickSolveTol, quickSeriesTol
      real*8  DEQP_ref

!
!---------------------------------------------------------------------72
!
!------- params_xtal for global newton iterations
!
      read(XTAL_I, *) maxIterN
      read(XTAL_I, *) tolerN, zeroTolerN, divTolerN
!
!------- params_xtal for global line search iterations
!
      read(XTAL_I, *) alwaysLineSearch
      read(XTAL_I, *) searchIters
      read(XTAL_I, *) maxStepSize, orthoToler
!
!------- params_xtal to control the increase/decrease of time step
!
      read(XTAL_I, *) maxIncrCuts
      read(XTAL_I, *) quickSolveTol
      read(XTAL_I, *) quickSeriesTol
      read(XTAL_I, *) DEQP_ref
!
!------- echo input data 
!
      write(XTAL_E, 1000) maxIterN, tolerN, zeroTolerN, divTolerN
      write(XTAL_E, 2000) alwaysLineSearch, searchIters, maxStepSize,  &
                         orthoToler
      write(XTAL_E, 3000) maxIncrCuts, quickSolveTol, quickSeriesTol,  &
                         DEQP_ref
!
!------- formats
! 
1000  format(/'*-----   Global Newton Iterations   -----*'/,  &
             7x,'  maxIterN    = ',i6/  &
             7x,'  tolerN      = ',e12.5/  &
             7x,'  zeroTolerN  = ',e12.5/  &
             7x,'  divTolerN   = ',e12.5)
2000  format(/'*-----   Global Line Search Iterations   -----*'/,  &
             7x,'  alwaysLS    = ',l6/  &
             7x,'  searchIters = ',i6/  &
             7x,'  maxStepSize = ',e12.5/  &
             7x,'  orthoToler  = ',e12.5)
3000  format(/'*-----   Reset Params for controlling dtime   -----*'/,  &
             7x,'  maxIncrCuts    = ',i6/  &
             7x,'  quickSolveTol  = ',i6/  &
             7x,'  quickSeriesTol = ',i6/  &
             7x,'  DEQP_ref       = ',e12.5)
          
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SetCrystalType(  &
         crystalID, numgrn, props, mprops, FILE_I, FILE_E, FILE_O  &
         )

      implicit none
      include 'params_xtal.inc'

      integer crystalID, numgrn, FILE_I, FILE_E, FILE_O
      integer mprops
      real*8  props(mprops)
!
!---------------------------------------------------------------------72
!
!------- crystal type and number of grains/IP
!
      read(FILE_I, *) crystalID, numgrn
!      crystalID = nint (props(1))
!      numgrn   = nint (props(2))
      if (crystalID .ne. kFCC .and.  &
          crystalID .ne. kBCC .and.  &
          crystalID .ne. kHCP)  &
         call RunTimeError(FILE_O, 'SetCrystalType: crystalID ?')
!
!------- echo input data
      write(FILE_E, 1000) crystalID, numgrn
!
!------- format
!
1000  format(/'*-----   Crystal Type and Number of Grns per IP -----*'/,  &
              7x,'  crystal type     = ',i5 /,  &
              7x,'  number grains/ip = ',i5)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SetCrystalElasticity(  &
         crystalID, fCeDev, fCeiDev, matProp, fCeDevVol, fCeVol,  &
         props, mprops, SXFILE_I, FILE_E, FILE_O  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer crystalID, SXFILE_I, FILE_E, FILE_O
      real*8  fCeDevVol(NVECS), fCeVol
      real*8  matProp(NPROPS)
      real*8  fCeDev(NVECS,NVECS), fCeiDev(NVECS,NVECS)

      integer mprops
      real*8  props(mprops)

      integer elastID, i, j
      real*8  eMod, eNu, gMod, eBulk
      real*8  c1, c2, c3, c4, c5, const

!
!---------------------------------------------------------------------72
!
!------- type of elasticity (in crystal coordinates)
!-------    elastID = kELAS_ISO : isotropic elasticity
!-------    elastID = kELAS_ANI : anisotropic elasticity
!
      read(SXFILE_I, *) elastID
!      elastID = nint (props(3))
      if (elastID .ne. kELAS_ISO .and.  &
          elastID .ne. kELAS_ANI)  &
         call RunTimeError(FILE_O, 'SetCrystalElasticity: elastID ?')
      write(FILE_E, 1000) elastID
!
!------- initialize deviatoric elastic stiffness, its inverse and
!-------  coupled deviatoric-volumetric term
!
      call SetTensor(fCeDev, pzero, NVECS*NVECS)
      call SetTensor(fCeiDev, pzero, NVECS*NVECS)
      call SetTensor(fCeDevVol, pzero, NVECS)
!
!------- read elastic constants, i.e.,
!-------     ISO   FCC  HCP
!-------     eMod  C11  C11
!-------     eMu   C12  C12
!-------           C44  C13
!-------                C33
!-------                C44
!------- and build deviatoric elastic stiffness, as well as
!------- coupled deviatoric-volumetric and volumetric terms.
!-------  For Isotropic and Cubic Crystala: 
!-------      fCeDevVol = 0, fCeVol = Bulk Modulues
!-------  For Hexagonal Closed Packed Crystals:
!-------      si_dev = fCeDev     *ei_dev + fCeDevVol*e_kk
!-------      press  = fCeDevVol^T*ei_dev + fCeVol   *e_kk
!-------   Here: 
!-------      si_dev = {(s11-s22)/sq2,sq32*s33,sq2*s12,sq2*s13,sq2*s23}
!-------      ei_dev = {(e11-e22)/sq2,sq32*e33,sq2*e12,sq2*e13,sq2*e23}
!

      if (elastID .eq. kELAS_ISO) then
!
!---------- isotropic elastic stiffness
         read(SXFILE_I, *)  eMod, eNu 
!         eMod = props(4)
!         eNu  = props(5)
         write(FILE_E, 2000) eMod, eNu
         matProp(1) = eMod / (1. + eNu) / 2.       ! gMod
         matProp(2) = eMod / (1. - 2. * eNu) / 3.  ! eBulk
         fCeDev(1,1) = 2.0*matProp(1)
         fCeDev(2,2) = 2.0*matProp(1)
         fCeDev(3,3) = 2.0*matProp(1)
         fCeDev(4,4) = 2.0*matProp(1)
         fCeDev(5,5) = 2.0*matProp(1)
         fCeVol      = matProp(2)
      elseif ((elastID .eq. kELAS_ANI .and. crystalID .eq. kFCC) .or.  &
              (elastID .eq. kELAS_ANI .and. crystalID .eq. kBCC)) then
!
!---------- anisotropic elastic stiffness, FCC or BCC
         read(SXFILE_I, *) c1, c2, c3
!         c1 = props(4)             ! C11
!         c2 = props(5)             ! C12
!         c3 = props(6)             ! C44
         write(FILE_E, 3000) c1, c2, c3
         matProp(1) = (2.0*(c1 - c2) + 6.0*c3) / 10.0
         matProp(2) = (c1 + 2.0*c2) / 3.0
         fCeDev(1,1) = c1 - c2
         fCeDev(2,2) = c1 - c2
         fCeDev(3,3) = 2.0*c3
         fCeDev(4,4) = 2.0*c3
         fCeDev(5,5) = 2.0*c3
         fCeVol      = (c1 + 2.0*c2) / 3.0
      elseif (elastID .eq. 2 .and. crystalID .eq. kHCP) then
!
!---------- anisotropic elastic stiffness, HCP
         read(SXFILE_I, *) c1, c2, c3, c4, c5
!         c1 = props(4)             ! C11
!         c2 = props(5)             ! C12
!         c3 = props(6)             ! C13
!         c4 = props(7)             ! C33
!         c5 = props(8)             ! C44
         write(FILE_E, 4000) c1, c2, c3, c4, c5
         const = c1 + c2 - 4.0*c3 + 2.0*c4
         matProp(1) = (6.0*(c1-c2) + const + 12.0*c5) / 30.0
         matProp(2) = ((c1 + c2)*c4 - 2.0*c3*c3) / const
         fCeDev(1,1) = c1 - c2
         fCeDev(2,2) = const / 3.0
         fCeDev(3,3) = c1 - c2
         fCeDev(4,4) = 2.0*c5
         fCeDev(5,5) = 2.0*c5
         fCeDevVol(2) = - dsqrt(2.d0/27.d0) * (c1 + c2 - c3 - c4)
         fCeVol       = (2.0*c1 + 2.0*c2 + 4.0*c3 + c4) / 9.0
      endif         
!
!------- effective shear and bulk modulus (to compute effective Ce)
!
      gMod  = matProp(1)
      eBulk = matProp(2)
!
!------- inverse of deviatoric elastic stiffness
!
      do i = 1, NVECS
         fCeiDev(i,i) = 1. / fCeDev(i,i)
      enddo
!
!------- echo some computed data
!
      write(FILE_E, 5000) matProp(1), matProp(2), fCeVol,  &
                          (fCeDevVol(i),i=1,NVECS)
      write(FILE_E, 6000) ((fCeDev(i,j), j=1,5), i=1,5)
!
!------- format
!
1000  format(/'*-----   Crystal Elasticity -----*'/,  &
              7x,'  elasticity type  = ',i5)
2000  format( 7x,'  young modulus    = ',e12.5/  &
              7x,'  poisson ratio    = ',e12.5)
3000  format( 7x,'  c_1              = ',e12.5/  &
              7x,'  c_2              = ',e12.5/  &
              7x,'  c_3              = ',e12.5)
4000  format( 7x,'  c_1              = ',e12.5/  &
              7x,'  c_2              = ',e12.5/  &
              7x,'  c_3              = ',e12.5/  &
              7x,'  c_4              = ',e12.5/  &
              7x,'  c_5              = ',e12.5)
5000  format( 7x,'  shear modulus    = ',e12.5/  &
              7x,'  bulk modulus     = ',e12.5/  &
              7x,'  fCeVol           = ',e12.5/  &
              7x,'  fCeDevVol(1...5) = ',5(e12.5,1x)/  &
              7x,'  CeDev : ')
6000  format((15x,5(e12.5,2x)))

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SetCrystalPlasticity(  &
         crystalID, xtalProp, crss0, props, mprops, mode,  &
         SXFILE_I, FILE_E  &
         )

      implicit none
      include 'params_xtal.inc'

      integer crystalID, mode, SXFILE_I, FILE_E
      real*8  crss0
      real*8  xtalProp(NPROPS)

      integer mprops
      real*8  props(mprops)

      integer i
!
!---------------------------------------------------------------------72
!
!------- slip system kinetic equation: power law (~thermal activation)
!-------   xtalProp(3) xm   : strain-rate sensitivity (m)
!-------   xtalProp(4) gam0 : reference shear rate (d_0)
!
      read(SXFILE_I, *)  xtalProp(3), xtalProp(4)
!      xtalProp(3) = props( 9)
!      xtalProp(4) = props(10)
!
!------ bounds for argument of power law
!
      !call BoundForArgPowLaw(xtalProp(3))
!
!------- slip system kinetic equation: drag-controlled plastic flow
!-------   xtalProp(11) bdrag : drag coefficient (B)
!
      read(SXFILE_I, *)  xtalProp(11)
!
!------ (rss/crss)_limit to switch from (powerLaw & drag) to pure drag
!------  xtalProp(12) = rss/crss_limit
!
      !call LimitRatioForDrag(xtalProp)
!
!------- slip system hardening law : Voce's type
!-------   xtalProp(5) h0    : initial work hardening rate (Theta_0)
!-------   xtalProp(6) tausi : initial slip system strength (tau_0)
!-------   xtalProp(7) taus0 : saturation threshold strength (tau_s0)
!-------   xtalProp(8) xms   : strain-rate sensitivity - state evol (m')
!-------   xtalProp(9) gamss0: ref deformation rate - state evol 
!
      read(SXFILE_I, *)  (xtalProp(i), i = 5,9)
!      xtalProp(5) = props(11)
!      xtalProp(6) = props(12)
!      xtalProp(7) = props(13)
!      xtalProp(8) = props(14)
!      xtalProp(9) = props(15)
!
!------- initial value of state variables: slip system hardness crss0
!
      read(SXFILE_I, *) crss0          ! kappa0
!      crss0 = props(16)
      xtalProp(10) = crss0
!
!------- echo date
!
      write(FILE_E, 1000) mode, xtalProp(3), xtalProp(4), xtalProp(11),  &
                          xtalProp(12)
      write(FILE_E, 2000) (xtalProp(i), i=5,9), crss0
!
!------- format
!
1000  format(/'*-----   Crystal Plasticity, mode', i3, '-----*'/,  &
              7x,'  Kinetic Equation : '/,  &
              7x,'  m                = ',e12.5/  &
              7x,'  gam0             = ',e12.5/  &
              7x,'  bdrag            = ',e12.5/  &
              7x,'  (rss/crss)_limit = ',e12.5)
2000  format( 7x,'  Hardening Law : '/,  &
              7x,'  h0               = ',e12.5/  &
              7x,'  tausi            = ',e12.5/  &
              7x,'  taus0            = ',e12.5/  &
              7x,'  xms              = ',e12.5/  &
              7x,'  gamss0           = ',e12.5/  &
              7x,'  crss0 (kappa0)   = ',e12.5)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SetCrystalSlipGeometry(  &
         crystalID, numslip, numvtx, tauSlip, scysFile, zBar0, &
         pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0, KBar0, &
         KpBar0, KqBar0, KpBar0Vec, KpBar0devVec, KqBar0Vec, KpKpTBar0, &
         FILE_E, filePath  &
         )
      
      implicit none
      include 'params_xtal.inc'

      character*80 filePath
      integer crystalID, numslip, numvtx, FILE_E
      real*8  tauSlip(MAX_SLIP)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), qBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0Vec(NVECS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      real*8  KBar0(DIMS, DIMS, MAX_SLIP)
      real*8  KpBar0(DIMS, DIMS, MAX_SLIP), KqBar0(DIMS, DIMS, MAX_SLIP)
      real*8  KpBar0Vec(NVECS, MAX_SLIP), KpBar0devVec(NVECS, MAX_SLIP), KqBar0Vec(DIMS, MAX_SLIP)
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      real*8  KpKpTBar0(NVECS, NVECS, MAX_SLIP)
      character*(*)  scysFile

      real*8  vecM(DIMS, MAX_SLIP), vecS(DIMS, MAX_SLIP)
!
!---------------------------------------------------------------------72
!
!------- Set up slip system vectors based on type of crystal
!
      if (crystalID .eq. kFCC)  &
          call SetSlipSystemFCC(numslip, vecM, vecS, scysFile, tauSlip,  &
                                FILE_E)

      if (crystalID .eq. kBCC)  &
          call SetSlipSystemBCC(numslip, vecM, vecS, scysFile, tauSlip,  &
                                FILE_E)
    

      if (crystalID .eq. kHCP)  &
          call SetSlipSystemHCP(numslip, numvtx, vecM, vecS, scysFile,  &
                                tauSlip, FILE_E, filePath)
!
!------- Set up Schmid tensors/vectors for each slip system
!
      call SetSlipSystemTensors(numslip, vecM, vecS, zBar0,  &
                         pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0, &
                         KBar0, KpBar0, KqBar0, KpBar0Vec, KpBar0devVec, KqBar0Vec, &
                         KpKpTBar0, FILE_E)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SetSlipSystemFCC(  &
         numslip, vecM, vecS, scysFile, tauSlip, FILE_E  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, FILE_E
      real*8  vecM(DIMS, MAX_SLIP), vecS(DIMS, MAX_SLIP)
      real*8  tauSlip(MAX_SLIP)
      character*(*)  scysFile

      integer is, i
      real*8 indexM(3,12), indexS(3,12)
      real*8 sDotm(12)
      real*8 InnerProductVec

      data indexM /1.,  1., -1.,  &
                   1.,  1., -1.,  &
                   1.,  1., -1.,  &
                   1., -1., -1.,  &
                   1., -1., -1.,  &
                   1., -1., -1.,  &
                   1., -1.,  1.,  &
                   1., -1.,  1.,  &
                   1., -1.,  1.,  &
                   1.,  1.,  1.,  &
                   1.,  1.,  1.,  &
                   1.,  1.,  1./
      data indexS /0.,  1.,  1.,  &
                   1.,  0.,  1.,  &
                   1., -1.,  0.,  &
                   0.,  1., -1.,  &
                   1.,  0.,  1.,  &
                   1.,  1.,  0.,  &
                   0.,  1.,  1.,  &
                   1.,  0., -1.,  &
                   1.,  1.,  0.,  &
                   0.,  1., -1.,  &
                   1.,  0., -1.,  &
                   1., -1.,  0./

!
!---------------------------------------------------------------------72
!
!------- set number of slip systems and slip system reference stress
!
      numslip = kSlipFCC
      call SetTensor(tauSlip, pone, numslip)
!
!------- slip system normals and slip directions: unit vectors
!
      do is = 1, numslip
         call UnitVector(indexS(1,is), vecS(1,is), DIMS)
         call UnitVector(indexM(1,is), vecM(1,is), DIMS)
      enddo
!
!------- file with RI vertex stresses
!
      scysFile = 'vert_fcc.028.00'     ! not used
!
!------- check normality of vecS and vecM
!
      do is = 1, numslip
         sDotm(is) = InnerProductVec(vecS(1,is), vecM(1,is), DIMS)
      enddo
!
!------- echo values
!
      write(FILE_E, 1000) numslip, scysFile
      do is = 1, numslip
         write(FILE_E, 2000) is, tauSlip(is), (vecS(i,is),i=1,3),  &
                             (vecM(i,is),i=1,3), sDotm(is)
      enddo
!
!------- format
!
1000  format(/'*-----   Slip Systems for FCC -----*'/,  &
              7x,'  number slip syst = ',i4/  &
              7x,'  vtx stress file  = ',a15/  &
              7x,'  SS#',5x,'tau',18x,'vecS',30x,'vecM',20x,'SdotM')
2000  format( 7x, i4, 4x, f5.2, 4x, 3f10.5, 4x, 3f10.5, 4x, f10.5)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SetSlipSystemBCC(  &
         numslip, vecM, vecS, scysFile, tauSlip, FILE_E  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
      
      integer numslip, FILE_E
      real*8  vecm(DIMS, MAX_SLIP), vecs(DIMS, MAX_SLIP)
      real*8  tauSlip(MAX_SLIP)
      character*(*)  scysFile

      integer is, i
      real*8 indexM(3,12), indexS(3,12)
      real*8 sDotm(12)
      real*8 InnerProductVec

      data indexM /0.,  1.,  1.,  &
                   1.,  0.,  1.,  &
                   1., -1.,  0.,  &
                   0.,  1., -1.,  &
                   1.,  0.,  1.,  &
                   1.,  1.,  0.,  &
                   0.,  1.,  1.,  &
                   1.,  0., -1.,  &
                   1.,  1.,  0.,  &
                   0.,  1., -1.,  &
                   1.,  0., -1.,  &
                   1., -1.,  0./
      data indexS /1.,  1., -1.,  &
                   1.,  1., -1.,  &
                   1.,  1., -1.,  &
                   1., -1., -1.,  &
                   1., -1., -1.,  &
                   1., -1., -1.,  &
                   1., -1.,  1.,  &
                   1., -1.,  1.,  &
                   1., -1.,  1.,  &
                   1.,  1.,  1.,  &
                   1.,  1.,  1.,  &
                   1.,  1.,  1./
!
!---------------------------------------------------------------------72
!
!------- set number of slip systems and slip system reference stress
!
      numslip = kSlipBCC
      call SetTensor(tauSlip, pone, numslip)
!
!------- slip system normals and slip directions: unit vectors
!
      do is = 1, numslip
         call UnitVector(indexS(1,is), vecS(1,is), DIMS)
         call UnitVector(indexM(1,is), vecM(1,is), DIMS)
      enddo
!
!------- file with RI vertex stresses
!
      scysFile = 'vert_bcc.028.00'     ! not used
!
!------- check normality of vecS and vecM
!
      do is = 1, numslip
         sDotm(is) = InnerProductVec(vecS(1,is), vecM(1,is), DIMS)
      enddo
!
!------- echo values
!
      write(FILE_E, 1000) numslip, scysFile
      do is = 1, numslip
         write(FILE_E, 2000) is, tauSlip(is), (vecS(i,is),i=1,3),  &
                             (vecM(i,is),i=1,3), sDotm(is)
      enddo
!
!------- format
!
1000  format(/'*-----   Slip Systems for BCC -----*'/,  &
              7x,'  number slip syst = ',i4/  &
              7x,'  vtx stress file  = ',a15/  &
              7x,'  SS#',5x,'tau',18x,'vecS',30x,'vecM',20x,'SdotM')
2000  format( 7x, i4, 4x, f5.2, 4x, 3f10.5, 4x, 3f10.5, 4x, f10.5)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SetSlipSystemHCP(  &
         numslip, numvtx, vecM, vecS, scysFile, tauSlip,  &
         FILE_E, filePath  &
         )

      implicit  none
      include   'params_xtal.inc'

      character*80 filePath
      integer   numvtx, numslip, FILE_E
      real*8    vecM(DIMS, MAX_SLIP), vecS(DIMS, MAX_SLIP)
      real*8    tauSlip(MAX_SLIP)

      integer   nmodesx, nmodes, kount, modex, nsmx, mode(5)
      integer   i, j, is, ix, io, nm
      real*8    rca, tau
      real*8    vecm4(48, 4), vecs4(48, 4), sDotm(48)
      character prosa(9)*8
      character*(*) scysFile

      real*8 InnerProductVec

      integer      lengthFile
      character*80 filename
!
!---------------------------------------------------------------------72
!
!     need input file 'slip_hcp.in'

      lengthFile = index(filePath,' ') - 1
      filename   = filePath(1:lengthFile)//'slip_hcp.in'

      io = 99
!      open(unit=io, file='slip_hcp.in', status='old')
      open(unit=io, file=filename, status='old')

      read(io,*) rca
      read(io,*) nmodesx
      read(io,*) nmodes
      read(io,*) (mode(i),i=1,nmodes)
      read(io,*) numvtx
      read(io,5) scysFile
5     format(3x,a15)
!      print *, scysFile
!      print *, '  nVertx for HCP = ', numvtx
!
      kount=1
      i=0
      do nm=1,nmodesx
         read(io,6) prosa
6        format(10a8)
         read(io,*) modex,nsmx,tau
         if(modex.ne.mode(kount)) then
            do ix=1,nsmx
               read(io,6) prosa
            enddo
         else
            kount=kount+1
            do is=1,nsmx
               i=i+1
               read(io,*) (vecm4(i,j),j=1,4),(vecs4(i,j),j=1,4)
               vecM(1,i)= vecm4(i,1)
               vecM(2,i)=(vecm4(i,1)+2.*vecm4(i,2))/sqrt(3.)
               vecM(3,i)= vecm4(i,4)/rca
               vecS(1,i)= 3./2.*vecs4(i,1)
               vecS(2,i)=(vecs4(i,1)/2.+vecs4(i,2))*sqrt(3.)
               vecS(3,i)= vecs4(i,4)*rca
               tauSlip(i) = tau
            enddo
         endif
      enddo

      close(io)
      numslip=i
!
!------- slip system normals and slip directions: unit vectors
!
      do is = 1, numslip
         call UnitVector(vecS(1,is), vecS(1,is), DIMS)
         call UnitVector(vecM(1,is), vecM(1,is), DIMS)
      enddo
!
!------- check normality of vecS and vecM
!
      do is = 1, numslip
         sDotm(is) = InnerProductVec(vecS(1,is), vecM(1,is), DIMS)
      enddo
!
!------- echo values
!
      write(FILE_E, 1000) numslip, scysFile
      do is = 1, numslip
         write(FILE_E, 2000) is, tauSlip(is), (vecS(i,is),i=1,3),  &
                             (vecM(i,is),i=1,3), sDotm(is)
      enddo
!
!------- format
!
1000  format(/'*-----   Slip Systems for HCP -----*'/,  &
              7x,'  number slip syst = ',i4/  &
              7x,'  vtx stress file  = ',a15/  &
              7x,'  SS#',5x,'tau',18x,'vecS',30x,'vecM',20x,'SdotM')
2000  format( 7x, i4, 4x, f5.2, 4x, 3f10.5, 4x, 3f10.5, 4x, f10.5)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SetSlipSystemTensors(  &
         numslip, vecM, vecS, zBar0, pBar0, qBar0, pBar0Vec, &
         qBar0Vec, ppTBar0, Kbar0, KpBar0, KqBar0, KpBar0Vec, KpBar0devVec, KqBar0Vec, &
         KpKpTBar0, FILE_E  &
         )

      implicit  none
      include  'params_xtal.inc'

      integer numslip, FILE_E
      real*8  vecM(DIMS, MAX_SLIP), vecS(DIMS, MAX_SLIP)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), pBar0Vec(NVECS, MAX_SLIP)
      real*8  qBar0(DIMS, DIMS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
      
      real*8  psai, pi
      real*8  KBar0(DIMS, DIMS, MAX_SLIP),KBar0_a(DIMS, DIMS, MAX_SLIP)
      real*8  vecXAI(DIMS, MAX_SLIP),KBar0_b(DIMS, DIMS, MAX_SLIP)
      real*8  KpBar0(DIMS, DIMS, MAX_SLIP), KpBar0Vec(NVECS, MAX_SLIP), &
              KpBar0dev(DIMS, DIMS, MAX_SLIP), KpBar0devVec(NVECS, MAX_SLIP)
      real*8  KqBar0(DIMS, DIMS, MAX_SLIP), KqBar0Vec(DIMS, MAX_SLIP)			  
      
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      real*8  KpKpTBar0(NVECS, NVECS, MAX_SLIP)

      integer is, j, k

      pi = dacos(-1.d0)
      psai = 30d0*pi/180d0 ! psai need to be calibrated for climb mechanism
      write(*,*) "psai (rad and degree)=", psai, psai*180d0/pi
!
!---------------------------------------------------------------------72
!
!------- Schmid Orientation Tensors/Vectors in Crystal Coordinates  
!
      do is = 1, numslip
!
!---------- Tensor zBar0
         call OuterProductVec(vecS(1,is), vecM(1,is), zBar0(1,1,is),  &
                              DIMS)
!
!---------- Symmetric and Skew parts of zBar0
         call SymmetrizeTensor(zBar0(1,1,is), pBar0(1,1,is), DIMS)
         call SkewSymmetrizeTensor(zBar0(1,1,is), qBar0(1,1,is), DIMS)
!
!---------- Vector form for pBar0 and qBar0
         call Mat3x3ToVec5x1Symm(pBar0(1,1,is), pBar0Vec(1,is), DIMS)
         call Mat3x3ToVec3x1Skew(qBar0(1,1,is), qbar0Vec(1,is), DIMS)
!
!---------- Form matrix {P}{P}^T
         call OuterProductVec(pBar0Vec(1,is), pBar0Vec(1,is),  &
                              ppTBar0(1,1,is), NVECS)
!
!---------- Cross Product vecXAI = vecS x vecM
         call CrossProductVec(vecS(1,is), vecM(1,is), vecXAI(1,is),  &
                              DIMS)
!               
!     Check cross produc vecXAI
!      do is=1, numslip
!            write(*,*) "Cross produc vecXAI"
!            write(*,*) vecXAI(1,is), vecXAI(2,is), vecXAI(3,is)
!      enddo
!---------- Tensor Product KBar0_a
         call OuterProductVec(vecS(1,is), vecXAI(1,is),  &
                              KBar0_a(1,1,is), DIMS)
!---------- Tensor Product KBar0_b
         call OuterProductVec(vecS(1,is), vecS(1,is),  &
                              KBar0_b(1,1,is), DIMS)
!
!------- Climb tensor KBar0 (paper: c = k + r = kd + kh + r
!                                        = s_i xai_j cos(psai) + s_i s_j sin(psai) )
!
         call AddTensors(dcos(psai),KBar0_a(1,1,is),  &
                         dsin(psai),KBar0_b(1,1,is), KBar0(1,1,is),  &
                         DIMS*DIMS)
!
!---------- Symmetric part KpBar0 of KBar0 (paper: kd+kh) and its deviatoric part (paper: kd)
         call SymmetrizeTensor(KBar0(1,1,is), KpBar0(1,1,is), DIMS)
         call DeviatoricTensor(KpBar0(1,1,is), KpBar0dev(1,1,is), DIMS)
		 call Mat3x3ToVec5x1Symm(KpBar0dev(1,1,is), KpBar0devVec(1,is), DIMS)
!
!---------- Antisymmetric part KqBar0 of KBar0 (paper: r)
         call SkewSymmetrizeTensor(KBar0(1,1,is), KqBar0(1,1,is), DIMS)
!
!---------- Vector form for KpBar0 and KqBar0
         call Mat3x3ToVec5x1Symm(KpBar0(1,1,is), KpBar0Vec(1,is), DIMS)
         call Mat3x3ToVec3x1Skew(KqBar0(1,1,is), Kqbar0Vec(1,is), DIMS)
!
!---------- Form matrix {Kp}{Kp}^T
         call OuterProductVec(KpBar0Vec(1,is), KpBar0Vec(1,is),  &
                              KpKpTBar0(1,1,is), NVECS)         
         
      enddo
!
!------- echo Schmid orientation tensors/vectors
!
      write(FILE_E, 1000)
      do is = 1, numslip
         write(FILE_E, 1500) is
         do j = 1, DIMS
            write(FILE_E, 2000) (zBar0(j,k,is),k=1,DIMS),  &
                                (pBar0(j,k,is),k=1,DIMS),  &
                                (qBar0(j,k,is),k=1,DIMS)
         enddo
      enddo

      write(FILE_E, 3000)
      do is = 1, numslip
         write(FILE_E, 4000) is, (KpBar0Vec(j,is),j=1,NVECS),  &
                                 (KpBar0devVec(j,is),j=1,NVECS),  &
                                 (KqBar0Vec(j,is),j=1,DIMS)
      enddo
!
!-------- echo climb tensors      
      write(FILE_E, 5000)
      do is = 1, numslip
         write(FILE_E, 5500) is
         do j = 1, DIMS
            write(FILE_E, 6000) (KBar0(j,k,is),k=1,DIMS),  &
                                (KpBar0(j,k,is),k=1,DIMS),  &
                                (KpBar0dev(j,k,is),k=1,DIMS),  &
                                (KqBar0(j,k,is),k=1,DIMS)
         enddo
      enddo      
!
!------- format
!
1000  format(/'*-----   Schmid Tensors -----*')
1500  format( 7x, ' SS # ', i2 /  &
             10x, ' zBar0, pBar0, qBar0: ') 
2000  format(10x, 3(3x, 3f9.5))
3000  format(/'*-----   Schmid Vectors -----*'/  &
              7x, ' SS#, KpBar0Vec, KpBar0devVec, KqBar0Vec: ') 
4000  format( 7x, i3, 3(3x, 5f9.5))
5000  format(/'*-----   Climb Tensors -----*')
5500  format( 7x, ' SS # ', i2 /  &
              7x, ' c_ij, k_ij, kd_ij, r_ij: ') 
6000  format(2x, 4(2x, 3f9.3))      

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SetCrystalLatticeOrient(  &
         numgrn, numqpt, numel, kODFout,  props, mprops,  &
         kODF, numor, seed, angles, FILE_I, FILE_E, FILE_O,  &
         filePath, FILE_TXT_OUT  &
         )

      implicit none
      include 'params_xtal.inc'

      integer numgrn, numqpt, numel, kODFout
      !real*8  gcrot0(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

      integer mprops
      real*8  props(mprops)

      character*80 filePath
      integer kODF, numor, seed, FILE_I, FILE_E, FILE_O
      integer FILE_TXT_OUT
      real*8  angles(DIMS, MAX_ORIEN)
!      real*8  euler(DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

      integer iikc, iidr, i, j
      real*8  pi, pi180, piby2, ph

      integer    length1, length2
      character  textureFile*20, filename*80

!
!---------------------------------------------------------------------72
!
      pi = 4.0 * datan(1.0d+00)
      pi180 = pi/180.
!
!------------------------------------- Data read from FILE_I
!
!------- read code to assign angles at each GR/IP/EL
!
      read(FILE_I, *) kODF
!      kODF = nint (props(17))
!
!------- read multiples of increment to output texture
!
      read(FILE_I, *) kODFout
!      kODFout = nint (props(18))
!
!------- read root name for input/output texture file
!
      read(FILE_I,'(a18)') textureFile
!      textureFile = 'texture'

      length1 = index(filePath,' ') - 1
      length2 = index(textureFile,' ') - 1

      filename = filePath(1:length1)//textureFile(1:length2)//'.txti'
      open(unit=XTAL_TXT_IN, file=filename, status='unknown',  &
                                         access='sequential')
      rewind(XTAL_TXT_IN)

      filename = filePath(1:length1)//textureFile(1:length2)//'.txto'
      open(unit=FILE_TXT_OUT, file=filename, status='unknown')
!
!------- echo data read from FILE_I
!
      write(FILE_E, 1000) kODF, kODFout, textureFile
!
!------------------------------------- Data read from XTAL_TXT_IN
!
!------- number of orientations in file
!
      read(XTAL_TXT_IN, *)  numor
      if (numor .lt. numgrn)  &
         call RunTimeError(FILE_O,  &
                             'SetCrystalLatticeOrient: numor < numgrn')
!
!------- read the flag for angle convention and set some constants
!-------   iikc = 0 : angles input in Kocks convention :  (psi,the,phi)
!                 1 : angles input in Canova convention : (ph,th,om)
!                 ph = 90 + phi; th = the; om = 90 - psi
!-------   iidr = 0 : angles input in degrees
!-------          1 : angles input in radians
      read(XTAL_TXT_IN, *) iikc, iidr
      piby2 = 90.
      if (iidr .eq. 1) piby2 = pi / 2.0
!
!-------- read Euler angles in file & storage them in : (Kocks, radians)
!
      do i = 1, numor
         read(XTAL_TXT_IN, *) (angles(j,i), j=1,DIMS)

         if (iikc .eq. 1) then
            ph = angles(1,i)
            angles(1,i) = piby2 - angles(3,i)
            angles(3,i) = ph - piby2
         endif

         if (iidr .eq. 0)  &
           call SetToScaledTensor(pi180, angles(1,i), angles(1,i), DIMS)
      enddo

      if (kODF .eq. kODF_from_abq) then
         angles(1,1) = props(3)
         angles(2,1) = props(4)
         angles(3,1) = props(5)

         if (iikc .eq. 1) then
            ph = angles(1,1)
            angles(1,1) = piby2 - angles(3,1)
            angles(3,1) = ph - piby2
         endif

         if (iidr .eq. 0)  &
           call SetToScaledTensor(pi180, angles(1,1), angles(1,1), DIMS)
      endif
!
!------- seed for a random angle distribution in polycrystal (needed?)
!
      read(XTAL_TXT_IN, *) seed
!
!------- close input file with texture info
!
      close(XTAL_TXT_IN)
!
!------------------------------------- Assign ODF to each GRs/IPs/ELs
!
!      call AssignCrystalODF(numel, numqpt)  ! No of args has changed!
!
!------- formats
!
1000  format(/'*-----   Lattice Orientation -----*'/  &
              7x,'  kODF             = ',i6/  &
              7x,'  kODFout          = ',i6/  &
              7x,'  root name for I/O txt file  = ',a18)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE VerticesSCYS_Cubic(  &
         !numvtx, kappa0, sigfs, scys_cub, FILE_E  &
         numvtx, slipr0, backst0, backst0c, sigfs, scys_cub, FILE_E  &
         )

      implicit none
      include  'params_xtal.inc'
      include  'numbers.inc'

      integer numvtx, FILE_E
      !real*8  kappa0, sigfs(5,MAX_VTX)
      real*8  slipr0, backst0, backst0c, sigfs(5,MAX_VTX)
      character*(*)  scys_cub

      integer nvtx, ivtx, j
      real*8  sigca(5,28), strVtx(3,3)

      data nvtx   /28/
      data sigca  /-1.73204, -1.00000,   .00000,   .00000,   .00000,  &
                    1.73204, -1.00000,   .00000,   .00000,   .00000,  &
                     .00000,  2.00000,   .00000,   .00000,   .00000,  &
                     .00000,   .00000,  1.73204,  1.73204,  1.73204,  &
                     .00000,   .00000, -1.73204,  1.73204,  1.73204,  &
                     .00000,   .00000,  1.73204, -1.73204,  1.73204,  &
                     .00000,   .00000,  1.73204,  1.73204, -1.73204,  &
                     .00000,   .00000,  3.46410,   .00000,   .00000,  &
                     .00000,   .00000,   .00000,  3.46410,   .00000,  &
                     .00000,   .00000,   .00000,   .00000,  3.46410,  &
                    -.86602,  -.50000,   .00000,  1.73204,  1.73204,  &
                    -.86602,  -.50000,   .00000, -1.73204, -1.73204,  &
                    -.86602,  -.50000,   .00000,  1.73204, -1.73204,  &
                    -.86602,  -.50000,   .00000, -1.73204,  1.73204,  &
                     .86602,  -.50000,  1.73204,   .00000,  1.73204,  &
                     .86602,  -.50000, -1.73204,   .00000, -1.73204,  &
                     .86602,  -.50000, -1.73204,   .00000,  1.73204,  &
                     .86602,  -.50000,  1.73204,   .00000, -1.73204,  &
                     .00000,  1.00000,  1.73204,  1.73204,   .00000,  &
                     .00000,  1.00000, -1.73204, -1.73204,   .00000,  &
                     .00000,  1.00000,  1.73204, -1.73204,   .00000,  &
                     .00000,  1.00000, -1.73204,  1.73204,   .00000,  &
                     .86602, -1.50000,  1.73204,   .00000,   .00000,  &
                     .86602, -1.50000, -1.73204,   .00000,   .00000,  &
                     .86602,  1.50000,   .00000,  1.73204,   .00000,  &
                     .86602,  1.50000,   .00000, -1.73204,   .00000,  &
                   -1.73204,   .00000,   .00000,   .00000,  1.73204,  &
                   -1.73204,   .00000,   .00000,   .00000, -1.73204/
!
!---------------------------------------------------------------------72
!
!------- sigca: vertices of single crystal yield surface
!------- {sigca}={-(11-22)/sqr2,sqr32*33,sqr2*32,sqr2*31,sqr2*21}

      call SetTensor(strVtx, pzero, DIMS)

      numvtx = nvtx
      do ivtx = 1, numvtx

!         strVtx(1,1) = -(sigca(1,ivtx) + sigca(2,ivtx) / sqr3) / sqr2
!         strVtx(2,2) =  (sigca(1,ivtx) - sigca(2,ivtx) / sqr3) / sqr2
!         strVtx(3,3) = sigca(2,ivtx) * sqr2 / sqr3
!         strVtx(2,1) = sigca(5,ivtx) / sqr2
!         strVtx(3,1) = sigca(4,ivtx) / sqr2
!         strVtx(3,2) = sigca(3,ivtx) / sqr2

         sigfs(1, ivtx) = -sigca(1, ivtx)
         sigfs(2, ivtx) =  sigca(2, ivtx)
         sigfs(3, ivtx) =  sigca(5, ivtx)
         sigfs(4, ivtx) =  sigca(4, ivtx)
         sigfs(5, ivtx) =  sigca(3, ivtx)

!         call mat3x3ToVec5x1Symm(strVtx, sigfs(1,ivtx), DIMS)
         !call SetToScaledTensor(kappa0, sigfs(1,ivtx), sigfs(1,ivtx), 5)
         call SetToScaledTensor(slipr0, sigfs(1,ivtx), sigfs(1,ivtx), 5)

      enddo
!
!------- echo vettices
!
      write(FILE_E, 1000) numvtx
      do ivtx = 1, numvtx
         write (FILE_E, 2000) ivtx, (sigfs(j,ivtx),j=1,5)
      enddo
!
!------- formats
!
1000  format(/'*-----   Vertices of SCYS - Cubic Crystal -----*'/  &
              7x,'  No of vertices   = ',i4/  &
              7x,'  ivtx     sigfs(1...5)')
2000  format( 7x, i5, 3x, 5(1x, e15.5))

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE VerticesSCYS_HCP(  &
         !numvtx, kappa0, sigfs, scys_hcp, FILE_E, filePath  &
         numvtx, slipr0, backst0, backst0c, sigfs, scys_hcp, FILE_E, filePath  &
         )

      implicit  none
      include   'params_xtal.inc'

      character*80 filePath
      integer   numvtx, FILE_E
      !real*8    kappa0, sigfs(5,MAX_VTX)
      real*8    slipr0, backst0, backst0c, sigfs(5,MAX_VTX)
      character*(*)      scys_hcp

      integer   io, ivtx, j

      integer      length1, length2
      character*80 filename
!
!---------------------------------------------------------------------72
!
!------- sigca: vertices of single crystal yield surface
!-------   {sigca}={(11-22)/sqr2,sqr32*33,sqr2*21,sqr2*31,sqr2*32}
!
      length1 = index(filePath,' ') - 1
      length2 = index(scys_hcp,' ') - 1

      filename = filePath(1:length1)//scys_hcp(1:length2)

      io = 99
!      open(io, file = scys_hcp, status='old')
      open(io, file = filename, status='old')

      do ivtx = 1, numvtx
         read(io,*) (sigfs(j,ivtx),j=1,5) 
         !call SetToScaledTensor(kappa0, sigfs(1,ivtx), sigfs(1,ivtx), 5)
         call SetToScaledTensor(slipr0, sigfs(1,ivtx), sigfs(1,ivtx), 5)
      enddo

      close(io)
!
!------- echo vettices
!
      write(FILE_E, 1000) numvtx
      do ivtx = 1, numvtx
         write (FILE_E, 2000) ivtx, (sigfs(j,ivtx),j=1,5)
      enddo
!
!------- formats
!
1000  format(/'*-----   Vertices of SCYS - HCP Crystal -----*'/  &
              7x,'  No of vertices   = ',i4/  &
              7x,'  ivtx   sigfs(1...5)')
2000  format( 7x, i5, 3x, 5f13.5)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE AssignCrystalODF(grainid, numel, numqpt, noel, npt,  &
         numgrn, numslip, kODF,  &
         numor, seed,  &
         gcrot0,  &
         angles,  &
!         euler,  &
         FILE_O, FILE_TXT_OUT  &
         )

      implicit none
      include  'params_xtal.inc'
!
      character*30 cmname
      integer numel, numqpt, noel, npt
      integer FILE_O, FILE_TXT_OUT

      integer numgrn, numslip, kODF
      integer numor, seed
      real*8  gcrot0(DIMS, DIMS, 1, NUMQPT_T, NUMEL_T)
      real*8  angles(DIMS, MAX_ORIEN)
!      real*8  euler(DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  euler(DIMS, 1, NUMQPT_T, NUMEL_T)
      real    ran1_nr

      integer ie, ip, ig, i, random, grainid
      real*8  pi, pi180

      integer seed_nr
      save    seed_nr
!
!---------------------------------------------------------------------72
!
      pi = 4.0 * datan(1.0d+00)
      pi180 = pi/180.

      if (noel .eq. 1 .and. npt .eq. 1) seed_nr = -seed
      if (noel .eq. 1 .and. npt .eq. 1) write(FILE_TXT_OUT, 2000)

!
!------- same ODF in all ELs and IPs
      if (kODF .eq. kODF_same_all) then
         do ie = 1, numel
            do ip = 1, numqpt
               do ig =1, numgrn
                  call EqualTensors(angles(1,grainid),  &
                                              euler(1,ig,ip,ie), DIMS)
                  !call EqualTensors(angles(1,ig),  &
                                              !euler(1,ig,ip,ie), DIMS)
               enddo
            enddo
         enddo
!
      else
         call RunTimeError(FILE_O, 'setCrystalLatticeOrient: Bad kODF')
      endif
!
!------- build rotation matrices C0: {x}_sm = [C0] {x}_cr 
!
      do ie = 1, numel
         do ip = 1, numqpt
            do ig =1, numgrn
               call AnglesToRotMatrix(euler(1,ig,ip,ie),  &
                                       gcrot0(1,1,ig,ip,ie), DIMS)
            enddo
         enddo
      enddo
!
!------- write initial assigned orientations
!-------  note that numel & numqpt are both dummy variables (=1)
!
!      write(*,*)
      do ie = 1, numel
         do ip = 1, numqpt
            do ig =1, numgrn
                write(FILE_TXT_OUT, 3000)  &
                       (euler(i,ig,ip,ie)/pi180,i=1,3), ig, npt, noel
            enddo
         enddo
      enddo
!
!------- formats
!
2000  format(/'*-----   Initial Assigned Orientations -----*'/  &
              4x,'ang1',4x,'ang2',4x,'ang3',9x,' igrn',7x,'intpt',  &
              6x,'ielem ')
3000  format( 3f8.2, 3(4x,i8))

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SetCrystalGeometry(  &
         crystalID, numslip, numvtx, scysFile, zBar0, pBar0, &
         qBar0, pBar0Vec, qBar0Vec, ppTBar0, KBar0, KpBar0, KqBar0, &
         KpBar0Vec, KpBar0devVec, KqBar0Vec, KpKpTBar0, xtalProp, matProp, &
         !kappa0, tauSlip, props, mprops, SXFILE_I, FILE_E  &
         hardmtx, &
         slipr0, backst0, backst0c, tauSlip, props, mprops, SXFILE_I, FILE_E  &
         )
      
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      character*(*)  scysFile
      integer crystalID, numslip, numvtx, SXFILE_I, FILE_E
      !real*8  kappa0(NKAPP)
      real*8  slipr0(NKAPP), backst0(NKAPP), backst0c(NKAPP)
      real*8  zBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0(DIMS, DIMS, MAX_SLIP), qBar0(DIMS, DIMS, MAX_SLIP)
      real*8  pBar0Vec(NVECS, MAX_SLIP), qBar0Vec(DIMS, MAX_SLIP)
	  real*8  KBar0(DIMS, DIMS, MAX_SLIP)      
	  real*8  KpBar0(DIMS, DIMS, MAX_SLIP), KqBar0(DIMS, DIMS, MAX_SLIP)
      real*8  KpBar0Vec(NVECS, MAX_SLIP), KpBar0devVec(NVECS, MAX_SLIP), KqBar0Vec(DIMS, MAX_SLIP)	        
      real*8  ppTBar0(NVECS, NVECS, MAX_SLIP)
      real*8  KpKpTBar0(NVECS, NVECS, MAX_SLIP)
      real*8  xtalProp(NPROPS), matProp(NPROPS, MAX_SLIP)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP), tauSlip(MAX_SLIP)

      integer mprops
      real*8  props(mprops)

      integer   nmodesx, nmodes, kount, modex, nsmx
      integer   i, j, is, js, nm, im, jm, isensex
      real*8    rca, crss0
      real*8    twshx, isectwx, thres1x, thres2x
      real*8    vecM(DIMS, MAX_SLIP), vecS(DIMS, MAX_SLIP)
      real*8    vecm4(48, 4), vecs4(48, 4), sDotm(48)
      character prosa*80

      real*8 InnerProductVec

      integer   NMFILE
      parameter (NMFILE=12)
      integer   mode(NMFILE), nsm(NMFILE)
      real*8    hselfx(NMFILE), hlatex(NMFILE, NMFILE)
!
!---------------------------------------------------------------------72
!
!---- zero out some arrays
! 
      call SetTensor(vecM,  pzero, DIMS*MAX_SLIP)
      call SetTensor(vecS,  pzero, DIMS*MAX_SLIP)
      call SetTensor(vecm4, pzero, 48*4)
      call SetTensor(vecs4, pzero, 48*4)
!
!---- read SXFILE_I to set slip / twinning data
!
      read(SXFILE_I,*) rca
      read(SXFILE_I,*) nmodesx               ! total # modes in file
      read(SXFILE_I,*) nmodes                ! # modes used (active)
      read(SXFILE_I,*) (mode(i),i=1,nmodes)  ! labels of active modes
      read(SXFILE_I,*) numvtx
      read(SXFILE_I,'(a15)') scysFile

      kount=1                         ! counter for active modes
      i=0                             ! counter for # slip/twin systems
      do nm=1,nmodesx

         read(SXFILE_I,'(a)') prosa
         read(SXFILE_I,*) modex,nsmx,isensex
         if(modex.ne.mode(kount)) then
!........... skip data
            read(SXFILE_I,'(a)') prosa
            read(SXFILE_I,'(a)') prosa
            read(SXFILE_I,'(a)') prosa
            read(SXFILE_I,'(a)') prosa
            read(SXFILE_I,'(a)') prosa
            read(SXFILE_I,'(a)') prosa
            do is=1,nsmx
               read(SXFILE_I,*)
            enddo
         else
!........... PlastProp: 4 lines; Twindata: 1 line; LatentHard: 1 line
            call SetCrystalPlasticity(crystalID, xtalProp, crss0,  &
                               props, mprops, modex, SXFILE_I, FILE_E)
            read(SXFILE_I,*) twshx, isectwx, thres1x, thres2x
            read(SXFILE_I,*) (hlatex(kount, jm), jm=1, nmodes)
            hselfx(kount) = 1.0
            nsm(kount) = nsmx

!........... indices for each slip/twin mode
            if (crystalID .eq. kFCC .or. crystalID .eq. kBCC) then
              do is=1,nsmx
                i=i+1
                !kappa0(i) = crss0
                slipr0(i) = crss0
                call EqualTensors(xtalProp, matProp(1,i), NPROPS)
                read(SXFILE_I,*) (vecM(j,i),j=1,3),(vecS(j,i),j=1,3)
                call UnitVector(vecS(1,i), vecS(1,i), DIMS)
                call UnitVector(vecM(1,i), vecM(1,i), DIMS)
              enddo
            endif

            if (crystalID .eq. kHCP) then
              do is=1,nsmx
                i=i+1
                !kappa0(i) = crss0
                slipr0(i) = crss0
                call EqualTensors(xtalProp, matProp(1,i), NPROPS)
                read(SXFILE_I,*) (vecm4(i,j),j=1,4),(vecs4(i,j),j=1,4)
                vecM(1,i)= vecm4(i,1)
                vecM(2,i)=(vecm4(i,1)+2.*vecm4(i,2))/sqrt(3.)
                vecM(3,i)= vecm4(i,4)/rca
                vecS(1,i)= 3./2.*vecs4(i,1)
                vecS(2,i)=(vecs4(i,1)/2.+vecs4(i,2))*sqrt(3.)
                vecS(3,i)= vecs4(i,4)*rca
                call UnitVector(vecS(1,i), vecS(1,i), DIMS)
                call UnitVector(vecM(1,i), vecM(1,i), DIMS)
              enddo
            endif
            kount=kount+1
         endif

      enddo

      numSlip=i
!
!------- ratio of slip system's kappa0: kappa0(is)/kappa0(1) 
!------- check normality of vecS and vecM
!
      do is = 1, numSlip
         !tauSlip(is) = kappa0(is)/kappa0(1)
         tauSlip(is) = slipr0(is)/slipr0(1)
         sDotm(is)   = InnerProductVec(vecS(1,is), vecM(1,is), DIMS)
      enddo
!
!------- set up latent hardening matrix
!
      i=0
      do im = 1, nmodes
        do is = 1, nsm(im)
          i=i+1
          j=0
          do jm = 1, nmodes
            do js = 1, nsm(jm)
              j=j+1
              hardmtx(i,j) = hlatex(im,jm)
            enddo
          enddo
          hardmtx(i,i)=hselfx(im)
        enddo
      enddo
!
!------- echo values
!
      write(FILE_E, 1000) crystalID, numslip, scysFile
      do is = 1, numslip
         write(FILE_E, 2000) is, tauSlip(is), (vecS(i,is),i=1,3),  &
                             (vecM(i,is),i=1,3), sDotm(is)
      enddo

      write(FILE_E,'(''*-----   Latent Hardening Matrix -----*'')')
      do is = 1, numSlip
         write(FILE_E, '(8x,24F5.1)') (hardmtx(is,js), js=1,numSlip)
      enddo
!
!------- Set up Schmid and climb tensors/vectors for each slip system
!
      call SetSlipSystemTensors(numslip, vecM, vecS, zBar0,  &
                       pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0,  & 
                       KBar0, KpBar0, KqBar0, KpBar0Vec, KpBar0devVec, KqBar0Vec,  &
                       KpKpTBar0, FILE_E)
!
!------- format
!
1000  format(/'*-----   Slip Systems for',i3,' (Crystal Type)-----*'/,  &
              7x,'  number slip syst = ',i4/  &
              7x,'  vtx stress file  = ',a15/  &
              7x,'  SS#',5x,'tau',18x,'vecS',30x,'vecM',20x,'SdotM')
2000  format( 7x, i4, 4x, f5.2, 4x, 3f10.5, 4x, 3f10.5, 4x, f10.5)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE stressCrystal(  &
         sigxx, sigyy, sigth, sigxy, epsxx, epsyy, epsth, epsxy,  &
         epsum, spinn, epsxz, epsyz, spinxz, spinyz, sigxz, sigyz,  &
         qptpo, qptp, dtime, time, ielem, incr, kstep, numel,  &
         numqpt, iprint, sigma, ddsdde, statusFlag, numIncrs, noel,  &
         npt, numel_aba, numqpt_aba, crystalvar, crystalvar_n, gcrot0, grainid  &
         )
!
      use     SlipData
      use     DataType
      implicit none
     
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer ielem, incr, numel, numqpt, iprint, statusFlag, numIncrs, kstep, grainid
      integer noel, npt, numel_aba, numqpt_aba
      real*8  dtime, time(2)
      real*8  sigxx(numqpt), sigyy(numqpt), sigth(numqpt), sigxy(numqpt)
      real*8  epsxx(numqpt), epsyy(numqpt), epsth(numqpt)
      real*8  epsxy(numqpt), epsum(numqpt), spinn(numqpt)
      real*8  qptpo(numqpt), qptp(numqpt)
      real*8  sigma(DIMS, DIMS, NUMQPT_T)
      real*8  ddsdde(DIMS2, DIMS2, NUMQPT_T)
!
      real*8  gcrot0(DIMS, DIMS, 1, NUMQPT_T, NUMEL_T)
!    
      integer iqpt
      real*8  thetao, theta, d_kk
      real*8  d_vec(NVECS), w_vec(3)

      real*8  epsxz(NUMQPT_T),  epsyz(NUMQPT_T)
      real*8  spinxz(NUMQPT_T), spinyz(NUMQPT_T)
      real*8  sigxz(NUMQPT_T),  sigyz(NUMQPT_T)
      type (xtalVars) crystalvar
      type (xtalVars_n) crystalvar_n
!
!---------------------------------------------------------------------72
!
!------ initialize vectors for symm/skew parts of velocity gradients
!
      call SetTensor(d_vec, pzero, NVECS)
      call SetTensor(w_vec, pzero, 3)
!
!------ loop over integration points
!
      do iqpt = 1, numqpt
!
!---------- recover symm/skew parts of velocity gradient at ielem/iqpt
!---------- eps_ij are deviatoric quantities
!
         d_vec(1) = (epsxx(iqpt) - epsyy(iqpt)) / sqr2 
         d_vec(2) = epsth(iqpt) * sqr32
         d_vec(3) = epsxy(iqpt) * sqr2
         if (NSHR .gt. 1) then
            d_vec(4) = epsxz(iqpt) * sqr2
            d_vec(5) = epsyz(iqpt) * sqr2
         endif
         d_kk = epsum(iqpt)

         w_vec(1) = -spinn(iqpt)
         if (NSHR .gt. 1) then
            w_vec(2) = -spinxz(iqpt)
            w_vec(3) = -spinyz(iqpt)
         endif
!
!---------- recover temperature at ielem/iqpt
!
         thetao = qptpo(iqpt)
         theta  = qptp(iqpt)
!
!---------- evolve state at ielem/iqpt
!
     !write(*,*) '6noel:', noel, '6grainid', grainid
         !write(*,*) 'drive in'
         call DriverCrystalEvolve(sigxx(iqpt), sigyy(iqpt), sigth(iqpt),  &
            sigxy(iqpt), sigxz(iqpt), sigyz(iqpt), d_vec, w_vec, d_kk,  &
            thetao, theta, dtime, time, iqpt, ielem, incr, kstep, numqpt,  &
            numel, iprint, sigma(1, 1, iqpt), ddsdde(1, 1, iqpt),  &
            statusFlag, numIncrs, noel, npt, numel_aba, numqpt_aba, crystalvar, crystalvar_n, gcrot0, grainid)
         if (statusFlag .eq. kGLOBAL_FAILED) return
         !write(*,*) 'drive out'

      enddo

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE DriverCrystalEvolve(  &
         sigxx, sigyy, sigth, sigxy, sigxz, sigyz, d_vec, w_vec, d_kk,  &
         thetao, theta, dtime, time, iqpt, ielem, incr, kstep, numqpt, numel,  &
         iprint, savg_ij, cavg_ijkl, statusFlag, numIncrs, noel, npt,  &
         numel_aba, numqpt_aba, crystalvar, crystalvar_n, gcrot0, grainid  &
         )
!
      use     WriteControl
!      use     InitialLocalRotMatrix
      use     IterData
      use     SlipData
      use     InitialSlipSys
      use     EProp
      use     XtalPar
      use     DataType
      implicit  none
      include  'params_xtal.inc'
      include  'numbers.inc'
!
!---- arguments
!
      integer iqpt, ielem, incr, numqpt, numel, iprint, statusFlag, kstep, grainid
      integer numIncrs, noel, npt, numel_aba, numqpt_aba 
      real*8  sigxx, sigyy, sigth, sigxy, sigxz, sigyz
      real*8  thetao, theta, d_kk, dtime, time(2)
      real*8  d_vec(NVECS), w_vec(DIMS)
      real*8  cavg_ijkl(DIMS2, DIMS2), savg_ij(DIMS, DIMS)
!
!
!      real*8  gstress  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
!      real*8  gestran  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
!      real*8  gslipr   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
!      real*8  gbackst   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
!      real*8  gstatev  (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
!      real*8  gstress_n(NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
!      real*8  gestran_n(NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
!
      real*8   gcrot0(DIMS, DIMS, 1, NUMQPT_T, NUMEL_T)
!
!      real*8  gslipr_n (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
!      real*8  gbackst_n (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
!      real*8  gstatev_n(NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
!      real*8  geqvalues(NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
!      real*8  ggamdot  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
!      real*8  gcrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
!      real*8  grrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
!      real*8  gcrot_n  (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
!      real*8  grrot_n  (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

!
!---- local variables
!
      integer iterCounterS, iterCounterN, ierr, igrn, nxtals_nc, i, j
      real*8  tmp, wpnorm_avg
      real*8  s_ij(DIMS,DIMS), ee_ij(DIMS,DIMS), tau(NVECS)

      real*8  stress(NVECS), estran(NVECS), slipr(NKAPP), backst(NKAPP), backstc(NKAPP), statev(NSTAV)
      real*8  stress_n(NVECS), estran_n(NVECS), slipr_n(NKAPP), backst_n(NKAPP), backstc_n(NKAPP)
      real*8  statev_n(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP), gamdotC(MAX_SLIP)
      real*8  crot(DIMS,DIMS), rrot(DIMS,DIMS), drot(DIMS,DIMS)
      real*8  crot_n(DIMS,DIMS), rrot_n(DIMS,DIMS), crot0(DIMS,DIMS)
      real*8  c_ijkl(DIMS2,DIMS2), sdev_ij(DIMS,DIMS)
      real*8  rss(MAX_SLIP), rssc(MAX_SLIP)
      type (xtalVars) crystalvar
      type (xtalVars_n) crystalvar_n
      real*8  sGamPPt(NVECS,NVECS), sGamKpKpt(NVECS,NVECS)
!
!---------------------------------------------------------------------72
!
      sGamPPt=0.d0
      sGamKpKpt=0.d0
!
!---------------------------------------------------------------------72
!
!------- initialize average stress (sij) and consistent moduli (cijkl)
!
      call SetTensor(savg_ij, pzero, DIMS*DIMS)
      call SetTensor(cavg_ijkl, pzero, DIMS2*DIMS2)
!
!------- initialize average value of wp_norm
!
      wpnorm_avg = pzero
!
!------- counter for crystals that do not converge - Large Scale Appls.
      nxtals_nc = 0
!
!------- loop over all grains in aggregate at iqpt/ielem
     !write(*,*) '7noel:', noel, '7grainid', grainid
!
      do igrn = 1, numgrn
!
!----------- zero local arrays
!
     !write(*,*) '7.1noel:', noel, '7.1grainid', grainid
         !write(*,*) 'drv 1'
         !call InitCrystalLocalArrays(tau, stress, estran, kappa, statev,  &
         call InitCrystalLocalArrays(tau, stress, estran, slipr, backst, backstc, statev,  &
            !stress_n, estran_n, kappa_n, statev_n, eqvalues, gamdot,  &
            stress_n, estran_n, slipr_n, backst_n, backstc_n, statev_n, eqvalues,  &
            gamdot, rss, gamdotC, rssc,  &
            crot, rrot, drot, crot_n, rrot_n, crot0, c_ijkl)
     !write(*,*) '8noel:', noel, '8grainid', grainid
!
!---------- fetch variables at ielem/iqpt/igrn
!
         !write(*,*) 'drv 2'
         !call FetchCrystalVariablesAtIP(gstress_n, gestran_n, gkappa_n,  &
         call FetchCrystalVariablesAtIP(  &
            crystalvar_n%gstress_n, crystalvar_n%gestran_n,  &
            crystalvar_n%gslipr_n, crystalvar_n%gbackst_n, crystalvar_n%gbackstc_n, &
            crystalvar_n%gstatev_n,  &
            crystalvar%geqvalues,  &
            crystalvar_n%gcrot_n, crystalvar_n%grrot_n,  &
            gcrot0,  &
            stress_n,  &
            estran_n, slipr_n, backst_n, backstc_n, statev_n, eqvalues, crot_n, rrot_n,  &
            crot0, igrn, iqpt, ielem)

     !write(*,*) '9noel:', noel, '9grainid', grainid
!
!---------- monitor integration/convergence of XTAL constitutive eqns
!
!        Check whether the input data is passed correctly, OK!
         !write(*,*) 'Check input for MonitorCrystalInt: KpBar0devVec, OK'
		 !do i = 1, numslip
         !	write(*,*) (KpBar0devVec(j,i),j=1,NVECS)
         !enddo
         call MonitorCrystalIntegration(d_vec, w_vec, d_kk, s_ij, ee_ij, tau,  &
            !stress, estran, kappa, statev, eqvalues, gamdot, rss,  &
            stress, estran, slipr, backst, backstc, statev, eqvalues,  &
            gamdot, rss, gamdotC, rssc,  &
            crot, rrot,  &
            !drot, stress_n, estran_n, kappa_n, statev_n, crot_n, rrot_n,  &
            drot, stress_n, estran_n, slipr_n, backst_n, backstc_n, statev_n, crot_n, rrot_n,  &
            crot0, fCeDev, fCeiDev, matProp, sigfs, tauSlip,  &
            zBar0, pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0,  &
            KBar0, KpBar0, KqBar0, KpBar0Vec, KpBar0devVec, KqBar0Vec, KpKpTBar0,  &
            numslip, numvtx, maxIterState, maxIterNewt,  &
            tolerState, tolerNewt, dtime, time(2), &
            fCeDevVol, fCeVol,  &
            theta, thetao, igrn, npt, noel, incr, kstep, iterCounterS,  &
            iterCounterN, ierr, c_ijkl, hardmtx, sGamPPt, sGamKpKpt, grainid)
     !write(*,*) '90noel:', noel, '90grainid', grainid
!
!---------- check convergence; if not then:
!              either reset crystal quantities or force time reduction
!
         !write(*,*) 'drv 4'
         if (ierr .ne. XTAL_CONVERGED) then
            call WriteMessage (XTAL_O,  &
                      'DriverXtalEvolve: Sub-Incrs failed!')
            call WriteMessage (XTAL_O,  &
                      'DriverXtalEvolve: nxtals_nc incremented by one')
            nxtals_nc = nxtals_nc + 1 
            write(XTAL_O, 1000) nxtals_nc, incr, igrn, npt, noel

            if (nxtals_nc .ge. NXTALS_NC_LIMIT) then      
              call WriteMessage (XTAL_O,  &
                      'DriverXtalEvolve: nxtals_nc > LIMIT')
              call WriteMessage(XTAL_O,  &
                      '... will force dtime reduction by')
              call WriteMessage(XTAL_O,  &
                      '... setting: status=kGLOBAL_FAILED')
              statusFlag = kGLOBAL_FAILED
              return
            endif

            call WriteMessage(XTAL_O, 'Resetting xtal quantities')
            !call ResetCrystalQnts(stress, estran, kappa, statev,  &
            !write(*,*) 'In'
!            
!           Check whether the input data is passed correctly, OK!
            !write(*,*) 'Check input for ResetCrystalQnts: KpBar0devVec, OK'
		   	!do i = 1, numslip
           	!	write(*,*) (KpBar0devVec(j,i),j=1,NVECS)
         	!enddo            
            call ResetCrystalQnts(stress, estran, slipr, backst, backstc, statev,  &
                   eqvalues, gamdot, gamdotC,  &
                   crot, rrot, s_ij, c_ijkl,  &
                   stress_n,  &
                   !estran_n, kappa_n, statev_n, crot_n, rrot_n,  &
                   estran_n, slipr_n, backst_n, backstc_n, statev_n, crot_n, rrot_n, matProp,  &
                   pBar0Vec, KpBar0Vec, KpBar0devVec,  &
                   fCeiDev, fCeDevVol, fCeVol, numslip, dtime,  &
                   sGamPPt, sGamKpKpt)
            !write(*,*) 'Out'

         endif
!
!---------- save computed variables at ielem/iqpt/igrn
!
         !write(*,*) 'drv 5'
         !call SaveCrystalVariablesAtIP(stress, estran, kappa, statev,  &
         call SaveCrystalVariablesAtIP(stress, estran, slipr, backst, backstc, statev,  &
            !gamdot, eqvalues, crot, rrot, gstress, gestran, gkappa,  &
            gamdot, gamdotC,  &
            eqvalues, crot, rrot,  &
            crystalvar%gstress, crystalvar%gestran,  &
            crystalvar%gslipr, crystalvar%gbackst, crystalvar%gbackstc, crystalvar%gstatev,  &
            crystalvar%geqvalues, crystalvar%ggamdot,  &
            crystalvar%ggamdot2,  &
            crystalvar%gcrot, crystalvar%grrot, igrn, iqpt, ielem)
!
!---------- output computed quantities at selected "noel,npt,igrn"
!
         !write(*,*) 'drv 6'
         if (iprint .eq. kPRINT_MODEL .and.  &
              noel  .eq. kPRINT_ELEM  .and.  &
               npt  .eq. kPRINT_QPT   .and.  &
               igrn .eq. kPRINT_GRN        ) then

            !call OutputQnts(d_vec, s_ij, ee_ij, kappa, gamdot, rss,  &
            call OutputQnts(d_vec, s_ij, ee_ij, slipr, backst, backstc, gamdot, rss,  &
                   gamdotC, rssc,  &
                   statev, eqvalues, crot, rrot, c_ijkl, d_kk, dtime,  &
                   time, numslip, incr, iterCounterS, iterCounterN,  &
                   noel, npt, igrn)

         endif
!
!---------- average values of sij and cijkl for aggregate

         !write(*,*) 'drv 7'
         tmp = 1.d0 / numgrn
         call AddScaledTensor(tmp, s_ij, savg_ij, DIMS*DIMS)
         call AddScaledTensor(tmp, c_ijkl, cavg_ijkl, DIMS2*DIMS2)
!
!---------- average value of wp_norm
!
         !write(*,*) 'drv 8'
         wpnorm_avg = wpnorm_avg + tmp*statev(kWPNORM)

      enddo
!
!------- move aggregate deviatoric stress to main-program arrays
!
      call DeviatoricTensor(savg_ij, sdev_ij, DIMS)
      sigxx = sdev_ij(1,1)
      sigyy = sdev_ij(2,2)
      sigth = sdev_ij(3,3)
      sigxy = sdev_ij(1,2)
      if (NSHR .gt. 1) then
         sigxz = sdev_ij(1,3)
         sigyz = sdev_ij(2,3)
      endif
!
!------- output aggregate qnts at selected "ielem,iqpt"
!-------  note that ielem=1,iqpt=1 always !!
!
!zx      if (iprint .eq. kPRINT_MODEL .and.  &
!zx           ielem .eq. kPRINT_ELEM  .and.  &
!zx           iqpt  .eq. kPRINT_QPT        ) then
!
!---------- write effective stress-strain curve
!
!zx         call WriteStressStrainCurve(sdev_ij, d_vec,  &
!zx                                     wpnorm_avg, dtime, time,  &
!zx                                     numgrn, incr, npt, noel, numel_aba, numqpt_aba)
!
!---------- write texture at specified increments
!
!         if ( (incr/kODFout*kODFout) .eq. incr .or.  &
!                             incr .eq. numIncrs) then
         if ( (kstep==TarStep1 .or. kstep==TarStep2  &
            .or. kstep==TarStep3) .and. ( abs((time(1)+dtime-TarTime2)) .le. 1.0 &
            .or. abs(time(1) +dtime - TarTime1) .le. 1.0e-3 .or. abs(time(1)+dtime-TarTime3) .le. 1.0e-6)) then
            !call WriteTexture(gcrot, gkappa, d_vec, numgrn, incr,  &
            call WriteTexture(crot, slipr, backst, backstc, d_vec, numgrn, kstep, incr,  &
                              iqpt, ielem, npt, noel)

         endif

!zx      endif

1000  format(3x, 'Grains did not converged: ', i3, ':', 3x,  &
                 'incr  #', i8, ';', 1x,  &
                 'grain #', i5, ';', 1x,  &
                 'iqpt  #', i3, ';', 1x,  &
                 'elem  #', i5)

      return

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE InitCrystalLocalArrays(  &
         !tau, stress, estran, kappa, statev, stress_n, estran_n,  &
         tau, stress, estran, slipr, backst, backstc, statev, stress_n, estran_n,  &
         !kappa_n, statev_n, eqvalues, gamdot, rss, crot, rrot, drot,  &
         slipr_n, backst_n, backstc_n, statev_n, eqvalues,  &
         gamdot, rss, gamdotC, rssc,  &
         crot, rrot, drot,  &
         crot_n, rrot_n, crot0, c_ijkl  &
         )

      implicit  none
      include  'params_xtal.inc'
      include  'numbers.inc'

      real*8  tau(NVECS)
      !real*8  stress(NVECS), estran(NVECS), kappa(NKAPP), statev(NSTAV)
      real*8  stress(NVECS), estran(NVECS), slipr(NKAPP), backst(NKAPP), &
      backstc(NKAPP), statev(NSTAV)
      !real*8  stress_n(NVECS), estran_n(NVECS), kappa_n(NKAPP)
      real*8  stress_n(NVECS), estran_n(NVECS), slipr_n(NKAPP),  &
      backst_n(NKAPP), backstc_n(NKAPP)
      real*8  statev_n(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP), gamdotC(MAX_SLIP)
      real*8  crot(DIMS,DIMS), rrot(DIMS,DIMS), drot(DIMS,DIMS)
      real*8  crot_n(DIMS,DIMS), rrot_n(DIMS,DIMS), crot0(DIMS,DIMS)
      real*8  c_ijkl(DIMS2, DIMS2)
      real*8  rss(MAX_SLIP), rssc(MAX_SLIP)
!
!---------------------------------------------------------------------72
!
      call SetTensor(tau,      pzero, NVECS)
      call SetTensor(stress,   pzero, NVECS)
      call SetTensor(estran,   pzero, NVECS)
      !call SetTensor(kappa,    pzero, NKAPP)
      call SetTensor(slipr,    pzero, NKAPP)
      call SetTensor(backst,    pzero, NKAPP)
      call SetTensor(backstc,    pzero, NKAPP)      
      call SetTensor(statev,   pzero, NSTAV)
      call SetTensor(stress_n, pzero, NVECS)
      call SetTensor(estran_n, pzero, NVECS)
      !call SetTensor(kappa_n,  pzero, NKAPP)
      call SetTensor(slipr_n,  pzero, NKAPP)
      call SetTensor(backst_n,  pzero, NKAPP)
      call SetTensor(backstc_n,  pzero, NKAPP)      
      call SetTensor(statev_n, pzero, NSTAV)
      call SetTensor(eqvalues, pzero, NEQVA)
      call SetTensor(gamdot,   pzero, MAX_SLIP)
      call SetTensor(gamdotC,  pzero, MAX_SLIP) ! climb rate
      call SetTensor(crot,     pzero, DIMS*DIMS)
      call SetTensor(rrot,     pzero, DIMS*DIMS)
      call SetTensor(drot,     pzero, DIMS*DIMS)
      call SetTensor(crot_n,   pzero, DIMS*DIMS)
      call SetTensor(rrot_n,   pzero, DIMS*DIMS)
      call SetTensor(crot0,    pzero, DIMS*DIMS)
      call SetTensor(c_ijkl,   pzero, DIMS2*DIMS2)

      call SetTensor(rss,      pzero, MAX_SLIP)
      call SetTensor(rssc,     pzero, MAX_SLIP) ! climb-induced plastic strain rate

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE FetchCrystalVariablesAtIP(  &
         !gstress_n, gestran_n, gkappa_n, gstatev_n, geqvalues,  &
         gstress_n, gestran_n, gslipr_n, gbackst_n, gbackstc_n, gstatev_n, geqvalues,  &
         !gcrot_n, grrot_n, gcrot0, stress_n, estran_n, kappa_n,  &
         gcrot_n, grrot_n, gcrot0, stress_n, estran_n, slipr_n, backst_n, backstc_n, &
         statev_n, eqvalues, crot_n, rrot_n, crot0, igrn, iqpt,  &
         ielem  &
         )

      implicit  none
      include  'params_xtal.inc'

      integer igrn, iqpt, ielem

      real*8  gstress_n (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran_n (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      !real*8  gkappa_n  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gslipr_n  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackst_n  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackstc_n  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gstatev_n (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues (NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot_n   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot_n   (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot0    (DIMS, DIMS, 1, NUMQPT_T, NUMEL_T)

      !real*8  stress_n(NVECS), estran_n(NVECS), kappa_n(NKAPP)
      real*8  stress_n(NVECS), estran_n(NVECS), slipr_n(NKAPP), backst_n(NKAPP)
      real*8  statev_n(NSTAV), eqvalues(NEQVA), backstc_n(NKAPP)
      real*8  crot_n(DIMS,DIMS), rrot_n(DIMS,DIMS), crot0(DIMS,DIMS)
!
!---------------------------------------------------------------------72
!
      call EqualTensors(gstress_n(1,igrn,iqpt,ielem), stress_n, NVECS)
      call EqualTensors(gestran_n(1,igrn,iqpt,ielem), estran_n, NVECS)
      !call EqualTensors(gkappa_n (1,igrn,iqpt,ielem), kappa_n,  NKAPP)
      call EqualTensors(gslipr_n (1,igrn,iqpt,ielem), slipr_n,  NKAPP)
      call EqualTensors(gbackst_n (1,igrn,iqpt,ielem), backst_n,  NKAPP)
      call EqualTensors(gbackstc_n (1,igrn,iqpt,ielem), backstc_n,  NKAPP)      
      call EqualTensors(gstatev_n(1,igrn,iqpt,ielem), statev_n, NSTAV)

      call EqualTensors(gcrot_n(1,1,igrn,iqpt,ielem), crot_n, DIMS*DIMS)
      call EqualTensors(grrot_n(1,1,igrn,iqpt,ielem), rrot_n, DIMS*DIMS)
      call EqualTensors(gcrot0(1,1,igrn,iqpt,ielem),  crot0,  DIMS*DIMS)

      eqvalues(kEQP_n)    = geqvalues(kEQP_n,   igrn, iqpt, ielem)
      eqvalues(kMISES_n)  = geqvalues(kMISES_n, igrn, iqpt, ielem)
      eqvalues(kSHRATE_n) = geqvalues(kSHRATE_n,igrn, iqpt, ielem)
      eqvalues(kGAMTOT_n) = geqvalues(kGAMTOT_n,igrn, iqpt, ielem)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SaveCrystalVariablesAtIP(  &
         !stress, estran, kappa, statev, gamdot, eqvalues, crot, rrot,  &
         stress, estran, slipr, backst, backstc, statev, gamdot,  &
         gamdotC,  &
         eqvalues, crot, rrot,  &
         !gstress, gestran, gkappa, gstatev, geqvalues, ggamdot, gcrot,  &
         gstress, gestran, gslipr, gbackst, gbackstc, gstatev, geqvalues, ggamdot,  &
         ggamdot2,  &
         gcrot, grrot, igrn, iqpt, ielem  &
         )

      implicit  none
      include  'params_xtal.inc'

      integer igrn, iqpt, ielem

      real*8  gstress  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gestran  (NVECS, MAX_GRN, NUMQPT_T, NUMEL_T)
      !real*8  gkappa   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gslipr   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackst   (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackstc  (NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)      
      real*8  gstatev  (NSTAV, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  geqvalues(NEQVA, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot  (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  ggamdot2 (MAX_SLIP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gcrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  grrot    (DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)

      !real*8  stress(NVECS), estran(NVECS), kappa(NKAPP)
      real*8  stress(NVECS), estran(NVECS), slipr(NKAPP), backst(NKAPP), backstc(NKAPP)
      real*8  statev(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP), gamdotC(MAX_SLIP)
      real*8  crot(DIMS,DIMS), rrot(DIMS,DIMS)
!
!---------------------------------------------------------------------72
!
      call EqualTensors(stress, gstress(1,igrn,iqpt,ielem), NVECS)
      call EqualTensors(estran, gestran(1,igrn,iqpt,ielem), NVECS)
      !call EqualTensors(kappa,  gkappa (1,igrn,iqpt,ielem), NKAPP)
      call EqualTensors(slipr,  gslipr (1,igrn,iqpt,ielem), NKAPP)
      call EqualTensors(backst,  gbackst (1,igrn,iqpt,ielem), NKAPP)
      call EqualTensors(backstc, gbackstc (1,igrn,iqpt,ielem), NKAPP)      
      call EqualTensors(statev, gstatev(1,igrn,iqpt,ielem), NSTAV)
      call EqualTensors(gamdot, ggamdot(1,igrn,iqpt,ielem), MAX_SLIP)
      call EqualTensors(gamdotC, ggamdot2(1,igrn,iqpt,ielem), MAX_SLIP)
      call EqualTensors(crot, gcrot(1,1,igrn,iqpt,ielem), DIMS*DIMS)
      call EqualTensors(rrot, grrot(1,1,igrn,iqpt,ielem), DIMS*DIMS)

      geqvalues(kEQP,   igrn,iqpt,ielem) = eqvalues(kEQP)
      geqvalues(kMISES, igrn,iqpt,ielem) = eqvalues(kMISES)
      geqvalues(kSHRATE,igrn,iqpt,ielem) = eqvalues(kSHRATE)
      geqvalues(kGAMTOT,igrn,iqpt,ielem) = eqvalues(kGAMTOT)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE OutputQnts(  &
         !d_vec, s_ij, ee_ij, kappa, gamdot, rss, statev, eqvalues, crot,  &
         d_vec, s_ij, ee_ij, slipr, backst, backstc, gamdot, rss,  &
         gamdotC, rssc,  &
         statev, eqvalues, crot,  &
         rrot, c_ijkl, d_kk, dtime, time, numslip, incr,  &
         iterCounterS, iterCounterN, ielem, iqpt, igrn  &
         )
!
      use     SlipOverStress
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, incr, iterCounterS, iterCounterN
      integer ielem, iqpt, igrn
      real*8  d_kk, dtime, time(2)
      real*8  d_vec(NVECS), s_ij(DIMS,DIMS), ee_ij(DIMS,DIMS)
      !real*8  kappa(NKAPP), gamdot(MAX_SLIP), statev(NSTAV)
      real*8  slipr(NKAPP), backst(NKAPP), gamdot(MAX_SLIP), gamdotC(MAX_SLIP)
      real*8  statev(NSTAV), backstc(NKAPP)
      real*8  eqvalues(NEQVA), crot(DIMS,DIMS), rrot(DIMS,DIMS)
      real*8  c_ijkl(DIMS2,DIMS2)
      real*8  rss(MAX_SLIP), rssc(MAX_SLIP)

      integer i, j
!      character*80 message
      real*8  strain_11, stress_11, d_eff, eqstran, deqstran
      real*8  d_ij(DIMS,DIMS)

      real*8  InnerProductVec

      save    strain_11, eqstran
      data    strain_11, eqstran /0.d0, 0.d0/
!
!----- NOTE: 'strain_11' makes sense only for coaxial loading
!
!---------------------------------------------------------------------72
!
!------- write headers in output files
!
      if (incr .eq. 1 .and. PrintFlag==1) then
         write(XTAL_STRESS_O,  800) ielem, iqpt, igrn
         write(XTAL_STRAIN_O,  900) ielem, iqpt, igrn
         write(XTAL_EFFSS_O,  1000) ielem, iqpt, igrn
         write(XTAL_TRUESS_O, 2000) ielem, iqpt, igrn
         write(XTAL_ITER_O,   3000) ielem, iqpt, igrn
      endif
!
!------- equivalent total strain and strain_11
!
      d_eff = dsqrt(2./3.*InnerProductVec(d_vec, d_vec, NVECS))
      deqstran = dtime * d_eff
      eqstran  = eqstran + deqstran

      call Vec5x1ToMat3x3Symm(d_vec, d_ij, DIMS)
      strain_11 = strain_11 + dtime*(d_ij(1,1) + pthird*d_kk)
      stress_11 = s_ij(1,1)
!
!------- write stress / consistent tangent
!
      if (incr .eq. 1 .and. PrintFlag==1) then
          write(XTAL_STRESS_O, 4000) incr, ((s_ij(i,j), j=i,3), i=1,3)
          write(XTAL_STRESS_O, 4500) ((c_ijkl(i,j), j=1,DIMS2), i=1,DIMS2)

    !      write(XTAL_STRESS_O, 4500) (rss(i)/(kappa(1)+overstress(i)),
    !     &                            i=1,numslip)
          !write(XTAL_STRESS_O, 4500) (rss(i)/(kappa(i)), i=1,numslip)
    !
    !------- write elastic strains / crot / rotation tensors; gammadot
    !
          write(XTAL_STRAIN_O, '(/i5)') incr
          do i = 1, 3
             write(XTAL_STRAIN_O, 5000) (ee_ij(i,j), j=1,3),  &
                                 (crot(i,j),j=1,3), (rrot(i,j), j=1,3)
          enddo

          write(XTAL_STRAIN_O, 5200) eqvalues(kEQP), statev(kWPNORM),  &
                                     int( statev(kSSACT) )
          write(XTAL_STRAIN_O, 5500) (gamdot(i), i=1,numslip)
    !
    !------- write effective quantities
    !
          write(XTAL_EFFSS_O, 5800) incr, int(statev(kSSACT)), dtime,  &
                          time+dtime, d_eff, deqstran,  &
                          statev(kWPNORM), eqvalues(kEQP),  &
                          (eqvalues(kEQP)-eqvalues(kEQP_n))/dtime,  &
                          eqvalues(kMISES)
    !
    !------- write true stress-strain curve (X-direction) (uniaxial/MPS)
    !
          write(XTAL_TRUESS_O, 6000) incr, time, strain_11, stress_11
    !
    !------- write iteration counters and escalar variables
    !
          write(XTAL_ITER_O, 7000) incr, iterCounterS, iterCounterN,  &
                            statev(kEVOL), statev(kPRES), statev(kDETVe)!,  &
                            !(kappa(i),i=1,numslip)
    !     &                  (kappa(1)+overstress(i),i=1,numslip)
      endif
!
!------- formats
!
 800  format(' INCR', 15x, 'CAUCHY STRESS' ,/,  &
             20x,'CONSISTENT TANGENT ',/, 20x, 'RSS/KAPPA ',/,  &
             ' (elem # ', i5, ',  qpt # ', i2, ',  grn # ', i5, ')'/) 
 900  format(' INCR',/,18x,' EE_ij',32x,'Crot_ij',30x,'Rrot_ij',/,  &
             19x, 'EQPS', 33x, 'WP_NORM', 30x, 'SS_ACTIV',/,  &
             48x, ' GAMMADOT(1,...,NUMSLIP)',/,38x,  &
             ' (elem # ', i5, ',  qpt # ', i2, ',  grn # ', i5, ')') 
1000  format(' INCR',3x,'SS_ACT',5x,'DTIME',8x,'TIME',8x,'D_EFF',10x,  &
            'DEQSTRAN',8x,'WP_NORM',7x,'EQP',10x,'DP_EFF',8x,'MISES',/,  &
             ' (elem # ', i5, ',  qpt # ', i2, ',  grn # ', i5, ')')
2000  format(' INCR',7x,'TIME',9x,'E_11',9x,'SIGMA_11',/,  &
             ' (elem # ', i5, ',  qpt # ', i2, ',  grn # ', i5, ')',/,  &
             ' (for uniaxial loading along X-direction / MPS runs)' ) 
3000  format(' INCR',1x,'IT-S',1x,'IT-N',6x,'EVOL', 8x,'PRESS',  &
             8x,'DETVe',6x,'KAPPA(1) ... KAPPA(NKAPP)'/,  &
             ' (elem # ', i5, ',  qpt # ', i2, ',  grn # ', i5, ')') 
4000  format(i8,2x,3(1x,e11.4)/22x,2(1x,e11.4)/34x,1(1x,e11.4)/)
4500  format(6(10x,6(1x,e11.4)/))
5000  format((3(2x,3(1x,e11.4))))
5200  format(/(15x,e11.4,27x,e11.4,28x,i5))
5500  format(/(22x,6(1x,e11.4)))
5800  format(i8,2x,i3,8(2x,e12.5))
6000  format(i8,8(2x,e12.5))
7000  format(i8,1x,i3,2x,i3,1x,3(2x,e11.4),18(1x,e11.4))

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE MonitorCrystalIntegration(  &
         !d_vec, w_vec, d_kk, s_ij, ee_ij, tau, stress, estran, kappa,  &
         d_vec, w_vec, d_kk, s_ij, ee_ij, tau,  &
         stress, estran, slipr, backst, backstc, statev, eqvalues,  &
         gamdot, rss, gamdotC, rssc,  &
         crot, rrot,  &
         drot, stress_n, estran_n,  &
         !kappa_n, statev_n, crot_n, rrot_n, crot0, fCeDev, fCeiDev,  &
         slipr_n, backst_n, backstc_n, statev_n, crot_n, rrot_n, crot0, fCeDev, fCeiDev,  &
         matProp, sigfs, tauSlip,  &
         zBar0, pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0,  &
         KBar0, KpBar0, KqBar0, KpBar0Vec, KpBar0devVec, KqBar0Vec, KpKpTBar0,  &
         numslip, numvtx, maxIterState, maxIterNewt,  &
         tolerState, tolerNewt, dtime, time, fCeDevVol, fCeVol,  &
         theta, thetao, igrn, iqpt,  &
         ielem, incr, kstep, iterCounterS, iterCounterN, ierr, cepmod, hardmtx, &
         sGamPPt, sGamKpKpt, grainid  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, numvtx, maxIterState, maxIterNewt, grainid
      integer igrn, iqpt, ielem, incr, iterCounterS, iterCounterN, ierr, kstep
      real*8  tolerState, tolerNewt, dtime, time, theta, thetao, d_kk
      real*8  fCeDevVol(NVECS), fCeVol
      real*8  d_vec(NVECS), w_vec(DIMS)
      real*8  s_ij(DIMS,DIMS), ee_ij(DIMS,DIMS), tau(NVECS)
      !real*8  stress(NVECS), estran(NVECS), kappa(NKAPP)
      real*8  stress(NVECS), estran(NVECS), slipr(NKAPP), backst(NKAPP), backstc(NKAPP)
      real*8  statev(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP), gamdotC(MAX_SLIP)
      real*8  crot(DIMS,DIMS), rrot(DIMS,DIMS), drot(DIMS,DIMS)
      !real*8  stress_n(NVECS), estran_n(NVECS), kappa_n(NKAPP)
      real*8  stress_n(NVECS), estran_n(NVECS), slipr_n(NKAPP), backst_n(NKAPP), backstc_n(NKAPP)
      real*8  statev_n(NSTAV), crot_n(DIMS,DIMS), rrot_n(DIMS,DIMS)
      real*8  crot0(DIMS,DIMS)
      real*8  fCeDev(NVECS,NVECS), fCeiDev(NVECS,NVECS)
      real*8  matProp(NPROPS, MAX_SLIP)
      real*8  sigfs(NVECS,MAX_VTX), tauSlip(MAX_SLIP)
      real*8  zBar0(DIMS,DIMS,MAX_SLIP)
      real*8  pBar0(DIMS,DIMS,MAX_SLIP), qBar0(DIMS,DIMS,MAX_SLIP)
      real*8  pBar0Vec(NVECS,MAX_SLIP), qBar0Vec(DIMS,MAX_SLIP)
      real*8  ppTBar0(NVECS,NVECS,MAX_SLIP)
      real*8  KBar0(DIMS,DIMS,MAX_SLIP)
      real*8  KpBar0(DIMS,DIMS,MAX_SLIP), KqBar0(DIMS,DIMS,MAX_SLIP)
      real*8  KpBar0Vec(NVECS,MAX_SLIP), KpBar0devVec(NVECS,MAX_SLIP), KqBar0Vec(DIMS,MAX_SLIP)
      real*8  KpKpTBar0(NVECS,NVECS,MAX_SLIP)
      real*8  cepmod(DIMS2,DIMS2)
      real*8  rss(MAX_SLIP), rssc(MAX_SLIP)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      integer subIncr, totSubIncrs
      real*8  dtime_tr, theta_tr, tmp, d_kk_n, d_kk_tr
      real*8  d_vec_tr(NVECS), w_vec_tr(DIMS)
      real*8  d_vec_n(NVECS), w_vec_n(DIMS), tau_save(NSTAV)
      real*8  sGamPPt(NVECS,NVECS), sGamKpKpt(NVECS,NVECS)
!
!---------------------------------------------------------------------72
!
!------- counters for sub-increments
!
      subIncr = 1
      totSubIncrs = 1
!  
!------- integrate crystal constitutive equations
!
      !write(*,*) '10noel:', ielem, '10grainid', grainid
      !write(*,*) 'moncryint 1'
      iterCounterS = 0
      call IntegrateCrystalEqns(d_vec, w_vec, d_kk, s_ij, ee_ij, tau,  &
         !stress, estran, kappa, statev, eqvalues, gamdot, rss,  &
         stress, estran, slipr, backst, backstc, statev, eqvalues,  &
         gamdot, rss, gamdotC, rssc,  &
         crot, rrot,  &
         !drot, stress_n, estran_n, kappa_n, statev_n, crot_n, rrot_n,  &
         drot, stress_n, estran_n, slipr_n, backst_n, backstc_n, statev_n, crot_n, rrot_n,  &
         crot0, fCeDev, fCeiDev, matProp, sigfs, tauSlip,  &
         zBar0, pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0,  &
         KBar0, KpBar0, KqBar0, KpBar0Vec, KpBar0devVec, KqBar0Vec, KpKpTBar0,  &
         numslip, numvtx, maxIterState, maxIterNewt,  &
         tolerState, tolerNewt, dtime,  time, &
         fCeDevVol, fCeVol,  &
         theta, thetao, iterCounterS, iterCounterN, incr, kstep, ierr,  &
         cepmod, subIncr, totSubIncrs, hardmtx, igrn, sGamPPt, sGamKpKpt)

      !write(*,*) '11noel:', ielem, '11grainid', grainid
!
!------- if converged -> return, else do -> subincrementation
!------- NOTE: here, subincrementation improves the initial guess 
!-------       for the solution variables 'stress'
!
      if (ierr .eq. XTAL_CONVERGED) return

      call WriteMessage  &
           (XTAL_O, 'MonitorCrystalIntegration: using Sub-Incrs!')
!
!------- compute velocity gradient at t_n
!
      !write(*,*) 'moncryint 2'
      call DeformationRate_n(d_vec_n, w_vec_n, d_kk_n, iqpt)
!     MAY NEED Vec6 to Vec5 routine.........
!
!------- loop for sub-incrementation procedure
!
      do while (.true.)
!
!---------- if not converged, increase # of sub-increments 
         if (ierr .ne. XTAL_CONVERGED) then
            subIncr = 2 * subIncr - 1
            totSubIncrs = 2 * totSubIncrs
            if (totSubIncrs .gt. 128) then
               call WriteMessage(XTAL_O,  &
                        'MonitorCrystalIntegration: totSubIncrs > 128')
               return
            endif
            call SetTensor(tau, pzero, NVECS)
            if (subIncr .gt. 1)  &
                 call EqualTensors(tau_save, tau, NVECS)
!
!---------- if converged, adjust subincrements
         else if (subIncr .lt. totSubIncrs) then
            if ( (subIncr/2*2) .eq. subIncr) then
               subIncr = subIncr / 2 + 1
               totSubIncrs = totSubIncrs / 2
            else
               subIncr = subIncr + 1
            endif
!
!---------- successful return for subincrementation
         else
            call WriteMessage(XTAL_O, 'Sub-Incrs successful !!!')
            return
         endif
!
!---------- report sub-incrs status
         write(XTAL_O, 1000) kstep, incr, igrn, iqpt, ielem, subIncr, totSubIncrs
!
!---------- trial time, velocity gradient and temperature (assumes 
!---------- linear variation of velocity gradient and temperature in dt)
         tmp = real(subIncr) / real(totSubIncrs)   
         dtime_tr = dtime * tmp
         theta_tr = (1.-tmp)*thetao + tmp*theta
         d_kk_tr  = (1.-tmp)*d_kk_n + tmp*d_kk

         call AddTensors((1.-tmp), d_vec_n, tmp, d_vec, d_vec_tr, NVECS)
         call AddTensors((1.-tmp), w_vec_n, tmp, w_vec, w_vec_tr, DIMS)
!
!---------- save current convergent solution before getting next solution
         if (subIncr .gt. 1 .and. ierr .eq. XTAL_CONVERGED)  &
                       call EqualTensors(tau, tau_save, NVECS)
!
!---------- integrate crystal constitutive equations
         !write(*,*) 'moncryint 3'
         iterCounterS = 0
         call IntegrateCrystalEqns(d_vec_tr, w_vec_tr, d_kk_tr, s_ij,  &
           !ee_ij, tau, stress, estran, kappa, statev, eqvalues, gamdot,  &
           ee_ij, tau, stress, estran, slipr, backst, backstc, statev, eqvalues,  &
           gamdot, rss, gamdotC, rssc,  &
           !crot, rrot, drot, stress_n, estran_n, kappa_n, statev_n,  &
           crot, rrot, drot, stress_n, estran_n, slipr_n, backst_n, backstc_n, statev_n,  &
           crot_n, rrot_n, crot0, fCeDev, fCeiDev, matProp, sigfs, tauSlip,  &
           zBar0, pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0,  &
           KBar0, KpBar0, KqBar0, KpBar0Vec, KpBar0devVec, KqBar0Vec, KpKpTBar0, &
           numslip, numvtx, maxIterState, maxIterNewt, tolerState,  &
           tolerNewt, dtime_tr, time, fCeDevVol, fCeVol,  &
           theta_tr, thetao, iterCounterS, iterCounterN, incr, kstep,  &
           ierr, cepmod, subIncr, totSubIncrs, hardmtx, igrn, sGamPPt, sGamKpKpt)

         !write(*,*) 'moncryint 4'
      enddo

1000  format(3x, 'kstep #',i5, ';', 1x,      & 
                 'kinc #',i5, ';', 1x,       &
                 'at grain #', i5, ';', 1x,  &
                 'iqpt #', i3, ';', 1x,  &
                 'elem #', i5, ';', 1x,  &
                 'subInc/totSubIncrs = ', i3, '/', i3)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE WriteTexture(  &
         !gcrot, gkappa, d_vec, numgrn, incr, iqpt, ielem, npt, noel  &
         gcrot, gslipr, gbackst, gbackstc, d_vec, numgrn, kstep, incr, iqpt, ielem, npt, noel  &
         )

      implicit none
      include  'params_xtal.inc'

      integer numgrn, incr, iqpt, ielem, npt, noel, kstep
      real*8  gcrot(DIMS, DIMS, MAX_GRN, NUMQPT_T, NUMEL_T)
      !real*8  gkappa(NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gslipr(NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackst(NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  gbackstc(NKAPP, MAX_GRN, NUMQPT_T, NUMEL_T)
      real*8  d_vec(NVECS)

      integer igrn, i
      real*8  pi, pi180, d_eff
      real*8  angle(DIMS)

      real*8  InnerProductVec
!
!---------------------------------------------------------------------72
!
      pi = 4.0 * datan(1.0d+00)
      pi180 = pi/180.
!
!-------  effective strain rate at the ip/elem
!
      d_eff = dsqrt(2./3.*InnerProductVec(d_vec, d_vec, NVECS))
!
!------- output heading
!
      if (noel .eq. 1 .and. npt .eq. 1)  &
                  write(XTAL_TXT_OUT, 1000) kstep, incr
!
!------- output orientation
!

      do igrn = 1, numgrn

         call RotMatrixToAngles(gcrot(1,1,igrn,iqpt,ielem), angle, DIMS)

         write(XTAL_TXT_OUT, 2000) (angle(i)/pi180, i=1,3), igrn,  &
                                !npt, noel, gkappa(1,igrn,iqpt,ielem), &
                                npt, noel, gslipr(1,igrn,iqpt,ielem), &
                                gbackst(1,igrn,iqpt,ielem), &
                                d_eff

      enddo
!
!------- formats
!
1000  format(/'*-----   Euler Angles at step'  i8, '    incr ', i8, ' -----*'/  &
              4x,'ang1',4x,'ang2',4x,'ang3',7x,' igrn',7x,'intpt',  &
              !3x,'ielem ',4x,'kappa(1)',4x,'d_eff')
              7x,'ielem ',4x,'slipr(1)',6x,'backst(1)',6x,'d_eff')
2000  format( 3f8.2, 3(3x,i8), 3x, e12.5, 4x, e12.5, 4x, e12.2)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE WriteStressStrainCurve(  &
         sdev_ij, d_vec, wpnorm_avg, dtime, time, numgrn,  &
         incr, npt, noel, numel_aba, numqpt_aba  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numgrn, incr, npt, noel
      real*8  dtime, time, wpnorm_avg
      real*8  sdev_ij(DIMS,DIMS), d_vec(NVECS)

      integer denom
      real*8  s_eff, d_eff, eqstran, deqstran
      real*8  sdev_vec(NVECS)
      real*8  d_eff_avg, s_eff_avg

      real*8  InnerProductVec

      save    eqstran
      data    eqstran /0.d0/

      save    d_eff_avg, s_eff_avg, deqstran

      integer numel_aba, numqpt_aba
!
!---------------------------------------------------------------------72
!
!------- write headers in output file
!
      if (noel .eq. 1 .and. npt .eq. 1) then
         s_eff_avg = 0.0d0
         d_eff_avg = 0.0d0
         deqstran  = 0.0d0
         if (incr .eq. 1 .and. PrintFlag==1)  &
              write(AGG_EFFSS_O, 800) numel_aba, numqpt_aba, numgrn
      endif
!
!------- effective stress
!
      call Mat3x3ToVec5x1Symm(sdev_ij, sdev_vec, DIMS)
      s_eff = dsqrt(3./2.*InnerProductVec(sdev_vec, sdev_vec, NVECS))

      s_eff_avg = s_eff_avg + s_eff
!
!-------  effective strain rate and total equivalent strain
!
      d_eff = dsqrt(2./3.*InnerProductVec(d_vec, d_vec, NVECS))
      d_eff_avg = d_eff_avg + d_eff

      deqstran = deqstran + dtime * d_eff
!
!------- write stress/strain curve
!
      if (noel .eq. numel_aba .and. npt .eq. numqpt_aba .and. PrintFlag==1) then
!         denom = numgrn * numel_aba * numqpt_aba
         denom = numel_aba * numqpt_aba
         s_eff_avg = s_eff_avg / denom
         d_eff_avg = d_eff_avg / denom
         deqstran = deqstran / denom
         eqstran  = eqstran + deqstran
         write(AGG_EFFSS_O, 1000) incr, dtime, time+dtime, d_eff_avg,  &
                 deqstran, eqstran, wpnorm_avg, s_eff_avg
      endif
!
!------- formats
!
800   format(' INCR',9x,'DTIME',9x,'TIME',9x,'D_EFF',8x,'DEQSTRAN',  &
             8x,'EQSTRAN',8x,'WP_NORM',6x,'S_EFF'/,  &
             ' (nelem: ', i5, ',  nqpts: ', i2, ',  ngrns: ', i5, ')') 
1000  format(i8,10(2x,e12.5))

      return
      END
!
!=====================================================================72
!
      SUBROUTINE ResetCrystalQnts(  &
         !stress, estran, kappa, statev, eqvalues, gamdot, crot, rrot,  &
         stress, estran, slipr, backst, backstc, statev,  &
         eqvalues, gamdot, gamdotC,  & 
         !s_ij, c_ijkl, stress_n, estran_n, kappa_n, statev_n, crot_n,  &
         crot, rrot, s_ij, c_ijkl,  &
         stress_n, estran_n, slipr_n, backst_n, backstc_n, &
         statev_n, crot_n, rrot_n, matProp,  &
         pBar0Vec, KpBar0Vec, KpBar0devVec,  & 
         fCeiDev, fCeDevVol, fCeVol, numslip, dtime,  &
         sGamPPt, sGamKpKpt  &
         )
!
      use     SlipOverStress
      implicit  none
      include  'params_xtal.inc'
      include  'numbers.inc'

      integer numslip
      real*8  dtime

      !real*8  stress(NVECS), estran(NVECS), kappa(NKAPP)
      real*8  stress(NVECS), estran(NVECS), slipr(NKAPP), backst(NKAPP), backstc(NKAPP)
      real*8  statev(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP), gamdotC(MAX_SLIP)
      real*8  crot(DIMS,DIMS), rrot(DIMS,DIMS), backst_n(NKAPP), backstc_n(NKAPP)
      real*8  s_ij(DIMS,DIMS), c_ijkl(DIMS2,DIMS2)

      !real*8  stress_n(NVECS), estran_n(NVECS), kappa_n(NKAPP)
      real*8  stress_n(NVECS), estran_n(NVECS), slipr_n(NKAPP)
      real*8  statev_n(NSTAV), crot_n(DIMS,DIMS), rrot_n(DIMS,DIMS)

      real*8  matProp(NPROPS, MAX_SLIP)
      real*8  pBar0Vec(NVECS,MAX_SLIP), KpBar0Vec(NVECS,MAX_SLIP), KpBar0devVec(NVECS,MAX_SLIP)
      real*8  fCeiDev(NVECS,NVECS), fCeDevVol(NVECS), fCeVol
!
!---- Local variables
!
      integer is
      real*8  qr5x5(NVECS,NVECS)
      real*8  pHatVec(NVECS,MAX_SLIP), ppTHat(NVECS,NVECS,MAX_SLIP)
      real*8  KpHatVec(NVECS,MAX_SLIP), KpHatdevVec(NVECS,MAX_SLIP), KpKpTHat(NVECS,NVECS,MAX_SLIP)      
      real*8  fCeiDevHat(NVECS,NVECS),fCeDevVolHat(NVECS)

      real*8  crss, bss, bssc
      real*8  tau(NVECS)
      real*8  rss(MAX_SLIP), dGamdTau(MAX_SLIP)
      real*8  rssc(MAX_SLIP), dGamCdTauC(MAX_SLIP)      

      real*8  InnerProductVec, SSKineticEqn2, SSKineticEqn3
!
!
      real*8  sGamPPt(NVECS,NVECS), sGamKpKpt(NVECS,NVECS)
!
!---------------------------------------------------------------------72
!
!---- Reset crystal quantities
!
      call EqualTensors(stress_n, stress, NVECS)
      call EqualTensors(estran_n, estran, NVECS)
      !call EqualTensors(kappa_n , kappa , NKAPP)
      call EqualTensors(slipr_n , slipr , NKAPP)
      call EqualTensors(backst_n , backst , NKAPP)
      call EqualTensors(backstc_n, backstc, NKAPP)
      call EqualTensors(statev_n, statev, NSTAV)
      call EqualTensors(crot_n  , crot  , DIMS*DIMS)
      call EqualTensors(rrot_n  , rrot  , DIMS*DIMS)
                                                                                    
      eqvalues(kEQP)    = eqvalues(kEQP_n)
      eqvalues(kMISES)  = eqvalues(kMISES_n)             
      eqvalues(kSHRATE) = eqvalues(kSHRATE_n)
      eqvalues(kGAMTOT) = eqvalues(kGAMTOT_n)

      call SetToScaledtensor(statev(kDETVe), stress, tau, NVECS)
!
!------- Cauchy stress (tensor form)
      call Vec5x1ToMat3x3Symm(stress, s_ij, DIMS)
      call AddScaledTensor(statev(kPRES), Ident2nd, s_ij, DIMS*DIMS)
!
!---- Compute Approximate Consistent Tangent
!
!------- Rotation tensor 
      call RotMat5x5ForSymm(crot, qr5x5, DIMS)
!
!------- Compute gamdot, gamdotC, sGamPPt and sGamKpKpt
      call SetTensor(sGamPPt, pzero, NVECS*NVECS)
      call SetTensor(sGamKpKpt, pzero, NVECS*NVECS)      
      do is = 1, numslip
!
!------- Rotate symmetric part of Schimdt tensor
         call MultAxu(qr5x5, pBar0Vec(1,is), pHatVec(1,is), NVECS)
         call OuterProductVec(pHatVec(1,is), pHatVec(1,is), ppTHat(1,1,is), NVECS)
		 !write(*,*) 'rescryqnt 1'
!
!------- Resolve shear stresses
         rss(is) = InnerProductVec(tau, pHatVec(1,is), NVECS)
!         crss = kappa(1) + overstress(is)
         !crss = kappa(is)
         crss = slipr(is)
         bss = backst(is)
         bssc = backstc(is)
!
!------- Shear strain rate
         gamdot(is) = SSKineticEqn2(rss(is), crss, bss, matProp(1,is), kGAMDOT)
!
!------- Derivative of shear strain rate
         dGamdTau(is) = SSKineticEqn2(rss(is), crss, bss, matProp(1,is),  &
                                         kdGAMdTAU)
         !write(*,*) 'rescryqnt 2'
         call AddScaledTensor(dtime*dGamdTau(is), ppTHat(1,1,is),  &
                              sGamPPt, NVECS*NVECS)
!
!		 
!
!------- Rotate symmetric part of CLIMB tensor
         call MultAxu(qr5x5, KpBar0Vec(1,is), KpHatVec(1,is), NVECS)
         call MultAxu(qr5x5, KpBar0devVec(1,is), KpHatdevVec(1,is), NVECS)
         call OuterProductVec(KpHatVec(1,is), KpHatVec(1,is),  &
                              KpKpTHat(1,1,is), NVECS)
!
!------- Resolve shear stresses of CLIMB
         rssc(is) = InnerProductVec(tau, KpHatdevVec(1,is), NVECS)
		 !write(*,*) "rssc=", rssc(is)
!
!------- Shear strain rate of CLIMB
         gamdotC(is) = SSKineticEqn3(rssc(is), bssc, matProp(1,is), kGAMDOTc)
!
!------- Derivative of shear strain rate of CLIMB
         dGamCdTauC(is) = SSKineticEqn3(rssc(is), bssc, matProp(1,is), kdGAMcdTAUc)
         !write(*,*) 'rescryqnt 2'
         call AddScaledTensor(dtime*dGamCdTauC(is), KpKpTHat(1,1,is),  &
                              sGamKpKpt, NVECS*NVECS)

      enddo
!
!------- Rotate elasticities
      call MultQAQT(qr5x5, fCeiDev, fCeiDevHat, NVECS)
      call MultAxu(qr5x5, fCeDevVol, fCeDevVolHat, NVECS)
      !write(*,*) 'rescryqnt 3'
!
!------- Approximate consistent tangent
      call PlasticModuli(c_ijkl, fCeiDevHat, fCeDevVolHat, fCeVol,  &
                         statev(kDETVe), sGamPPt+sGamKpKpt)
      !write(*,*) 'rescryqnt 4'

      return
      end
!
!=====================================================================72
!
!=====================================================================72
!
      SUBROUTINE IntegrateCrystalEqns(  &
         d_vec, w_vec, d_kk, s_ij, ee_ij, tau, stress,  &
         !estran, kappa, statev, eqvalues, gamdot, rss, crot, rrot, drot,  &
         estran, slipr, backst, backstc, statev, eqvalues,  &
         gamdot, rss, gamdotC, rssc,  &
         crot, rrot, drot,  &
         !stress_n, estran_n, kappa_n, statev_n, crot_n, rrot_n,  &
         stress_n, estran_n, slipr_n, backst_n, backstc_n, statev_n, crot_n, rrot_n,  &
         crot0, fCeDev, fCeiDev, matProp, sigfs, tauSlip,  &
         zBar0, pBar0, qBar0, pBar0Vec, qBar0Vec, ppTBar0,  &
         KBar0, KpBar0, KqBar0, KpBar0Vec, KpBar0devVec, KqBar0Vec, KpKpTBar0,  &
         numslip, numvtx, maxIterState, maxIterNewt,  &
         tolerState, tolerNewt, dtime, time, &
         fCeDevVol, fCeVol,  &
         theta, thetao, iterCounterS, iterCounterN, globalIncr, kstep, ierr,  &
         cepmod, subIncr, totSubIncrs, hardmtx, igrn, sGamPPt, sGamKpKpt)

      implicit none
      include 'params_xtal.inc' 
      include 'numbers.inc'

      integer numslip, numvtx, maxIterState, maxIterNewt, globalIncr, kstep
      integer iterCounterS, iterCounterN, ierr, subIncr, totSubIncrs
      real*8  tolerState, tolerNewt, dtime, time, theta, thetao, d_kk
      real*8  fCeDevVol(NVECS), fCeVol
      real*8  d_vec(NVECS), w_vec(DIMS)
      real*8  s_ij(DIMS,DIMS), ee_ij(DIMS,DIMS), tau(NVECS)
      !real*8  stress(NVECS), estran(NVECS), kappa(NKAPP)
      real*8  stress(NVECS), estran(NVECS), slipr(NKAPP), backst(NKAPP), backstc(NKAPP)
      real*8  statev(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP), gamdotC(MAX_SLIP)
      real*8  crot(DIMS,DIMS), rrot(DIMS,DIMS), drot(DIMS,DIMS)
      !real*8  stress_n(NVECS), estran_n(NVECS), kappa_n(NKAPP)
      real*8  stress_n(NVECS), estran_n(NVECS), slipr_n(NKAPP)
      real*8  statev_n(NSTAV), crot_n(DIMS,DIMS), rrot_n(DIMS,DIMS)
      real*8  crot0(DIMS,DIMS), backst_n(NKAPP), backstc_n(NKAPP)
      real*8  fCeDev(NVECS,NVECS), fCeiDev(NVECS,NVECS)
      real*8  matProp(NPROPS,MAX_SLIP)
      real*8  sigfs(NVECS,MAX_VTX), tauSlip(MAX_SLIP)
      real*8  zBar0(DIMS,DIMS,MAX_SLIP)
      real*8  pBar0(DIMS,DIMS,MAX_SLIP), qBar0(DIMS,DIMS,MAX_SLIP)
      real*8  pBar0Vec(NVECS,MAX_SLIP), qBar0Vec(DIMS,MAX_SLIP)
      real*8  KBar0(DIMS,DIMS,MAX_SLIP)
      real*8  KpBar0(DIMS,DIMS,MAX_SLIP), KqBar0(DIMS,DIMS,MAX_SLIP)
      real*8  KpBar0Vec(NVECS,MAX_SLIP), KqBar0Vec(DIMS,MAX_SLIP)
      real*8  KpBar0devVec(NVECS,MAX_SLIP)
      real*8  ppTBar0(NVECS,NVECS,MAX_SLIP), KpKpTBar0(NVECS,NVECS,MAX_SLIP)
      real*8  cepmod(DIMS2,DIMS2)
      real*8  rss(MAX_SLIP), rssc(MAX_SLIP)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      integer is, iterState, igrn, i, j
      logical converged
      !real*8  e_kk, norm_tau0, norm_kapp0, norm_tau, norm_kapp, epsdot
      real*8  e_kk, norm_tau0, norm_slipr0, norm_backst0, norm_tau, &
      norm_slipr, norm_backst, epsdot, norm_backst0c, norm_backstc
      real*8  d_vec_lat(NVECS), tau_lat(NVECS)
      real*8  qr5x5(NVECS,NVECS), qr3x3(DIMS,DIMS)
      real*8  pHatVec(NVECS,MAX_SLIP), qHatVec(DIMS,MAX_SLIP)
      real*8  KpHatVec(NVECS,MAX_SLIP), KqHatVec(DIMS,MAX_SLIP)
      real*8  KpHatdevVec(NVECS,MAX_SLIP)      
      real*8  ppTHat(NVECS,NVECS,MAX_SLIP), KpKpTHat(NVECS,NVECS,MAX_SLIP)
      real*8  fCeiDevHat(NVECS,NVECS),fCeDevVolHat(NVECS)

      logical ConvergeState
      real*8  InnerProductVec
      real*8  sGamPPt(NVECS,NVECS), sGamKpKpt(NVECS,NVECS)
!
!----------------------------------------------------------------------
!------------------------- DEVIATORIC RESPONSE
!
!------- effective deformation rate 
!      write(*,*) 'flag1'
!
      epsdot = dsqrt(2./3.*InnerProductVec(d_vec, d_vec, NVECS))
      epsdot = max(epsdot, TINY)
!
!
      call RotMat5x5ForSymm(crot_n, qr5x5, DIMS)
      call RotMat3x3ForSkew(crot_n, qr3x3, DIMS)
!
!------- rotate Scmidt tensors (vectors) to Btilde_n configuration
!
      do is = 1, numslip
         call MultAxu(qr5x5, pBar0Vec(1,is), pHatVec(1,is), NVECS)
         call MultAxu(qr3x3, qBar0Vec(1,is), qHatVec(1,is), DIMS)
      enddo
!
!------- rotate climb tensors (vectors) to Btilde_n configuration
!
	  !write(*,*) 'Check input: KqBarVec, OK'
	  !do i = 1, numslip
      	!write(*,*) (KpBar0Vec(j,i),j=1,NVECS)
        !write(*,*) (KqBar0Vec(j,i),j=1,DIMS)
      !enddo
      do is = 1, numslip
         call MultAxu(qr5x5, KpBar0Vec(1,is), KpHatVec(1,is), NVECS)
         call MultAxu(qr5x5, KpBar0devVec(1,is), KpHatdevVec(1,is), NVECS)
         call MultAxu(qr3x3, KqBar0Vec(1,is), KqHatVec(1,is), DIMS)
      enddo
      !write(*,*) 'Check input in ICEqns: KpBar0devVec, OK'
	  !do i = 1, numslip
      	!write(*,*) (KpBar0devVec(j,i),j=1,NVECS)
      !enddo
!
!
!------- initialze slip hardening at time t
!
      !call EqualTensors(kappa_n, kappa, NKAPP)
      call EqualTensors(slipr_n, slipr, NKAPP)
      call EqualTensors(backst_n, backst, NKAPP)
      call EqualTensors(backstc_n, backstc, NKAPP)
!
!------- initial estimate for the stress
!
      if (subIncr .eq. 1) then

         if (globalIncr .eq. 1 .and. time == 0.0) then
               !write (*,*) 'Initial guess of stress set to zero'
!
!------------- initial guess from viscoplastic solution
!---------------- transform D from sample to crystal coordinates
            !call MultATxu(qr5x5, d_vec, d_vec_lat, NVECS)
!
!---------------- viscoplastic solution in crystal coordinates
            !call StressSolveViscoPlastic(tau_lat, d_vec_lat, kappa,  &
            !write(*,*) 'intcryeqn 1'
            !call StressSolveViscoPlastic(tau_lat, d_vec_lat, slipr, backst, &
                    !pBar0Vec, ppTBar0, matProp, sigfs, tauSlip, dtime,  &
                    !theta, thetao, epsdot, tolerNewt, maxIterNewt,  &
                    !numslip, numvtx)
!
!---------------- transform computed stresses from crystal to sample axis
            !call MultAxu(qr5x5, tau_lat, tau, NVECS)

            !New Initialization - Robby
            !write(*,*) 'intcryeqn 1'
            call SetTensor(tau, pzero, NVECS)
         else
!
!------------- initial guess from previous solution
            call EqualTensors(stress_n, tau, NVECS)
         endif

      endif

!
!------- initial estimate for hardness and rotation tensor
!
      !call HardeningRotationSolve(kappa, kappa_n, eqvalues, rrot,  &
      !write(*,*) 'intcryeqn 2'
      call HardeningRotationSolve(slipr, backst, backstc, slipr_n,  &
      		  backst_n, backstc_n, eqvalues, rrot,  &
              rrot_n, crot, crot0, drot, tau, w_vec,  &
              gamdot, rss, gamdotC, rssc,  &
              pHatVec, qHatVec, KpHatdevVec, KqHatVec,  &
              matProp, hardmtx, epsdot, dtime, numslip,  &
              kHARD_EXPL)
!
!------- initial valuies for the norm of stress and kappa
!
      norm_tau0  = dsqrt(InnerProductVec(tau, tau, NVECS))
      !norm_kapp0 = dsqrt(InnerProductVec(kappa, kappa, NKAPP))
      norm_slipr0 = dsqrt(InnerProductVec(slipr, slipr, NKAPP))
      norm_backst0 = dsqrt(InnerProductVec(backst, backst, NKAPP))
      norm_backst0c = dsqrt(InnerProductVec(backstc, backstc, NKAPP))
!
!------- initialized global flag to monitor Newton/State convergence
!
      ierr = XTAL_CONVERGED
!
!------- predictor volumetric elastic strain in Btilde
!
      e_kk = statev_n(kEVOL) + dtime * d_kk
!
!------- iterate for the material state
!
      iterState    = 0
      iterCounterN = 0
      converged = .false.
      !write(*,*) 'intcryeqn 3'

      do while(iterState .lt. maxIterState  .and. .not.converged)

         iterState = iterState + 1
!
!---------- rotate Scmidt tensors to Btilde configuration
!
         call RotateSlipGeometry(crot, pBar0Vec, qBar0Vec,  &
         					     qr5x5, qr3x3,  &
                                 pHatVec, qHatVec, ppTHat,  &
                                 numslip)
!
!---------- rotate CLIMB tensors to Btilde configuration
!
		 !write(*,*) 'Check input: KpHatdevVec, OK'
		 !do i = 1, numslip
      	 !	write(*,*) (KpHatdevVec(j,i),j=1,NVECS)
         !enddo
         call RotateSlipGeometry2(crot, KpBar0Vec, KpBar0devVec, KqBar0Vec,  &
         						  qr5x5, qr3x3,  &
         						  KpHatVec, KpHatdevVec, KqHatVec, KpKpTHat, &
         						  numslip)
	     !write(*,*) 'Check output: KpHatdevVec, OK'
		 !do i = 1, numslip
      	  ! write(*,*) (KpHatdevVec(j,i),j=1,NVECS)
         !enddo
!         
!
!---------- rotate deviatoric & dev/volumetric elasticities to Btilde:
!
         call MultQAQT(qr5x5, fCeiDev, fCeiDevHat, NVECS)
         call MultAxu(qr5x5, fCeDevVol, fCeDevVolHat, NVECS)
!
!---------- solve for the crystal stresses
!
         !write(*,*) 'intcryeqn 4'
         !call StressSolveDeviatoric(tau, d_vec, estran, estran_n, kappa,  &
         call StressSolveDeviatoric(tau, d_vec, estran, estran_n, slipr,  &
         	     backst, backstc, &
                 gamdot, rss, gamdotC, rssc,  &
                 drot, fCeiDevHat, fCeDevVolHat,  &
                 pHatVec, ppTHat,  &
                 KpHatVec, KpHatdevVec, KpKpTHat, &
                 matProp, e_kk, dtime, tolerNewt, numslip,  &
                 maxIterNewt, iterCounterN, ierr, igrn, globalIncr, kstep,  &
                 sGamPPt, sGamKpKpt)
         if (ierr .ne. XTAL_CONVERGED) return
         norm_tau = dsqrt(InnerProductVec(tau, tau, NVECS))
!
!---------- solve for hardness and rotation tensor
!
         !write(*,*) 'intcryeqn 5'
         !call HardeningRotationSolve(kappa, kappa_n, eqvalues, rrot,  &
         call HardeningRotationSolve(slipr, backst, backstc, slipr_n,  &
         		 backst_n, backstc_n, eqvalues, rrot,  &
         		 rrot_n, crot, crot0, drot, tau, w_vec,  &
                 gamdot, rss, gamdotC, rssc,  &
                 pHatVec, qHatVec, KpHatdevVec, KqHatVec,  &
                 matProp, hardmtx, epsdot, dtime, numslip,  &
!     &           kHARD_MIDP)
!                 kHARD_ANAL)
                 kHARD_EXPL)
         !norm_kapp = dsqrt(InnerProductVec(kappa, kappa, NKAPP))
         norm_slipr = dsqrt(InnerProductVec(slipr, slipr, NKAPP))
         norm_backst = dsqrt(InnerProductVec(backst, backst, NKAPP))
         norm_backstc = dsqrt(InnerProductVec(backstc, backstc, NKAPP))         
!
!---------- check convergence
!
         !converged = ConvergeState(norm_tau, norm_kapp, norm_tau0,  &
                                   !norm_kapp0, tolerState)
         converged = ConvergeState(norm_tau, norm_slipr,  &
                                   norm_backst, norm_backstc,  &
                                   norm_tau0, norm_slipr0,  &
                                   norm_backst0, norm_backst0c, &
                                   tolerState)
         !write(*,*) 'intcryeqn 6'

      enddo
      !write(*,*) "rss(is)=", (rss(is), is=1,numslip)
      !write(*,*) "rssc(is)=", (rssc(is), is=1,numslip)
!
!------- keep track of state iteration and check number of iterations
!
      iterCounterS = iterState
      if (iterState .ge. maxIterState) then
         call WriteWarning(XTAL_O,  &
                  'StressSolveElastoViscoPlastic: iters > maxIters')
         ierr = XTAL_MAX_ITERS_HIT
         return
      endif
!
!------- continue sub-incrementation process if needed
!
      if (subIncr .lt. totSubIncrs) return
!
!------------------------ VOLUMETRIC RESPONSE
!
!------- integrate crystal volumetric  response
!
!      call RotMat5x5ForSymm(crot, qr5x5, DIMS)
!      call MultAxu(qr5x5, fCeDevVol, fCeDevVolHat, NVECS)
      call StressSolveVolumetric(statev, fCeDevVolHat, fCeVol, estran,  &
                                 e_kk, dtime)
!
!----------------------- CAUCHY STRESS - CONSISTENT TANGENT
!
!------- compute Cauchy stress s_ij (note that "tau" is deviatoric)
!-------  s_ij = stress + statev(kPRES)*I  (i.e., stress is deviatoric)
!-------  ee_ij = estran + statev(KEVOL)*I  (i.e., estran is deviatoric)
!
      call UpdateStress(s_ij, stress, tau, ee_ij, estran, statev,  &
                        eqvalues)
!
!------- compute norm of spin, slip system activity, and eqp
!
      call UpdateDeformQnts(gamdot, gamdotC,  &
                            eqvalues, statev, crot,  &
                            zBar0, KBar0,  &
                            dtime, numslip)

!      call RotMat5x5ForSymm(crot, qr5x5, DIMS)
!      call MultQAQT(qr5x5, fCeiDev, fCeiDevHat, NVECS)
!      call MultAxu(qr5x5, fCeDevVol, fCeDevVolHat, NVECS)

      call PlasticModuli(cepmod, fCeiDevHat, fCeDevVolHat, fCeVol,  &
                         statev(kDETVe), sGamPPt+sGamKpKpt)

      !write(*,*) 'flag2'

      return
      END

!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE RotateSlipGeometry(  &
         crot, pBar0Vec, qBar0Vec, qr5x5, qr3x3, pHatVec, qHatVec,  &
         ppTHat, numslip  &
         )

      implicit none
      include 'params_xtal.inc'

      integer numslip
      real*8  crot(DIMS,DIMS)
      real*8  pBar0Vec(NVECS,MAX_SLIP), qBar0Vec(DIMS,MAX_SLIP)
      real*8  qr5x5(NVECS,NVECS), qr3x3(DIMS,DIMS)
      real*8  pHatVec(NVECS,MAX_SLIP), qHatVec(DIMS,MAX_SLIP)
      real*8  ppTHat(NVECS,NVECS,MAX_SLIP)

      integer is
!
!---------------------------------------------------------------------72
!
!------- rotation matrices for symm and skew quantities
!
      call RotMat5x5ForSymm(crot, qr5x5, DIMS)
      call RotMat3x3ForSkew(crot, qr3x3, DIMS)
!
!------- rotate Scmidt tensors (vectors) to Btilde configuration
!
      do is = 1, numslip
         call MultAxu(qr5x5, pBar0Vec(1,is), pHatVec(1,is), NVECS)
         call MultAxu(qr3x3, qBar0Vec(1,is), qHatVec(1,is), DIMS)
         call OuterProductVec(pHatVec(1,is), pHatVec(1,is),  &
                              ppTHat(1,1,is), NVECS)
      enddo

      return
      END
!
!=====================================================================72
!
!=====================================================================72
!
      SUBROUTINE RotateSlipGeometry2(  &
         crot, KpBar0Vec, KpBar0devVec, KqBar0Vec, qr5x5, qr3x3, KpHatVec,  &
         KpHatdevVec, KqHatVec, KpKpTHat, numslip)

      implicit none
      include 'params_xtal.inc'

      integer numslip
      real*8  crot(DIMS,DIMS)
      real*8  KpBar0Vec(NVECS,MAX_SLIP), KqBar0Vec(DIMS,MAX_SLIP)
      real*8  KpBar0devVec(NVECS,MAX_SLIP)
      real*8  qr5x5(NVECS,NVECS), qr3x3(DIMS,DIMS)
      real*8  KpHatVec(NVECS,MAX_SLIP), KqHatVec(DIMS,MAX_SLIP)
      real*8  KpHatdevVec(NVECS,MAX_SLIP)
      real*8  KpKpTHat(NVECS,NVECS,MAX_SLIP)

      integer is
!
!---------------------------------------------------------------------72
!
!------- rotation matrices for symm and skew quantities
!
      call RotMat5x5ForSymm(crot, qr5x5, DIMS)
      call RotMat3x3ForSkew(crot, qr3x3, DIMS)
!
!------- rotate Climb tensors (vectors) to Btilde configuration
!
      do is = 1, numslip
         call MultAxu(qr5x5, KpBar0Vec(1,is), KpHatVec(1,is), NVECS)
         call MultAxu(qr5x5, KpBar0devVec(1,is), KpHatdevVec(1,is), NVECS)
         call MultAxu(qr3x3, KqBar0Vec(1,is), KqHatVec(1,is), DIMS)
         call OuterProductVec(KpHatVec(1,is), KpHatVec(1,is),  &
                              KpKpTHat(1,1,is), NVECS)
      enddo

      return
      END
!
!=====================================================================72
!
!=====================================================================72
!
      SUBROUTINE HardeningRotationSolve(  &
         !kappa, kappa_n, eqvalues, rrot, rrot_n, crot, crot0, drot,  &
         slipr, backst, backstc, slipr_n, backst_n, backstc_n, eqvalues, rrot,  &
         rrot_n, crot, crot0, drot, tau, w_vec,  &
         gamdot, rss, gamdotC, rssc,  &
         pHatVec, qHatVec, KpHatdevVec, KqHatVec,  &
         matProp, hardmtx, epsdot, dtime, numslip, kInteg_Hard  &
         )
!
      use     SlipOverStress
      implicit none
      include 'params_xtal.inc'

      integer numslip, kInteg_Hard
      real*8  epsdot, dtime
      !real*8  kappa(NKAPP), kappa_n(NKAPP), eqvalues(NEQVA)
      real*8  slipr(NKAPP), backst(NKAPP), slipr_n(NKAPP), backst_n(NKAPP), eqvalues(NEQVA)
      real*8  backstc(NKAPP), backstc_n(NKAPP)
      real*8  rrot(DIMS,DIMS), rrot_n(DIMS,DIMS), crot(DIMS,DIMS)
      real*8  crot0(DIMS,DIMS), drot(DIMS,DIMS), tau(NVECS), w_vec(DIMS)
      real*8  gamdot(MAX_SLIP), gamdotC(MAX_SLIP), matProp(NPROPS,MAX_SLIP)
      real*8  pHatVec(NVECS,MAX_SLIP), qHatVec(DIMS,MAX_SLIP)
      real*8  KpHatdevVec(NVECS,MAX_SLIP), KqHatVec(DIMS,MAX_SLIP)
      real*8  rss(MAX_SLIP), rssc(MAX_SLIP)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      integer is, i, j
      !real*8  InnerProductVec, SSKineticEqn
      real*8  InnerProductVec, SSKineticEqn2, SSKineticEqn3

      real*8  crss, bss, bssc
!
!---------------------------------------------------------------------72
!
!------- slip quantities: resolve shear stress and shear strain rate
!        glide + climb
!
      do is = 1, numslip
!      
         rss(is) = InnerProductVec(tau, pHatVec(1,is), NVECS)
!         crss = kappa(1) + overstress(is)
         !crss = kappa(is)
         crss = slipr(is)
         bss = backst(is)
         gamdot(is) = SSKineticEqn2(rss(is),crss,bss,matProp(1,is),kGAMDOT)
!
!         
         rssc(is) = InnerProductVec(tau, KpHatdevVec(1,is), NVECS)
         bssc = backstc(is)
         gamdotC(is) = SSKineticEqn3(rssc(is), bssc, matProp(1,is), kGAMDOTc)
      enddo
!
!------- solve for hardness (evolution equations of glide and climb rates)
      call IntegrateHardening2(slipr, backst, backstc, slipr_n, backst_n, backstc_n,  &
      						   eqvalues, gamdot, gamdotC, matProp,  &
      						   hardmtx, epsdot, dtime, numslip, kInteg_Hard)
   
!
!------------------------------------------------------------------------      
!
!------- solve for rotation: dotR = (W - Wp)*R
!------- update slip system orientation matrix: C = R * C0
!
!        Check whether the input data is passed correctly, OK!
      !write(*,*) 'Check input for IntegrateRotation: KqHatVec, OK'
	  !do i = 1, numslip
      !	write(*,*) (KqHatVec(j,i),j=1,DIMS)
      !enddo            

      call IntegrateRotation(w_vec, rrot_n, crot0,  &
                             qHatVec, gamdot,  &
                             KqHatVec, gamdotC,  &
                             rrot, crot, drot, dtime, numslip)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE IntegrateHardening2(  &
         slipr, backst, backstc, slipr_n, backst_n, backstc_n,  &
         eqvalues, gamdot, gamdotC, matProp, &
         hardmtx, epsdot, dtime, numslip, kInteg_Code )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, kInteg_Code
      real*8  epsdot, dtime
      real*8  slipr(numslip), backst(numslip), backstc(numslip)
      real*8  slipr_n(numslip), backst_n(numslip), backstc_n(numslip)
      real*8  eqvalues(NEQVA), matProp(NPROPS, MAX_SLIP)
      real*8  gamdot(MAX_SLIP), gamdotC(MAX_SLIP)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      integer is, ik, jk
      real*8  h_0, tausi, taus0, xms, gamss0
      real*8  shr_min, shr_max, kappa_sat
      real*8  c, g_n, g_s, g
      real*8  dkappa, gamtot_n, delgam, fac
      real*8  kTHETA
      data    kTHETA /1.0d0/
!	  Parameters related to glide slip evolution
      real*8  S0, fc, hB, hS, dD, mu, mu0, mu0p, lambda, rD, SHat, SBar, staticHard
!	  Parameters related to climb slip evolution
      real*8  hC, backstINF, hC2
      
!
!---------------------------------------------------------------------72
!
!------- total shear strain rate: Sum(|gamdot(i)|), i=1,numslip
!

      S0 =143.41e6
      fc = 0.42d0
      hB = 1.043e8
      hS =3.9773e8
      dD =5.07362e3
      mu = 79.35d9
      mu0 = 265.3d9
      mu0p = 31.74d9
      lambda = 0.86d0
      SBar=18.03e6
      SHat=SBar+hS/dD
      
!------ climb parameters
	  hC = 0.5d0
	  backstINF = 300d6
	  hC2 = 32d0

      if ((hS/dD) .gt. (S0-SBar) )     then

             write(*,*) 'make sure hS/dD < S0-SBar'
             stop
      endif

      eqvalues(kSHRATE) = pzero
      do is = 1, numslip
         eqvalues(kSHRATE) = eqvalues(kSHRATE) + dabs(gamdot(is))
      enddo

      shr_min = 1.0d-6 * epsdot
      shr_max = 1.0d+1 * epsdot

      if (eqvalues(kSHRATE) .le. shr_min) eqvalues(kSHRATE) = shr_min
      if (eqvalues(kSHRATE) .ge. shr_max) eqvalues(kSHRATE) = shr_max
!
!------- accumulated shear strain: gamtot
!
      eqvalues(kGAMTOT) = eqvalues(kGAMTOT_n) + eqvalues(kSHRATE)*dtime
!
!------- integration of Voce's hardening law (one hardness/slip system)
!
      if (kInteg_Code .eq. kHARD_EXPL) then  
!
!---------- explicit update
!         do ik = 1, NKAPP
         do ik = 1, numslip
         
         	staticHard=0.d0
!         	if (abs(gamdot(ik)).le. 1.0d-13)  staticHard=1.d0
!          	staticHard=( atan((1.0d-8- eqvalues(kSHRATE))*1.0d10)/ atan (1.d0)/2+1.d0)/2.d0
         	if (eqvalues(kSHRATE) .le. 1.01d-6)  staticHard=1.d0
         	
         	!============================= Back stress evolution of glide ======================================
            ! S_alpha=... Eq (20) in MSMS paper (Xiang & Caglar, 2016)
            slipr(ik) = slipr_n(ik) + &
                    	dtime*(hS-dD*(slipr_n(ik)-SBar))*dabs(gamdot(ik)) &
                    	- 0.015d0*dtime*(slipr_n(ik)-S0)*staticHard
            
            ! D_alpha=... Eq (22) in MSMS paper (Xiang & Caglar, 2016)
            rD = (hB*mu0/slipr_n(ik))/(mu0p/(fc*lambda)-mu)
            
            ! B_alpha=... Eq (21) in MSMS paper (Xiang & Caglar, 2016)
            backst(ik) = backst_n(ik) + &
                  	     dtime*( hB*gamdot(ik)-rD*backst_n(ik)*dabs(gamdot(ik)) )
         end do	
 		         	
         do ik = 1, numslip 	
         	!============================= Back stress evolution of climb ======================================
         	! Bc_alpha=... Eq (14) Staroselsky and Caeenti, IJSS 2011
            backstc(ik) = backstc_n(ik) + dtime*( hC*backstINF*gamdotC(ik)  &
                                                - hC2*backstc_n(ik)*dabs(gamdotC(ik)) )
         end do

      else
!
!------- wrong code number
         call RunTimeError(XTAL_O,  &
                            'IntegrateHardening: Wrong kInteg_Code!')

      endif

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE IntegrateHardening(  &
         kappa, kappa_n, eqvalues, gamdot, matProp, hardmtx, epsdot,  &
         dtime, numslip, kInteg_Code  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, kInteg_Code
      real*8  epsdot, dtime
      real*8  kappa(NKAPP), kappa_n(NKAPP), eqvalues(NEQVA)
      real*8  gamdot(MAX_SLIP), matProp(NPROPS, MAX_SLIP)
      real*8  hardmtx(MAX_SLIP, MAX_SLIP)

      integer is, ik, jk
      real*8  h_0, tausi, taus0, xms, gamss0
      real*8  shr_min, shr_max, kappa_sat
      real*8  c, g_n, g_s, g
      real*8  dkappa, gamtot_n, delgam, fac
      real*8  kTHETA
      data    kTHETA /1.0d0/
!
!---------------------------------------------------------------------72
!
!------- total shear strain rate: Sum(|gamdot(i)|), i=1,numslip
!
      eqvalues(kSHRATE) = pzero
      do is = 1, numslip
         eqvalues(kSHRATE) = eqvalues(kSHRATE) + dabs(gamdot(is))
      enddo

      shr_min = 1.0d-6 * epsdot
      shr_max = 1.0d+1 * epsdot

      if (eqvalues(kSHRATE) .le. shr_min) eqvalues(kSHRATE) = shr_min
      if (eqvalues(kSHRATE) .ge. shr_max) eqvalues(kSHRATE) = shr_max
!
!------- accumulated shear strain: gamtot
!
      eqvalues(kGAMTOT) = eqvalues(kGAMTOT_n) + eqvalues(kSHRATE)*dtime
!
!------- integration of Voce's hardening law (one hardness/slip system)
!
      if (kInteg_Code .eq. kHARD_EXPL) then  
!
!---------- explicit update
!         do ik = 1, NKAPP
         do ik = 1, numslip
            h_0    = matProp(5,ik)
            tausi  = matProp(6,ik)
            taus0  = matProp(7,ik)
            xms    = matProp(8,ik)
            gamss0 = matProp(9,ik)

            kappa_sat = taus0 * ((eqvalues(kSHRATE) / gamss0)**xms)

            c = dtime*h_0
            g_s = kappa_sat - tausi
            g_n = kappa_n(ik) - tausi

            if ( (g_n/g_s) .le. 1.0 ) then
               g = g_n + c*(1.0-g_n/g_s)*eqvalues(kSHRATE)
            else
               g = g_n
            endif
            kappa(ik) = g + tausi
         enddo

      else if (kInteg_Code .eq. kHARD_MIDP) then
!
!---------- generalized mid-point rule 
!         do ik = 1, NKAPP
         do ik = 1, numslip
            h_0    = matProp(5,ik)
            tausi  = matProp(6,ik)
            taus0  = matProp(7,ik)
            xms    = matProp(8,ik)
            gamss0 = matProp(9,ik)

            kappa_sat = taus0 * ((eqvalues(kSHRATE) / gamss0)**xms)

            c = dtime*h_0
            g_s = kappa_sat - tausi
            g_n = kappa_n(ik) - tausi

            if ( (g_n/g_s) .le. 1.0 ) then
               g = g_n + c*(  &
                      (1.0-kTHETA)*(1.0-g_n/g_s)*eqvalues(kSHRATE_n)  &
                      +   (kTHETA)*eqvalues(kSHRATE)  &
                           )
               g = g / (1.0 + c*kTHETA*eqvalues(kSHRATE)/g_s)
            else
               g = g_n
            endif
            kappa(ik) = g + tausi
         enddo

      else if (kInteg_Code .eq. kHARD_ANAL) then

         gamtot_n = eqvalues(kGAMTOT_n)
         delgam   = eqvalues(kSHRATE) * dtime

         do ik = 1, numslip

            h_0    = matProp(5,ik)
            tausi  = matProp(6,ik)
            taus0  = matProp(7,ik)
            xms    = matProp(8,ik)
            gamss0 = matProp(9,ik)

            kappa_sat = taus0 * ((eqvalues(kSHRATE) / gamss0)**xms)
            g_s = kappa_sat - tausi
            fac = dabs(h_0/g_s)

            dkappa = 0.0
            do jk = 1, numslip
              dkappa = dkappa +  &
                          hardmtx(ik,jk)*dabs(gamdot(jk))*dtime/delgam
            enddo
            dkappa = dkappa * g_s * dexp(-gamtot_n*fac) *  &
                                              (1.0 - dexp(-delgam*fac))

            kappa(ik) = kappa_n(ik) + dkappa

         enddo
         
      else
!
!------- wrong code number
         call RunTimeError(XTAL_O,  &
                            'IntegrateHardening: Wrong kInteg_Code!')

      endif

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE IntegrateRotation(w_vec, rrot_n, crot0,  &
         qHatVec, gamdot,  &
         KqHatVec, gamdotC,  &         
         rrot, crot, drot, dtime, numslip  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip
      real*8  dtime
      real*8  w_vec(DIMS), rrot_n(DIMS,DIMS), crot0(DIMS,DIMS)
      real*8  qHatVec(DIMS,MAX_SLIP), gamdot(MAX_SLIP)
      real*8  KqHatVec(DIMS,MAX_SLIP), gamdotC(MAX_SLIP)
      real*8  rrot(DIMS,DIMS), crot(DIMS,DIMS), drot(DIMS,DIMS)

      integer is
      real*8  wpHat(DIMS), spin(DIMS)
!
!---------------------------------------------------------------------72
!
!------- plastic spin from slip system activity (glide + climb)
!
      call SetTensor(wpHat, pzero, DIMS)
      do is = 1, numslip
         call AddScaledTensor(gamdot(is), qHatVec(1,is), wpHat, DIMS)
         call AddScaledTensor(gamdotC(is), KqHatVec(1,is), wpHat, DIMS)
      enddo
!
!------- net spin: (W - Wp)
!
      call SetTensor(spin, pzero, DIMS)
      
      call AddTensors(pone, w_vec, -pone, wpHat, spin, DIMS)
!
!------- incremental rot: dR = exp(spin*dtime)
!
      call IncrementalRotTensor(drot, spin, dtime)
!
!------- rotation tensor: R = dR*R_n
!
      call MultAxB(drot, rrot_n, rrot, DIMS)
!
!------- slip system orientation matrix crot: C = R*C_0
!
      call MultAxB(rrot, crot0, crot, DIMS)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE IncrementalRotTensor(  &
         drot, spin, dtime  &
         )

      implicit none

      real*8  dtime
      real*8  drot(3, 3), spin(3)

      real*8  th1, th2, th3, tau, taua
!
!---------------------------------------------------------------------72
!
      th1 = spin(1) * dtime
      th2 = spin(2) * dtime
      th3 = spin(3) * dtime

      tau = sqrt(th1 * th1 + th2 * th2 + th3 * th3)
      if (tau .ne. 0.0) then
         taua = tan(tau / 2.0) / tau
         th1 = taua * th1
         th2 = taua * th2
         th3 = taua * th3
         tau = taua * tau
      endif

      tau = 2.0 / (1.0 + tau * tau)

      drot(1, 1) = 1.0 - tau * (th1 * th1 + th2 * th2)
      drot(1, 2) = - tau * (th1 + th2 * th3)
      drot(1, 3) = tau * ( - th2 + th1 * th3)
      drot(2, 1) = tau * (th1 - th2 * th3)
      drot(2, 2) = 1.0 - tau * (th1 * th1 + th3 * th3)
      drot(2, 3) = - tau * (th3 + th1 * th2)
      drot(3, 1) = tau * (th2 + th1 * th3)
      drot(3, 2) = tau * (th3 - th1 * th2)
      drot(3, 3) = 1.0 - tau * (th2 * th2 + th3 * th3)
 
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      logical FUNCTION ConvergeState(  &
         norm_a, norm_b, norm_c, norm_d, norm_a0, norm_b0, norm_c0, norm_d0, toler  &
         )

      implicit none

      real*8  norm_a, norm_b, norm_c, norm_d,  &
      	      norm_a0, norm_b0, norm_c0, norm_d0, toler
!
!---------------------------------------------------------------------72
!
      ConvergeState = (dabs(norm_a - norm_a0) .lt. max(toler*norm_a0,1d-16)) .and.  &
                      (dabs(norm_b - norm_b0) .lt. max(toler*norm_b0,1d-16)) .and.  &
                      (dabs(norm_c - norm_c0) .lt. max(toler*norm_c0,1d-16)) .and.  &
                      (dabs(norm_d - norm_d0) .lt. max(toler*norm_d0,1d-16))

      if (.not.ConvergeState) then
          norm_a0 = norm_a
          norm_b0 = norm_b
          norm_c0 = norm_c
          norm_d0 = norm_d
      endif

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE UpdateStress(  &
         s_ij, stress, tau, ee_ij, estran, statev, eqvalues  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      real*8  s_ij(DIMS,DIMS), stress(NVECS), tau(NVECS)
      real*8  ee_ij(DIMS,DIMS), estran(NVECS), statev(NSTAV)
      real*8  eqvalues(NEQVA)

      real*8  ekk3
      real*8  Ve_ij(DIMS, DIMS)

      real*8  DetOfMat3x3, InnerProductVec
!
!---------------------------------------------------------------------72
!
!------- Elastic strain (ee) and V tensors (V = ee + I)
!
      call Vec5x1ToMat3x3Symm(estran, ee_ij, DIMS)
      ekk3 = statev(kEVOL)/3.0
      call AddScaledTensor(ekk3, Ident2nd, ee_ij, DIMS*DIMS)

      call AddTensors(pone, Ident2nd, pone, ee_ij, Ve_ij, DIMS*DIMS)
!	  det(I + ee) in Eq. 3.46
      statev(kDETVe) = DetOfMat3x3(Ve_ij)
!
!------- deviatoric (vector) and volumetric (pressure) Cauchy stress
!
      call SetToScaledtensor(1./statev(kDETVe), tau, stress, NVECS)
      statev(kPRES) = statev(kPRES)/statev(kDETVe)
!
!------- VonMises stress
!
      eqvalues(kMISES) = dsqrt(3./2.*InnerProductVec(stress, stress,  &
                                                      NVECS))
!
!------- Cauchy stress (tensor form)
!
      call Vec5x1ToMat3x3Symm(stress, s_ij, DIMS)
      call AddScaledTensor(statev(kPRES), Ident2nd, s_ij, DIMS*DIMS)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE PlasticModuli(  &
         fCep, fCeiDevHat, fCeDevVolHat, fCeVol, jac_e, lhs5x5  &
         )
!
      use     TransData
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      real*8  fCeVol, jac_e
      real*8  fCep(DIMS2,DIMS2)
      real*8  fCeiDevHat(NVECS,NVECS), fCeDevVolHat(NVECS)

      real*8  det
      real*8  fMat5x5(NVECS, NVECS), fVect5(NVECS)
      real*8  fOuter6x6(DIMS2, DIMS2), fVect6(DIMS2)
      real*8  fMat5x6(NVECS, DIMS2), fOuter5x6(NVECS, DIMS2)
      real*8  fCepDev(NVECS, DIMS2), fCepVol(DIMS2)

      real*8  lhs5x5(NVECS, NVECS)
!
!
!---------------------------------------------------------------------72
!
!---- Deviatoric moduli: fCepDev_5x6
!
!-------------- ([Ce_d]^-1 + [S])^-1
      !write(*,*) 'PlasMod 1'
      call AddTensors(pone, fCeiDevHat, pone, lhs5x5, fMat5x5,  &
                      NVECS*NVECS) 
      call MatrixInverse(fMat5x5, NVECS, NVECS, det)
      !write(*,*) 'PlasMod 2'
!
!-------------- ([T][Pdev] + [Ce_d]^-1 {H} {1}^T)
      call MultAxu(fCeiDevHat, fCeDevVolHat, fVect5, NVECS)
      call OuterProductVec_G(fvect5, Ident, fOuter5x6, NVECS, DIMS2)
      call Addtensors(pone, fMatTId5x6, pone, fOuter5x6, fMat5x6,  &
                      NVECS*DIMS2)
      !write(*,*) 'PlasMod 3'
!
!-------------- (([Ce_d]^-1 + [S])^-1) ([T][Pdev] + [Ce_d]^-1 {H} {1}^T)
      call MultAxB_G(fMat5x5, fMat5x6, fCepDev, NVECS, NVECS, DIMS2)
      !write(*,*) 'PlasMod 4'
!
!---- Volumetric moduli: fCepVol_6
!
!-------------- {H}^T [T] [Pdev] + M {1}^T
      call MultATxu_G(fMatTId5x6, fCeDevVolHat, fVect6, NVECS, DIMS2)
      call AddTensors(pone, fVect6, fCeVol, Ident, fCepVol, DIMS2)
      !write(*,*) 'PlasMod 5'
!
!-------------- {H}^T [S] [Cep_d]
      call MultAxB_G(lhs5x5, fCepDev, fMat5x6, NVECS, NVECS, DIMS2)
      call MultATxu_G(fMat5x6, fCeDevVolHat, fVect6, NVECS, DIMS2)
      !write(*,*) 'PlasMod 6'
!
!-------------- {H}^T [T] [Pdev] + M {1}^T - {H}^T [S] [Cep_d]
      call AddScaledTensor(-pone, fVect6, fCepVol, DIMS2)
      !write(*,*) 'PlasMod 7'
!
!---- Algorithmic Moduli: fCep_6x6
!
!-------------- [J] [Cep_d]
      call MultAxB_G(fDevMat6x5, fCepDev, fCep, DIMS2, NVECS, DIMS2)
      !write(*,*) 'PlasMod 8'
!
!-------------- [J] [Cep_d] + {1} {Cep_v}^T -> [Cep]: {s} = [Cep] {d}
      call OuterProductVec(Ident, fCepVol, fOuter6x6, DIMS2)
      call AddScaledTensor(pone, fOuter6x6, fCep, DIMS2*DIMS2)
      call SetToScaledTensor(pone/jac_e, fCep, fCep, DIMS2*DIMS2)
      !write(*,*) 'PlasMod 9'
!
!---- Change fCep to be consistent with abaqus representation of {d}
!
      call SetToScaledTensor(phalf, fCep(1,4), fCep(1,4), DIMS2*DIMS2/2)
      !write(*,*) 'PlasMod 10'

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE UpdateDeformQnts(  &
         gamdot, gamdotC,  &
         eqvalues, statev, crot,  &
         zBar0, KBar0,  &
         dtime, numslip  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip
      real*8  dtime
      real*8  statev(NSTAV), eqvalues(NEQVA), gamdot(MAX_SLIP), gamdotC(MAX_SLIP)
      real*8  crot(DIMS,DIMS), zBar0(DIMS,DIMS,MAX_SLIP)
	  real*8  KBar0(DIMS,DIMS,MAX_SLIP)
      integer is
      real*8  eqpdot, maxGamdot, ratio
      real*8  wp_vec(DIMS), matx_1(DIMS,DIMS)
      real*8  lp_ij(DIMS,DIMS), dp_ij(DIMS,DIMS), wp_ij(DIMS, DIMS)

      real*8  InnerProductVec, MaxAbsValueOfVec
!
!---------------------------------------------------------------------72
!
!------- Plastic velocity gradient lp in crystal axis
!	     glide + climb
!
      call setTensor(lp_ij, pzero, DIMS*DIMS)
      
      do is = 1, numslip
         call AddScaledTensor(gamdot(is), zBar0(1,1,is), lp_ij, DIMS*DIMS)
         call AddScaledTensor(gamdotC(is), KBar0(1,1,is), lp_ij, DIMS*DIMS)
      enddo
!
!------- Symm and Skew parts of lp: lp =  dp + wp
!
!      call EqualTensors(lp_ij, matx_1, DIMS*DIMS)
!      call MultQAQT(crot, matx_1, lp_ij, DIMS)

      call SymmetrizeTensor(lp_ij, dp_ij, DIMS)
      call SkewSymmetrizeTensor(lp_ij, wp_ij, DIMS)
!
!------- Equivalent plastic strain included glide + climb
!
      eqpdot = dsqrt(2./3.*InnerProductVec(dp_ij, dp_ij, DIMS*DIMS))
      eqvalues(kEQP) = eqvalues(kEQP_n) + dtime*eqpdot
!
!------- Norm of plastic spin
!
      call Mat3x3ToVec3x1Skew(wp_ij, wp_vec, DIMS)
      statev(kWPNORM) = dsqrt(InnerProductVec(wp_vec, wp_vec, DIMS))
!
!------- Slip system activity ???
!
!	  Find max(gamdot(is))
      maxGamdot = MaxAbsValueOfVec(gamdot, numslip)
      if (maxGamdot .lt. TINY) maxGamdot = TINY

      statev(kSSACT) = 0.0
      do is = 1, numslip
         ratio = dabs(gamdot(is))/maxGamdot
         if (ratio .ge. 0.05d0) statev(kSSACT) = statev(kSSACT) + pone
      enddo
        
      return
      END
!
!=====================================================================72
!
!		NOT USED in this UMAT
!=====================================================================72
!
      SUBROUTINE StressSolveViscoPlastic(  &
         !tau_lat, d_vec_lat, kappa, pBar0Vec, ppTBar0, matProp, sigfs,  &
         tau_lat, d_vec_lat, slipr, backst, pBar0Vec, ppTBar0, KpBar0Vec, KpBar0devVec, &
         KpKpTBar0, matProp, sigfs,  &
         tauSlip, dtime, theta, thetao, epsdot, tolerNewt, maxIterNewt,  &
         numslip, numvtx  &
         )

      implicit none
      include 'params_xtal.inc'

      integer maxIterNewt, numslip, numvtx
      real*8  dtime, theta, thetao, epsdot, tolerNewt
      !real*8  tau_lat(NVECS), d_vec_lat(NVECS), kappa(NKAPP)
      real*8  tau_lat(NVECS), d_vec_lat(NVECS), slipr(NKAPP), backst(NKAPP)
      real*8  pBar0Vec(NVECS,MAX_SLIP), ppTBar0(NVECS,NVECS,MAX_SLIP)
      real*8  KpBar0Vec(NVECS,MAX_SLIP), KpBar0devVec(NVECS, MAX_SLIP), KpKpTBar0(NVECS,NVECS,MAX_SLIP)
      real*8  sigfs(NVECS,MAX_VTX), tauSlip(MAX_SLIP)
      real*8  matProp(NPROPS,MAX_SLIP)

      integer iv, ierr
      real*8  rescale
      real*8  sdotd(MAX_VTX)
      real*8  InnerProductVec
!
!---------------------------------------------------------------------72
!
!------- solve for stresses in lattice reference frame
!
!------- use effective deformation rate to normalized (scale-down) D
!
!      epsdot = dsqrt(2./3.*InnerProductVec(d_vec_lat, d_vec_lat, NVECS))
      call SetToScaledTensor(1./epsdot, d_vec_lat, d_vec_lat, NVECS)
!
!------- compute work to select guessed values from vertex stresses
!
      do iv = 1, numvtx
         sdotd(iv) = InnerProductVec(sigfs(1,iv), d_vec_lat, NVECS) 
      enddo
!
!------- initialized local flag to monitor convergence of Newton iters
!
      ierr = XTAL_CONVERGED
!
!------- solve fot the stresses. If needed, will try all vertices
!
      do iv = 1, numvtx
!
!---------- compute and scaled initial guess
         !call InitialGuessStress(tau_lat, sdotd, kappa, sigfs, pBar0Vec,  &
         call InitialGuessStress(tau_lat, sdotd, slipr, backst, sigfs, pBar0Vec,  &
                                 matProp, numslip, numvtx)
!
!---------- compute stresses
         !call StressViscoPlastic(tau_lat, d_vec_lat, kappa, pBar0Vec,  &
         !write(*,*) 'strsolvispla 1'
         call StressViscoPlastic(tau_lat, d_vec_lat, slipr, backst, pBar0Vec,  &
                                 ppTBar0, KpBar0Vec, KpBar0devVec, KpKpTBar0, matProp, tolerNewt,  &
                                 maxIterNewt, numslip, ierr)
         !write(*,*) 'strsolvispla 2'
! 
!---------- if not converged, try another initial guess
         if (ierr .ne. XTAL_CONVERGED) then
            call WriteMessage(XTAL_O,  &
               'StressSolveViscoPlastic: did not converged !!!')
            call WriteMessage(XTAL_O,  &
               'StressSolveViscoPlastic: will try next vertex stresses')
            ierr = XTAL_CONVERGED
! 
!---------- or else converged; then, rescale the stresses, and return
         else
            rescale = epsdot**matProp(3,1)   ! epsdot**m
            call SetToScaledTensor(rescale, tau_lat, tau_lat, NVECS)
            return
         endif

      enddo
!
!------- if gets here, then the stress solution did not converged
!-------   for all trials. 
!
      call RunTimeError(XTAL_O,  &
             'StressSolveViscoPlastic: Newton method did not converged')

      return
      END
!
!=====================================================================72
!
!	NOT USED
!=====================================================================72
!
      SUBROUTINE InitialGuessStress(  &
         !s, sdotd, kappa, sigfs, pBar0Vec, matProp, numslip, numvtx  &
         s, sdotd, slipr, backst, sigfs, pBar0Vec, matProp, numslip, numvtx  &
         )
!
      use     SlipOverStress
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, numvtx
      !real*8  s(NVECS), sdotd(MAX_VTX), kappa(NKAPP)
      real*8  s(NVECS), sdotd(MAX_VTX), slipr(NKAPP), backst(NKAPP)
      real*8  sigfs(NVECS, MAX_VTX), pBar0Vec(NVECS, MAX_SLIP)
      real*8  matProp(NPROPS,MAX_SLIP)

      integer kmax, is
      real*8  signmax, factor
      real*8  rssk(MAX_SLIP)

      integer IndexMaxAbsValueOfVec
      real*8  MaxAbsValueOfVec, InnerProductVec, SignOf

      real*8  crss, bss
!
!---------------------------------------------------------------------72
!
!------- locate sign of max value of sdotd (plastic work)
!
      kmax = IndexMaxAbsValueOfVec(sdotd, numvtx)
      signmax = SignOf(sdotd(kmax))
      sdotd(kmax) = pzero        ! avoids 'kmax' vertex being re-used
!
!------- asign correct sign to selected vertex stress (initial guess)
!
      call SetToScaledTensor(signmax, sigfs(1,kmax), s, NVECS)
!
!------ compute factor to scaled down initial guess
!------   factor = |rss|_max / kappa * gam0**m
!
      do is = 1, numslip
!         crss = kappa(1) + overstress(is)
         !crss = kappa(is)
         crss = slipr(is)
         bss = backst(is)
         !rssk(is) = InnerProductVec(pBar0Vec(1,is), s, NVECS)/crss
         rssk(is) = InnerProductVec(pBar0Vec(1,is), s, NVECS)
      enddo
      factor = 1.0
      !factor = MaxAbsValueOfVec(rssk, numslip) *  &
                                            !matProp(4,1)**matProp(3,1)
!
!------- scaled initial guess for stress
!
      call SetToScaledTensor(1./factor, s, s, NVECS)

      return
      END
!
!=====================================================================72
!
!	NOt USED
!=====================================================================72
!
      SUBROUTINE StressViscoPlastic(  &
         !s, d_vec_lat, kappa, pBar0Vec, ppTBar0, matProp, toler,  &
         s, d_vec_lat, slipr, backst, pBar0Vec, ppTBar0, KpBar0Vec, KpBar0devVec,  &
         KpKpTBar0, matProp, toler,  &
         maxIters, numslip, ierr  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer maxIters, numslip, ierr
      real*8  toler
      !real*8  s(NVECS), d_vec_lat(NVECS), kappa(NKAPP)
      real*8  s(NVECS), d_vec_lat(NVECS), slipr(NKAPP), backst(NKAPP)
      real*8  pBar0Vec(NVECS,MAX_SLIP), ppTBar0(NVECS,NVECS,MAX_SLIP)
      real*8  KpBar0Vec(NVECS,MAX_SLIP), KpBar0devVec(NVECS, MAX_SLIP), KpKpTBar0(NVECS,NVECS,MAX_SLIP)
      real*8  matProp(NPROPS,MAX_SLIP)

      integer iters
      real*8  rhs_norm_0, rhs_norm, search
      real*8  delts(NVECS), s0(NVECS)
      real*8  rhs(NVECS), rss(MAX_SLIP), lhs(NVECS,NVECS)

      logical Converged
!      real*8  InnerProductVec
      real*8  NormOfVector
!
!---------------------------------------------------------------------72
!
!------- compute initial residual and its norm
!
      !call FormResidual(rhs, s, rss, d_vec_lat, kappa, pBar0Vec,  &
     !write(*,*) 'strvispla 1'
      call FormResidual(rhs, s, rss, d_vec_lat, slipr, backst, pBar0Vec,  &
                        matProp, numslip)
!      rhs_norm_0 = dsqrt(InnerProductVec(rhs, rhs, NVECS))
      rhs_norm_0 = NormOfVector(rhs, NVECS)
!
!------- initialize flags and start iterations
!
      iters = 0
      do while(.not.Converged(rhs, toler) .and. iters .lt. maxIters)
         iters = iters + 1
         call EqualTensors(s, s0, NVECS)
!
!---------- compute local jacobian
         !call FormJacobian(lhs, rss, ppTBar0, kappa, matProp, numslip)
         !write(*,*) 'strvispla 2'
         call FormJacobian(lhs, rss, ppTBar0, KpKpTBar0, slipr, backst, matProp, numslip)
!
!---------- solve for the increment of stress
!---------- use the lower half of the symmetric matrix lhs
         !write(*,*) 'strvispla 2.1'
         write(*,*) 'lhs',lhs
         write(*,*) 'rhs',rhs
         call lsolve(lhs, rhs, 1, ierr, NVECS)
         !write(*,*) 'strvispla 2.2'
         if (ierr .eq. XTAL_SING_JACOBIAN) return
         call EqualTensors(rhs, delts, NVECS)
!
         !write(*,*) 'strvispla 2.3'
!---------- update variables
         search = pone
         call AddTensors(pone, s0, search, delts, s, NVECS)
!
!---------- compute new residual and its norm
         !call FormResidual(rhs, s, rss, d_vec_lat, kappa, pBar0Vec,  &
         !write(*,*) 'strvispla 3'
         call FormResidual(rhs, s, rss, d_vec_lat, slipr, backst, pBar0Vec,  &
                           matProp, numslip)
!         rhs_norm = dsqrt(InnerProductVec(rhs, rhs, NVECS))
         rhs_norm = NormOfVector(rhs, NVECS)
!
!---------- simple line search
         do while ( rhs_norm .gt. rhs_norm_0 ) 

            search = search*0.5d0
            if (search .lt. TINY) then
               call WriteWarning(XTAL_O,  &
                       'StressViscoPlastic: LS Failed, search < TINY')
               ierr = XTAL_LS_FAILED
               return
            endif
            call AddTensors(pone, s0, search, delts, s, NVECS)

            !call FormResidual(rhs, s, rss, d_vec_lat, kappa, pBar0Vec,  &
            !write(*,*) 'strvispla 4'
            call FormResidual(rhs, s, rss, d_vec_lat, slipr, backst, pBar0Vec,  &
                              matProp, numslip)
!            rhs_norm = dsqrt(InnerProductVec(rhs, rhs, NVECS))
            rhs_norm = NormOfVector(rhs, NVECS)

         enddo
!
!---------- update norm of residual
         rhs_norm_0 = rhs_norm

      enddo
!
!------- check convergence on iters < maxIters
!
      if (iters .ge. maxIters) then
        call WriteMessage(XTAL_O,'StressViscoPlastic: iters > maxIters')
        ierr = XTAL_MAX_ITERS_HIT
        return
      endif
!
!------- keep track of Newton's iterations
!
      write(XTAL_O, *) 'iterCounter (viscoplastic) = ', iters

      return
      END
!
!=====================================================================72
!
!	NOT USED
!=====================================================================72
!
      SUBROUTINE FormResidual(  &
         !rhs, s, rss, d_vec_lat, kappa, pBar0Vec, matProp, numslip  &
         rhs, s, rss, d_vec_lat, slipr, backst, pBar0Vec, matProp, numslip  &
         )
!
      use     SlipOverStress
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip
      real*8  rhs(NVECS), s(NVECS), rss(MAX_SLIP)
      !real*8  d_vec_lat(NVECS), kappa(NKAPP)
      real*8  d_vec_lat(NVECS), slipr(NKAPP), backst(NKAPP)
      real*8  pBar0Vec(NVECS,MAX_SLIP), matProp(NPROPS,MAX_SLIP)

      integer is
      real*8  gamdot(MAX_SLIP), d_curr(NVECS)
      !real*8  InnerProductVec, SSKineticEqn
      real*8  InnerProductVec, SSKineticEqn2

      real*8  crss, bss
!
!---------------------------------------------------------------------72
!
!------- resolve shear stresses and shear strain rates
!
      do is = 1, numslip
         rss(is) = InnerProductVec(s, pBar0Vec(1,is), NVECS)
!         crss = kappa(1) + overstress(is)
         !crss = kappa(is)
         crss = slipr(is)
         bss = backst(is)
         gamdot(is) = SSKineticEqn2(rss(is),crss,bss,matProp(1,is),kGAMDOT)
      enddo
!
!------- rate of deformation due to slip system activity
!
      call SetTensor(d_curr, pzero, NVECS)
      do is = 1, numslip
         call AddScaledTensor(gamdot(is), pBar0Vec(1,is), d_curr, NVECS)
      enddo
!
!------- residual
!
      call SetTensor(rhs, pzero, NVECS)
      call AddTensors(+pone, d_vec_lat, -pone, d_curr, rhs, NVECS)

      return
      END
!
!=====================================================================72
!
!	NOT USED
!=====================================================================72
!
      SUBROUTINE FormJacobian(  &
         !lhs, rss, ppTBar0, kappa, matProp, numslip  &
         lhs, rss, ppTBar0, KpKpTBar0, slipr, backst, matProp, numslip  &
         )
!
      use     SlipOverStress
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip
      !real*8  lhs(NVECS,NVECS), rss(MAX_SLIP), kappa(NKAPP)
      real*8  lhs(NVECS,NVECS), rss(MAX_SLIP), slipr(NKAPP), backst(NKAPP)
      real*8  lhs1(NVECS,NVECS), lhs2(NVECS,NVECS)
      real*8  ppTBar0(NVECS,NVECS,MAX_SLIP), matProp(NPROPS,MAX_SLIP)
      real*8  KpKpTBar0(NVECS,NVECS,MAX_SLIP)

      integer is
      !real*8  SSKineticEqn
      real*8  SSKineticEqn2
      real*8  dGamdTau(MAX_SLIP)

      real*8  crss, bss
!
!---------------------------------------------------------------------72
!
!------- d(gamdot)/d(tau)
!
      do is = 1, numslip
!         crss = kappa(1) + overstress(is)
         !crss = kappa(is)
         crss = slipr(is)
         bss = backst(is)
         dGamdTau(is) = SSKineticEqn2(rss(is),crss,bss,matProp(1,is),  &
                                     kdGAMdTAU)
      enddo
!
!------- jacobian
!
      call SetTensor(lhs, pzero, NVECS*NVECS)
      do is = 1, numslip
         call AddScaledTensor(dGamdTau(is), ppTBar0(1,1,is), lhs,  &
                              NVECS*NVECS)
      enddo

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      logical FUNCTION Converged(  &
         res, toler  &
         )

      implicit none
      include 'params_xtal.inc'

      real*8  toler
      real*8  res(NVECS)

      integer i
!
!---------------------------------------------------------------------72
!     
!------- check convergence on residual
!
      Converged = ( dabs(res(1)) .lt. toler )
      do i = 2, NVECS
         Converged = ( (dabs(res(i)) .lt. toler) .and. Converged)
      enddo

      return
      END
!     
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE StressSolveDeviatoric(  &
         !tau, d_vec, estran, estran_n, kappa,  &
         tau, d_vec, estran, estran_n, slipr, backst, backstc, &
         gamdot, rss, gamdotC, rssc,  &
         drot, fCeiDevHat, fCeDevVolHat,  &
         pHatVec, ppTHat,  &
         KpHatVec, KpHatdevVec, KpKpTHat,  &
         matProp, e_kk, dtime, tolerNewt, numslip,  &
         maxIterNewt, iterCounterN, ierr, igrn, globalIncr, kstep,  &
         sGamPPt, sGamKpKpt  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, maxIterNewt, iterCounterN, ierr, igrn, globalIncr, kstep
      real*8  e_kk, dtime, tolerNewt
      real*8  tau(NVECS), d_vec(NVECS), estran(NVECS), estran_n(NVECS)
      !real*8  kappa(NKAPP), gamdot(MAX_SLIP), drot(DIMS,DIMS)
      real*8  slipr(NKAPP), backst(NKAPP), backstc(NKAPP), gamdot(MAX_SLIP), gamdotC(MAX_SLIP), drot(DIMS,DIMS)
      real*8  fCeiDevHat(NVECS,NVECS), fCeDevVolHat(NVECS)
      real*8  pHatVec(NVECS, MAX_SLIP), ppTHat(NVECS,NVECS, MAX_SLIP)
      real*8  KpHatVec(NVECS, MAX_SLIP), KpHatdevVec(NVECS, MAX_SLIP), KpKpTHat(NVECS, NVECS, MAX_SLIP)
      real*8  matProp(NPROPS, MAX_SLIP)
      real*8  rss(MAX_SLIP), rssc(MAX_SLIP)

      integer subIncr, totSubIncrs, is
      real*8  xm_o(MAX_SLIP), tmp
      real*8  tau_o(NVECS), tau_save(NVECS)
      real*8  sGamPPt(NVECS,NVECS), sGamKpKpt(NVECS, NVECS)
!
!---------------------------------------------------------------------72
!
!------- counters for sub-increments / backup some variables
!
      subIncr = 1
      totSubIncrs = 1
      do is = 1, numslip
        matProp(3,is)= 0.5148d-18
        xm_o(is) = matProp(3,is)
      enddo
      call EqualTensors(tau, tau_o, NVECS)
!
!------- integrate crystal constitutive equations
!
      !write(*,*) 'strsoldev 1'
      !call SolveDeviatoric(tau, d_vec, estran, estran_n, kappa,  &
      call SolveDeviatoric(tau, d_vec, estran, estran_n, slipr, backst, backstc, &
                 gamdot, rss, gamdotC, rssc,  &
                 drot, fCeiDevHat, fCeDevVolHat,  &
                 pHatVec, ppTHat,  &
                 KpHatVec, KpHatdevVec, KpKpTHat,  &
                 matProp, e_kk, dtime, tolerNewt, numslip,  &
                 maxIterNewt, iterCounterN, ierr, igrn, globalIncr, kstep,  &
                 sGamPPt, sGamKpKpt )
      !write(*,*) 'strsoldev 2'
!
!------- if converged -> return, else do -> continuation method in "m"
!------- NOTE: here, continuation improves the initial guess 
!-------       for the solution variables 'tau'
!
      if (ierr .eq. XTAL_CONVERGED) return
!      write(XTAL_O, '(15x, A10, 1x, i3, A10, 1x, i3, A10, 1x, i3)'), 'kstep=', kstep, 'incr=', globalIncr, 'ingn=', igrn
      call WriteMessage  &
           (XTAL_O, 'DriverStressSolveDev: using Contin. for F0!')

!
!------- loop for sub-incrementation procedure
!
      do while (.true.)

         !write(*,*) 'Checker'
!
!---------- if not converged, increase # of sub-increments 
         if (ierr .ne. XTAL_CONVERGED) then
            !write(*,*) 'Whack'
            subIncr = 2 * subIncr - 1
            totSubIncrs = 2 * totSubIncrs
            if (totSubIncrs .gt. 16) then
               call WriteMessage(XTAL_O,  &
                        'DriverStressSolveDev: totSubIncrs > 16')
               do is = 1, numslip
                  matProp(3,is) = xm_o(is)
                  !call BoundForArgPowLaw(matProp(3,is))
               enddo
               return
            endif
            if (subIncr .eq. 1)  &
                 call EqualTensors(tau_o, tau, NVECS)
            if (subIncr .gt. 1)  &
                 call EqualTensors(tau_save, tau, NVECS)
!
!---------- if converged, adjust subincrements
         else if (subIncr .lt. totSubIncrs) then
            if ( (subIncr/2*2) .eq. subIncr) then
               subIncr = subIncr / 2 + 1
               totSubIncrs = totSubIncrs / 2
            else
               subIncr = subIncr + 1
            endif
!
!---------- successful return for continuation method
         else
            call WriteMessage(XTAL_O, 'Contin. for F0 successful !!!')
            return
         endif
!
!---------- trial rate sensitivity coefficient m
         tmp = real(subIncr) / real(totSubIncrs)   
         do is = 1, numslip
            !write(*,*) 'Whack',tmp
            matProp(3,is) = xm_o(is) / tmp
            !write(*,*) 'Jimmy'
            !call BoundForArgPowLaw(matProp(3,is))
         enddo
!
!---------- save current convergent solution before getting next solution
         if (subIncr .gt. 1 .and. ierr .eq. XTAL_CONVERGED)  &
                       call EqualTensors(tau, tau_save, NVECS)
!
!---------- report sub-incrs status
         write(XTAL_O, 1000) subIncr, totSubIncrs, matProp(3,1)
!
!---------- integrate crystal constitutive equations
!        iterCounterN = 0
         ierr = XTAL_CONVERGED
         !call SolveDeviatoric(tau, d_vec, estran, estran_n, kappa,  &
         !write(*,*) 'In'
         call SolveDeviatoric(tau, d_vec, estran, estran_n, slipr, backst, backstc, &
                 gamdot, rss, gamdotC, rssc,  &
                 drot, fCeiDevHat, fCeDevVolHat,  &
                 pHatVec, ppTHat, KpHatVec, KpHatdevVec, KpKpTHat,  &
                 matProp, e_kk, dtime, tolerNewt, numslip,  &
                 maxIterNewt, iterCounterN, ierr, igrn, globalIncr, kstep,  &
                 sGamPPt, sGamKpKpt)

      enddo

1000  format(3x, 'subInc/totSubIncrs = ', i2, '/', i3, 5x,  &
            'm(1)=', e12.5)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SolveDeviatoric(  &
         !tau, d_vec, estran, estran_n, kappa, gamdot, rss, drot,  &
         tau, d_vec, estran, estran_n, slipr, backst, backstc,  &
         gamdot, rss, gamdotC, rssc,  &
         drot, fCeiDevHat, fCeDevVolHat,  &
         pHatVec, ppTHat, KpHatVec, KpHatdevVec, KpKpTHat,  &
         matProp, e_kk, dtime, toler, numslip,  &
         maxIters, iterCounterN, ierr, igrn, globalIncr, kstep,  &
         sGamPPt, sGamKpKpt  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip, maxIters, iterCounterN, ierr, igrn, globalIncr, kstep
      real*8  e_kk, dtime, toler
      real*8  tau(NVECS), d_vec(NVECS), estran(NVECS), estran_n(NVECS)
      !real*8  kappa(NKAPP), gamdot(MAX_SLIP), drot(DIMS,DIMS)
      real*8  slipr(NKAPP), backst(NKAPP), backstc(NKAPP), gamdot(MAX_SLIP), gamdotC(MAX_SLIP)
      real*8  fCeiDevHat(NVECS,NVECS), fCeDevVolHat(NVECS), drot(DIMS,DIMS)
      real*8  pHatVec(NVECS, MAX_SLIP), ppTHat(NVECS,NVECS, MAX_SLIP)
      real*8  KpHatVec(NVECS, MAX_SLIP), KpHatdevVec(NVECS, MAX_SLIP), KpKpTHat(NVECS,NVECS, MAX_SLIP)
      real*8  matProp(NPROPS, MAX_SLIP)
      real*8  rss(MAX_SLIP), rssc(MAX_SLIP)

      integer iters, i, j
      real*8  rhs_norm_0, rhs_norm, search
      real*8  estranHat(NVECS), estranStar(NVECS)
      real*8  deltau(NVECS), tau0(NVECS)
      real*8  rhs(NVECS), lhs(NVECS,NVECS)
      real*8  dqr5x5(NVECS,NVECS)

      logical Converged
!      real*8  InnerProductVec
      real*8  NormOfVector, tol
      integer numnum, lwa, iflag, info
      real*8, allocatable :: wa(:)
      real*8  sGamPPt(NVECS,NVECS), sGamKpKpt(NVECS,NVECS)
!
!---------------------------------------------------------------------72
!
!------- predictor elastic strain in Btilde: dR*e_n*dR^T + D*dt 
!
      numnum = NVECS 
      lwa = (numnum*(3*numnum+13))/2
      allocate(wa(lwa))
      tol = toler/100.0
      
      call RotMat5x5ForSymm(drot, dqr5x5, DIMS)
      
      call MultAxu(dqr5x5, estran_n, estranHat, NVECS)
      
!------- Eq (3.31), page 29 (3. Numerical Implementation)      
      call AddTensors(pone, estranHat, dtime, d_vec, estranStar, NVECS)
      
      !write(*,*) 'tau:',tau
      call hybrd1(fcn,numnum,tau,rhs,tol,info,wa,lwa)
      !write(*,*) 'hybrd_rhs',info,rhs
      deallocate(wa)
!
!------- compute initial residual   
!
      !write(*,*) 'soldev 1'
      !write(*,*) 'Check KpHatdevVec for CompResidual, OK'
	  !do i = 1, numslip
      	!write(*,*) (KpHatdevVec(j,i),j=1,NVECS)
      !enddo
      call ComputeResidual(rhs, tau, estran, estranStar, fCeiDevHat,  &
                           !fCeDevVolHat, rss, gamdot, kappa, pHatVec,  &
                           fCeDevVolHat,  &
                           rss, gamdot, rssc, gamdotC,  &
                           slipr, backst, backstc, &
                           pHatVec, KpHatVec, KpHatdevVec,  &
                           matProp, e_kk, dtime, numslip)
!      rhs_norm_0 = dsqrt(InnerProductVec(rhs, rhs, NVECS))
      rhs_norm_0 = NormOfVector(rhs, NVECS)

!
!------- initialize flags and start Newton iterations for stress
!
      !write(*,*) 'soldev 2'
      iters = 0
      do while(.not.Converged(rhs, toler) .and. iters .lt. maxIters)
         !call xit()
         iters = iters + 1
         call EqualTensors(tau, tau0, NVECS)
!
!---------- compute local jacobian
         !call ComputeJacobian(lhs, rss, kappa, ppTHat, fCeiDevHat,  &
         !write(*,*) 'soldev 3'
         call ComputeJacobian(lhs, rss, rssc,  &
                              slipr, backst, backstc,  &
                              ppTHat, KpKpTHat,  &
                              fCeiDevHat, matProp, dtime, numslip,  &
                              sGamPPt, sGamKpKpt)
!
!---------- solve for the increment of stress
         call lsolve(lhs, rhs, 1, ierr, NVECS)
         if (ierr .eq. XTAL_SING_JACOBIAN) then
            write(XTAL_O, '(10x, A10, 1x, i3, A10, 1x, i3, A10, 1x, i3)'), &
                        'kstep=', kstep, 'incr=', globalIncr, 'ingn=', igrn
            call WriteWarning(XTAL_O,  &
                   'StressSolveDeviatoric: Jacobian is singular')
            return
         endif

         call EqualTensors(rhs, deltau, NVECS)

!
!---------- update stresses
         search = pone
         call AddTensors(pone, tau0, search, deltau, tau, NVECS)
!
!---------- compute new residual
         !write(*,*) 'soldev 4'
         call ComputeResidual(rhs, tau, estran, estranStar, fCeiDevHat,  &
                              !fCeDevVolHat, rss, gamdot, kappa, pHatVec,  &
                              fCeDevVolHat,  &
                              rss, gamdot, rssc, gamdotC,  &
                              slipr, backst, backstc, &
                              pHatVec, KpHatVec, KpHatdevVec,  &
                              matProp, e_kk, dtime, numslip)
!         rhs_norm = dsqrt(InnerProductVec(rhs, rhs, NVECS))
         rhs_norm = NormOfVector(rhs, NVECS)
         !write(*,*) 'soldev 5'
!
!---------- simple line search
         do while ( rhs_norm .gt. rhs_norm_0 )

            search = search*0.5d0
            if (search .lt. TINY) then
               call WriteWarning(XTAL_O,  &
                      'StressSolveDeviatoric: LS Failed, search < TINY')
                      !write(*,*) 'Crap'
               ierr = XTAL_LS_FAILED
               return
            endif
            call AddTensors(pone, tau0, search, deltau, tau, NVECS)

            call ComputeResidual(rhs, tau, estran, estranStar,  &
                                 fCeiDevHat, fCeDevVolHat,  &
                                 rss, gamdot, rssc, gamdotC,  &
                                 slipr, backst, backstc, &
                                 pHatVec, KpHatVec, KpHatdevVec,  &
                                 matProp, e_kk, dtime, numslip)

            rhs_norm = NormOfVector(rhs, NVECS)

         enddo
         !write(*,*) 'soldev 9'
         !write(*,*) 'rhs',rhs
         !write(*,*) 'tol',toler
!
!---------- update norm of residual
         !write(*,*) 'norm:', rhs_norm
         !write(*,*) 'toler:', toler
         !write(*,*) 'rhs:', rhs
         rhs_norm_0 = rhs_norm
      enddo
      !write(*,*) 'soldev 10'
!
!------- keep track of max number of newton iterations
!
      iterCounterN = max(iterCounterN, iters)
      !write(*,*) 'soldev 11'
!
!------- check convergence on iters < maxIters
!
      if (iters .ge. maxIters) then
         call WriteWarning(XTAL_O,  &
                      'StressSolveDeviatoric: iters > maxIters')
         ierr = XTAL_MAX_ITERS_HIT
         !write(*,*) 'soldev max iters', maxiters
         return
      endif
      !write(*,*) 'soldev 12'

      return

      CONTAINS

          SUBROUTINE fcn(n,tautau,rhsrhs,iflag)

              integer :: n, iflag
              real*8 :: tautau(n), rhsrhs(n)

              call ComputeResidual(rhsrhs, tautau, estran, estranStar, fCeiDevHat,  &
                  				   fCeDevVolHat,  &
			                       rss, gamdot, rssc, gamdotC,  &
                 				   slipr, backst, backstc, &
                 				   pHatVec, KpHatVec, KpHatdevVec,  &
                 				   matProp, e_kk, dtime, numslip)

              RETURN

          END SUBROUTINE FCN

      END
!

!=====================================================================72
!
      SUBROUTINE ComputeResidual(  &
         rhs, tau, estran, estranStar, fCeiDevHat, fCeDevVolHat,  &
         !gamdot, kappa, pHatVec, matProp, e_kk, dtime, numslip  &
         rss, gamdot, rssc, gamdotC,  &
         slipr, backst, backstc, &
         pHatVec, KpHatVec, KpHatdevVec,  &
         matProp, e_kk, dtime, numslip)
!
      use     SlipOverStress
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip
      real*8  e_kk, dtime
      real*8  rhs(NVECS), tau(NVECS), estran(NVECS), estranStar(NVECS)
      real*8  fCeiDevHat(NVECS,NVECS), fCeDevVolHat(NVECS)
      !real*8  rss(MAX_SLIP), gamdot(MAX_SLIP), kappa(NKAPP)
      real*8  rss(MAX_SLIP), gamdot(MAX_SLIP), slipr(NKAPP), backst(NKAPP), backstc(NKAPP)
      real*8  rssc(MAX_SLIP), gamdotC(MAX_SLIP)
      real*8  pHatVec(NVECS,MAX_SLIP), KpHatVec(NVECS,MAX_SLIP), KpHatdevVec(NVECS,MAX_SLIP)
      real*8  matProp(NPROPS,MAX_SLIP)

      integer is
      real*8  tmp(NVECS)
      !real*8  InnerProductVec, SSKineticEqn
      real*8  InnerProductVec, SSKineticEqn2, SSKineticEqn3

      real*8  crss, bss, bssc
!
!---------------------------------------------------------------------72
!
!------- elastic strain: -[Cei]({tau}-ekk*{fCeDevVol})+{e*} = -{e}+{e*}
!		 Page 30 (3. Numerical implementation, Ma & Marin 2006)
!
      call AddTensors(pone, tau, -e_kk, fCeDevVolHat, tmp, NVECS)
      
      call MultAxu(fCeiDevHat, tmp, estran, NVECS)
      
      call AddTensors(-pone, estran, +pone, estranStar, rhs, NVECS)
!
!------- resolve shear stresses and shear strain rates
!
      do is = 1, numslip
         rss(is) = InnerProductVec(tau, pHatVec(1,is), NVECS)
!         crss = kappa(1) + overstress(is)
         !crss = kappa(is)
         crss = slipr(is)
         bss = backst(is)
         gamdot(is) = SSKineticEqn2(rss(is), crss, bss, matProp(1,is),  &
                                   kGAMDOT)
         !if (dabs(gamdot(is)) > 1e-12) then
            !write(*,*) 'numslip',is
            !write(*,*) 'rss',rss(is)
            !write(*,*) 'slipr',slipr(is)
            !write(*,*) 'backst',backst(is)
            !write(*,*) 'gamdot',gamdot(is)
         !end if
!------- Climb quantities
		 bssc = backstc(is)
	     rssc(is) = InnerProductVec(tau, KpHatdevVec(1,is), NVECS)
	     gamdotC(is) = SSKineticEqn3(rssc(is), bssc, matProp(1,is), kGAMDOTc)
      enddo

!
!------- residual
!
      do is = 1, numslip
         call AddScaledTensor(-dtime*gamdot(is), pHatVec(1,is), rhs, NVECS)
         call AddScaledTensor(-dtime*gamdotC(is), KpHatVec(1,is), rhs, NVECS)
      enddo

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE ComputeJacobian(  &
         !lhs, rss, kappa, ppTHat, fCeiDevHat, matProp, dtime, numslip  &
         lhs, rss, rssc,  &
         slipr, backst, backstc, &
         ppTHat, KpKpTHat,  &
         fCeiDevHat, matProp, dtime, numslip,  &
         sGamPPt, sGamKpKpt  &
         )
!
      use     SlipOverStress
      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer numslip
      real*8  dtime
      !real*8  lhs(NVECS,NVECS), rss(MAX_SLIP), kappa(NKAPP)
      real*8  lhs(NVECS,NVECS), rss(MAX_SLIP), rssc(MAX_SLIP)
      real*8  slipr(NKAPP), backst(NKAPP), backstc(NKAPP)
      real*8  ppTHat(NVECS,NVECS,MAX_SLIP), KpKpTHat(NVECS,NVECS,MAX_SLIP)
      real*8  fCeiDevHat(NVECS,NVECS)
      real*8  matProp(NPROPS,MAX_SLIP)

      integer is
      real*8  dGamdTau(MAX_SLIP), dGamCdTauC(MAX_SLIP)
      !real*8  SSKineticEqn
      real*8  SSKineticEqn2, SSKineticEqn3

      real*8  crss, bss, bssc
      real*8  sGamPPt(NVECS,NVECS), sGamKpKpt(NVECS,NVECS)

!
!---------------------------------------------------------------------72
!
      do is = 1, numslip
!      
!------- d(gamdot)/d(tau)
!
         !crss = kappa(is)
         crss = slipr(is)
         bss = backst(is)
!         crss = kappa(1) + overstress(is)
         dGamdTau(is) = SSKineticEqn2(rss(is), crss, bss, matProp(1,is),  &
                                     kdGAMdTAU)
         !if (dabs(dGamdTau(is)) > 1e-12) then
            !write(*,*) 'numslip',is
            !write(*,*) 'rss',rss(is)
            !write(*,*) 'slipr',slipr(is)
            !write(*,*) 'backst',backst(is)
            !write(*,*) 'dGamdTau',dGamdTau(is)
         !end if
!         
!------- d(gamdotC)/d(tauC)
!
		 bssc = backstc(is)
         dGamCdTauC(is) = SSKineticEqn3(rssc(is), bssc, matProp(1,is), kdGAMcdTAUc)
   
      enddo
!
!------- jacobian
!
      call SetTensor(sGamPPt, pzero, NVECS*NVECS)
      call SetTensor(sGamKpKpt, pzero, NVECS*NVECS)

      do is = 1, numslip
         call AddScaledTensor(dtime*dGamdTau(is), ppTHat(1,1,is),  &
                              sGamPPt, NVECS*NVECS)
         call AddScaledTensor(dtime*dGamCdTauC(is), KpKpTHat(1,1,is),  &
                              sGamKpKpt, NVECS*NVECS)
      enddo
      call AddTensors(pone, fCeiDevHat, pone, sGamPPt+sGamKpKpt, lhs, NVECS*NVECS)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE StressSolveVolumetric(  &
         statev, fCeDevVolHat, fCeVol, estran, e_kk, dtime  &
         )

      implicit none
      include 'params_xtal.inc'

      real*8  fCeVol, e_kk, dtime
      real*8  fCeDevVolHat(NVECS), estran(NVECS), statev(NSTAV)

      real*8  InnerProductVec
!
!---------------------------------------------------------------------72
!
      statev(kEVOL) = e_kk
      statev(kPRES) = InnerProductVec(fCeDevVolHat, estran, NVECS)  &
                      + fCeVol * e_kk
      
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SetTensor(  &
         tensor, value, n  &
         )
      implicit none

      integer n
      real*8  value
      real*8  tensor(n)

      integer i
!
!---------------------------------------------------------------------72
!
      do i = 1, n
         tensor(i) = value
      enddo

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SymmetrizeTensor(  &
         tensor, symmTensor, n  &
         )

      implicit none
      include 'params_xtal.inc'

      integer  n
      real*8   tensor(n*n), symmTensor(n*n)
!
!---------------------------------------------------------------------72
!
      if ((n .ne. 2) .and. (n .ne. 3))  &
         call RunTimeError(XTAL_O, 'SkewSymmetrizeTensor: n =! 2 or 3')

      if (n .eq. 2) then
         symmTensor(1) = tensor(1)
         symmTensor(4) = tensor(4)
         symmTensor(3) = 0.5*(tensor(3) + tensor(2))
         symmTensor(2) = symmTensor(3)
      else   ! n = 3
         symmTensor(1) = tensor(1)                     ! 11
         symmTensor(5) = tensor(5)                     ! 22
         symmTensor(9) = tensor(9)                     ! 33
         symmTensor(4) = 0.5*(tensor(4) + tensor(2))   ! 12
         symmTensor(7) = 0.5*(tensor(7) + tensor(3))   ! 13
         symmTensor(8) = 0.5*(tensor(8) + tensor(6))   ! 23
         symmTensor(2) = symmTensor(4)                 ! 21
         symmTensor(3) = symmTensor(7)                 ! 31
         symmTensor(6) = symmTensor(8)                 ! 32
      endif

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SkewSymmetrizeTensor(  &
         tensor, skewTensor, n  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer  n
      real*8   tensor(n*n), skewTensor(n*n)
!
!---------------------------------------------------------------------72
!
      if ((n .ne. 2) .and. (n .ne. 3))  &
         call RunTimeError(XTAL_O, 'SkewSymmetrizeTensor: n =! 2 or 3')

      call SetTensor(skewTensor, pzero, n*n)

      if (n .eq. 2) then
         skewTensor(3) = 0.5*(tensor(3) - tensor(2))
         skewTensor(2) = - skewTensor(3)
      else   ! n = 3
         skewTensor(4) = 0.5*(tensor(4) - tensor(2))
         skewTensor(7) = 0.5*(tensor(7) - tensor(3))
         skewTensor(8) = 0.5*(tensor(8) - tensor(6))
         skewTensor(2) = - skewTensor(4)
         skewTensor(3) = - skewTensor(7)
         skewTensor(6) = - skewTensor(8)
      endif

      return
      END
!

!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE Mat3x3ToVec5x1Symm(  &
         matrix, vector, n  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
     
      integer  n
      real*8   matrix(n*n), vector(5)
!
!---------------------------------------------------------------------72
!
!------- matrix is deviatoric
!
      if (n .ne. 3)  &
         call RunTimeError(XTAL_O, 'Mat3x3ToVect5x1Symm: n =! 3')

      vector(1) = (matrix(1) - matrix(5)) / sqr2
      vector(2) = matrix(9) * sqr32
      vector(3) = matrix(2) * sqr2
      vector(4) = matrix(3) * sqr2
      vector(5) = matrix(6) * sqr2
 
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE Vec5x1ToMat3x3Symm(  &
         vector, matrix, n  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'
     
      integer n
      real*8  matrix(n*n), vector(5)
!
!---------------------------------------------------------------------72
!
!------- matrix is deviatoric
!
      if (n .ne. 3)  &
         call RunTimeError(XTAL_O, 'Vect5x1ToMat3x3Symm: n =! 3')

      matrix(1) = 0.5 * (sqr2 * vector(1) - sqr23 * vector(2))
      matrix(5) =-0.5 * (sqr2 * vector(1) + sqr23 * vector(2))
      matrix(9) = vector(2) * sqr23
      matrix(2) = vector(3) / sqr2
      matrix(3) = vector(4) / sqr2
      matrix(6) = vector(5) / sqr2
      matrix(4) = matrix(2)
      matrix(7) = matrix(3)
      matrix(8) = matrix(6)
 
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE Mat3x3ToVec3x1Skew(  &
         matrix, vector, n  &
         )

      implicit none
      include 'params_xtal.inc'
     
      integer n
      real*8  matrix(n*n), vector(3)
!
!---------------------------------------------------------------------72
!
      if (n .ne. 3)  &
         call RunTimeError(XTAL_O, 'Vect5x1ToMat3x3Skew: n =! 3')
!
!------- uses tensor components below main diagonal
!
      vector(1) = matrix(2)
      vector(2) = matrix(3)
      vector(3) = matrix(6)
 
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE Vec3x1ToMat3x3Skew(  &
         vector, matrix, n  &
         )

      implicit none
      include 'params_xtal.inc'
     
      integer  n
      real*8   matrix(n*n), vector(3)
!
!---------------------------------------------------------------------72
!
      if (n .ne. 3)  &
         call RunTimeError(XTAL_O, 'Vect5x1ToMat3x3Skew: n =! 3')

      matrix(1) = 0.
      matrix(5) = 0.
      matrix(9) = 0.
      matrix(2) = vector(1)
      matrix(3) = vector(2)
      matrix(6) = vector(3)
      matrix(4) = - vector(1)
      matrix(7) = - vector(2)
      matrix(8) = - vector(3)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE RotMat5x5ForSymm(  &
         cmatrix, qr5x5, n  &
         )

      implicit 	none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer  n
      real*8   cmatrix(n*n), qr5x5(5, 5)

      real*8   c11, c21, c31, c12, c22, c32, c13, c23, c33
!
!---------------------------------------------------------------------72
!
!----- construct 5x5 rotation matrix for symm 2nd order tensors: 
!         [A_sm]=[C][A_lat][C]' <=>  {A_sm} = [qr5x5]{A_lat}
!----- with: {A}={()/sqr2,sqr32*(),sqr2*(),sqr2*(),sqr2*()}

      if (n .ne. 3)  &
         call RunTimeError(XTAL_O, 'RotMat5x5ForSymm: n =! 3')

      c11 = cmatrix(1)
      c21 = cmatrix(2)
      c31 = cmatrix(3)
      c12 = cmatrix(4)
      c22 = cmatrix(5)
      c32 = cmatrix(6)
      c13 = cmatrix(7)
      c23 = cmatrix(8)
      c33 = cmatrix(9)

      qr5x5(1, 1)  =  0.5d0 * (c11 * c11 - c12 * c12 -  &
                                                 c21 * c21 + c22 * c22)
      qr5x5(1, 2)  =  sqr3 / 2.d0 * (c13 * c13 - c23 * c23)
      qr5x5(1, 3)  =  c11 * c12 - c21 * c22
      qr5x5(1, 4)  =  c11 * c13 - c21 * c23
      qr5x5(1, 5)  =  c12 * c13 - c22 * c23
      qr5x5(2, 1)  =  sqr3 / 2.d0 * (c31 * c31 - c32 * c32)
      qr5x5(2, 2)  =  1.5d0 * c33 * c33 - 0.5d0
      qr5x5(2, 3)  =  sqr3 * c31 * c32
      qr5x5(2, 4)  =  sqr3 * c31 * c33
      qr5x5(2, 5)  =  sqr3 * c32 * c33
      qr5x5(3, 1)  =  c11 * c21 - c12 * c22
      qr5x5(3, 2)  =  sqr3 * c13 * c23
      qr5x5(3, 3)  =  c11 * c22 + c12 * c21
      qr5x5(3, 4)  =  c11 * c23 + c13 * c21
      qr5x5(3, 5)  =  c12 * c23 + c13 * c22
      qr5x5(4, 1)  =  c11 * c31 - c12 * c32
      qr5x5(4, 2)  =  sqr3 * c13 * c33
      qr5x5(4, 3)  =  c11 * c32 + c12 * c31
      qr5x5(4, 4)  =  c11 * c33 + c13 * c31
      qr5x5(4, 5)  =  c12 * c33 + c13 * c32
      qr5x5(5, 1)  =  c21 * c31 - c22 * c32
      qr5x5(5, 2)  =  sqr3 * c23 * c33
      qr5x5(5, 3)  =  c21 * c32 + c22 * c31
      qr5x5(5, 4)  =  c21 * c33 + c23 * c31
      qr5x5(5, 5)  =  c22 * c33 + c23 * c32

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE RotMat3x3ForSkew(  &
         cmatrix, qr3x3, n  &
         )

      implicit  none
      include 'params_xtal.inc'

      integer  n
      real*8   cmatrix(n*n), qr3x3(3, 3)

      real*8   c11, c21, c31, c12, c22, c32, c13, c23, c33
!
!---------------------------------------------------------------------72
!
!---- Construct 3x3 rotation matrix for skew 2nd order tensors
!        [W_sm]=[C][W_lat][C]'  <=>  {W_sm} = [qr3x3]{W_lat}

      if (n .ne. 3)  &
         call RunTimeError(XTAL_O, 'RotMat3x3ForSkew: n =! 3')

      c11 = cmatrix(1)
      c21 = cmatrix(2)
      c31 = cmatrix(3)
      c12 = cmatrix(4)
      c22 = cmatrix(5)
      c32 = cmatrix(6)
      c13 = cmatrix(7)
      c23 = cmatrix(8)
      c33 = cmatrix(9)

      qr3x3(1, 1) = c22 * c11 - c21 * c12
      qr3x3(1, 2) = c23 * c11 - c21 * c13
      qr3x3(1, 3) = c23 * c12 - c22 * c13
      qr3x3(2, 1) = c32 * c11 - c31 * c12
      qr3x3(2, 2) = c33 * c11 - c31 * c13
      qr3x3(2, 3) = c33 * c12 - c32 * c13
      qr3x3(3, 1) = c32 * c21 - c31 * c22
      qr3x3(3, 2) = c33 * c21 - c31 * c23
      qr3x3(3, 3) = c33 * c22 - c32 * c23

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE MultAxu(  &
         A, u, v, n  &
         )

      implicit none
      include 'numbers.inc'

      integer n
      real*8  A(n,n), u(n), v(n)

      integer i, j
!
!---------------------------------------------------------------------72
!
      call SetTensor(v, pzero, n)

      do i = 1, n
         do j = 1, n
            v(i) = v(i) + A(i,j) * u(j)
         end do
      end do

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE MultATxu(  &
         A, u, v, n  &
         )

      implicit none
      include 'numbers.inc'

      integer n
      real*8  A(n,n), u(n), v(n)

      integer i, j
!
!---------------------------------------------------------------------72
!
      call SetTensor(v, pzero, n)

      do i = 1, n
         do j = 1, n
            v(i) = v(i) + A(j,i) * u(j)
         end do
      end do

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE MultAxB(  &
         A, B, C, n  &
         )

      implicit none
      include 'numbers.inc'
     
      integer n
      real*8  A(n, n), B(n, n), C(n, n)

      integer i, j, k
!
!---------------------------------------------------------------------72
!
      call SetTensor(C, pzero, n*n)

      do i = 1, n 
         do j = 1, n
            do k = 1, n
               C(i,j) = C(i,j) + A(i,k) * B(k,j)
            enddo
         enddo
      enddo
   
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE MultATxB(  &
         A, B, C, n  &
         )

      implicit none
      include 'numbers.inc'
     
      integer n
      real*8  A(n, n), B(n, n), C(n, n)

      integer i, j, k
!
!---------------------------------------------------------------------72
!
      call SetTensor(C, pzero, n*n)

      do i = 1, n 
         do j = 1, n
            do k = 1, n
               C(i,j) = C(i,j) + A(k,i) * B(k,j)
            enddo
         enddo
      enddo
   
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE MultAxBT(  &
         A, B, C, n  &
         )

      implicit none
      include 'numbers.inc'
     
      integer n
      real*8  A(n, n), B(n, n), C(n, n)

      integer i, j, k
!
!---------------------------------------------------------------------72
!
      call SetTensor(C, pzero, n*n)

      do i = 1, n 
         do j = 1, n
            do k = 1, n
               C(i,j) = C(i,j) + A(i,k) * B(j,k)
            enddo
         enddo
      enddo
   
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE DeviatoricTensor(  &
         tensor, tensorDev, n  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer n
      real*8  tensor(n*n), tensorDev(n*n)

      real*8  TraceOfTensor, tensHyd
!
!---------------------------------------------------------------------72
!
      if (n .ne. 3)  &
         call RunTimeError(XTAL_O, 'DeviatoricTensor: n =! 3')

      tensHyd = TraceOfTensor(tensor, n) / pthree

      tensorDev(1) = tensor(1) - tensHyd
      tensorDev(5) = tensor(5) - tensHyd
      tensorDev(9) = tensor(9) - tensHyd

      tensorDev(2) = tensor(2)
      tensorDev(3) = tensor(3)
      tensorDev(6) = tensor(6)

      tensorDev(4) = tensor(4)
      tensorDev(7) = tensor(7)
      tensorDev(8) = tensor(8)
      
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      real*8 FUNCTION TraceOfTensor(  &
         tensor, n  &
         )

      implicit none
      include 'params_xtal.inc'

      integer n
      real*8  tensor(n*n)
!
!---------------------------------------------------------------------72
!
      if (n .ne. 3)  &
         call RunTimeError(XTAL_O, 'TraceOfTensor: n =! 3')

      TraceOfTensor = tensor(1) + tensor(5) + tensor(9)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      real*8 FUNCTION ScalarProduct(  &
         tensA, tensB, n  &
         )

      implicit none
      include 'params_xtal.inc'
      include 'numbers.inc'

      integer n
      real*8  tensA(n*n), tensB(n*n)
!
!---------------------------------------------------------------------72
!
      if (n .ne. 3)  &
         call RunTimeError(XTAL_O, 'ScalarProduct: n =! 3')

      ScalarProduct = tensA(1)*tensB(1) + tensA(2)*tensB(2) +  &
                      tensA(3)*tensB(3) + tensA(4)*tensB(4) +  &
                      tensA(5)*tensB(5) + tensA(6)*tensB(6) +  &
                      tensA(7)*tensB(7) + tensA(8)*tensB(8) +  &
                      tensA(9)*tensB(9)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE OuterProductVec(  &
         vecU, vecV, outer, n  &
         )

      implicit none
      include 'numbers.inc'

      integer n
      real*8  vecU(n), vecV(n)
      real*8  outer(n, n)
     
      integer i, j
!
!----------------- dyadic product-------------------------------------72
!
      do i = 1, n
         do j = 1, n
            outer(i,j) = vecU(i) * vecV(j)
         enddo
      enddo

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      real*8 FUNCTION InnerProductVec(  &
         vecU, vecV, n  &
         )

      implicit none
      include 'numbers.inc'

      integer n
      real*8  vecU(n), vecV(n)

      integer i
!
!---------------------------------------------------------------------72
!
!---- compute inner (dot) product of two vectors: product = {u}^T*{v}
!
      InnerProductVec = pzero
      do i = 1, n
         InnerproductVec = InnerProductVec + vecU(i) * vecV(i)
      enddo

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE CrossProductVec(  &
         vecU, vecV, cross, n  &
         )

      implicit none
      include 'numbers.inc'

      integer n
      real*8  vecU(n), vecV(n)
      real*8  cross(n)
!
!----------------- cross product in 3D -------------------------------72
!
      cross(1) = vecU(2) * vecV(3) - vecU(3) * vecV(2)
      cross(2) = vecU(3) * vecV(1) - vecU(1) * vecV(3)
  	cross(3) = vecU(1) * vecV(2) - vecU(2) * vecV(1)
      
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE MultQAQT(  &
         Q, A, B, n  &
         )

      implicit none

      integer n
      real*8  Q(n,n), A(n,n), B(n,n)

      real*8  tmp(n,n)
!
!---------------------------------------------------------------------72
!
      call MultAxB(Q, A, tmp, n)       ! tmp = Q*A
      call MultAxBT(tmp, Q, B, n)      ! B = tmp*QT

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE MultQTAQ(  &
         Q, A, B, n  &
         )

      implicit none

      integer n
      real*8  Q(n,n), A(n,n), B(n,n)

      real*8  tmp(n,n)
!
!---------------------------------------------------------------------72
!
      call MultATxB(Q, A, tmp, n)      ! tmp = QT*A
      call MultAxB(tmp, Q, B, n)       ! B = tmp*Q

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE EqualTensors(  &
         tensA, tensB, n  &
         )

      implicit none

      integer n
      real*8  tensA(n), tensB(n)

      integer i
!
!---------------------------------------------------------------------72
!
      do i = 1, n
         tensB(i) = tensA(i)
      enddo

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE AddTensors(  &
         coef_A, tensA, coef_B, tensB, tensC, n  &
         )

      implicit none

      integer n
      real*8  coef_A, coef_B
      real*8  tensA(n), tensB(n), tensC(n)

      integer i
!
!---------------------------------------------------------------------72
!
       do i = 1, n
          tensC(i) = coef_A * tensA(i) + coef_B * tensB(i) 
       enddo

       return
       END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE SetToScaledTensor(  &
         fac_A, tensA, tensB, n  &
         )

      implicit none

      integer n
      real*8  fac_A
      real*8  tensA(n), tensB(n)

      integer i
!
!---------------------------------------------------------------------72
!
      do i = 1, n
          tensB(i) = fac_A * tensA(i)
      enddo

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE AddScaledTensor(  &
         fac_A, tensA, tensB, n  &
         )

      implicit none

      integer n
      real*8  fac_A
      real*8  tensA(n), tensB(n)

      integer i
!
!---------------------------------------------------------------------72
!
      do i = 1, n
          tensB(i) = tensB(i) + fac_A * tensA(i)
      enddo

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      real*8 FUNCTION DetOfMat3x3(  &
         tensor  &
         )

      implicit none
      
      real*8  tensor(3, 3)

      real*8  det11, det12, det13, det21, det22, det23
!
!---------------------------------------------------------------------72
!
!  Determinant of a second order tensor (3x3 matrix)
 
      det11 = tensor(1, 1) * tensor(2, 2) * tensor(3, 3)
      det12 = tensor(1, 2) * tensor(2, 3) * tensor(3, 1)
      det13 = tensor(1, 3) * tensor(2, 1) * tensor(3, 2)
      det21 = tensor(1, 1) * tensor(2, 3) * tensor(3, 2)
      det22 = tensor(1, 2) * tensor(2, 1) * tensor(3, 3)
      det23 = tensor(1, 3) * tensor(2, 2) * tensor(3, 1)
      
      DetOfMat3x3 = det11 + det12 + det13 - det21 - det22 - det23

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE InverseOfMat3x3(  &
         tensor, tensor_inv, tensor_det  &
         )

      implicit 	none

      real*8  tensor_det
      real*8  tensor(3, 3), tensor_inv(3, 3)

      real*8  tinv11, tinv12, tinv13
      real*8  tinv21, tinv22, tinv23
      real*8  tinv31, tinv32, tinv33
     
      real*8  DetOfMat3x3
!
!---------------------------------------------------------------------72
!
      tensor_det = DetOfMat3x3(tensor)

      tinv11 = tensor(2, 2) * tensor(3, 3) - tensor(2, 3) * tensor(3, 2)
      tinv12 = tensor(1, 3) * tensor(3, 2) - tensor(1, 2) * tensor(3, 3)
      tinv13 = tensor(1, 2) * tensor(2, 3) - tensor(1, 3) * tensor(2, 2)
      tinv21 = tensor(2, 3) * tensor(3, 1) - tensor(2, 1) * tensor(3, 3)
      tinv22 = tensor(1, 1) * tensor(3, 3) - tensor(1, 3) * tensor(3, 1)
      tinv23 = tensor(1, 3) * tensor(2, 1) - tensor(1, 1) * tensor(2, 3)
      tinv31 = tensor(2, 1) * tensor(3, 2) - tensor(2, 2) * tensor(3, 1)
      tinv32 = tensor(3, 1) * tensor(1, 2) - tensor(1, 1) * tensor(3, 2)
      tinv33 = tensor(1, 1) * tensor(2, 2) - tensor(1, 2) * tensor(2, 1)

      tensor_inv(1, 1) = tinv11 / tensor_det
      tensor_inv(1, 2) = tinv12 / tensor_det
      tensor_inv(1, 3) = tinv13 / tensor_det
      tensor_inv(2, 1) = tinv21 / tensor_det
      tensor_inv(2, 2) = tinv22 / tensor_det
      tensor_inv(2, 3) = tinv23 / tensor_det
      tensor_inv(3, 1) = tinv31 / tensor_det
      tensor_inv(3, 2) = tinv32 / tensor_det
      tensor_inv(3, 3) = tinv33 / tensor_det

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      real*8 FUNCTION SignOf(  &
         value  &
         )

      implicit none

      real*8  value
!
!---------------------------------------------------------------------72
!
!------- compute the sign of a number
!
      SignOf = 1.0
      if (value .lt. 0.0) SignOf = -1.0

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      real*8 FUNCTION Power(  &
         x, y  &
         )

      implicit none

      real*8  x, y
!
!---------------------------------------------------------------------72
!
!---- evaluates  x^y
!
      if (x .eq. 0.0) then
         if (y .gt. 0.0) then
            Power = 0.d0
         elseif (y .lt. 0.0) then
            Power = 1.d+300
         else
            Power = 1.d0
         endif
      else
         Power = y * log10(dabs(x))
         if (Power .gt. 300.0) then
            Power = 1.d+300
         else
            Power = 10.d0 ** Power
         endif
         if (x .lt. 0.0) Power = -Power
      endif

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE UnitVector(  &
         vector, unitV, n  &
         )
      implicit none
      include 'params_xtal.inc'

      integer n
      real*8  vector(n), unitV(n)

      real*8  magnitude, InnerProductVec
!
!---------------------------------------------------------------------72
!
      magnitude = dsqrt(InnerProductVec(vector, vector, n))
      if (magnitude .le. 0.d0)  &
              call RunTimeError(XTAL_O, 'UnitVector: magnitude <= 0.0')
      call SetToScaledTensor(1./magnitude, vector, unitV, n)

      return
      END
!
!=====================================================================72
!
!=====================================================================72
!
      integer FUNCTION IndexMaxAbsValueOfVec(  &
         vector, n  &
         )

      implicit none

      integer n
      real*8  vector(n)

      integer indexVal, i
      real*8  maxValue
!
!---------------------------------------------------------------------72
!
      indexVal = 1
      maxValue = dabs(vector(1))
      do i = 2, n
         if (dabs(vector(i)) .gt. maxValue) then
            indexVal = i
            maxValue = dabs(vector(i))
         endif
      enddo

      IndexMaxAbsValueOfVec = indexVal

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      real*8 FUNCTION MaxAbsValueOfVec(  &
         vector, n  &
         )

      implicit none

      integer n
      real*8  vector(n)

      integer i
      real*8  maxValue
!
!---------------------------------------------------------------------72
!
      maxValue = dabs(vector(1))
      do i = 2, n
         if (dabs(vector(i)) .gt. maxValue) maxValue = dabs(vector(i))
      enddo

      MaxAbsValueOfVec = maxValue

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE MatrixInverse(  &
         a, n, np, det  &
         )
!
!---- routine borrowed from Hammid Youssef
!---- input:  a   : nonsym matrix to invert
!----         n   : order of matrix to invert
!----         np  : dimension of matrix in calling routine
!---- output: a   : matrix inverse
!----         det : determinant of matrix
!---- local:  k   : VECTEUR DE TRAVAIL ENTIER 
!
      IMPLICIT DOUBLE  PRECISION (A-H,O-Z)
      DIMENSION A(NP,NP),K(6)
      DATA ZERO,UN,EPS/0.D0,1.D0,1.D-13/
!
!---------------------------------------------------------------------72
!
!-------  Initialization
!
      DET=UN
      DO I=1,N
         K(I)=I
      END DO
!
!-------  inversion de la matrix A
!
      DO II=1,N
!
!-------  search for non-zero pivot in column II
!
         DO I=II,N
            XPIV=A(I,II)
            IF(DABS(XPIV).GT.EPS) GO TO 10
         END DO
         DET=ZERO
         GOTO 1000
!
!-------  exchange rows II and I
!
 10      DET=DET*XPIV
         IF(I.EQ.II) GO TO 20
         I1=K(II)
         K(II)=K(I)
         K(I)=I1
         DO J=1,N
            C=A(I,J)
            A(I,J)=A(II,J)
            A(II,J)=C
         END DO
         DET=-DET
!
!-------  normalize the row of the pivot
!
 20      C=UN/XPIV
         A(II,II)=UN
         DO J=1,N
            A(II,J)=A(II,J)*C
         END DO
!
!-------  Elimination
!
         DO I=1,N
            IF(I.EQ.II) GO TO 30
            C=A(I,II)
            A(I,II)=ZERO
            DO J=1,N
               A(I,J)=A(I,J)-C*A(II,J)
            END DO
 30      END DO
      END DO
!
!-------  re-order columns of inverse
!
      DO J=1,N
!-------  search J1 until K(J1)=J
         DO J1=J,N
            JJ=K(J1)
            IF(JJ.EQ.J) GO TO 100
         END DO
100      IF(J.EQ.J1) GO TO 110
!-------  exchange columns J and J1
         K(J1)=K(J)
         DO I=1,N
            C=A(I,J)
            A(I,J)=A(I,J1)
            A(I,J1)=C
         END DO
110   END DO

1000  CONTINUE

      return
      END
!
!=====================================================================72
!
!     
!=====================================================================72
!     
      SUBROUTINE lsolve(a, x, nrhs, mystatus, n)
!     
!---- routine borrowed from Paul Dawson (Cornell University)
!---- Solve symmetric positive definite 5x5 linear systems. By calling 
!----  with the identity as right hand sides, you can use this routine
!----  to invert the matrix.
!
!     a    -- "n x n" matrix                                   (in/out)     
!     x    -- "n x nrhs" matrix, the array of right hand sides (in/out)
!     nrhs -- number of right hand sides                       (in)
!     mystatus -- return status                                (out)
!     n    -- size of system (must be equal to 5)
!     
      implicit none
      include 'params_xtal.inc'

      integer nrhs, mystatus, n
      real*8  a(n,n), x(n,*)
     
      integer i
      real*8  v1, v2, v3, v4, sum, not_zero, amax
!     
!---------------------------------------------------------------------72
!     
!------- check size of matrix
!
      if (n .ne. 5)  &
         call RunTimeError(XTAL_O, 'lsolve : n .ne. 5')
!     
!------- set the scale according to the largest entry of a.
!     
      amax = 0.0d0
      amax = max(a(1,1), a(2,1), a(3,1), a(4,1), a(5,1))
      amax = max(a(2,2), a(3,2), a(4,2), a(5,2), amax)
      amax = max(a(3,3), a(4,3), a(5,3), amax)
      amax = max(a(4,4), a(5,4), amax)
      amax = max(a(5,5), amax)
     
      not_zero = amax * TINY
!     
!------- A = LDL'.
!
!------- At each step, check the size of each pivot
!---------- j = 1.
      if (a(1,1) .LT. not_zero) then
        mystatus = XTAL_SING_JACOBIAN
        return
      endif
     
      sum = 1.0/a(1,1)

      a(2,1) = a(2,1) * sum
      a(3,1) = a(3,1) * sum
      a(4,1) = a(4,1) * sum
      a(5,1) = a(5,1) * sum
!
!---------- j = 2.
      v1 = a(2,1)*a(1,1)
      a(2,2) = a(2,2) - a(2,1)*v1
      
      a(3,2) = a(3,2)-a(3,1)*v1
      a(4,2) = a(4,2)-a(4,1)*v1
      a(5,2) = a(5,2)-a(5,1)*v1

      if (a(2,2) .LT. not_zero) then
        mystatus = XTAL_SING_JACOBIAN
        return
      endif

      sum = 1.0/a(2,2)

      a(3,2) = a(3,2) * sum
      a(4,2) = a(4,2) * sum
      a(5,2) = a(5,2) * sum
!
!---------- j = 3.
      v1 = a(3,1)*a(1,1)
      v2 = a(3,2)*a(2,2)
      sum = a(3,1)*v1 + a(3,2)*v2

      a(3,3) = a(3,3) - sum
      
      if (a(3,3) .LT. not_zero) then
        mystatus = XTAL_SING_JACOBIAN
        return
      endif

      sum = 1.0/a(3,3)

      a(4,3) = a(4,3)-a(4,1)*v1-a(4,2)*v2
      a(5,3) = a(5,3)-a(5,1)*v1-a(5,2)*v2
      a(4,3) = a(4,3) * sum
      a(5,3) = a(5,3) * sum
!
!---------- j = 4.
      v1 = a(4,1)*a(1,1)
      v2 = a(4,2)*a(2,2)
      v3 = a(4,3)*a(3,3)
      sum = a(4,1)*v1+a(4,2)*v2+a(4,3)*v3

      a(4,4) = a(4,4) - sum

      if (a(4,4) .LT. not_zero) then
        mystatus = XTAL_SING_JACOBIAN
        return
      endif

      a(5,4) = a(5,4)-a(5,1)*v1-a(5,2)*v2-a(5,3)*v3
      a(5,4) = a(5,4)/a(4,4)
!
!---------- j = 5.
      v1 = a(5,1)*a(1,1)
      v2 = a(5,2)*a(2,2)
      v3 = a(5,3)*a(3,3)
      v4 = a(5,4)*a(4,4)
      sum = a(5,1)*v1+a(5,2)*v2+a(5,3)*v3+a(5,4)*v4

      a(5,5) = a(5,5) - sum

      if (a(5,5) .LT. not_zero) then
        mystatus = XTAL_SING_JACOBIAN
        return
      endif
!     
!------- solve for RHS
!
      do i=1, nrhs
!
!---------- Ly=b. 
        x(2, i) = x(2, i)  &
           - a(2,1)*x(1, i)
        x(3, i) = x(3, i)  &
           - a(3,1)*x(1, i) - a(3,2)*x(2, i)
        x(4, i) = x(4, i)  &
           - a(4,1)*x(1, i) - a(4,2)*x(2, i) - a(4,3)*x(3, i)
        x(5, i) = x(5, i)  &
           - a(5,1)*x(1, i) - a(5,2)*x(2, i) - a(5,3)*x(3, i)  &
           - a(5,4)*x(4, i)
!
!---------- Dz=y.
        x(1, i) = x(1, i)/a(1,1)
        x(2, i) = x(2, i)/a(2,2)
        x(3, i) = x(3, i)/a(3,3)
        x(4, i) = x(4, i)/a(4,4)
        x(5, i) = x(5, i)/a(5,5)
!
!---------- L'x=z.
        x(4, i) = x(4, i)  &
           - a(5,4)*x(5, i)
        x(3, i) = x(3, i)  &
           - a(4,3)*x(4, i) - a(5,3)*x(5, i)
        x(2, i) = x(2, i)  &
           - a(3,2)*x(3, i) - a(4,2)*x(4, i) - a(5,2)*x(5, i)
        x(1, i) = x(1, i)  &
           - a(2,1)*x(2, i) - a(3,1)*x(3, i) - a(4,1)*x(4, i)  &
           - a(5,1)*x(5, i)

      enddo

      return
      END
!     
!=====================================================================72
!     
!
!=====================================================================72
!
      SUBROUTINE Vec5x1ToVec6x1(  &
         vec5x1, vec6x1  &
         )

      implicit none
      include 'numbers.inc'

      real*8  vec5x1(5), vec6x1(6)
!
!---------------------------------------------------------------------72
!
      vec6x1(1) = 0.5 * (sqr2*vec5x1(1) - sqr23*vec5x1(2)) 
      vec6x1(2) =-0.5 * (sqr2*vec5x1(1) + sqr23*vec5x1(2)) 
      vec6x1(3) = vec5x1(2) * sqr23
      vec6x1(4) = vec5x1(3) / sqr2
      vec6x1(5) = vec5x1(4) / sqr2
      vec6x1(6) = vec5x1(5) / sqr2

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE MultAxB_G(  &
         A, B, C, m1, m2, m3  &
         )

      implicit none
      include 'numbers.inc'
     
      integer m1, m2, m3
      real*8  A(m1, m2), B(m2, m3), C(m1, m3)

      integer i, j, k
!
!---------------------------------------------------------------------72
!
!---- Dimensions: A is m1 x m2; B is m2 x m3; C is m1 x m3 
!
      call SetTensor(C, pzero, m1*m3)

      do i = 1, m1 
         do j = 1, m3
            do k = 1, m2
               C(i,j) = C(i,j) + A(i,k) * B(k,j)
            enddo
         enddo
      enddo
   
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE MultAxu_G(  &
         A, u, v, m1, m2  &
         )

      implicit none
      include 'numbers.inc'

      integer m1, m2
      real*8  A(m1,m2), u(m2), v(m1)

      integer i, j
!
!---------------------------------------------------------------------72
!
!---- Dimensions: A is m1 x m2; u is m2; v is m1
!
      call SetTensor(v, pzero, m1)

      do i = 1, m1
         do j = 1, m2
            v(i) = v(i) + A(i,j) * u(j)
         end do
      end do

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE MultATxu_G(  &
         A, u, v, m1, m2  &
         )

      implicit none
      include 'numbers.inc'

      integer m1, m2
      real*8  A(m1,m2), u(m1), v(m2)

      integer i, j
!
!---------------------------------------------------------------------72
!
!---- Dimensions: A is m2 x m1; u is m1; v is m2
!
      call SetTensor(v, pzero, m2)

      do i = 1, m2
         do j = 1, m1
            v(i) = v(i) + A(j,i) * u(j)
         end do
      end do

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE OuterProductVec_G(  &
         vecU, vecV, outer, m1, m2  &
         )

      implicit none
      include 'numbers.inc'

      integer m1, m2
      real*8  vecU(m1), vecV(m2)
      real*8  outer(m1, m2)
     
      integer i, j
!
!---------------------------------------------------------------------72
!
      do i = 1, m1
         do j = 1, m2
            outer(i,j) = vecU(i) * vecV(j)
         enddo
      enddo

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      real*8 FUNCTION NormOfVector(  &
         vec,  n  &
         )

      implicit none
      include 'numbers.inc'

      integer n
      real*8  vec(n)

      real*8  max_vec
      real*8  MaxAbsValueOfVec, InnerProductVec
!
!---------------------------------------------------------------------72
!
      max_vec = MaxAbsValueOfVec(vec, n)
      if (max_vec .gt. pzero)  &
          call SetToScaledTensor(pone/max_vec, vec, vec, n)

      NormOfVector = dsqrt(InnerProductVec(vec, vec, n))

      if (max_vec .gt. pzero) then
         NormOfVector = max_vec*NormOfVector
         call SetToScaledTensor(max_vec, vec, vec, n)
      endif

      return
      END
!
!=====================================================================72
!
!
!
!=====================================================================72
!
      SUBROUTINE AnglesToRotMatrix(  &
         angle, crot, n  &
         )

      implicit none

      integer n
      real*8  angle(n)
      real*8  crot(n,n)

      real*8  sps, cps, sth, cth, sph, cph
!
!---------------------------------------------------------------------72
!
!------- Construct [C] matrix from euler angles (in radians)
!-------     {a}_sm = [C] {a}_xtal
!
      sps = dsin(angle(1))
      cps = dcos(angle(1))
      sth = dsin(angle(2))
      cth = dcos(angle(2))
      sph = dsin(angle(3))
      cph = dcos(angle(3))

      crot(1,1) = -sps * sph - cps * cph * cth
      crot(2,1) =  cps * sph - sps * cph * cth
      crot(3,1) =  cph * sth
      crot(1,2) =  cph * sps - sph * cps * cth
      crot(2,2) = -cps * cph - sps * sph * cth
      crot(3,2) =  sph * sth
      crot(1,3) =  cps * sth
      crot(2,3) =  sps * sth
      crot(3,3) =  cth

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE RotMatrixToAngles(  &
         crot, angle, n  &
         )

      implicit none

      integer n
      real*8  angle(n)
      real*8  crot(n,n)
      
      real*8  sth
!
!---------------------------------------------------------------------72
!
!------- compute euler angles from [C] (in radians)
!
      angle(2) = acos(crot(3,3))
      if (dabs(crot(3,3)) .ne. 1.0) then
         sth = dsin(angle(2))
         angle(1) = datan2(crot(2,3)/sth, crot(1,3)/sth)
         angle(3) = datan2(crot(3,2)/sth, crot(3,1)/sth)
      else
         angle(1) = 0.
         angle(3) = datan2(-crot(1,2), -crot(1,1))
      endif

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
!=====================================================================72
!
      recursive FUNCTION SSKineticEqn2(  &
         rss, crss, bss, matProp, kflag  &
         ) result (sskine)

      implicit none
      include 'params_xtal.inc'

      integer kflag, flaginner, flagoutter
      real*8  rss, crss, bss, rssc
      real*8  matProp(NPROPS)

      real*8 F0, p, q, gam0, tau0, theta, kappa, mu, mu0
      real*8 argc
      real*8 SignOf, Power
      real*8 arg, p1, p2, eps
      real*8 sskine, SSKineticEqn, SmoothMB

!
!---------------------------------------------------------------------72
! 
      flaginner=1
      flagoutter=2
!------- recover power law material constants
!
      F0 =  matProp(3)
      p = 0.18d0
      q = 1.63d0
      gam0 = 1.44d-3
      tau0 = 268.2d06
      theta = 1223.d0
      kappa = 1.3806488d-23
      mu = 79.35d9
      mu0 = 265.3d9

      if (kflag .eq. kGAMDOT) then 
        arg = (DABS(rss-bss)-crss*(mu/mu0))/(tau0*mu/mu0)

!     Evaluated the inner Macaulay bracket
        arg=SmoothMB(arg, flaginner, p, q)
        arg = 1.d0-arg

!     Evaluated the outter Macaulay bracket
        arg=SmoothMB(arg, flagoutter, p, q)

!	  Flow rule included dislocation glide gamma
        sskine = gam0*DEXP(-arg*F0/(kappa*theta))*SignOf(rss-bss)
        
      else if (kflag .eq. kdGAMdTAU) then
        eps = epsilon(0.0d0)
        eps = max(rss*eps,eps)
        p1 = SSKineticEqn2(rss+0.5*eps, crss, bss, matProp, kGAMDOT)
        p2 = SSKineticEqn2(rss-0.5*eps, crss, bss, matProp, kGAMDOT)
        sskine = (p1-p2)/eps
      end if

      return
      END
!
!=====================================================================72
!
!	Formulation of gamdotC(is)
!
!=====================================================================72
!
      recursive FUNCTION SSKineticEqn3(  &
         rssc, bssc, matProp, kflag  &
         ) result (sskineC)

      implicit none
      include 'params_xtal.inc'

      integer kflag
      real*8  rssc, bssc
      real*8  matProp(NPROPS)

      real*8 F0, gam0, theta, kappa
      real*8 pc, tau0c, argc, tauch
      real*8 SignOf, Power
      real*8 p1, p2, eps
      real*8 sskineC, SSKineticEqn, SmoothMB

!
!---------------------------------------------------------------------72
!      
!------- Flow rule parameters need to be calibrated for climb mechanism
!
      F0 =  matProp(3)
      gam0 = 1.44d-3
      theta = 1223.d0
      kappa = 1.3806488d-23
      
!------- Parameters need to be calibrated      
!        creep exponent accepted for thermal creep (Nabarro, Phil. Mag. 16, 1967, p231)
	  pc = 3d0
      tau0c = 0.0073d06 ! (Pa)

      if (kflag .eq. kGAMDOTc) then 
        argc = (DABS(rssc-bssc)/tau0c)
        argc = argc**pc

!	  Flow rule relevant to dislocation climp gamdotC
        sskineC = gam0*DEXP(-F0/(kappa*theta))*argc*SignOf(rssc-bssc)
                
      else if (kflag .eq. kdGAMcdTAUc) then
        eps = epsilon(0.0d0)
        eps = max(rssc*eps,eps)
        p1 = SSKineticEqn3(rssc+0.5*eps, bssc, matProp, kGAMDOTc)
        p2 = SSKineticEqn3(rssc-0.5*eps, bssc, matProp, kGAMDOTc)
        sskineC = (p1-p2)/eps
      end if

      return
      END
!
!=====================================================================72
!
     real*8 FUNCTION SmoothMB(fun, flag, p, q)
!
     implicit none
     
     integer flag
     real*8 fun
     real*8 arclen, tmp, p, q, parpow

!--- This function is using a arc that are tangenc with the original Macaulay brackets
!    at point (-c, 0) and (sqrt(2)/2*c, sqrt(2)/2*c) to substitute the original function 
!    between this two points to smoothen the original functions
        
!    the inner Macaulay bracket = 0~O(10e9), while the outer =0~O(10)
!    different c is chosen 
     
     arclen=0.d0
     select case (flag)

     case(1)        ! inner MB  selected
        arclen=0.001d0
        parpow=p
        !write(*,*) 'Inner ', fun
     case(2)
        arclen=0.001d0   ! outter MB selected
        parpow=q
        !write(*,*) 'Outter ', fun
     case default 
         write(*,*) 'Error flag used when calling function SmoothMB'
         stop
     end select     

     if (fun .le. -arclen) then
            fun=0.d0
            SmoothMB=0.d0
     elseif ((fun .gt. -arclen)  .and. (fun .le. sqrt(2.d0)/2.d0*arclen)) then
            tmp=tan(22.5d0/180.d0*3.141592653d0)  
            fun=arclen/tmp-sqrt(arclen*arclen*(1.d0/tmp/tmp-1.d0)-fun*fun-2.d0*fun*arclen)
            if (fun .lt. epsilon(0.d0)) then 
                SmoothMB=0.d0
            else
                SmoothMB=fun**parpow
            endif
     else
            SmoothMB=fun**parpow
     endif

     return
     end FUNCTION SmoothMB
      
!=====================================================================72
!
!
!
!=====================================================================72
!
!=====================================================================72
!
      SUBROUTINE RunTimeError(  &
         io, message  &
         )

      implicit none
      
      character*(*) message
      integer io
!
!---------------------------------------------------------------------72
!
      write(io, 1000) message
      write(*, 1000) message
!      call closfl( )
      call CrystalCloseIOFiles( )
      call CrystalCloseIOFiles_2( )
      stop

1000  format(/,'***ERROR Message: '/, 3x, a)

      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE WriteWarning(  &
         io, message  &
         )

      implicit none
      
      character*(*) message
      integer io
!
!---------------------------------------------------------------------72
!
      write(io, 1000) message

1000  format(/,'***WARNING Message: '/, 3x, a)

      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE WriteMessage(  &
         io, message  &
         )

      implicit none
      
      character*(*) message
      integer io
!
!---------------------------------------------------------------------72
!
      write(io, 1000) message

1000  format('***Message: ', a)

      return
      END
!
!=====================================================================72
!


SUBROUTINE uel (rhs, amatrx, svars, energy, ndofel, nrhs,  &
    nsvars, props, nprops, coords, mcrd, nnode, u, du, v, a, jtype,  &
    time, dtime, kstep, kinc, jelem, params, ndload, jdltyp,  &
    adlmag, predef, npredf, lflags, mlvarx, ddlmag, mdload, pnewdt,  &
    jprops, njprop, period)

 INCLUDE 'ABA_PARAM.INC'


     DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*), &
      SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL), &
      DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*), &
      JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*), &
      PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)


      DIMENSION ds1(3),ds2(3),dn(3),trac(mcrd,nrhs),  &
        trac_jacob(mcrd,mcrd),r(mcrd,mcrd),coord_l(mcrd,nnode),  &
        gp_coord(2),sf(3),b(mcrd,ndofel),co_de_m(3,3),  &
        b_t(ndofel,mcrd), transformation_m(ndofel,ndofel),  &
        transformation_m_t(ndofel,ndofel),temp1(mcrd,ndofel)

       DIMENSION stiff_l(ndofel,ndofel),temp2(ndofel,ndofel),  &
         stiff_g(ndofel,ndofel),residual_l(ndofel,nrhs),  &
         residual_g(ndofel,nrhs),ajacob_m(2,3),delu_loc_gp(mcrd), &
         co_de(mcrd,nnode),SVARSi(9)

    REAL*8 :: alpha,sk_n,delc,f_m,f_n,c_r,c_p,c_C,T_c,Y_0,AA,   &
     & sk_c 
    character*100 filename,eleid,prefix,suffix,gpn 


	n_GP = 3
	gp_w = 1d0/3d0

!=======================================================================
!	Initialize Matrices and Vectors
!======================================================================= 

	! shear and normal local opening displacments
	CALL k_vector_zero(ds1, 3)
	CALL k_vector_zero(ds2, 3)
	CALL k_vector_zero(dn, 3)

	CALL k_matrix_zero(trac, mcrd, nrhs)
	CALL k_matrix_zero(trac_jacob, mcrd, mcrd)
	CALL k_matrix_zero(r, mcrd, mcrd)
	CALL k_matrix_zero(coord_l, mcrd, nnode)
	CALL k_vector_zero(gp_coord, 2)

	! Shape functions
	CALL k_vector_zero(sf, 3)

	CALL k_matrix_zero(transformation_m,ndofel,ndofel)
	CALL k_matrix_zero(transformation_m_t,ndofel,ndofel)

	CALL k_matrix_zero(b,mcrd,ndofel)
	CALL k_matrix_zero(b_t,ndofel,mcrd)
	CALL k_matrix_zero(temp1,mcrd,ndofel)
	CALL k_matrix_zero(stiff_l,ndofel,ndofel)
	CALL k_matrix_zero(temp2,ndofel,ndofel)
	CALL k_matrix_zero(stiff_g,ndofel,ndofel)
	CALL k_matrix_zero(residual_l,ndofel,nrhs)
	CALL k_matrix_zero(residual_g,ndofel,nrhs)
	CALL k_matrix_zero(ajacob_m,2,3)
	CALL k_matrix_zero(rhs,ndofel,nrhs)
	CALL k_matrix_zero(amatrx,ndofel,ndofel)
	CALL k_matrix_zero(co_de,mcrd,nnode)

	a_jacob = 0.d0
	Mflag = 0

! open(unit=2000,file='~/1820/Damage.dat',status='unknown',access='sequential')


      
!======================================================================= 
!	Do local computations
!======================================================================= 

	do i = 1, mcrd
  		do j = 1, nnode
    		co_de(i,j)=coords(i,j)+u(3.0*(j-1.0)+i)
  		enddo
  	enddo

!	convert to local coordinates
	CALL k_local_coordinates(co_de, r, coord_l, transformation_m,  &
	    transformation_m_t, a_jacob, ajacob_m, coords, u, ndofel, nnode, mcrd)

!=======================================================================    
! 	Compute shear and normal local opening displacments
!======================================================================= 

	do i = 1, 3
  		ds1(i)=coord_l(1,i+3)-coord_l(1,i)
		ds2(i)=coord_l(2,i+3)-coord_l(2,i)
  		dn(i) =coord_l(3,i+3)-coord_l(3,i)
	enddo

!=======================================================================    
!	Do Calculations at Gauss Points
!======================================================================= 

	DO i = 1, n_GP
  
!		determine the values of the shape function at each GP
  		CALL k_shape_fun(i, sf)

!		initialize 1 normal & 2 shear separation
  		CALL k_vector_zero(delu_loc_gp, mcrd)
  
! 		determine 2 shear and 1 normal opening displamenets at each GP
  		do j = 1, 3
    		delu_loc_gp(1) = delu_loc_gp(1) + ds1(j)*sf(j)
		    delu_loc_gp(2) = delu_loc_gp(2) + ds2(j)*sf(j)
		    delu_loc_gp(3) = delu_loc_gp(3) + dn(j)*sf(j)
  		enddo
  
  
! 		initialize the variables Svars
! 		each GP has 9 vars => 27 vars/3 GPs
  		Svarsi(1) = Svars(9*(i-1)+1)	! D_f
		Svarsi(2) = Svars(9*(i-1)+2) 	! D_c
  		Svarsi(3) = Svars(9*(i-1)+3) 	! Y
  		Svarsi(4) = Svars(9*(i-1)+4) 	! U_t
	    Svarsi(5) = Svars(9*(i-1)+5) 	! U_s
  		Svarsi(6) = Svars(9*(i-1)+6) 	! U_n
  		Svarsi(7) = Svars(9*(i-1)+7) 	! T_t
  		Svarsi(8) = Svars(9*(i-1)+8) 	! T_s
  		Svarsi(9) = Svars(9*(i-1)+9) 	! T_n

		!if ((JELEM .eq. 1598) .and. (i .eq. 1)) then
		!	write(*,*) 'CHZ #',JELEM,'Gause point', i
		 !   write(*,*) 'Un^(t+dt)=', delu_loc_gp(3)
		  !  write(*,*) 'Un^(t)=', svarsi(6)		    
		!endif
  
  
! 		determine traction vector and tangent modulus matrix
   		call k_cohesive_law(Trac, Trac_Jacob, delu_loc_gp, mcrd, nrhs, &
     & 			Svarsi, PROPS, DTIME, Mflag)
     
!		Print out results of one element in files stored in the folder "508"
  		if ((JELEM .eq. 150000)) then
  		!if ( (JELEM .eq. 1640) .or. (JELEM .eq. 1613) .or. (JELEM .eq. 1657) .or. (JELEM .eq. 1641) ) then
  		!if ((JELEM .eq. 15116) .or. (JELEM .eq. 15022) .or. (JELEM .eq. 15657) &
      !& .or. (JELEM .eq. 15115) .or. (JELEM .eq. 15656) .or. (JELEM .eq. 15053) .or. (JELEM .eq. 15071) ) then
  		!if ((JELEM .eq. 131861) .or. (JELEM .eq. 131860) .or. (JELEM .eq. 115536) .or. (JELEM .eq. 115535) .or. (JELEM .eq. 115534) &
  	   !& .or. (JELEM .eq. 131870)) then
      		if ((Svarsi(1) .ge. 0d0) .or. (Svarsi(2) .ge. 0d0)) then
  				write(JELEM*10+1000+i,'(I10,11e14.6)') &
      			!write(2000,'(2I10,10e14.6)') KINC,JELEM,TIME(2),      &
     &				KINC,TIME(2), &			           ! Kinc Time	! C1 C2
     & 				Svarsi(1), Svarsi(2), Svarsi(3), & ! Df Dc Y	! C3 C4 C5
     &				Svarsi(6), Svarsi(5), Svarsi(4), & ! Un Us Ut	! C6 C7 C8
     !& 				Svarsi(9), Svarsi(8), Svarsi(7), & ! Tn Ts Tt
     &		  sqrt( ((Svarsi(9)+abs(Svarsi(9)))/2d0)**2d0 + 2d0*Svarsi(8)**2d0 + 2d0*Svarsi(7)**2d0	)! Tdet
  			endif
                endif

  		if (Mflag .eq. 1) then 
  			write(*,*) JELEM
    		write(*,*) 'new step size is requested'
    		PNEWDT=0.5d0 ! restart
    		return        
  		endif                                         
                                                                        
!		store the variables
 		Svars(9*(i-1)+1) = Svarsi(1) 
    	Svars(9*(i-1)+2) = Svarsi(2) 
    	Svars(9*(i-1)+3) = Svarsi(3) 
    	Svars(9*(i-1)+4) = Svarsi(4) 
    	Svars(9*(i-1)+5) = Svarsi(5) 
    	Svars(9*(i-1)+6) = Svarsi(6) 
    	Svars(9*(i-1)+7) = Svarsi(7)
    	Svars(9*(i-1)+8) = Svarsi(8) 
    	Svars(9*(i-1)+9) = Svarsi(9) 
  
! 		determine B matrix and its transpose
  		CALL k_bmatrix(sf, b, mcrd, ndofel)
 		CALL k_matrix_transpose(b, b_t, mcrd, ndofel)
  
! 		Compute the stiffness matrix
! 		Local Stiffness = B_t * Trac_Jacob * B
   		CALL k_matrix_multiply(trac_jacob, b, temp1, mcrd, mcrd, ndofel)
  		CALL k_matrix_multiply(b_t, temp1, stiff_l, ndofel, mcrd, ndofel)
  
! 		Compute Global stiffness matrix
! 		Global_K = Transpose(T) * K * T
   		CALL k_matrix_multiply(transformation_m_t, stiff_l,  &
      		temp2, ndofel, ndofel, ndofel)
  		CALL k_matrix_multiply(temp2, transformation_m, stiff_g,  &
  		    ndofel, ndofel, ndofel)
  
! 		Multiply Jacobian with the Global stiffness and add contribution
! 		from each Gauss Point
   		a_mult = a_jacob*gp_w
  		CALL k_matrix_plus_scalar(amatrx, stiff_g, a_mult,  &
  			ndofel, ndofel)
  
! 		Compute the global residual vector
! 		Local_residual = B_t * Trac
! 		Global_residual = Transpose(T) * Local_residual

  		CALL k_matrix_multiply(b_t, trac, residual_l, ndofel,  &
	  		mcrd, nrhs)
  		CALL k_matrix_multiply(transformation_m_t, residual_l,  &
      			residual_g, ndofel, ndofel, nrhs)
  
! 		Multiply the Global residual by the Jacobian and add the
! 		contribution from each point
   		CALL k_matrix_plus_scalar(rhs, residual_g, a_mult,  &
   			ndofel, nrhs)

	END DO

RETURN
END SUBROUTINE uel
!=====================================================================


!=====================================================================
! subroutine: Determine the global displacement-separation (B) matrix
!=====================================================================
SUBROUTINE k_bmatrix(sf,b,mcrd,ndofel)
INCLUDE 'ABA_PARAM.INC'

	dimension sf(3),B(mcrd,ndofel)

	b(1,1) =  sf(1)
	b(1,4) =  sf(2)
	b(1,7) =  sf(3)
	b(1,10)= -sf(1)
	b(1,13)= -sf(2)
	b(1,16)= -sf(3)
	b(2,2) =  sf(1)
	b(2,5) =  sf(2)
	b(2,8) =  sf(3)
	b(2,11)= -sf(1)
	b(2,14)= -sf(2)
	b(2,17)= -sf(3)
	b(3,3) =  sf(1)
	b(3,6) =  sf(2)
	b(3,9) =  sf(3)
	b(3,12)= -sf(1)
	b(3,15)= -sf(2)
	b(3,18)= -sf(3)

	RETURN

END SUBROUTINE k_bmatrix
!=====================================================================

!=====================================================================
! subroutine: Cohesive law
!=====================================================================

subroutine k_cohesive_law(T,T_d,delu,mcrd,nrhs,Svarsi,PROPS,h,Mflag) 
                                                                        
      INCLUDE 'ABA_PARAM.INC' 
!       IMPLICIT REAL*8(A-H,O-Z)                                        
                                                                        
      dimension T(mcrd,nrhs),T_d(mcrd,mcrd),delu(mcrd),SVARSi(9),       &
     & U_v(3),T_v(3),dEdD(mcrd),dDfdU(mcrd),dDcdU(mcrd),dDdU(mcrd),     &
     & D_dot(mcrd,mcrd),Em(mcrd,mcrd),Res(2),Res_D(2,2),PROPS(14)       &
     & ,pivot(2),Res_Dv(2,2),rhs(2),Smatrix(2,2)                                            
                                                                        
       real*8 Df,Dc,Y,alpha,sk_n,sk_c,delc,f_m,f_n,c_r,       &
     &  c_p,c_C,T_c,Y_0,AA,As,h,T_d,T,delu,U_v,T_v,SVARSi,Res_Dv,       &
     &  ResDdet                                                         
                                                                        
       real*8 popn,popt,Y_cu,dY,Df_int,RDF,RDF_D,Df_cu,       &
     &  Dc_int,Dc_cu,RDc,RDc_D,D_cu,Un_inr,Us_inr,U_det,Y_crt,          &
     &  dYcrtdYcu,dEdD,dDfdU,dDcdU,dDdU,D_dot,Em,error,Res,Res_D        
                                                                        
       real*8 T_r                                            
	   real*8 residue_f, residue_c, jacob_f, jacob_c                      
	   real*8 jacob_fc, jacob_cf
	   integer ncount                                          
                                                                        

!=======================================================================
!	cohesive zone law parameters                                          
!=======================================================================
                                                                        
      alpha= PROPS(1) ! ratio k_n = alpha*k_t
      sk_n = PROPS(2)  ! initial normal stiffness k_n
      delc = PROPS(3)  ! critical opening
      f_m  = PROPS(4)   ! exponent (1 - w_f)^m of fatigue damage evolution
      f_n  = PROPS(5)   ! exponent in fatigue traction evolution
      c_r  = PROPS(6)   ! exponent in creep traction evolution
      c_p  = PROPS(7)   ! exponent (1 - w_cf)^p of creep damage evolution
      c_C  = PROPS(8)   ! An adjustable parameter in MPa
      T_c  = PROPS(9)   ! traction threshold for creep damage
      Y_0  = PROPS(10)  ! cohesive zone threshold 
      AA   = PROPS(11)  ! Parameter A controls the damage rate
      sk_c = PROPS(12) 	! high valued penalization coefficient k_c=10*k_n
      As   = PROPS(13)  ! A* not using
      T0   = PROPS(14)  ! traction threshold for fatigue damage  
     
  
	  !write(*,*) 'Y0,T_c',PROPS(10),PROPS(9) 
 
!=======================================================================
!	variables - used for recording last step values
!=======================================================================

      Df 	 = Svarsi(1) ! Df_cu scalar damage w_f variable for fatigue
      Dc 	 = Svarsi(2) ! Dc_cu scalar damage variable w_c for creep
      Y  	 = Svarsi(3) ! thermodynamic force associated with damage variable
      
      U_v(1) = Svarsi(4) ! Sliding separation 1
      U_v(2) = Svarsi(5) ! Sliding separation 2
      U_v(3) = Svarsi(6) ! Normal separation
      
      T_v(1) = Svarsi(7) ! Tangential traction 1
      T_v(2) = Svarsi(8) ! Tangential traction 2
      T_v(3) = Svarsi(9) ! Normal traction

!=======================================================================
!	Opening displacement in normal and shear
!=======================================================================
      popn = delu(3) 
      popt = sqrt( delu(1)**2d0 + delu(2)**2d0 ) 
      
      	                                                                                
!	Initialize                                                                         
      call k_matrix_zero(T_d, mcrd, mcrd) 
      call k_matrix_zero(Em, mcrd, mcrd) 

!=======================================================================
!   Thermaldynamic force for damage parameter Df_dot not zero
!   only when thermaldyanmic force is increasing
!	Eq. 5 in Bouvard et al 2009 ???
!=======================================================================
                                                                        
       Y_cu = sk_n/2d0/delc*( ((popn+abs(popn))/2d0)**2d0 + alpha*popt**2d0)
       
       !write(*,*) 'Y=', Y_cu
       
       dY = Y_cu - Y 
       
!=======================================================================
!   Df_cu is calculated using Newton Raphson
!=======================================================================
       
       ! Assume the initial value             
       Df_int = Df 
       Dc_int = Dc 
        
       ! Consider thermodynamic force increment > 0 or not
       If (dY .GE. 0d0) then 
                                                                        
        error=1d0 
                                                                        
        ncount=0 
                                                                        
        do while (error .GT. 1d-05)
        
        	Nflag=0    
                                                                      
        	RDc=residue_c(Dc_int,Dc,Df_int,Df,T_c,T_v,delu,c_C,c_r,c_p,     &
     & 			sk_n,delc,h,alpha,sk_c,Nflag)
     
        	RDc_D=jacob_c(Dc_int,Dc,Df_int,                                 &
     & 			Df,T_c,T_v,delu,c_r,c_p,c_C,sk_n,delc,h,alpha,sk_c)              
            
            !write(*,*) 'predict Dc',RDc, RDc_D 
                                                                        
        	Dc_cu = Dc_int - RDc/RDc_D
        	
        	RDc=residue_c(Dc_cu,Dc,Df_int,Df,T_c,T_v,delu,c_C,c_r,c_p,      &
     & 			sk_n,delc,h,alpha,sk_c,Nflag)                                          
        	!write(*,*) 'check Dc error, current Dc', abs(RDc), Dc_cu 
        	
        	error=abs(RDc) 
        
        	if (Nflag==1) then                                                                
        		ncount = 50
        	endif

        	ncount = ncount + 1 
        	
        	if (ncount .GT. 50) then 
        	!if ( (ncount .GT. 50) .and. (Dc_cu .le. 1d0) ) then 
        		write(*,*) 'An error has occured' 
        		Mflag = 1
        		return 
        	endif 
        	
        	if (Dc_cu .gt. 1d0) then
            	Dc_cu = 1d0
            endif
                    	    
            Dc_int = Dc_cu                                                
                                                                        
        enddo 
                                                                        
                                                                        
       else 
                                                                        
        ncount = 0 
        error = 1d0 
         
        do while (error .GT. 1D-5)
        !do while ( (error .GT. 1D-5) .and. (Dc_cu .le. 1d0) ) ! creep damage evolution
        
        	Nflag = 0
        	
        	RDc=residue_c(Dc_int,Dc,Df_int,Df,T_c,T_v,delu,c_C,c_r,c_p,     &
     & 			sk_n,delc,h,alpha,sk_c,Nflag)
     
        	RDc_D=jacob_c(Dc_int,Dc,Df_int,                                 &
     & 			Df,T_c,T_v,delu,c_r,c_p,c_C,sk_n,delc,h,alpha,sk_c)              
            
            !write(*,*) 'predict Dc',RDc, RDc_D 
                                                                        
        	Dc_cu = Dc_int - RDc/RDc_D
        	
        	RDc=residue_c(Dc_cu,Dc,Df_int,Df,T_c,T_v,delu,c_C,c_r,c_p,      &
     & 			sk_n,delc,h,alpha,sk_c,Nflag)                                          
        	!write(*,*) 'check Dc error, current Dc', abs(RDc), Dc_cu 
        	
        	error=abs(RDc) 
        
        	if (Nflag==1) then                                                                
        		ncount = 50
        	endif

        	ncount = ncount + 1 
        	
        	if (ncount .GT. 50) then 
        	!if ( (ncount .GT. 50) .and. (Dc_cu .le. 1d0) ) then 
        		write(*,*) 'An error has occured' 
        		Mflag = 1
        		return 
        	endif 
        	
        	if (Dc_cu .gt. 1d0) then
            	Dc_cu = 1d0
            endif
            
        	Dc_int = Dc_cu 
        	
        enddo 
                                                                        
        !write(*,*) 'final Dc',Dc_cu
                                                                               
        Dc_cu = Dc_int
                                                                        
      endif ! End of consideration for thermodynamic force increment
      
      ! Still keep the input Df at the beginning
      Df_cu = Df_int
                                        
      !write(*,*) 'Df_cu, Dc_cu', Df_cu, Dc_cu                              

!     damage parameter for current time step                            
      D_cu = Dc_cu !+ Df_cu
!      pause
                                                                  
!=========================================================================
!     local traction in the cohesive zone elements at the current time s
!     hard code for three dimension
                                                                        
      T(3,1) = (1d0-D_cu)*sk_n*(popn + abs(popn))/2d0/delc  & 
     &        + sk_c*(popn - abs(popn))/2d0/delc
      T(1,1) = alpha*(1d0 - D_cu)*sk_n*delu(1)/delc                           
      T(2,1) = alpha*(1d0 - D_cu)*sk_n*delu(2)/delc  
      
!=========================================================================      

    
!      write(*,*) 'Damage parameter'
!      write(*,*) Df_cu,Dc_cu,D_cu
!      write(*,*) 'Traction'
!      write(*,*) T(3,1),T(1,1),T(2,1)
!      write(*,*) 'Sep'
!      write(*,*) popn,delu(1),delu(2)
                
                                                                        
      dEdD(1)=-alpha*sk_n/delc*delu(1)                                  
      dEdD(2)=-alpha*sk_n/delc*delu(2)                                  
      
      if (popn .ge. 0d0) then                                               
      	dEdD(3)=-sk_n/delc*delu(3)                                        
      else                                                            
      	dEdD(3)=0
      endif                                                                     
                                                                        
      Tn_cu=T(3,1)     
      Tt_cu=T(1,1)
      Ts_cu=T(2,1)
      
!=====================================================================================      
      
      ! || T^{t+dt} + T^{t} ||
      T_det = sqrt((Tn_cu+T_v(3))**2d0+(Tt_cu+T_v(1))**2d0/alpha+(Ts_cu    &
     &      + T_v(2))**2d0/alpha)
     
      !T_crt = T_det/2d0/(1d0-(Df_cu+Df)/2d0)-T0
      T_crt = T_det/2d0/(1d0/2d0)-T0
      T_n = ((T_crt+abs(T_crt))/2d0)**(f_n)
      T_r = (((T_det/2d0-T_c)+abs(T_det/2d0-T_c))/2.0D0/c_C)**c_r
      
      ! Updated - Tung
      !dTrdTdet = c_r/((2D0*c_C)**c_r)/2d0
      !dTrdTdet = dTrdTdet*( (T_det - 2d0*T_c)  &
     !&                 + abs(T_det - 2d0*T_c) )**(c_r-1d0)
     
      if (T_det .eq. 0d0) then
      	T_detinv = 0d0
      else
      	T_detinv = 1d0/T_det
      endif
      
      Un_inr = (delu(3)-U_v(3))/delc
      Ut_inr = (delu(1)-U_v(1))/delc
      Us_inr = (delu(2)-U_v(2))/delc
      
      ! || U^{t+dt}-U^{t} ||
      U_det = sqrt(((Un_inr+abs(Un_inr))/2d0)**2d0+alpha*(Ut_inr)**2d0   &
     & 		+ alpha*(Us_inr)**2d0)
      
      if (U_det .eq. 0d0) then
      	U_detinv = 0d0
      else
      	U_detinv = 1d0/U_det
      endif
      
      ! Y^{t+dt}
      Y_cu = sk_n/2d0/delc*(((delu(3)+abs(delu(3)))/2d0)**2d0+alpha*delu(1)**2d0+alpha*delu(2)**2d0)
      
      ! Y^{t}
      Y_old = sk_n/2d0/delc*(((U_v(3)+abs(U_v(3)))/2d0)**2d0+alpha*U_v(1)**2d0+alpha*U_v(2)**2d0)         
      
      Y_crt = sqrt(Y_cu/2d0+Y_old/2d0)-sqrt(Y_0) 
      
	  ! <sqrt(Y^{t+dt}/2 + Y^{t}/2) - sqrt(Y0)>^n
      Y_crtn =( Y_crt/2d0 + abs(Y_crt)/2d0 )**f_n                                                
      
      if ( (Y_cu .ne. 0d0) .and. (Y_old .ne. 0d0) ) then
      	Ydet_inv = 1d0/sqrt(Y_cu/2d0+Y_old/2d0)
      else
      	Ydet_inv = 0d0
      endif
      
      
	  !write(*,*) Ydet_inv,Y_crt,Udet
      
      D_old = Dc !+ Df
      
      if ( 1d0-(D_cu+D_old)/2d0 .LE. 0d0) then
        Para1 = 0d0
      	Para2 = 0d0
      	Para3 = 0d0
      else
      	Para1 = AA*f_n*(1d0-(D_cu+D_old)/2d0)**(f_m)*((Y_crt+abs(Y_crt))/2d0) &
     & 	    **(f_n-1d0)*U_det*sk_n/delc/4d0*Ydet_inv
        Para2 = AA*(1d0-(D_cu+D_old)/2d0)**(f_m)*((Y_crt+abs(Y_crt))/2d0)**(f_n)
        Para3 = AA*f_m/2d0*(1d0-(D_cu+D_old)/2d0)**(f_m-1d0)*((Y_crt+abs(Y_crt))/2d0) &
     & 	    **(f_n)*U_det
      endif   
      
      !**********************************************************
      !	Tung's modification
      !**********************************************************    
      !if ( (T_det - 2d0*T_c) .le. 0d0 ) then
      !	T_r = 0d0
      !	dTrdTdet = 0d0
      !	Para4 = 0d0
      !	Para5 = 0d0
      !else
      !	dTrdTdet = c_r/((2D0*c_C)**c_r)/2d0
      !	dTrdTdet = dTrdTdet*( (T_det - 2d0*T_c)  &
     !& 	                 + abs(T_det - 2d0*T_c) )**(c_r-1d0)
     !	Para4 = h*(1d0-D_cu/2d0-Df/2d0-Dc/2d0)**(-c_p)*dTrdTdet*T_detinv
     !	Para5 =-h*c_p/2d0*(1d0-D_cu/2d0-Df/2d0-Dc/2d0)**(-c_p-1d0)*T_r
      !endif
      !**********************************************************
      
      
      !Para4=h*(1d0-D_cu/2d0-Df/2d0-Dc/2d0)**(-c_p)*dTrdTdet*T_detinv
      Para4=h*(1d0-D_cu/2d0-Dc/2d0)**(-c_p)*dTrdTdet*T_detinv
      
      !Para5=-h*c_p/2d0*(1d0-D_cu/2d0-Df/2d0-Dc/2d0)**(-c_p-1d0)*T_r
      Para5=-h*c_p/2d0*(1d0-D_cu/2d0-Dc/2d0)**(-c_p-1d0)*T_r
                                            
     
!===========dDfdU(3) dDcdU(3)==================================================
      
      if (delu(3) .ge. 0d0) then        
 
      	Rhs(1)=Para1*(abs(delu(3))+delu(3))/2d0+Para2*(Un_inr+abs(Un_inr))/2d0/delc*U_detinv
 
      	Rhs(2)=Para4*(Tn_cu+T_v(3))*sk_n*(1d0-D_cu)/delc
       
      	Smatrix(2,1)=Para5+Para4*sk_n/delc*(delu(3)*(Tn_cu+T_v(3))+delu(1)*(Tt_cu+T_v(1))+delu(2)*(Ts_cu+T_v(2)))
      
      	Smatrix(2,2)=1+Para5+Para4*sk_n/delc*(delu(3)*(Tn_cu+T_v(3))+delu(1)*(Tt_cu+T_v(1))+delu(2)*(Ts_cu+T_v(2)))
 
      else
 
      	Rhs(1)=Para1*(abs(delu(3))+delu(3))/2d0+Para2*(Un_inr+abs(Un_inr))/2d0/delc*U_detinv
      	Rhs(2)=Para4*(Tn_cu+T_v(3))*sk_c/delc
      
        Smatrix(2,1)=Para5+Para4*sk_n/delc*(delu(1)*(Tt_cu+T_v(1))+delu(2)*(Ts_cu+T_v(2)))
      	Smatrix(2,2)=1+Para5+Para4*sk_n/delc*(delu(1)*(Tt_cu+T_v(1))+delu(2)*(Ts_cu+T_v(2)))
      
      endif          
           
      Smatrix(1,1)=1.0+Para3
      Smatrix(1,2)=Para3
     
            
      call DGESV(2,1,Smatrix,2,pivot,Rhs,2,INFO)
        
      dDfdU(3)=Rhs(1)
      dDcdU(3)=Rhs(2)
      
!      write(*,*) 'dDfdU(3),dDcdU(3)',dDfdU(3),dDcdU(3)

      
!===========dDfdU(1) dDcdU(1)==================================================      
      
            
      Rhs(1)=Para1*alpha*delu(1)+Para2*alpha*Ut_inr/delc*U_detinv
      
      Rhs(2)=Para4*(Tt_cu+T_v(1))*sk_n*(1d0-D_cu)/delc
      
      call DGESV(2,1,Smatrix,2,pivot,Rhs,2,INFO)
        
      dDfdU(1)=Rhs(1)
      dDcdU(1)=Rhs(2)
      
!      write(*,*) 'dDfdU(1),dDcdU(1)',dDfdU(1),dDcdU(1)
      
!===========dDfdU(2) dDcdU(2)==================================================      
      
      
            
      Rhs(1)=Para1*alpha*delu(2)+Para2*alpha*Us_inr/delc*U_detinv
      
      Rhs(2)=Para4*(Ts_cu+T_v(2))*sk_n*(1d0-D_cu)/delc
  
      call DGESV(2,1,Smatrix,2,pivot,Rhs,2,INFO)
        
      dDfdU(2)=Rhs(1)
      dDcdU(2)=Rhs(2)   
               
! numerically calculating dDcdU      
      if (Un_inr .eq. 0d0) then
      	dDcdU(3)=0
      else      
      	!dDcdU(3)=(Dc_cu-Dc)/(delu(3)-U_v(3))
        !dDcdU(3) = -h*(1d0-Df_cu/2d0-Dc_cu/2d0-Df/2d0-Dc/2d0)**(-c_p)*dTrdTdet &
     !&                   *sk_n*(Tn_cu+T_v(3))*(1d0-Df_cu-Dc_cu)*T_detinv/delc
     	dDcdU(3) = -h*(1d0-Dc_cu/2d0-Dc/2d0)**(-c_p)*dTrdTdet &
     &                   *sk_n*(Tn_cu+T_v(3))*(1d0-Dc_cu)*T_detinv/delc
      endif
      
      if (Ut_inr .eq. 0d0) then
      	dDcdU(1)=0
      else
      	!dDcdU(1)=(Dc_cu-Dc)/(delu(1)-U_v(1))
        !dDcdU(1) = -h*(1d0-Df_cu/2d0-Dc_cu/2d0-Df/2d0-Dc/2d0)**(-c_p)*dTrdTdet &
     !&                   *sk_n*(Tt_cu+T_v(1))*(1d0-Df_cu-Dc_cu)*T_detinv/delc
     	dDcdU(1) = -h*(1d0-Dc_cu/2d0-Dc/2d0)**(-c_p)*dTrdTdet &
     &                   *sk_n*(Tt_cu+T_v(1))*(1d0-Dc_cu)*T_detinv/delc	
      endif
      
      if (Us_inr .eq. 0d0) then
      	dDcdU(2)=0
      else
      	!dDcdU(2)=(Dc_cu-Dc)/(delu(2)-U_v(2))
        !dDcdU(2) = -h*(1d0-Df_cu/2d0-Dc_cu/2d0-Df/2d0-Dc/2d0)**(-c_p)*dTrdTdet &
     !&                   *sk_n*(Ts_cu+T_v(2))*(1d0-Df_cu-Dc_cu)*T_detinv/delc
        dDcdU(2) = -h*(1d0-Dc_cu/2d0-Dc/2d0)**(-c_p)*dTrdTdet &
     &                   *sk_n*(Ts_cu+T_v(2))*(1d0-Dc_cu)*T_detinv/delc

      endif
	  
	  !**********************************************************
      !	Tung's modification
      !**********************************************************      
      !if ( (T_det .eq. 0d0) .or. (T_det - 2d0*T_c) .le. 0d0 ) then
!	  	dDcdU(3) = 0d0
!		dDcdU(1) = 0d0
!		dDcdU(2) = 0d0
 !     else
  !    	T_detinv = 1d0/T_det
   !   	dTrdTdet = c_r/((2D0*c_C)**c_r)/2d0
    !  	dTrdTdet = dTrdTdet*( (T_det - 2d0*T_c)  &
     !&           + abs(T_det - 2d0*T_c) )**(c_r-1d0)     
      	
      	!dDcdU(3) = -h*(1d0-Df_cu/2d0-Dc_cu/2d0-Df/2d0-Dc/2d0)**(-c_p)*dTrdTdet &
     !&			 *sk_n*(Tn_cu+T_v(3))*(1d0-Df_cu-Dc_cu)*T_detinv/delc
       	!dDcdU(1) = -h*(1d0-Df_cu/2d0-Dc_cu/2d0-Df/2d0-Dc/2d0)**(-c_p)*dTrdTdet &
     !&			 *sk_n*(Tt_cu+T_v(1))*(1d0-Df_cu-Dc_cu)*T_detinv/delc
       	!dDcdU(2) = -h*(1d0-Df_cu/2d0-Dc_cu/2d0-Df/2d0-Dc/2d0)**(-c_p)*dTrdTdet &
     !&			 *sk_n*(Ts_cu+T_v(2))*(1d0-Df_cu-Dc_cu)*T_detinv/delc
      !endif
      !**********************************************************    
      
                                                                  
      dDdU(1) = dDcdU(1) !+ dDfdU(1)
      dDdU(2) = dDcdU(2) !+ dDfdU(2)
      dDdU(3) = dDcdU(3) !+ dDfdU(3)
                                                                        
                                                                  
      Em(1,1)=alpha*sk_n*(1.0-D_cu)/delc                                  
      Em(2,2)=alpha*sk_n*(1.0-D_cu)/delc                                  
                                                                        
      if (popn .ge. 0D0) then                                           
      	Em(3,3)=(1d0-D_cu)*sk_n/delc                                        
      else                                                              
      	Em(3,3)=sk_c/delc                                                 
      endif                                                             
                                                                    
      ! ~Dyadic product                                                                  
      call k_vector_multiply(dEdD,dDdU,D_dot,3,3)                       
      
      ! T_d = D_dot + Em
      call k_matrix_add(D_dot,Em,T_d,3,3)                               
      
!      write(*,*) 'jacobin'
!      write(*,*) T_d(1,1),T_d(2,2),T_d(3,3)
                                                                   
      Svarsi(1)=Df_cu                                                   
      Svarsi(2)=Dc_cu                                                   
      Svarsi(3)=Y_cu  
      Svarsi(4)=delu(1) 
      Svarsi(5)=delu(2) 
      Svarsi(6)=delu(3) 
      Svarsi(7)=T(1,1) 
      Svarsi(8)=T(2,1) 
      Svarsi(9)=T(3,1)  
      
!      do i=1,3
!      do j=1,3
      
!      T_d(i,j)=T_d(i,j)*1000000000d0
      
!      enddo
!      enddo   
      
      
!      T(3,1)=T(3,1)*1000000d0                                               
!      T(1,1)=T(1,1)*1000000d0                          
!      T(2,1)=T(2,1)*1000000d0                                              
                                                                        
                                                                        
      return                                                            
      END                                           
!=====================================================================  


!=====================================================================
!	Function: fatigue residual function OK
! 	This function is for the residue calculation of the Df
!=====================================================================
      REAL*8 FUNCTION residue_f(Df_int,Df,Dc_int,Dc,Y_0,delu,U_v,       &
     & AA,f_n,f_m,alpha,delc,sk_n,Nflag)                                      

      INCLUDE 'ABA_PARAM.INC' 
                                                                        
      Dimension delu(3), U_v(3)
      
      real*8 Df_int,Dc_int,Df,Dc,Y_0,AA,delu,U_v,f_n,sk_n,             &
     & f_m,alpha,Y_crt,Y_crt_n,Un_inr,Ut_inr,Us_inr,U_det,delc
                                                              
      Un_inr = ( delu(3) - U_v(3) )/delc ! normal separation increment
      Ut_inr = ( delu(1) - U_v(1) )/delc ! shear separation increment
      Us_inr = ( delu(2) - U_v(2) )/delc ! shear separation increment
                                                                             
      ! || U^{t+dt} - U^{t} || is calculated based on 1 normal and 2 shear
      U_det = sqrt( ((Un_inr+abs(Un_inr))/2d0)**2d0 + alpha*(Ut_inr)**2d0  &
     & 		+ alpha*(Us_inr)**2d0 )
                                 
	  ! Thermodynamic force                                                                       
      Y_cu = sk_n/2d0/delc*( ((delu(3)+abs(delu(3)))/2d0)**2d0    &
     &			+ alpha*delu(1)**2d0 + alpha*delu(2)**2d0) 
      Y_old= sk_n/2d0/delc*( ((U_v(3)+abs(U_v(3)))/2d0)**2d0      &
     &		    + alpha*U_v(1)**2d0 + alpha*U_v(2)**2d0) 
      
      Y_crt = sqrt( Y_cu/2d0 + Y_old/2d0 ) - sqrt(Y_0) 
      Y_crtn = ( Y_crt/2d0 + abs(Y_crt)/2d0 )**f_n
		
      if ( (1d0-Df_int/2d0-Df/2d0-Dc_int/2d0-Dc/2d0) .LE. 0d0 ) then
      	residue_f = Df_int - Df
      else
      	residue_f = Df_int - Df - AA*( (1d0 - Df_int/2d0 - Df/2d0 - Dc_int/2d0   &
     & 					   - Dc/2d0)**f_m )*Y_crtn*U_det
     	!Nflag=1
      	!return
      endif
                                                                        
      RETURN 
      END                                           

!=====================================================================
!	Function: creep residual function - check norm (T_det-2d0*T_c)???
! 	this subroutine is for the residue calculation of the Dc
!=====================================================================
      
      REAL*8 FUNCTION residue_c(Dc_int,Dc,Df_int,Df,T_c,T_v,            &
     & delu,cC,cr,cp,sk_n,delc,t_h,alpha,sk_c,Nflag)                          

      INCLUDE 'ABA_PARAM.INC' 
                                                                        
!      IMPLICIT REAL*8(A-H,O-Z)

      Dimension delu(3),T_v(3) 
      
      real*8 Dc_int,Dc,Df_int,Df,T_c,delu,cr,cp,cC,sk_n,      &
     & delc,alpha,t_h,Tn_cu,Tt_cu,Ts_cu,T_r,T_det,sk_c                  
                                                                        
	  ! compute normal traction
      if (delu(3) .ge. 0d0) then
      	! tension Un > 0
      	!Tn_cu=sk_n*(1d0-Df_int-Dc_int)*delu(3)/delc
      	Tn_cu=sk_n*(1d0-Dc_int)*delu(3)/delc
      else 
      	! compression Un < 0
      	Tn_cu=sk_c*delu(3)/delc 
      endif
      
      !Tt_cu = alpha*sk_n*(1d0-Df_int-Dc_int)*delu(1)/delc
      Tt_cu = alpha*sk_n*(1d0-Dc_int)*delu(1)/delc 

      !Ts_cu = alpha*sk_n*(1d0-Df_int-Dc_int)*delu(2)/delc
      Ts_cu = alpha*sk_n*(1d0-Dc_int)*delu(2)/delc
	  
	  ! || T^{t+dt} + T^{t} ||
      T_det = sqrt( ((Tn_cu+T_v(3))/2d0+abs(Tn_cu+T_v(3))/2d0)**2d0  &
     & + (1d0/alpha)*(Tt_cu+T_v(1))**2d0	&
     & + (1d0/alpha)*(Ts_cu+T_v(2))**2d0 )
      
      T_r = ( ( (T_det/2d0-T_c) + abs(T_det/2d0-T_c) )/2d0/cC )**cr 
      
      !if ((1d0-Dc_int/2d0-Dc/2d0-Df_int/2d0-Df/2d0) .gt. 0d0) then
      	!residue_c = Dc_int - Dc - t_h*((1d0-Dc_int/2d0-Dc/2d0-Df_int/2d0-Df/2d0)**(-cp))*T_r
      if ((1d0-Dc_int/2d0-Dc/2d0) .gt. 0d0) then
      	residue_c = Dc_int - Dc - t_h*((1d0-Dc_int/2d0-Dc/2d0)**(-cp))*T_r
      else
      	Nflag=1
      	return
      endif
      !-------------------------------------------------------------------------------------------                                
	                                                                          
      RETURN 
      END                                           
                                                                        
!=====================================================================
!	Function: df_1/dD_f OK
!=====================================================================
	  
	  REAL*8 FUNCTION jacob_f(Df_int,Df,Dc_int,Dc,Y_0,delu,U_v,AA,      &
     & 					f_n,f_m,alpha,delc,sk_n)                                         
! this subroutine is for the jacob calculation of the Df, need to verify
! stability as there is machale bracket                                 

      INCLUDE 'ABA_PARAM.INC' 
      
      Dimension delu(3),U_v(3) 
      REAL*8 Df_int,Df,Dc_int,Dc,AA,delu,U_v,f_n,f_m,sk_n,            &
     & alpha,Un_inr,Ut_inr,Us_inr,U_det,delc
                                                                        
      Un_inr = (delu(3)-U_v(3))/delc 
      Ut_inr = (delu(1)-U_v(1))/delc 
      Us_inr = (delu(2)-U_v(2))/delc 
      
      ! || U^{t+dt} - U^{t} || is calculated based on 1 normal and 2 shear
      U_det = sqrt(((Un_inr+abs(Un_inr))/2d0)**2d0+alpha*(Ut_inr)**2d0          &
     & 		+ alpha*(Us_inr)**2d0)                                              
      
      ! || Y^{t+dt} ||
      Y_cu  = sk_n/2d0/delc*(((delu(3)+abs(delu(3)))/2d0)**2d0+alpha*delu(1)**2d0+alpha*delu(2)**2d0)
      
      ! || Y^{t} ||
      Y_old = sk_n/2d0/delc*(((U_v(3)+abs(U_v(3)))/2d0)**2d0+alpha*U_v(1)**2d0+alpha*U_v(2)**2d0)
      
      Y_crt  = sqrt(Y_cu/2d0+Y_old/2d0)-sqrt(Y_0)  
      
      Y_crtn = (Y_crt/2d0+abs(Y_crt)/2d0)**f_n
      
  	  !write(*,*) 'vars in jacob_f:', 1d0-Df_int/2d0-Df/2d0-Dc_int/2d0-Dc/2d0
  	  
  	  if ((1d0-Df_int/2d0-Df/2d0-Dc_int/2d0-Dc/2d0) .LE. 0d0) then
      	  jacob_f = 1d0
      else  	  
  	   	  jacob_f = 1d0 + AA*f_m/2d0*((1d0-Df_int/2d0-Df/2d0-Dc_int/2d0-Dc/2d0)**(f_m-1d0))*Y_crtn*U_det
      !else
      	  !Nflag=1
      	  !return
  	  endif
                                                                        
      RETURN 
      END                                           

!=====================================================================
!	Function: df_2/dD_f and df2_/dD_c need to be checked again ???
!=====================================================================
                                                                        
      REAL*8 FUNCTION jacob_c(Dc_int,Dc,Df_int,Df,T_c,T_v,delu,         &
     & cr,cp,cC,sk_n,delc,t_h,alpha,sk_c)                               
! this subroutine is for the residue calculation of the Dc
             
      INCLUDE 'ABA_PARAM.INC' 
      
      Dimension delu(3),T_v(3) 
      
      REAL*8 Dc_int,Dc,Df_int,Df,T_c,T_v,delu,cr,  &
     & cp,cC,sk_n,delc,alpha,Un,Ut,Us,U_det,h,Tt_cu,        &
     & Ts_cu, T_det,T_r,dTrdTdet,dTdetdDc,t_h,sk_c                      

	  if (delu(3) .ge. 0d0 ) then 
	  	! If Tensile Un >= 0
      	!Tn_cu = sk_n*(1d0 - Df_int - Dc_int)*delu(3)/delc
      	Tn_cu = sk_n*(1d0 - Dc_int)*delu(3)/delc
      else
      	! If compression Un < 0
      	Tn_cu = sk_c*delu(3)/delc
      endif                                                                         
                                                                        
      !Tt_cu = alpha*sk_n*(1d0 - Df_int - Dc_int)*delu(1)/delc
      Tt_cu = alpha*sk_n*(1d0 - Dc_int)*delu(1)/delc 
      !Ts_cu = alpha*sk_n*(1d0 - Df_int - Dc_int)*delu(2)/delc
      Ts_cu = alpha*sk_n*(1d0 - Dc_int)*delu(2)/delc
      
      ! || T^{t+dt} + T^{t} ||
      T_det = sqrt( ((Tn_cu + T_v(3))/2d0+abs(Tn_cu+T_v(3))/2d0)**2d0 &
     & 			   + (Tt_cu + T_v(1))**2d0/alpha                      &
     &             + (Ts_cu + T_v(2))**2d0/alpha )
      
      T_r = ( ( (T_det/2d0-T_c) + abs(T_det/2d0-T_c) )/2.0D0/cC )**cr 

      !--------------------------------------	                                                                        
      ! Yumeng's derivative
      !--------------------------------------      
      !dTrdTdet = cr/2.0D0/(cC**cr)*(((T_det/2d0-T_c)+abs(T_det/2d0-T_c))/2.0D0)**(cr-1d0)
      
      ! Tung's 
      dTrdTdet = cr/((2D0*cC)**cr)/2d0
      dTrdTdet = dTrdTdet*( (T_det - 2d0*T_c)+abs(T_det - 2d0*T_c) )**(cr-1d0)
      
      if (T_det .ne. 0d0) then 
      	if (delu(3) .ge. 0D0) then 
      		dTdetdDc = ( (Tn_cu+T_v(3))*(-sk_n*delu(3)/delc)  &
     & 	           + (Tt_cu+T_v(1))*(-sk_n*delu(1)/delc)  &
     &                + (Ts_cu+T_v(2))*(-sk_n*delu(2)/delc) )/T_det
      	else 
      		dTdetdDc = ((Tt_cu+T_v(1))*(-sk_n*delu(1)/delc)  &
     & 			  + (Ts_cu+T_v(2))*(-sk_n*delu(2)/delc))/T_det 
      	endif       
      else                                                                         
      	dTdetdDc=0d0                                                                         
      endif
      
      !if ((1d0-Df_int/2d0-Df/2d0-Dc_int/2d0-Dc/2d0) .LE. 0d0) then
      if ((1d0-Dc_int/2d0-Dc/2d0) .LE. 0d0) then
        jacob_c = 1d0
      else
        !jacob_c = 1d0 - t_h*cp/2d0*( (1d0-Dc_int/2d0-Dc/2d0-Df_int/2d0-Df/2d0)**(-cp-1d0) )*T_r &
     !&		        - t_h*((1d0-Dc_int/2d0-Dc/2d0-Df_int/2d0-Df/2d0)**(-cp))*dTrdTdet*dTdetdDc
     	jacob_c = 1d0 - t_h*cp/2d0*( (1d0-Dc_int/2d0-Dc/2d0)**(-cp-1d0) )*T_r &
     &		        - t_h*((1d0-Dc_int/2d0-Dc/2d0)**(-cp))*dTrdTdet*dTdetdDc
      endif
 
      RETURN 
      END                                           
                   

!=====================================================================
!	subroutine: K in local coordinates
!=====================================================================
SUBROUTINE k_local_coordinates(co_de,r,coord_l,transformation_m,  &
    transformation_m_t,a_jacob,ajacob_m,coords,u,ndofel,nnode, mcrd)
 INCLUDE 'ABA_PARAM.INC'


    dimension R(mcrd,mcrd),coord_l(mcrd,nnode),aJacob_M(2,3),&
     & Transformation_M(ndofel,ndofel),coords(mcrd,nnode),&
     & Transformation_M_T(ndofel,ndofel),u(ndofel),&
     & co_de(mcrd,nnode), co_de_m(3,3),SFD(2,4)

CALL k_matrix_zero(co_de_m,3,4)

DO i = 1, 3
  co_de_m(i,1)=(co_de(i,1)+co_de(i,4))*0.5
  co_de_m(i,2)=(co_de(i,2)+co_de(i,5))*0.5
  co_de_m(i,3)=(co_de(i,3)+co_de(i,6))*0.5
END DO

sfd(1,1) =-1
sfd(1,2) = 1
sfd(1,3) = 0
sfd(2,1) =-1
sfd(2,2) = 0
sfd(2,3) = 1

DO i = 1,2
  DO j = 1,3
    DO k =1, 3
      ajacob_m(i,j) = ajacob_m(i,j) + sfd(i,k)*co_de_m(j,k)
    END DO
  END DO
END DO

dum1 = ajacob_m(1,2)*ajacob_m(2,3) - ajacob_m(1,3)*ajacob_m(2,2)
dum2 = ajacob_m(1,3)*ajacob_m(2,1) - ajacob_m(1,1)*ajacob_m(2,3)
dum3 = ajacob_m(1,1)*ajacob_m(2,2) - ajacob_m(1,2)*ajacob_m(2,1)

a_jacob = SQRT(dum1**2 + dum2**2 + dum3**2)/2.0D0
rn1 = SQRT(dum1**2 + dum2**2 + dum3**2)

r(3,1) = dum1/rn1
r(3,2) = dum2/rn1
r(3,3) = dum3/rn1

alen=SQRT(ajacob_m(1,1)**2.0 + ajacob_m(1,2)**2.0 + ajacob_m(1,3)**2.0)
r(1,1)=ajacob_m(1,1)/alen
r(1,2)=ajacob_m(1,2)/alen
r(1,3)=ajacob_m(1,3)/alen

r(2,1)=r(3,2)*r(1,3)-r(3,3)*r(1,2)
r(2,2)=r(3,3)*r(1,1)-r(3,1)*r(1,3)
r(2,3)=r(3,1)*r(1,2)-r(3,2)*r(1,1)

!=====================================================================
num=nnode

DO i = 1, num
  dum=3.0*(i-1.0)
  transformation_m(dum+1,dum+1)=r(1,1)
  transformation_m(dum+1,dum+2)=r(1,2)
  transformation_m(dum+1,dum+3)=r(1,3)
  transformation_m(dum+2,dum+1)=r(2,1)
  transformation_m(dum+2,dum+2)=r(2,2)
  transformation_m(dum+2,dum+3)=r(2,3)
  transformation_m(dum+3,dum+1)=r(3,1)
  transformation_m(dum+3,dum+2)=r(3,2)
  transformation_m(dum+3,dum+3)=r(3,3)
END DO

CALL k_matrix_transpose(transformation_m,transformation_m_t, ndofel,ndofel)

DO i = 1, nnode
  coord_l(1,i)=(r(1,1)*co_de(1,i)+r(1,2)*co_de(2,i) +r(1,3)*co_de(3,i))
  coord_l(2,i)=(r(2,1)*co_de(1,i)+r(2,2)*co_de(2,i) +r(2,3)*co_de(3,i))
  coord_l(3,i)=(r(3,1)*co_de(1,i)+r(3,2)*co_de(2,i) +r(3,3)*co_de(3,i))
END DO

RETURN
END SUBROUTINE k_local_coordinates

!=====================================================================
!	Subroutine: shape function
!=====================================================================
	SUBROUTINE k_shape_fun(i,sf)
	INCLUDE 'ABA_PARAM.INC'
	
	DIMENSION  sf(3),gp_coord(2)

IF (i == 1) THEN
  gp_coord(1)= 1.0D0/6.0D0
  gp_coord(2)= 2.0D0/3.0D0
ELSE IF (i == 2) THEN
  gp_coord(1)= 2.0D0/3.0D0
  gp_coord(2)= 1.0D0/6.0D0
ELSE IF (i == 3) THEN
  gp_coord(1)= 1.0D0/6.0D0
  gp_coord(2)= 1.0D0/6.0D0
END IF

sf(1)= 1.0 - gp_coord(1) - gp_coord(2)
sf(2)= gp_coord(1)
sf(3)= gp_coord(2)


	RETURN
	END SUBROUTINE k_shape_fun

!=====================================================================
!	Function:
!=====================================================================
!     A and B are vector, C=A*transpose(B)                              
	subroutine k_vector_multiply(A,B,C,l,m) 
     INCLUDE 'ABA_PARAM.INC' 
                                                                      
!      IMPLICIT REAL*8(A-H,O-Z)                                         
      dimension A(l),B(m),C(l,m) 
!                                                                       
      call k_matrix_zero(C,l,m) 
!                                                                       
      do i = 1, l 
         do j = 1, m 
               C(i,j)=A(i)*B(j) 
         end do 
      end do 
!                                                                       
      return 
      END   

!=====================================================================
!	Function:
!=====================================================================
!     A and B are matrix, C=A+b                                         
      subroutine k_matrix_add(A,B,C,l,m) 
      INCLUDE 'ABA_PARAM.INC' 
                                                                        
                                         
      dimension A(l,m),B(l,m),C(l,m) 
                                                                        
      call k_matrix_zero(C,l,m) 
                                                                        
      do i = 1, l 
         do j = 1, m 
               C(i,j)=A(i,j)+B(i,j) 
         end do 
      end do 
!                                                                       
      return 
      END
       
!=====================================================================
!	subroutine
!=====================================================================
SUBROUTINE k_matrix_multiply(a,b,c,l,n,m)
INCLUDE 'ABA_PARAM.INC'


dimension A(l,n),B(n,m),C(l,m)

CALL k_matrix_zero(c,l,m)

DO i = 1, l
  DO j = 1, m
    DO k = 1, n
      c(i,j)=c(i,j)+a(i,k)*b(k,j)
    END DO
  END DO
END DO

RETURN
END SUBROUTINE k_matrix_multiply

!=====================================================================
!	subroutine
!=====================================================================
SUBROUTINE k_matrix_plus_scalar(a,b,c,n,m)
INCLUDE 'ABA_PARAM.INC'


 dimension A(n,m),B(n,m) 


DO i = 1, n
  DO j = 1, m
    a(i,j)=a(i,j)+c*b(i,j)
  END DO
END DO

RETURN
END SUBROUTINE k_matrix_plus_scalar

!=====================================================================
!	subroutine
!=====================================================================
SUBROUTINE k_matrix_transpose(a,b,n,m)
INCLUDE 'ABA_PARAM.INC'


dimension A(n,m),B(m,n)


DO i = 1, n
  DO j = 1, m
    b(j,i)=a(i,j)
  END DO
END DO

RETURN
END SUBROUTINE k_matrix_transpose

!=====================================================================
!	subroutine
!=====================================================================
SUBROUTINE k_matrix_zero(a,n,m)
INCLUDE 'ABA_PARAM.INC'


 dimension A(n,m) 


DO i = 1, n
  DO j = 1, m
    a(i,j)=0.d0
  END DO
END DO

RETURN
END SUBROUTINE k_matrix_zero

!=====================================================================
!	subroutine
!=====================================================================
SUBROUTINE k_vector_zero(a,n)
INCLUDE 'ABA_PARAM.INC'

  dimension A(n) 


DO i = 1, n
  a(i)=0.d0
END DO

RETURN
END SUBROUTINE k_vector_zero


!=====================================================================
!	subroutine
!=====================================================================
SUBROUTINE k_mac(pm,a,b)
INCLUDE 'ABA_PARAM.INC'

IF ((a-b) >= 0.0) THEN
  pm=a-b
ELSE IF ((a-b) < 0.0) THEN
  pm=0.d0
END IF

RETURN
END SUBROUTINE k_mac
!=====================================================================
!=================================END=================================
!=====================================================================
