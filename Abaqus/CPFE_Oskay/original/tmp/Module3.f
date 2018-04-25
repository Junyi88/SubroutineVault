module WorkDir
!
!----Define the file path and root name

	character*80 filepath
	parameter (filePath='../')

	character*80 fileroot
	parameter (fileroot='test')
end
!----------------------------------------------------------------------
!----------------------------------------------------------------------
module DummyVar
       integer numel, numqpt, iprint      
       parameter (numel=1, numqpt=1, iprint=1 )

end
!----------------------------------------------------------------------
module SlipData
	integer numslip, numgrn, numvtx
        real*8  sigfs(5, 241)
 
end
!----------------------------------------------------------------------
module OriData
      integer kODF
      real*8  angles(3, 250)
!      real*8  euler(3, 2000, 1, 1)
      integer numor, seed
end 
!----------------------------------------------------------------------
module XtalPar
      real*8  matProp(20, 48)
      real*8  hardmtx(48, 48)
 
end
!----------------------------------------------------------------------
module SlipOverStress
        real*8  overstress(48)
end
!----------------------------------------------------------------------
module  IterData
        integer maxIterstate, MaxIterNewt
        real*8  tolerState, tolerNewt 
end
module WriteControl
	integer kODFout
end
!----------------------------------------------------------------------
module InitialState
       integer NKAPP1
       parameter (NKAPP1=12)
	   real*8 slipr0(NKAPP1)
       real*8 backst0(NKAPP1)
       real*8 backst0c(NKAPP1)
end
!----------------------------------------------------------------------
!module InitialLocalRotMatrix
	!real*8  globalcrot0   (3, 3, 2000, 1, 1) 
!end
!----------------------------------------------------------------------
module InitialSlipSys
	real*8  zBar0(3, 3, 12), KBar0(3, 3, 12)
	real*8  pBar0(3, 3, 12), qBar0(3, 3, 12)
	real*8  pBar0Vec(5, 12), qBar0Vec(3, 12)
	real*8  ppTBar0(5, 5, 12)
	real*8  KpBar0(3, 3, 12), KqBar0(3, 3, 12)
	real*8  KpBar0Vec(5, 12), KqBar0Vec(3, 12)
	real*8  KpBar0devVec(5, 12)
	real*8  KpKpTBar0(5, 5, 12)
	real*8  tauSlip(12)
end
!----------------------------------------------------------------------
module  EProp
	real*8  fCeDev(5,5)
	real*8  fCeiDev(5,5)
	real*8  fCeDevVol(5)
	real*8  fCeVol
end
!----------------------------------------------------------------------
module TransData
	real*8  fDevMat5x6(5, 6), fDevMat6x5(6, 5)
	real*8  fMatTId5x6(5, 6)
end
!----------------------------------------------------------------------
module DataType
	type xtalVars
	      real*8  gstress    (5, 1, 1, 1)
	      real*8  gestran    (5, 1, 1, 1)
	      real*8  gslipr     (12, 1, 1, 1)
	      real*8  gbackst    (12, 1, 1, 1)
	      real*8  gbackstc   (12, 1, 1, 1)
	      real*8  gstatev    (5, 1, 1, 1)
	      real*8  geqvalues  (8, 1, 1, 1)
	      real*8  ggamdot    (12, 1, 1, 1)
	      real*8  ggamdot2    (12, 1, 1, 1)
	      real*8  gcrot      (3, 3, 1, 1, 1)
	      real*8  grrot      (3, 3, 1, 1, 1)        
        end type xtalVars
!----------------------------------------------------------------------
	type xtalVars_n
	      real*8  gstress_n    (5, 1, 1, 1)
	      real*8  gestran_n    (5, 1, 1, 1)
	      real*8  gslipr_n     (12, 1, 1, 1)
	      real*8  gbackst_n    (12, 1, 1, 1)
	      real*8  gbackstc_n    (12, 1, 1, 1)
	      real*8  gstatev_n    (5, 1, 1, 1)
	      real*8  gcrot_n      (3, 3, 1, 1, 1)
	      real*8  grrot_n      (3, 3, 1, 1, 1)        
        end  type xtalVars_n
end    module DataType
