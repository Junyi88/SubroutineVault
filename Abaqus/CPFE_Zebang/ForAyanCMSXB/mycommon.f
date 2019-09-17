!Use this to avoid editing multiple common blocks independently
		integer, parameter :: nElements = 203723
		integer, parameter :: nintpts = 8
		integer, parameter :: nsdv  = 102
          real*8 :: kFp, kgausscoords, kcurlFp
          
      COMMON/UMPS/kFp(nElements, nintpts, 9),          
     +    kgausscoords(nElements,nintpts,3),
     +    kcurlFp(nElements, nintpts, 9)


