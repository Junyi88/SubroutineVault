c -----------------
c     SDV (1->9): ROTATION MATRIX
c     SDV (10->18): FP [10,11,12 | 13,14,15 | 16,17,18]
c     SDV (19->27): curlFP
c     SDV(28->30) :: COORDS

c     SDV(31->42) :: G
c     SDV(43->54) :: gamma

c          STATEV(55) = FaiValue
c        STATEV(56) = ITRNUM   

C XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c     SDV(28->45) :: Rho SSD (ISLIP) 
c     SDV(46->63) :: Rho CSD (ISLIP)

c     SDV(64->81) :: Rho GND S (SLIP)
c     SDV(82->99) :: Rho GND ET (SLIP)
c     SDV(100->117) :: Rho GND EN (SLIP)

c     SDV(118->135) :: Gamma Slip cumulative
c     SDV(136->138) :: COORDS

c     SDV(139) :: GAMMA SLIP ALL
c     SDV(140) :: RHO SSD ALL
c     SDV(141) :: RHO CSD ALL
c     SDV(142) :: RHO GND ALL

c ---------------------------------
c SDV(163) :: Cumulative gamma slip