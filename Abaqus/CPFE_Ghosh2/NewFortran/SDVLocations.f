C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c PROPS  Definitions
c
c PROPS(1->9) :: Orientation Matrix
c PROPS(1->27) :: Rho_SSD Initial
c PROPS(28->48) :: Compliance Matrix in Voigt Notation

c PROPS(49) :: TYPE (NOT USED)

c Note that these properties were precalculated so they don't have to be repeatedly calculated

c GetRhoPFMGND
c PROPS(50) :: (1) :: c_10 * kB * Theta / (G * b^2)

c GetTauSlips
c PROPS(51) :: (1) :: c_3 * G * b
c PROPS(52) :: (2) :: c_4 * kB * Theta/ (b^2)
c PROPS(53) :: (3) :: c_1 * (Theta ^ (c_2))

c GetCSDHTauC
c PROPS(54) :: (1) :: b / Gamma_111
c PROPS(55) :: (2) :: G * b^3 / (4 * pi)
c PROPS(56) :: (3) :: G * b^2 / (2 * pi * Gamma_111)
c ** PROPS(57) :: (4) ::  xi * G * b = xi_0 * exp(A/(Theta-Theta_c)) * G * b
c ** PROPS(58) :: (5) :: tau_cc
c PROPS(59) :: (6) :: C_H
c PROPS(60) :: (7) ::  h
c PROPS(61) :: (8) :: k_1
c PROPS(62) :: (9) :: k_2
c PROPS(63) :: (10) :: (1/sqrt(3)) - Gamma_010 / Gamma_111
c ** PROPS(64) :: (11) :: b / B
c ** PROPS(65) :: (12) :: rho_0
c PROPS(66) :: (13) :: kB * Theta

c GetGammaDot
c PROPS(67) :: (1) :: exp(-Q / (kB * Theta))
c PROPS(68) :: (2) :: p
c PROPS(69) :: (3) :: b

c GetRhoSSDEvolve
c PROPS(70) :: (1) :: c_5 / b
c PROPS(71) :: (2) :: (c_6 / b) * (sqrt(3) * G * b)/ (16 * (1-nu))
c PROPS(72) :: (3) :: c_7
c PROPS(73) :: (4) :: c_8 * ( (D_0 b^3) / (kB * Theta) ) * exp(- Q_Bulk / (kB * Theta)))
c PROPS(74) :: (5) :: c_9
c PROPS(75) :: (6) :: gamma_dot_ref

c No Longer In Use
c PROPS(76) :: (1) :: c_11
c PROPS(77) :: (2) :: c_12
c PROPS(78) :: (3) :: c_44

c No Longer In Use
c PROPS(79) :: (1) :: rho_ssd

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c SDV  Definitions
c
C SDV(1->54) :: Slip direction ((ISLIP-1)*3+COMPONENT)
C SDV(55->108) :: Slip normals ((ISLIP-1)*3+COMPONENT)+OFFSET
C SDV(109->126) :: Rho SSD (ISLIP)
C SDV(127->144) :: Rho CSD (ISLIP)

C SDV(145->162) :: Gamma Slip cumulative
c SDV(163) :: Cumulative gamma slip
c SDV(164->184) :: Compliance Matrix

c SDV(185->220) ::  Slip direction PE ((ISLIP-1)*3+COMPONENT)+OFFSET
c SDV(221->256) ::  Slip normals PE ((ISLIP-1)*3+COMPONENT)+OFFSET
c SDV(257->292) ::  Slip direction SE ((ISLIP-1)*3+COMPONENT)+OFFSET
c SDV(293->328) ::  Slip normals SE ((ISLIP-1)*3+COMPONENT)+OFFSET
c SDV(329->364) ::  Slip direction CB ((ISLIP-1)*3+COMPONENT)+OFFSET
c SDV(365->400) ::  Slip normals CB ((ISLIP-1)*3+COMPONENT)+OFFSET

c SDV(401->409) :: Plastic Deformation Tensor
c SDV(410->427) :: Rho GND S (SLIP)
c SDV(428->429) :: Nothing Made a mistake here previously
c SDV(430->447) :: Rho GND ET (SLIP)
c SDV(448->465) :: Rho GND EN (SLIP)

c *********************************************
