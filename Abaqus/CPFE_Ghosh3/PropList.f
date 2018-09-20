c -----------------
c     PROP(1-9): ROTATION MATRIX
c     PROP(10-12): C11, C12, C44

c ---- GetRhoPFMGND *** NOTE DIFFERENT***-----
c     PROPS(13) :: (1) :: c_10 * kB * Theta / (G * b^3)

c ---- GetTauSlips -----
c     PROPS(14) :: (1) :: c_3 * G * b
c     PROPS(15) :: (2) :: c_4 * kB * Theta/ (b^2)

c ---- GetCSDHTauC -----
c     PROPS(16) :: (1) :: b / Gamma_111
c     PROPS(17) :: (2) :: G * b^3 / (4 * pi)
c     PROPS(18) :: (3) :: G * b^2 / (2 * pi * Gamma_111)
c     ** PROPS(19) :: (4) ::  xi * G * b = xi_0 * exp(A/(Theta-Theta_c)) * G * b
c     ** PROPS(20) :: (5) :: tau_cc
c     PROPS(21) :: (6) :: C_H
c     PROPS(22) :: (7) ::  h
c     PROPS(23) :: (8) :: k_1
c     PROPS(24) :: (9) :: k_2
c     PROPS(25) :: (10) :: (1/sqrt(3)) - Gamma_010 / Gamma_111
c     ** PROPS(26) :: (11) :: b / B
c     ** PROPS(27) :: (12) :: rho_0
c     PROPS(28) :: (13) :: kB * Theta 

c ---- Calculate SLup Rate **** NEW **** -----
c     PROPS(29) :: (1) :: [c_10 * kB * Theta / (G * b^2) ] * [c_1*(THETA^c_2)*exp(-Q/kB*THETA)]
c     PROPS(30) :: (2) :: p