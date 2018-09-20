      subroutine PrepNewtonIter(
     1 DFGRD1,DFGRD0, DTIME,
     1 StiffR, STRESS,
     1 Ltot, DTotStran,
     1 StressV, StressVMat,
     1 StressTrial, StressTrialMat
     1 )
     
      IMPLICIT NONE
      INTEGER :: I,J,K, NN

      REAL*8, INTENT(IN) :: DFGRD1(3,3),DFGRD0(3,3), DTIME 
      REAL*8, INTENT(IN) :: StiffR(6,6), STRESS(6)

      REAL*8, INTENT(OUT) :: StressV(6), StressVMat(3,3)
      REAL*8, INTENT(OUT) :: StressTrial(6), StressTrialMat(3,3)
      REAL*8, INTENT(OUT) :: Ltot(3,3), DTotStran(6)
      
      INTEGER :: info
      REAL*8 :: F(3,3), Fdot(3,3), invF(3,3) 
      REAL*8 :: TempStrain(3,3), Spin(3,3)     
        
C ==== Calculate Strains =================     
      F = DFGRD0   
      Fdot = (DFGRD1-DFGRD0)/DTIME 
      CALL lapinverse(F,3,info,invF) 
      Ltot = matmul(Fdot,invF)  
       
      TempStrain=(Ltot+transpose(Ltot))*0.5*DTIME        
      Spin=(Ltot-transpose(Ltot))*0.5 
      CALL kmatvec6(TempStrain,DTotStran)
      DTotStran(4:6) = 2.0*DTotStran(4:6)     
     
C ==== Calculate Trial Stress =================  
      StressV = STRESS
      StressTrial = StressV + matmul(StiffR,DTotStran)
      CALL kvecmat6(StressTrial,StressTrialMat)        
      CALL kvecmat6(StressV,StressVMat)
      
      StressTrialMat = StressTrialMat + 
     2 (matmul(Spin,StressVMat)- matmul(StressVMat,Spin))*DTIME        
         
      return
      end 