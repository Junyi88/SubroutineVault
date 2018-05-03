      subroutine GetDSTRESS(DStress,GammaDot,DStrain,Stress,dTIME, 
     1 FCC_Mu,FCC_Ohm,Cubic_Mu,Cubic_Ohm,	   
     2 CinS)

C Subroutine to calculate forest parallel and mobile dislocations
      
      implicit none
      
      real*8,intent(in) :: dTIME 
      real*8,intent(in) :: GammaDot(18),DStrain(6),Stress(6)
      real*8,intent(in) :: FCC_Mu(3,3,18),FCC_Ohm(3,3,18)
      real*8,intent(in) :: Cubic_Mu(3,3,18),Cubic_Ohm(3,3,18)
      real*8,intent(out) :: DStress(6)

      real*8,intent(in) :: CinS(3)
      
      real*8:: HYDROSTRAIN, DGA(18), DUM1
      integer ISLIPS

      HYDROSTRAIN=DSTRAIN(1)+DSTRAIN(2)+DSTRAIN(3)
	  
      DO ISLIPS=1,18
	        DGA(ISLIPS)=GammaDot(ISLIPS)*DTIME
      END DO
      DO ISLIPS=1,6
	        DStress(ISLIPS)=0.0
      END DO
	  
      print *, 'CINS=', CinS
C ------------------------------------------------------	
      DSTRESS(1)=CINS(1)*DSTRAIN(1)+CINS(2)*(DSTRAIN(2)+DSTRAIN(3))


C ------------------------------------------------------	
      DSTRESS(2)=CINS(1)*DSTRAIN(2)+CINS(2)*(DSTRAIN(1)+DSTRAIN(3))


C ------------------------------------------------------	
      DSTRESS(3)=CINS(1)*DSTRAIN(3)+CINS(2)*(DSTRAIN(1)+DSTRAIN(2))

	  
C ===============================================================================
C ------------------------------------------------------	
      DSTRESS(4)=CINS(3)*DSTRAIN(4)-HYDROSTRAIN*STRESS(4)
	  
 
C ------------------------------------------------------	
      DSTRESS(5)=CINS(3)*DSTRAIN(5)-HYDROSTRAIN*STRESS(5)
	  
   
C ------------------------------------------------------	
      DSTRESS(6)=CINS(3)*DSTRAIN(6)-HYDROSTRAIN*STRESS(6)
	  
      return
      end subroutine GetDSTRESS