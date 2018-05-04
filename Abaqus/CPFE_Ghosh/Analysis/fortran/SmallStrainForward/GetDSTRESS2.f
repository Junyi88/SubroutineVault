      subroutine GetDSTRESS(DStress,GammaDot,DStrain,Stress,dTIME, 
     1 FCC_Mu,FCC_Ohm,Cubic_Mu,Cubic_Ohm,	   
     2 CinS)

C Subroutine to calculate forest parallel and mobile dislocations
      
      implicit none
      
      real*8,intent(in) :: dTIME 
      real*8,intent(in) :: GammaDot(18),DStrain(6),Stress(6)
      real*8,intent(in) :: FCC_Mu(3,3,12),FCC_Ohm(3,3,12)
      real*8,intent(in) :: Cubic_Mu(3,3,6),Cubic_Ohm(3,3,6)
      real*8,intent(out) :: DStress(6)

      real*8,intent(in) :: CinS(3)
      
      real*8:: HYDROSTRAIN, DGA(18), DUM1
      integer ISLIPS

      HYDROSTRAIN=DSTRAIN(1)+DSTRAIN(2)+DSTRAIN(3)
	  
	  
      DO ISLIPS=1,18
	        DGA(ISLIPS)=GammaDot(ISLIPS)*DTIME
      END DO
C XXX	  

      DO ISLIPS=1,6
	        DStress(ISLIPS)=0.0
      END DO
C ------------------------------------------------------	
      DSTRESS(1)=CINS(1)*DSTRAIN(1)+CINS(2)*(DSTRAIN(2)+DSTRAIN(3))
	  
      DO ISLIPS=1,12
       DUM1=CINS(1)*FCC_Mu(1,1,ISLIPS)+
     1 CINS(2)*(FCC_Mu(2,2,ISLIPS)+FCC_Mu(3,3,ISLIPS))
	 
       DSTRESS(1)=DSTRESS(1)-DUM1*DGA(ISLIPS)
	 
      END DO
      DO ISLIPS=1,6
       DUM1=CINS(1)*CUBIC_Mu(1,1,ISLIPS)+
     1 CINS(2)*(CUBIC_Mu(2,2,ISLIPS)+CUBIC_Mu(3,3,ISLIPS))
	 
       DSTRESS(1)=DSTRESS(1)-DUM1*DGA(ISLIPS+12)
   
      END DO

C ------------------------------------------------------	
      DSTRESS(2)=CINS(1)*DSTRAIN(2)+CINS(2)*(DSTRAIN(1)+DSTRAIN(3))
      DSTRESS(2)=DSTRESS(2)-HYDROSTRAIN*STRESS(2)
	  
      DO ISLIPS=1,12
       DUM1=CINS(1)*FCC_Mu(2,2,ISLIPS)+
     1 CINS(2)*(FCC_Mu(1,1,ISLIPS)+FCC_Mu(3,3,ISLIPS))
	 
       DSTRESS(2)=DSTRESS(2)-DUM1*DGA(ISLIPS)
 
      END DO
      DO ISLIPS=1,6
       DUM1=CINS(1)*CUBIC_Mu(2,2,ISLIPS)+
     1 CINS(2)*(CUBIC_Mu(1,1,ISLIPS)+CUBIC_Mu(3,3,ISLIPS))
	 
       DSTRESS(2)=DSTRESS(2)-DUM1*DGA(ISLIPS+12)
  
      END DO	  

C ------------------------------------------------------	
      DSTRESS(3)=CINS(1)*DSTRAIN(3)+CINS(2)*(DSTRAIN(1)+DSTRAIN(2))
      DSTRESS(3)=DSTRESS(3)-HYDROSTRAIN*STRESS(3)
	  
      DO ISLIPS=1,12
       DUM1=CINS(1)*FCC_Mu(3,3,ISLIPS)+
     1 CINS(2)*(FCC_Mu(1,1,ISLIPS)+FCC_Mu(2,2,ISLIPS))
	 
       DSTRESS(3)=DSTRESS(3)-DUM1*DGA(ISLIPS)
  
      END DO
      DO ISLIPS=1,6
       DUM1=CINS(1)*CUBIC_Mu(3,3,ISLIPS)+
     1 CINS(2)*(CUBIC_Mu(1,1,ISLIPS)+CUBIC_Mu(2,2,ISLIPS))
	 
       DSTRESS(3)=DSTRESS(3)-DUM1*DGA(ISLIPS+12)
   
      END DO		  
	  
	  
C ===============================================================================
C ------------------------------------------------------	
      DSTRESS(4)=CINS(3)*DSTRAIN(4)-HYDROSTRAIN*STRESS(4)
	  
      DO ISLIPS=1,12
       DUM1=CINS(3)*(FCC_Mu(1,2,ISLIPS)+FCC_Mu(2,1,ISLIPS))
       DSTRESS(4)=DSTRESS(4)-DUM1*DGA(ISLIPS)
   
      END DO
      DO ISLIPS=1,6
       DUM1=CINS(3)*(CUBIC_Mu(1,2,ISLIPS)+CUBIC_Mu(2,1,ISLIPS))
       DSTRESS(4)=DSTRESS(4)-DUM1*DGA(ISLIPS+12)
		  
      END DO	

C ------------------------------------------------------	
      DSTRESS(5)=CINS(3)*DSTRAIN(5)-HYDROSTRAIN*STRESS(5)
	  
      DO ISLIPS=1,12
       DUM1=CINS(3)*(FCC_Mu(1,3,ISLIPS)+FCC_Mu(3,1,ISLIPS))
       DSTRESS(5)=DSTRESS(5)-DUM1*DGA(ISLIPS)
   
      END DO
      DO ISLIPS=1,6
       DUM1=CINS(3)*(CUBIC_Mu(1,3,ISLIPS)+CUBIC_Mu(3,1,ISLIPS))
       DSTRESS(4)=DSTRESS(4)-DUM1*DGA(ISLIPS+12)
	  
      END DO	

C ------------------------------------------------------	
      DSTRESS(6)=CINS(3)*DSTRAIN(6)-HYDROSTRAIN*STRESS(6)
	  
      DO ISLIPS=1,12
       DUM1=CINS(3)*(FCC_Mu(2,3,ISLIPS)+FCC_Mu(3,2,ISLIPS))
       DSTRESS(6)=DSTRESS(6)-DUM1*DGA(ISLIPS)
 
      END DO
      DO ISLIPS=1,6
       DUM1=CINS(3)*(CUBIC_Mu(2,3,ISLIPS)+CUBIC_Mu(3,2,ISLIPS))
       DSTRESS(6)=DSTRESS(6)-DUM1*DGA(ISLIPS+12)
	  
      END DO		  
	  
      return
      end subroutine GetDSTRESS