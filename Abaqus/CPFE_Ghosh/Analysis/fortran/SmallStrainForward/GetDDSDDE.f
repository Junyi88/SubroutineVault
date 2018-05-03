      subroutine GetDDSDDE(DDSDDE,Stress,	   
     2 CinS)

C Subroutine to calculate forest parallel and mobile dislocations
      
      implicit none
      

      real*8,intent(in) :: Stress(6)
      real*8,intent(inout) :: DDSDDE(6,6)
      real*8:: CinS(3)

	  
C ------------------------------------------------------	
      DDSDDE(1,1)=CINS(1)-STRESS(1)
      DDSDDE(1,2)=CINS(2)-STRESS(1)
      DDSDDE(1,3)=CINS(2)-STRESS(1)

      DDSDDE(2,1)=CINS(2)-STRESS(2)
      DDSDDE(2,2)=CINS(1)-STRESS(2)
      DDSDDE(2,3)=CINS(2)-STRESS(2)

      DDSDDE(3,1)=CINS(2)-STRESS(3)
      DDSDDE(3,2)=CINS(2)-STRESS(3)
      DDSDDE(3,3)=CINS(1)-STRESS(3)
	  
C ------------------------------------------------------	
      DDSDDE(4,1)=-STRESS(4)
      DDSDDE(4,2)=-STRESS(4)
      DDSDDE(4,3)=-STRESS(4)
      DDSDDE(4,4)=CINS(3)

      DDSDDE(5,1)=-STRESS(5)
      DDSDDE(5,2)=-STRESS(5)
      DDSDDE(5,3)=-STRESS(5)
      DDSDDE(5,5)=CINS(3)

      DDSDDE(6,1)=-STRESS(6)
      DDSDDE(6,2)=-STRESS(6)
      DDSDDE(6,3)=-STRESS(6)
      DDSDDE(6,6)=CINS(3)
	  
      return
      end subroutine GetDDSDDE