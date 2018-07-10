      subroutine GetDDSDDE(DDSDDE,Stress,	   
     2 CinS)

C Subroutine to calculate forest parallel and mobile dislocations
      
      implicit none
      

      real*8,intent(in) :: Stress(6)
      real*8,intent(inout) :: DDSDDE(6,6)
      real*8:: CinS(21)
 	  
C ------------------------------------------------------	
      DDSDDE(1,1)=CINS(1) !-STRESS(1)
      DDSDDE(1,2)=CINS(7) !-STRESS(1)
      DDSDDE(1,3)=CINS(8) !-STRESS(1)
      DDSDDE(1,4)=CINS(9)
      DDSDDE(1,5)=CINS(10)
      DDSDDE(1,6)=CINS(11)

C ------------------------------------------------------	
      DDSDDE(2,1)=CINS(7) !-STRESS(2)
      DDSDDE(2,2)=CINS(2) !-STRESS(2)
      DDSDDE(2,3)=CINS(12) !-STRESS(2)
      DDSDDE(2,4)=CINS(13)
      DDSDDE(2,5)=CINS(14)
      DDSDDE(2,6)=CINS(15)

C ------------------------------------------------------	
      DDSDDE(3,1)=CINS(8) !-STRESS(3)
      DDSDDE(3,2)=CINS(12) !-STRESS(3)
      DDSDDE(3,3)=CINS(3) !-STRESS(3)
      DDSDDE(3,4)=CINS(16)
      DDSDDE(3,5)=CINS(17)
      DDSDDE(3,6)=CINS(18)

C ------------------------------------------------------	
      DDSDDE(4,1)=CINS(9) !-STRESS(4)
      DDSDDE(4,2)=CINS(13) !-STRESS(4)
      DDSDDE(4,3)=CINS(16) !-STRESS(4)
      DDSDDE(4,4)=CINS(4)
      DDSDDE(4,5)=CINS(19)
      DDSDDE(4,6)=CINS(20)

C ------------------------------------------------------	
      DDSDDE(5,1)=CINS(10) !-STRESS(5)
      DDSDDE(5,2)=CINS(14) !-STRESS(5)
      DDSDDE(5,3)=CINS(17) !-STRESS(5)
      DDSDDE(5,4)=CINS(19)
      DDSDDE(5,5)=CINS(5)
      DDSDDE(5,6)=CINS(21)

C ------------------------------------------------------	
      DDSDDE(6,1)=CINS(11) !-STRESS(6)
      DDSDDE(6,2)=CINS(15) !-STRESS(6)
      DDSDDE(6,3)=CINS(18) !-STRESS(6)
      DDSDDE(6,4)=CINS(20)
      DDSDDE(6,5)=CINS(21)
      DDSDDE(6,6)=CINS(6)      
	  
      return
      end subroutine GetDDSDDE