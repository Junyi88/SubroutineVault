!-------------------------------------------------
!start in xy plane; 2d used (-1,-1,-1), -ve z-plane as start: causes no change at all!
!     Natural coordinates  
      xnat(4,:) = (/-1.,-1.,-1./)
      xnat(8,:) = (/1.,-1.,-1./)
      xnat(7,:) = (/1.,1.,-1./)
      xnat(3,:) = (/-1.,1.,-1./)
      xnat(2,:) = (/-1.,-1.,1./)
      xnat(6,:) = (/1.,-1.,1./)
      xnat(5,:) = (/1.,1.,1./)
      xnat(1,:) = (/-1.,1.,1./)
      xnat(19,:) = (/0.,-1.,-1./)
      xnat(15,:) = (/1.,0.,-1./)
      xnat(20,:) = (/0.,1.,-1./)
      xnat(11,:) = (/-1.,0.,-1./)
      xnat(18,:) = (/0.,-1.,1./)
      xnat(13,:) = (/1.,0.,1./)
      xnat(17,:) = (/0.,1.,1./)
      xnat(9,:) = (/-1.,0.,1./)
      xnat(10,:) = (/-1.,-1.,0./)
      xnat(14,:) = (/1.,-1.,0./)
      xnat(16,:) = (/1.,1.,0./)
      xnat(12,:) = (/-1.,1.,0./)  
!=================================================      
!      do i = 1,8
!      gauss(i,:)=(/xnat(i,1),xnat(i,2),xnat(i,3)/)*xgauss
!      end do 

!!LAST BEST       
!      gauss(1,:) = (/-xgauss,xgauss,-xgauss/)  
!      gauss(2,:) = (/xgauss,xgauss,-xgauss/) 
!      gauss(3,:) = (/-xgauss,-xgauss,-xgauss/)  
!      gauss(4,:) = (/xgauss,-xgauss,-xgauss/) 
!      gauss(5,:) = (/-xgauss,xgauss,xgauss/)
!      gauss(6,:) = (/xgauss,xgauss,xgauss/)
!      gauss(7,:) = (/-xgauss,-xgauss,xgauss/) 
!      gauss(8,:) = (/xgauss,-xgauss,xgauss/)
            
!!LAST BEST: tied to the corresponding nodes to avoid changing them with node changes!     
!      gauss(1,:) = (/xnat(4,1),xnat(4,2),xnat(4,3)/)*xgauss 
!      gauss(2,:) = (/xnat(3,1),xnat(3,2),xnat(3,3)/)*xgauss 
!      gauss(3,:) = (/xnat(1,1),xnat(1,2),xnat(1,3)/)*xgauss  
!      gauss(4,:) = (/xnat(2,1),xnat(2,2),xnat(2,3)/)*xgauss 
!      gauss(5,:) = (/xnat(8,1),xnat(8,2),xnat(8,3)/)*xgauss 
!      gauss(6,:) = (/xnat(7,1),xnat(7,2),xnat(7,3)/)*xgauss 
!      gauss(7,:) = (/xnat(5,1),xnat(5,2),xnat(5,3)/)*xgauss 
!      gauss(8,:) = (/xnat(6,1),xnat(6,2),xnat(6,3)/)*xgauss      

!================================================= 
!Switch 1&2, 5&6: No!
!Reflect: This gives exactly the results before. It gets one of the faces in logical order!      
      gauss(1,:) = (/xnat(1,1),xnat(1,2),xnat(1,3)/)*xgauss  
      gauss(2,:) = (/xnat(2,1),xnat(2,2),xnat(2,3)/)*xgauss
      gauss(3,:) = (/xnat(3,1),xnat(3,2),xnat(3,3)/)*xgauss
      gauss(4,:) = (/xnat(4,1),xnat(4,2),xnat(4,3)/)*xgauss      
      gauss(5,:) = (/xnat(5,1),xnat(5,2),xnat(5,3)/)*xgauss 
      gauss(6,:) = (/xnat(6,1),xnat(6,2),xnat(6,3)/)*xgauss 
      gauss(7,:) = (/xnat(7,1),xnat(7,2),xnat(7,3)/)*xgauss  
      gauss(8,:) = (/xnat(8,1),xnat(8,2),xnat(8,3)/)*xgauss    
              
!================================================= 
        
