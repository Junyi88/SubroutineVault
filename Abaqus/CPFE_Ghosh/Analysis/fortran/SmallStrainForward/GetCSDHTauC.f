      subroutine GetCSDHTauC(TAUPE,TAUSE,TAUCB,
     1 H, RhoCSD, TAUC, 	   
     2 CinS)

C Subroutine to calculate forest parallel and mobile dislocations
      
      implicit none
      
      real*8,intent(in) :: TAUPE(12),TAUSE(12),TAUCB(12)
      real*8,intent(out) ::  H(12), RhoCSD(12), TAUC(18)

      real*8,intent(in) :: CinS(11)

      integer ISLIPS

C ------------------------------------------------------	
      DO ISLIPS=1,12
       TAUPE(ISLIPS)=TAUPE(ISLIPS)*CinS(1)
       TAUSE(ISLIPS)=TAUSE(ISLIPS)*CinS(1)
       TAUCB(ISLIPS)=TAUCB(ISLIPS)*CinS(1)
      END DO

C ------------------------------------------------------	
      DO ISLIPS=1,12
       H(ISLIPS)=CinS(6)*(CinS(7)+
     1	Cins(8)*(TAUPE(ISLIPS)-CinS(ISLIPS)*TAUSE(ISLIPS))+
     2	sqrt((CinS(10)+abs(TauCB(ISLIPS)))*Cins(11)))
      END DO	  
	  
C ------------------------------------------------------	
      DO ISLIPS=1,12
       RhoCSD(ISLIPS)=CinS(12)*exp(-H(ISLIPS)/CinS(13))
      END DO	  	  

C ------------------------------------------------------	
      DO ISLIPS=1,12
       TauC(ISLIPS)=CinS(4)*sqrt(RhoCSD(ISLIPS))
      END DO		  
	  DO ISLIPS=13,18
       TauC(ISLIPS)=CinS(5)
      END DO	
	  
      return
      end subroutine GetCSDHTauC