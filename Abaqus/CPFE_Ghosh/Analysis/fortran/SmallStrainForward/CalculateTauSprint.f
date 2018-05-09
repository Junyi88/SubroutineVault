      subroutine CalculateTauS(STRESS, TAU, TAUPE, TAUSE, TAUCB,
     +  FCC_N,FCC_S,
     +  FCC_NPE,FCC_SPE,
     +  FCC_NSE,FCC_SSE,
     +  FCC_NCB,FCC_SCB,
     +  CUBIC_N,CUBIC_S)

C Subroutine Calculating All Values of ResolveShearStress
      
      implicit none
      real*8,intent(in) :: STRESS(6)
      real*8,intent(in) :: FCC_N(12,3),FCC_S(12,3)
      real*8,intent(in) :: FCC_NPE(12,3),FCC_SPE(12,3)  
      real*8,intent(in) :: FCC_NSE(12,3),FCC_SSE(12,3)	  
      real*8,intent(in) :: FCC_NCB(12,3),FCC_SCB(12,3)	 
      real*8,intent(in) :: CUBIC_N(6,3),CUBIC_S(6,3)	
 	  
      real*8,intent(out) :: TAU(18), TAUPE(12), 
     1 TAUSE(12), TAUCB(12)
      integer ISLIPS, ISLIPSX
	  
      print *, '=xxx===================='
      DO ISLIPS=1,12
	     
	     print *, 'S=', FCC_S(ISLIPS,:)
	     print *, 'N=', FCC_N(ISLIPS,:)

      END DO	  
	  
C --------------------------------------------------------------------------
      DO ISLIPS=1,18    
       TAU(ISLIPS)=0.0
      END DO 
      DO ISLIPS=1,12    
       TAUPE(ISLIPS)=0.0 
       TAUSE(ISLIPS)=0.0
       TAUCB(ISLIPS)=0.0
      END DO 
C --------------------------------------------------------------------------
 
C Calculate Tau
      DO ISLIPS=1,12     
       TAU(ISLIPS)=TAU(ISLIPS)+
     1 (FCC_N(ISLIPS,1)*FCC_S(ISLIPS,1))*STRESS(1)+
     2 (FCC_N(ISLIPS,2)*FCC_S(ISLIPS,2))*STRESS(2)+
     3 (FCC_N(ISLIPS,3)*FCC_S(ISLIPS,3))*STRESS(3)+
     4 (FCC_N(ISLIPS,1)*FCC_S(ISLIPS,2)+
     +   FCC_N(ISLIPS,2)*FCC_S(ISLIPS,1))*STRESS(4)+
     5 (FCC_N(ISLIPS,1)*FCC_S(ISLIPS,3)+
     +   FCC_N(ISLIPS,3)*FCC_S(ISLIPS,1))*STRESS(5)+
     6 (FCC_N(ISLIPS,2)*FCC_S(ISLIPS,3)+
     +   FCC_N(ISLIPS,3)*FCC_S(ISLIPS,2))*STRESS(6)
      END DO    
      DO ISLIPS=1,6     
       ISLIPSX=ISLIPS+12
       TAU(ISLIPSX)=TAU(ISLIPSX)+
     1 (CUBIC_N(ISLIPS,1)*CUBIC_S(ISLIPS,1))*STRESS(1)+
     2 (CUBIC_N(ISLIPS,2)*CUBIC_S(ISLIPS,2))*STRESS(2)+
     3 (CUBIC_N(ISLIPS,3)*CUBIC_S(ISLIPS,3))*STRESS(3)+
     4 (CUBIC_N(ISLIPS,1)*CUBIC_S(ISLIPS,2)+
     +   CUBIC_N(ISLIPS,2)*CUBIC_S(ISLIPS,1))*STRESS(4)+
     5 (CUBIC_N(ISLIPS,1)*CUBIC_S(ISLIPS,3)+
     +   CUBIC_N(ISLIPS,3)*CUBIC_S(ISLIPS,1))*STRESS(5)+
     6 (CUBIC_N(ISLIPS,2)*CUBIC_S(ISLIPS,3)+
     +   CUBIC_N(ISLIPS,3)*CUBIC_S(ISLIPS,2))*STRESS(6)
      END DO    
	  
C --------------------------------------------------------------------------
C Calculate TauPE	  
      DO ISLIPS=1,12     
       TAUPE(ISLIPS)=TAUPE(ISLIPS)+
     1 (FCC_NPE(ISLIPS,1)*FCC_SPE(ISLIPS,1))*STRESS(1)+
     2 (FCC_NPE(ISLIPS,2)*FCC_SPE(ISLIPS,2))*STRESS(2)+
     3 (FCC_NPE(ISLIPS,3)*FCC_SPE(ISLIPS,3))*STRESS(3)+
     4 (FCC_NPE(ISLIPS,1)*FCC_SPE(ISLIPS,2)+
     +   FCC_NPE(ISLIPS,2)*FCC_SPE(ISLIPS,1))*STRESS(4)+
     5 (FCC_NPE(ISLIPS,1)*FCC_SPE(ISLIPS,3)+
     +   FCC_N(ISLIPS,3)*FCC_SPE(ISLIPS,1))*STRESS(5)+
     6 (FCC_NPE(ISLIPS,2)*FCC_SPE(ISLIPS,3)+
     +   FCC_NPE(ISLIPS,3)*FCC_SPE(ISLIPS,2))*STRESS(6)
      END DO    
	  
C --------------------------------------------------------------------------
C Calculate TauSE
      DO ISLIPS=1,12     
       TAUSE(ISLIPS)=TAUSE(ISLIPS)+
     1 (FCC_NSE(ISLIPS,1)*FCC_SSE(ISLIPS,1))*STRESS(1)+
     2 (FCC_NSE(ISLIPS,2)*FCC_SSE(ISLIPS,2))*STRESS(2)+
     3 (FCC_NSE(ISLIPS,3)*FCC_SSE(ISLIPS,3))*STRESS(3)+
     4 (FCC_NSE(ISLIPS,1)*FCC_SSE(ISLIPS,2)+
     +   FCC_NSE(ISLIPS,2)*FCC_SSE(ISLIPS,1))*STRESS(4)+
     5 (FCC_NSE(ISLIPS,1)*FCC_SSE(ISLIPS,3)+
     +   FCC_NSE(ISLIPS,3)*FCC_SSE(ISLIPS,1))*STRESS(5)+
     6 (FCC_NSE(ISLIPS,2)*FCC_SSE(ISLIPS,3)+
     +   FCC_NSE(ISLIPS,3)*FCC_SSE(ISLIPS,2))*STRESS(6)
      END DO   	  

C --------------------------------------------------------------------------
C Calculate TauCB
      DO ISLIPS=1,12     
       TAUCB(ISLIPS)=TAUCB(ISLIPS)+
     1 (FCC_NCB(ISLIPS,1)*FCC_SCB(ISLIPS,1))*STRESS(1)+
     2 (FCC_NCB(ISLIPS,2)*FCC_SCB(ISLIPS,2))*STRESS(2)+
     3 (FCC_NCB(ISLIPS,3)*FCC_SCB(ISLIPS,3))*STRESS(3)+
     4 (FCC_NCB(ISLIPS,1)*FCC_SCB(ISLIPS,2)+
     +   FCC_NCB(ISLIPS,2)*FCC_SCB(ISLIPS,1))*STRESS(4)+
     5 (FCC_NCB(ISLIPS,1)*FCC_SCB(ISLIPS,3)+
     +   FCC_NCB(ISLIPS,3)*FCC_SCB(ISLIPS,1))*STRESS(5)+
     6 (FCC_NCB(ISLIPS,2)*FCC_SCB(ISLIPS,3)+
     +   FCC_NCB(ISLIPS,3)*FCC_SCB(ISLIPS,2))*STRESS(6)
      END DO   
	  
      return
      end subroutine CalculateTauS