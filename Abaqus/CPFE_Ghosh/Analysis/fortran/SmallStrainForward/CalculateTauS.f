      subroutine CalculateTauS(STRESS, TAU, TAUPE, TAUSE, TAUCB,
     +  FCC_N,FCC_S,
     +  FCC_NPE,FCC_SPE,
     +  FCC_NSE,FCC_SSE,
     +  FCC_NCB,FCC_SCB,
     +  CUBIC_N,CUBIC_S)

C Subroutine Calculating All Values of ResolveShearStress
      
      implicit none
      real*8,intent(in) :: STRESS(6)
      real*8,intent(in) :: FCC_N(3,12),FCC_S(3,12)
      real*8,intent(in) :: FCC_NPE(3,12),FCC_SPE(3,12)  
      real*8,intent(in) :: FCC_NSE(3,12),FCC_SSE(3,12)	  
      real*8,intent(in) :: FCC_NCB(3,12),FCC_SCB(3,12)	 
      real*8,intent(in) :: CUBIC_N(3,6),CUBIC_S(3,6)	
 	  
      real*8,intent(out) :: TAU(18), TAUPE(12), 
     1 TAUSE(12), TAUCB(12)
      integer ISLIPS, ISLIPSX
	  
 
	  
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
     1 (FCC_N(1,ISLIPS)*FCC_S(1,ISLIPS))*STRESS(1)+
     2 (FCC_N(2,ISLIPS)*FCC_S(2,ISLIPS))*STRESS(2)+
     3 (FCC_N(3,ISLIPS)*FCC_S(3,ISLIPS))*STRESS(3)+
     4 (FCC_N(1,ISLIPS)*FCC_S(2,ISLIPS)+
     +   FCC_N(2,ISLIPS)*FCC_S(1,ISLIPS))*STRESS(4)+
     5 (FCC_N(1,ISLIPS)*FCC_S(3,ISLIPS)+
     +   FCC_N(3,ISLIPS)*FCC_S(1,ISLIPS))*STRESS(5)+
     6 (FCC_N(2,ISLIPS)*FCC_S(3,ISLIPS)+
     +   FCC_N(3,ISLIPS)*FCC_S(2,ISLIPS))*STRESS(6)
      END DO    
      DO ISLIPS=1,6     
       ISLIPSX=ISLIPS+12
       TAU(ISLIPSX)=TAU(ISLIPSX)+
     1 (CUBIC_N(1,ISLIPS)*CUBIC_S(1,ISLIPS))*STRESS(1)+
     2 (CUBIC_N(2,ISLIPS)*CUBIC_S(2,ISLIPS))*STRESS(2)+
     3 (CUBIC_N(3,ISLIPS)*CUBIC_S(3,ISLIPS))*STRESS(3)+
     4 (CUBIC_N(1,ISLIPS)*CUBIC_S(2,ISLIPS)+
     +   CUBIC_N(2,ISLIPS)*CUBIC_S(1,ISLIPS))*STRESS(4)+
     5 (CUBIC_N(1,ISLIPS)*CUBIC_S(3,ISLIPS)+
     +   CUBIC_N(3,ISLIPS)*CUBIC_S(1,ISLIPS))*STRESS(5)+
     6 (CUBIC_N(2,ISLIPS)*CUBIC_S(3,ISLIPS)+
     +   CUBIC_N(3,ISLIPS)*CUBIC_S(2,ISLIPS))*STRESS(6)
      END DO    
	  
C --------------------------------------------------------------------------
C Calculate TauPE	  
      DO ISLIPS=1,12     
       TAUPE(ISLIPS)=TAUPE(ISLIPS)+
     1 (FCC_NPE(1,ISLIPS)*FCC_SPE(1,ISLIPS))*STRESS(1)+
     2 (FCC_NPE(2,ISLIPS)*FCC_SPE(2,ISLIPS))*STRESS(2)+
     3 (FCC_NPE(3,ISLIPS)*FCC_SPE(3,ISLIPS))*STRESS(3)+
     4 (FCC_NPE(1,ISLIPS)*FCC_SPE(2,ISLIPS)+
     +   FCC_NPE(2,ISLIPS)*FCC_SPE(1,ISLIPS))*STRESS(4)+
     5 (FCC_NPE(1,ISLIPS)*FCC_SPE(3,ISLIPS)+
     +   FCC_N(3,ISLIPS)*FCC_SPE(1,ISLIPS))*STRESS(5)+
     6 (FCC_NPE(2,ISLIPS)*FCC_SPE(3,ISLIPS)+
     +   FCC_NPE(3,ISLIPS)*FCC_SPE(2,ISLIPS))*STRESS(6)
      END DO    
	  
C --------------------------------------------------------------------------
C Calculate TauSE
      DO ISLIPS=1,12     
       TAUSE(ISLIPS)=TAUSE(ISLIPS)+
     1 (FCC_NSE(1,ISLIPS)*FCC_SSE(1,ISLIPS))*STRESS(1)+
     2 (FCC_NSE(2,ISLIPS)*FCC_SSE(2,ISLIPS))*STRESS(2)+
     3 (FCC_NSE(3,ISLIPS)*FCC_SSE(3,ISLIPS))*STRESS(3)+
     4 (FCC_NSE(1,ISLIPS)*FCC_SSE(2,ISLIPS)+
     +   FCC_NSE(2,ISLIPS)*FCC_SSE(1,ISLIPS))*STRESS(4)+
     5 (FCC_NSE(1,ISLIPS)*FCC_SSE(3,ISLIPS)+
     +   FCC_NSE(3,ISLIPS)*FCC_SSE(1,ISLIPS))*STRESS(5)+
     6 (FCC_NSE(2,ISLIPS)*FCC_SSE(3,ISLIPS)+
     +   FCC_NSE(3,ISLIPS)*FCC_SSE(2,ISLIPS))*STRESS(6)
      END DO   	  

C --------------------------------------------------------------------------
C Calculate TauCB
      DO ISLIPS=1,12     
       TAUCB(ISLIPS)=TAUCB(ISLIPS)+
     1 (FCC_NCB(1,ISLIPS)*FCC_SCB(1,ISLIPS))*STRESS(1)+
     2 (FCC_NCB(2,ISLIPS)*FCC_SCB(2,ISLIPS))*STRESS(2)+
     3 (FCC_NCB(3,ISLIPS)*FCC_SCB(3,ISLIPS))*STRESS(3)+
     4 (FCC_NCB(1,ISLIPS)*FCC_SCB(2,ISLIPS)+
     +   FCC_NCB(2,ISLIPS)*FCC_SCB(1,ISLIPS))*STRESS(4)+
     5 (FCC_NCB(1,ISLIPS)*FCC_SCB(3,ISLIPS)+
     +   FCC_NCB(3,ISLIPS)*FCC_SCB(1,ISLIPS))*STRESS(5)+
     6 (FCC_NCB(2,ISLIPS)*FCC_SCB(3,ISLIPS)+
     +   FCC_NCB(3,ISLIPS)*FCC_SCB(2,ISLIPS))*STRESS(6)
      END DO   
	  
      return
      end subroutine CalculateTauS