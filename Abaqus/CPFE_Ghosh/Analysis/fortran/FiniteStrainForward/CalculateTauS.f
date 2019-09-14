      subroutine CalculateTauS(STRESS, TAU, TAUPE, TAUSE, TAUCB,
     +  SLIP_S, SLIP_N,
     +  SLIP_SPE, SLIP_NPE,
     +  SLIP_SSE, SLIP_NSE,
     +  SLIP_SCB, SLIP_NCB)	 

C Subroutine Calculating All Values of ResolveShearStress
      
      implicit none
      real*8,intent(in) :: STRESS(6)
      real*8,intent(in) :: SLIP_S(54), SLIP_N(54)
      real*8,intent(in) :: SLIP_SPE(36), SLIP_NPE(36)
      real*8,intent(in) :: SLIP_SSE(36), SLIP_NSE(36)	  
      real*8,intent(in) :: SLIP_SCB(36), SLIP_NCB(36)	 
	 	  
      real*8,intent(out) :: TAU(18), TAUPE(12), 
     1 TAUSE(12), TAUCB(12)
      integer ISLIPS, ISLIPSX, I_N, I_S, I_OFF
	  
 
	  
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
      DO ISLIPS=1,18     
       I_OFF=3*(ISLIPS-1)
       TAU(ISLIPS)=TAU(ISLIPS)+
     1 (SLIP_N(1+I_OFF)*SLIP_S(1+I_OFF))*STRESS(1)+
     2 (SLIP_N(2+I_OFF)*SLIP_S(2+I_OFF))*STRESS(2)+
     3 (SLIP_N(3+I_OFF)*SLIP_S(3+I_OFF))*STRESS(3)+
     4 (SLIP_N(1+I_OFF)*SLIP_S(2+I_OFF)+
     +   SLIP_N(2+I_OFF)*SLIP_S(1+I_OFF))*STRESS(4)+
     5 (SLIP_N(1+I_OFF)*SLIP_S(3+I_OFF)+
     +   SLIP_N(3+I_OFF)*SLIP_S(1+I_OFF))*STRESS(5)+
     6 (SLIP_N(2+I_OFF)*SLIP_S(3+I_OFF)+
     +   SLIP_N(3+I_OFF)*SLIP_S(2+I_OFF))*STRESS(6)
      END DO    

	  
C --------------------------------------------------------------------------
C Calculate TauPE	  
      DO ISLIPS=1,12
       I_OFF=3*(ISLIPS-1)	  
       TAUPE(ISLIPS)=TAUPE(ISLIPS)+
     1 (SLIP_NPE(1+I_OFF)*SLIP_SPE(1+I_OFF))*STRESS(1)+
     2 (SLIP_NPE(2+I_OFF)*SLIP_SPE(2+I_OFF))*STRESS(2)+
     3 (SLIP_NPE(3+I_OFF)*SLIP_SPE(3+I_OFF))*STRESS(3)+
     4 (SLIP_NPE(1+I_OFF)*SLIP_SPE(2+I_OFF)+
     +   SLIP_NPE(2+I_OFF)*SLIP_SPE(1+I_OFF))*STRESS(4)+
     5 (SLIP_NPE(1+I_OFF)*SLIP_SPE(3+I_OFF)+
     +   SLIP_N(3+I_OFF)*SLIP_SPE(1+I_OFF))*STRESS(5)+
     6 (SLIP_NPE(2+I_OFF)*SLIP_SPE(3+I_OFF)+
     +   SLIP_NPE(3+I_OFF)*SLIP_SPE(2+I_OFF))*STRESS(6)
      END DO    
	  
C --------------------------------------------------------------------------
C Calculate TauSE
      DO ISLIPS=1,12     
       I_OFF=3*(ISLIPS-1)	 
       TAUSE(ISLIPS)=TAUSE(ISLIPS)+
     1 (SLIP_NSE(1+I_OFF)*SLIP_SSE(1+I_OFF))*STRESS(1)+
     2 (SLIP_NSE(2+I_OFF)*SLIP_SSE(2+I_OFF))*STRESS(2)+
     3 (SLIP_NSE(3+I_OFF)*SLIP_SSE(3+I_OFF))*STRESS(3)+
     4 (SLIP_NSE(1+I_OFF)*SLIP_SSE(2+I_OFF)+
     +   SLIP_NSE(2+I_OFF)*SLIP_SSE(1+I_OFF))*STRESS(4)+
     5 (SLIP_NSE(1+I_OFF)*SLIP_SSE(3+I_OFF)+
     +   SLIP_NSE(3+I_OFF)*SLIP_SSE(1+I_OFF))*STRESS(5)+
     6 (SLIP_NSE(2+I_OFF)*SLIP_SSE(3+I_OFF)+
     +   SLIP_NSE(3+I_OFF)*SLIP_SSE(2+I_OFF))*STRESS(6)
      END DO   	  

C --------------------------------------------------------------------------
C Calculate TauCB
      DO ISLIPS=1,12    
       I_OFF=3*(ISLIPS-1)	 	  
       TAUCB(ISLIPS)=TAUCB(ISLIPS)+
     1 (SLIP_NCB(1+I_OFF)*SLIP_SCB(1+I_OFF))*STRESS(1)+
     2 (SLIP_NCB(2+I_OFF)*SLIP_SCB(2+I_OFF))*STRESS(2)+
     3 (SLIP_NCB(3+I_OFF)*SLIP_SCB(3+I_OFF))*STRESS(3)+
     4 (SLIP_NCB(1+I_OFF)*SLIP_SCB(2+I_OFF)+
     +   SLIP_NCB(2+I_OFF)*SLIP_SCB(1+I_OFF))*STRESS(4)+
     5 (SLIP_NCB(1+I_OFF)*SLIP_SCB(3+I_OFF)+
     +   SLIP_NCB(3+I_OFF)*SLIP_SCB(1+I_OFF))*STRESS(5)+
     6 (SLIP_NCB(2+I_OFF)*SLIP_SCB(3+I_OFF)+
     +   SLIP_NCB(3+I_OFF)*SLIP_SCB(2+I_OFF))*STRESS(6)
      END DO   
	  
      return
      end subroutine CalculateTauS