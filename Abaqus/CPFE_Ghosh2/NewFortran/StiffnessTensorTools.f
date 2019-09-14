      subroutine ROTATE_COMTEN(QROT,COM_TEN0,COM_TEN1)

      include 'ComplianceTensorMaps.f' 
	  
      real*8, Intent(in) :: QROT(3,3)
      real*8, Intent(in)  :: COM_TEN0(21) 
      real*8, Intent(inout)  :: COM_TEN1(21) 
      Integer :: nLIST	  
      Integer :: nI, nJ, nK, nL
      Integer :: nP, nQ, nR, nS

      DO nLIST=1,21
       COM_TEN1(nLIST)=0.0 
      END DO	
	 
      DO nLIST=1,21
       DO nP=1,3	  
       DO nQ=1,3	  
       DO nR=1,3	
       DO nS=1,3	
       nI=MLIST2FULL(1,nList) 
       nJ=MLIST2FULL(2,nList) 
       nK=MLIST2FULL(3,nList) 
       nL=MLIST2FULL(4,nList) 
	   
       COM_TEN1(nList)=COM_TEN1(nList)+
     1 QROT(nI,nP)*QROT(nJ,nQ)*
     1 QROT(nK,nR)*QROT(nL,nS)*
     1 COM_TEN0(MFULL2LIST(nP,nQ,nR,nS))
       END DO	
       END DO	
       END DO	
       END DO		   
      END DO	 
	 
      return
      end subroutine ROTATE_COMTEN
	  
c ------------------------------------------------	  
      subroutine ROTATE_Vec(QROT,VEC0,VEC1)
  
      include 'ComplianceTensorMaps.f' 
      real*8, Intent(in) :: QROT(3,3)
      real*8, Intent(in) :: VEC0(3) 
      real*8, Intent(inout) :: VEC1(3) 	 
      Integer :: I, J

      DO I=1,3
        VEC1(I) =0.0 
      END DO	

      DO I=1,3
      DO J=1,3
        VEC1(I) =VEC1(I)+ QROT(I,J)*VEC0(J)  
      END DO
      END DO
	  
	
      return
      end subroutine ROTATE_Vec

c ------------------------------------------------	  
      subroutine Get_TfromSN(S0,N0,T0)
       
      Integer :: I, J , K

      real*8, Intent(in) :: S0(54), N0(54) 
      real*8, Intent(inout) :: T0(54) 	 
      
      DO I=1,18
        J=1+(I-1)*3
        K=I*3
        call CROSS_Vec(S0(J:K),N0(J:K),T0(J:K))
      END DO	       
    
	   
      return
      end subroutine Get_TfromSN
	  
c ------------------------------------------------	  
      subroutine CROSS_Vec(S0,N0,T0)
  

      real*8, Intent(in) :: S0(3), N0(3) 
      real*8, Intent(inout) :: T0(3) 	 
      
       T0(1)=S0(2)*N0(3)-S0(3)*N0(2)
       T0(2)=S0(3)*N0(1)-S0(1)*N0(3)
       T0(3)=S0(1)*N0(2)-S0(2)*N0(1)	   
      return
      end subroutine CROSS_Vec