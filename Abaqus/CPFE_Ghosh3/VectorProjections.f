      subroutine VectorProjections(CosProjection,SinProjection,vA,vB)

C Subroutine for vector projection
C This is projection of vB onto vA for the cosine and sine values
      
      implicit none
      real*8,intent(in) :: vA(3),vB(3)
      real*8,intent(out) :: CosProjection, SinProjection
      real*8 :: MagA
	  
      MagA=vA(1)*vA(1)+vA(2)*vA(2)+vA(3)*vA(3)
      
      IF (MagA.GE.1.0e-12) THEN
      MagA=sqrt(MagA)
      CosProjection=vB(1)*vA(1)+vB(2)*vA(2)+vB(3)*vA(3)
      CosProjection=abs(CosProjection)/MagA
      SinProjection=sqrt(1.0-CosProjection*CosProjection)	
      
      IF ((CosProjection*CosProjection).GT.1.0) THEN
      CosProjection=1.0
      SinProjection=0.0
      ENDIF
      

      ELSE
      CosProjection=1.0
      SinProjection=0.0
      ENDIF 
	  	  
      return
      end subroutine VectorProjections