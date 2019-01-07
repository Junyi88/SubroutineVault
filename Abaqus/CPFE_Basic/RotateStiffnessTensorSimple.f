      subroutine RotateStiffnessTensorSimple(ROTM,
     1 Stiff0, StiffR)

      IMPLICIT NONE     

      INTEGER :: I,J,K
      INTEGER, PARAMETER :: NDOF = 6
      INTEGER, PARAMETER :: NDIR = 3
	  
	  
      REAL*8, INTENT(IN) :: ROTM(NDIR,NDIR)
      REAL*8, INTENT(IN) :: Stiff0(3)
      REAL*8, INTENT(OUT) :: StiffR(NDOF,NDOF)

      REAL*8 :: StiffT(NDOF,NDOF)
      REAL*8 :: tSig(NDOF,NDOF), tStr(NDOF,NDOF)
      REAL*8 :: tSiginv(NDOF,NDOF), prod6(NDOF,NDOF)  
      INTEGER :: info2
c ------------------------ 
      
      DO I = 1, NDOF
        DO J = 1, NDOF
            StiffT(I,J) = 0.0
        END DO
      END DO
      
      DO I = 1, NDIR
        StiffT(I,I) = Stiff0(1)
      END DO
      
      StiffT(1,2) = Stiff0(2)
      StiffT(2,1) = Stiff0(2) 
      StiffT(1,3) = Stiff0(2)
      StiffT(3,1) = Stiff0(2) 
      StiffT(2,3) = Stiff0(2)
      StiffT(3,2) = Stiff0(2) 
      
      DO I = 4, NDOF
        StiffT(I,I) = Stiff0(3)
      END DO
      
      CALL rotord4sig(ROTM,tSig)
      CALL rotord4str(ROTM,tStr)
      CALL lapinverse(tSig,NDOF,info2,tSiginv)
      prod6 = matmul(tSiginv,StiffT)      
      StiffR = matmul(prod6,tStr)	  
	  
      return
      end 