      subroutine RotateStiffnessTensor(ROTM,
     1 Stiff0, StiffR)

      IMPLICIT NONE     

      INTEGER :: I,J,K
      INTEGER, PARAMETER :: NDOF = 6
      INTEGER, PARAMETER :: NDIR = 3
	  
	  
      REAL*8, INTENT(IN) :: ROTM(NDIR,NDIR)
      REAL*8, INTENT(IN) :: Stiff0(NDOF,NDOF)
      REAL*8, INTENT(OUT) :: StiffR(NDOF,NDOF)

      REAL*8 :: StiffT(NDOF,NDOF)
      REAL*8 :: tSig(NDOF,NDOF), tStr(NDOF,NDOF)
      REAL*8 :: tSiginv(NDOF,NDOF), prod6(NDOF,NDOF)  
      INTEGER :: info2
c ------------------------ 
 
      StiffT=RESHAPE(Stiff0, (/NDOF, NDOF/))
      CALL rotord4sig(ROTM,tSig)
      CALL rotord4str(ROTM,tStr)
      CALL lapinverse(tSig,NDOF,info2,tSiginv)
      prod6 = matmul(tSiginv,StiffT)      
      StiffR = matmul(prod6,tStr)	  
	  
      return
      end 