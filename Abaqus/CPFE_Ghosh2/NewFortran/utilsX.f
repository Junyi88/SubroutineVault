      SUBROUTINE KMLT54991(DMIN,DVIN,DVOUT)
      INCLUDE 'aba_param.inc'
      PARAMETER (M=54,N=9)
      DIMENSION DMIN(M,N),DVIN(N),DVOUT(M)

      DO i=1,54
        x = 0.
        DO j=1,9
         x = x + DMIN(i,j)*DVIN(j)
        END DO
        
        DVOUT(i) = x
      END DO

      RETURN
      END