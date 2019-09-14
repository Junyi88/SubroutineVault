      SUBROUTINE UXFEMNONLOCALWEIGHT(WEIGHT, JELNO, NPT, COORDS, 
     1 CRACKTIPCOORD, NNCRD, RADIUS, KSTEP, KINC, TIME)
C

      INCLUDE 'ABA_PARAM.INC'

C
      DIMENSION TIME(2), COORDS(NNCRD), CRACKTIPCOORD(NNCRD)
      
      WEIGHT = 0
      
      DO K=1, NNCRD
          N = (COORDS(K) - CRACKTIPCOORD(K))
      	WEIGHT = WEIGHT + (N*N)
      END DO
      
      WEIGHT = SQRT(WEIGHT)
      
      WEIGHT = 1/(WEIGHT*WEIGHT)
      
      write(*,*) "Weight calculated: ",WEIGHT
      
      RETURN
      END