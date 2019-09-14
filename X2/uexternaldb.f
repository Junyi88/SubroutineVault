      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
!
      INCLUDE 'ABA_PARAM.INC'
!
      DIMENSION TIME(2)
!
   
      IF (LOP == 0 ) THEN         
          call MutexInit( 2 )
          call MutexInit( 3 )
          call MutexInit( 4 )
          call MutexInit( 5 )
          call MutexInit( 6 )
          call MutexInit( 7 )
      END IF
      
      
      
      
  
      RETURN
      END
	  
