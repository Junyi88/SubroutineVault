
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
!
      INCLUDE 'ABA_PARAM.INC'
!
      DIMENSION TIME(2)
!
  
c       write(6,*) "MWTF---------------------------"
      IF (LOP == 0 ) THEN                   
          call MutexInit( 1 )      ! initialize Mutex #1
          call MutexInit( 2 )
          call MutexInit( 3 )
          call MutexInit( 4 )
          call MutexInit( 5 )
c       write(6,*) "MUTEXINIT---------------------------"
      END IF
    


      RETURN
      END