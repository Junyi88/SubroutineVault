      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
!
      INCLUDE 'ABA_PARAM.INC'
!
      DIMENSION TIME(2)
!
      INTEGER, PARAMETER  :: K=3,M=3,N=8,nnodes=8,knsdv=18
      INTEGER, INTENT(IN) :: NOEL, NPT , I 
      include 'UserParameters.f'       
      REAL*8 :: kgausscoords, kFp, kcurlFp      
      COMMON/UMPS/kgausscoords(TOTALELEMENTNUM,8,3),
     1 kFp(TOTALELEMENTNUM,8, 9),
     1 kcurlFp(TOTALELEMENTNUM, 8, 9)  
	 
	 
c       write(6,*) "MWTF---------------------------"
      IF (LOP == 0 ) THEN                   
          call MutexInit( 1 )      ! initialize Mutex #1
          call MutexInit( 2 )
          call MutexInit( 3 )
          call MutexInit( 4 )
          call MutexInit( 5 )
          call MutexInit( 6 )
          call MutexInit( 7 )
c       write(6,*) "MUTEXINIT---------------------------"
      END IF
    
c       write(6,*) "MWTF---------------------------"
      IF (LOP == 1 ) THEN                   
          call PerformCurl(NOEL,NPT) 
      END IF
	  
c       write(6,*) "MWTF---------------------------"
      IF (LOP == 4 ) THEN                   

          call MutexLock( 7 )
          open (200, file = './FPTotal.txt', status = 'old')
          do NOEL=1,TOTALELEMENTNUM
		    do NPT = 1,8
			 do I = 1,9
             read(200,*) kFp(NOEL, NPT, I)
		    end do
		    end do
          end do
          close(200)

          END IF
          PRINT *, 'READ DONE = ', GND(2,1)	  
          call MutexUnLock( 7 )    
      END IF
	  
c       write(6,*) "MWTF---------------------------"	  
      IF (LRESTART.GT.0 ) THEN                   

          call MutexLock( 7 )
          open (200, file = './FPTotal.txt', status = 'NEW', 
     1 ACCESS='APPEND', ACTION='WRITE')
          do NOEL=1,TOTALELEMENTNUM
		    do NPT = 1,8
			 do I = 1,9
             write(200,*) kFp(NOEL, NPT, I)
		    end do
		    end do
          end do
          close(200)

          END IF  
          call MutexUnLock( 7 )    
      END IF

	  
      RETURN
      END
	  
      include 'PerformCurl.f'