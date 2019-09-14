      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
!
      INCLUDE 'ABA_PARAM.INC'
!
      DIMENSION TIME(2)
!
      INTEGER, PARAMETER  :: K=3,M=3,N=8,nnodes=8,knsdv=18
      INTEGER :: NOEL, NPT , I 
      include 'UserParameters.f'       
      REAL*8 :: kgausscoords, kFp, kcurlFp      
      COMMON/UMPS/kgausscoords(TOTALELEMENTNUM,8,3),
     1 kFp(TOTALELEMENTNUM,8, 9),
     1 kcurlFp(TOTALELEMENTNUM, 8, 9)  
	 
      write(6,*) "A0"	 
c       write(6,*) "MWTF---------------------------"
      IF (LOP == 0 ) THEN         
c      include 'HeaderUMAT.f'   
C           open (201, file = '/home/jl1908/projects/scrap/X/LogLog.txt', 
C      1 status = 'NEW', 
C      1 ACCESS='APPEND', ACTION='WRITE') 	  
          call MutexInit( 1 )      ! initialize Mutex #1
          call MutexInit( 2 )
          call MutexInit( 3 )
          call MutexInit( 4 )
          call MutexInit( 5 )
          call MutexInit( 6 )
          call MutexInit( 7 )
c       write(6,*) "MUTEXINIT---------------------------"
      END IF
      write(6,*) "A1"	    
c       write(6,*) "MWTF---------------------------"
      IF (LOP == 1 ) THEN                   
c          call PerformCurl(NOEL,NPT) 
      END IF
      write(6,*) "A2"		  
c       write(6,*) "MWTF---------------------------"
      IF (LOP == 4 ) THEN                   

          call MutexLock( 7 )
C           open (200, file = '/home/jl1908/projects/scrap/X/FPTotal.txt',
C      1		  status = 'old')
C           do NOEL=1,TOTALELEMENTNUM
C 		    do NPT = 1,8
C 			 do I = 1,9
C              read(200,*) kFp(NOEL, NPT, I)
C 		    end do
C 		    end do
C           end do
C           close(200)
C           close(201)	 
          call MutexUnLock( 7 )   
          END IF
	  
 
      write(6,*) "A3"	
	  
c       write(6,*) "MWTF---------------------------"	  
      IF (LRESTART.GT.0 ) THEN                   

          call MutexLock( 7 )
C           open (200, file = '/home/jl1908/projects/scrap/X/FPTotal.txt', 
C      1 status = 'NEW', 
C      1 ACCESS='APPEND', ACTION='WRITE')
C           do NOEL=1,TOTALELEMENTNUM
C 		    do NPT = 1,8
C 		    do I = 1,9
C              write(200,*) kFp(NOEL, NPT, I)
C 		    end do
C 		    end do
C           end do
C           close(200)


          call MutexUnLock( 7 )    
      END IF

      write(6,*) "A4"	  
      RETURN
      END
	  
      include 'PerformCurl.f'