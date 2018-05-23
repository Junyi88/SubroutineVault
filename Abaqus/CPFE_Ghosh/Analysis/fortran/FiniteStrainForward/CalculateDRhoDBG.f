      subroutine CalculateDRhoDBG(Fp,dtime,gammadot,SLIP_S,SLIP_N,
     1 SLIP_T,
     1 dRhoS,dRhoET,dRhoEN,Burgers,gausscoords,noel,npt,DBGOUT)
      
	  
      implicit none
      real*8,intent(in) :: dTIME, BURGERS 
      real*8,intent(in) :: GammaDot(18),FP(9)
      real*8,intent(in) :: SLIP_S(54),SLIP_N(54)
      real*8,intent(in) ::  SLIP_T(54)	  
      real*8,intent(out) :: dRhoS(18),dRhoET(18),dRhoEN(18)
      integer, intent(in) :: noel, npt
      real*8,intent(out) :: DBGOUT(150)
	  
      real*8 :: CURLCOM(54), SVARS(6), DGA(18),SVARSFULL(48)
      integer ISLIPS, IVAL, ISYS, I, J, K, L, ICOR
      real*8 :: BurgersI
      real*8 :: gausscoords(3,8)

      integer, parameter :: TOTALELEMENTNUM=1000
      real*8,parameter  :: xgauss = 0.577350269189626
      integer :: kint
	  
      real*8:: xnat8(8,3),xnat(20,3),gauss(8,3)
      real*8 :: kgausscoords, kFp, kcurlFp
c XDANGER
      COMMON/UMPS/kgausscoords(TOTALELEMENTNUM,8,3),
     1 kFp(TOTALELEMENTNUM,8, 3),
     1 kcurlFp(TOTALELEMENTNUM, 8, 3)	  
	  
      INCLUDE 'kgauss.f'     	  
      xnat8 = xnat(1:8,:) 
      

	  
	  BurgersI=1.0/Burgers
c-------------------------
      DO ISLIPS=1,18
	        DGA(ISLIPS)=GammaDot(ISLIPS)*DTIME
			dRhoS(ISLIPS)=0.0
			dRhoET(ISLIPS)=0.0
			dRhoEN(ISLIPS)=0.0
      END DO
      DO ISLIPS=1,54
			CURLCOM(ISLIPS)=0.0
      END DO
c-------------------------
 
      DO ISLIPS=1,18
	    ICOR=(ISLIPS-1)*3
        DO I=1,3
		SVARS(I)=0.0
        END DO
		
        DO I=1,3
        DO J=1,3
		SVARS(I)=SVARS(I)+DGA(ISLIPS)*FP(3*(I-1)+J)*SLIP_N(ICOR+J)
        END DO
        END DO
c-----------		
c-----------		
       call MutexLock( 1 )      ! lock Mutex #1 
      DO i=1,3                                                      
          kFp(noel,npt,i)= SVARS(I)
      END DO
      call MutexUnlock( 1 )      ! lock Mutex #1 
         DO kint =1,8
             
             DO i=1,3         
                 gausscoords(i,kint) = kgausscoords(noel,kint,i) 
c                 IF (ISLIPS.EQ.1) THEN 
c                 DBGOUT(I+(kint-1)*3)=gausscoords(i,kint)
c                 END IF
             END DO
         
             DO i=1,3          
                 SVARSFULL(i + 6*(kint-1)) = kFp(noel,kint,i)   	 
             END DO
         END DO
		 
		 
         DO kint =1,8	
             DO i=1,3 		 
                 IF (ISLIPS.EQ.1) THEN 
                 DBGOUT(I+(kint-1)*3)=gausscoords(i,kint)
                 DBGOUT(24+I+(kint-1)*3)=SVARSFULL(i + 6*(kint-1))
                 DBGOUT(48+I+(kint-1)*3)=gauss(kint,i)
                 DBGOUT(72+I+(kint-1)*3)=xnat8(kint,i)			 
                 END IF			
             END DO
         END DO			
	 
c-------------------------	 
       call VectorCurl(SVARSFULL,xnat8,gauss,gausscoords) 
	   
         DO kint =1,8	
             DO i=1,6 		 
                 IF (ISLIPS.EQ.1) THEN 
                 DBGOUT(96+I+(kint-1)*6)=SVARSFULL(i + 6*(kint-1))		 
                 END IF			
             END DO
         END DO
c---------------
      call MutexLock( 3 )      ! lock Mutex #1 
      DO kint =1, 8
          DO i=1, 3
              kcurlFp(noel,kint,i) =SVARSFULL(3+i + 6*(kint-1))
c                 IF (ISLIPS.EQ.1) THEN 
c                 DBGOUT(I+(kint-1)*3)=SVARSFULL(3+i + 6*(kint-1))
c                 END IF	
          END DO
      END DO
      call MutexUnlock( 3 )      ! lock Mutex #1 
      DO i=1,3
          svars(I+3) = kcurlfp(noel,npt,i)
          IF (ISLIPS.EQ.1) THEN 
          DBGOUT(144+I)=kcurlfp(noel,npt,i)
          END IF
      END DO		  
c---------------	  
	  
        DO I=1,3
		CURLCOM(ICOR+I)=svars(I+3)
          IF (ISLIPS.EQ.1) THEN 
          DBGOUT(147+I)=CURLCOM(ICOR+I)
          END IF
        END DO
	  
	  
	  
c---------------
      END DO
c-------------------------	 		
c-------------------------	  
      DO ISLIPS=1,18
        ICOR=(ISLIPS-1)*3
        DO I=1,3
c         dRhoS(ISLIPS)=dRhoS(ISLIPS)+CURLCOM(ICOR+I)*SLIP_S(ICOR+I)*BurgersI
c         dRhoET(ISLIPS)=dRhoET(ISLIPS)+CURLCOM(ICOR+I)*SLIP_T(ICOR+I)*BurgersI
c         dRhoEN(ISLIPS)=dRhoEN(ISLIPS)+CURLCOM(ICOR+I)*SLIP_N(ICOR+I)*BurgersI
		 
         dRhoS(ISLIPS)=dRhoS(ISLIPS)+CURLCOM(ICOR+I)*BurgersI
     1  *SLIP_S(ICOR+I)
         dRhoET(ISLIPS)=dRhoET(ISLIPS)+CURLCOM(ICOR+I)*BurgersI
     1  *SLIP_T(ICOR+I)
         dRhoEN(ISLIPS)=dRhoEN(ISLIPS)+CURLCOM(ICOR+I)*BurgersI
     1  *SLIP_N(ICOR+I)
        END DO
    
      END DO			
 
      DO I=1,54
       DBGOUT(I)=SLIP_T(I)

      END DO		  
	  
      return
      end subroutine CalculateDRhodbg