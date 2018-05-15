      subroutine CalculateDRho(Fp,dtime,gammadot,SLIP_S,SLIP_N,SLIP_T,
     1 dRhoS,dRhoET,dRhoEN,Burgers,gausscoords,noel,npt)
	  
      real*8,intent(in) :: dTIME, BURGERS 
      real*8,intent(in) :: GammaDot(18)
      real*8,intent(in) :: SLIP_S(54),SLIP_N(54),SLIP_T(54)	  
      real*8,intent(out) :: dRhoS(18),dRhoET(18),dRhoEN(18)
	  integer, intent(in) :: noel, npt
	  
      real*8 :: CURLCOM(54), SVARS(6), DGA(18),SVARSFULL(48)
      integer ISLIPS, IVAL, ISYS, I, J, K, L, ICOR
      real*8 :: BurgersI
      real*8 :: gausscoords(3,8)
	 
      real*8:: xnat8(8,3),xnat(20,3),gauss(8,3)
	  real*8 :: kgausscoords, kFp, kcurlFp
c XDANGER
      COMMON/UMPS/kgausscoords(1,8,3),kFp(1,8, 3),
     1 kcurlFp(1, 8, 3)	  	  
	  
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
		SVARS(I)=SVARS(I)+DGA(ISLIPS)*FP(3*(I-1)+J)*SLIPN(ICOR+J)
        END DO
        END DO
c-----------		
       call MutexLock( 1 )      ! lock Mutex #1 
      DO i=1,3                                                      
          kFp(noel,npt,i)= SVARS(I)
      END DO
      call MutexUnlock( 1 )      ! lock Mutex #1 
         DO kint =1,8
             
             DO i=1,3         
                 gausscoords(i,kint) = kgausscoords(noel,kint,i)                          
             END DO
         
             DO i=1,3          
                 SVARSFULL(i + 6*(kint-1)) = kFp(noel,kint,i)         
             END DO
         END DO
c-------------------------	 		
		
        call VectorCurl(svars,xnat8,gauss,gausscoords) 
c---------------
      call MutexLock( 3 )      ! lock Mutex #1 
      DO kint =1, 8
          DO i=1, 3
              kcurlFp(noel,kint,i) =SVARSFULL(3+i + 6*(kint-1))
          END DO
      END DO
      call MutexUnlock( 3 )      ! lock Mutex #1 
      DO i=1,3
          svars(I+3) = kcurlfp(noel,npt,i)
      END DO	  

c-------------
        DO I=1,3
		CURLCOM(ICOR+I)=svars(I+3)
        END DO

		      
      END DO
c-------------------------	  
      DO ISLIPS=1,18
	    ICOR=(ISLIPS-1)*3
        DO I=1,3
         dRhoS(ICOR)=dRhoS(ICOR)+CURLCOM(ICOR+I)*SLIP_S(ICOR+I)*BurgersI
		 dRhoET(ICOR)=dRhoET(ICOR)+CURLCOM(ICOR+I)*SLIP_T(ICOR+I)*BurgersI
		 dRhoEN(ICOR)=dRhoEN(ICOR)+CURLCOM(ICOR+I)*SLIP_N(ICOR+I)*BurgersI
        END DO
    
      END DO	  
	  
	  
      return
      end subroutine CalculateDRho