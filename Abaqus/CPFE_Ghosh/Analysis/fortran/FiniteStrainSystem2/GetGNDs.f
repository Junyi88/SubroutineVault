      SUBROUTINE GetGNDs(kCurlFP,burger,SLIP_S,SLIP_N,SLIP_T,
     1 dRhoS,dRhoET,dRhoEN)

      implicit none
      integer,parameter::MZ=9,NZ=54  
	  
      real*8,intent(in)::kCurlFP(9),burger
      real*8,intent(in)::SLIP_S(54),SLIP_N(54),SLIP_T(54)
      real*8,intent(out):: dRhoS(18),dRhoET(18),dRhoEN(18)
	  
      integer:: I,J,ISLIP,ICOR,ICORSTART,ICOREND
      integer:: NCOL,MX,NX, ialloc
      real*8:: curlfp(3,3),gv(9),burgX(3),sslip(3),nnorm(3),tnorm(3),
     1 screw(3,3),edgeN(3,3),edgeT(3,3),scol(9),encol(9),etcol(9),
     1 rhoOutput(54)	 
	  
      DOUBLE PRECISION,dimension(:,:),allocatable :: V,U,Sinv,A,AM,
     + Ainv,diagm,tempA, UT	
      REAL*8,dimension(:),allocatable :: S   	 
c===========================================	  

    !=========================================================================
    ! *** DETERMINE DENSITY OF GNDs
    ! Update P.Ashton November 2015
    !=========================================================================  

      DO i=1,3
        DO j=1,3 
         curlfp(i,j) = kCurlFP(i+(j-1)*3)
        END DO
      END DO

      IF(MAXVAL(ABS(curlfp)) <= 1.0e-8) THEN 
c        rhoGND = 0. !;gndab  = 0.;gndapr = 0.;gndapy = 0. 
        drhos =0.0
        drhoen = 0.0
        drhoet = 0.0    
      ELSE
      
    
      ALLOCATE(A(mz,nz),AM(mz,nz),V(nz,nz),S(mz),U(mz,mz),Sinv(nz,mz),
     + Ainv(nz,mz),tempA(nz,mz),UT(mz,mz),STAT=ialloc)     
      A=0.;AM=0.;V=0.;S=0.;U=0.;Sinv=0.;Ainv=0.;tempA=0.;UT=0.      
      gv = reshape(curlfp,(/9/))             

        DO I=1,18 ! 1-12 for fcc type      
            ICORStart=(I-1)*3
            ICOREND=I*3	
            burgX = SLIP_S(ICORStart:ICOREND)!slip direction
            burgX = burgX*burger       
			burgX = burgX
            sslip = SLIP_S(ICORStart:ICOREND)
            nnorm = SLIP_N(ICORStart:ICOREND)
            tnorm = SLIP_T(ICORStart:ICOREND)      
            CALL DyadicProd(sslip,burgX,screw)
            CALL DyadicProd(nnorm,burgX,edgeN)
            CALL DyadicProd(tnorm,burgX,edgeT)    
            CALL Convert2Col(screw,scol)
            CALL Convert2Col(edgeN,encol)
            CALL Convert2Col(edgeT,etcol)     
            DO J=1,9
                A(J,I)    = scol(J)
                A(J,I+18) = encol(J)
                A(J,I+36) = etcol(J)
            END DO                
        END DO   ! Generate current A-matrix loop
  ! ***********************************************************************    
      ! Matrix inversion by singular value decomposition: A+ = [V][S+][U*]
  ! ***********************************************************************   
        CALL SVD(A,U,S,V,MZ,NZ)       
        DO i = 1, ubound(S,1)
            IF (S(i) > 1e-6) THEN
                Sinv(i,i)= 1.0 / S(i)
            END IF
        END DO           
        UT=transpose(U)     
        NCOL=NZ;MX=NZ;NX=MZ 
        CALL KMLTM (V,Sinv,tempA,NCOL,MX,NX)     
        NCOL=MZ;MX=NZ;NX=MZ 
        CALL KMLTM (tempA,UT,Ainv,NCOL,MX,NX) 
        CALL KMLT54991(Ainv,gv,rhoOutput)

        DO i = 1, 18
        drhos(i) = rhoOutput(i)
        drhoen(i) = rhoOutput(i+18)
        drhoet(i) = rhoOutput(i+36)
        END DO  
   
    
      END IF ! if curl(Fp) < 1e-8

      RETURN
      end subroutine GetGNDs