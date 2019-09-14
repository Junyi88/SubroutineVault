      SUBROUTINE UMAT(stress,statev,ddsdde,sse,spd,scd,
     1 rpl, ddsddt, drplde, drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
	 
      include 'aba_param.inc'
c
      CHARACTER*8 CMNAME
      EXTERNAL F

      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     2 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     3 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
      
      PARAMETER (ND=50)
      INTEGER I,J, NELAS, NCOUNT, NN , P, Q
      INTEGER, PARAMETER:: NELASMAP(6,6)=reshape([
     1  1, 7, 8, 9, 10, 11, 
     2  7, 2, 12, 13, 14, 15, 
     3  8, 12, 3, 16, 17, 18, 
     4  9, 13, 16, 4, 19, 20, 
     5  10, 14, 17, 19, 5, 21, 
     6  11, 15, 18, 20, 10, 6 
     1  ], [6,6])
	 
      INTEGER, PARAMETER:: VOIGT(3,3)=reshape([
     1  1, 4, 5, 
     2  4, 2, 6, 
     3  5, 6, 3
     1  ], [3,3])
	 
      INTEGER, PARAMETER:: TOFULL(2,6)=reshape([
     1  1, 1, 
     1  2, 2,  
     1  3, 3, 
     1  1, 2, 
     1  1, 3, 
     1  2, 3
     1  ], [2,6])
c ------------------------------------------------
      Real*8:: DStress(6), OHM(3,3), LL(3,3)
	  
	  
	  
c ------------------------------------------------	  
C
C     CALCULATE VELOCITY GRADIENT FROM DEFORMATION GRADIENT.
C     REFERENCE: Li & al. Acta Mater. 52 (2004) 4859-4875
C     
      Real*8:: FTINV(3,3),STRATE(3,3),VELGRD(3,3),AUX1(3,3),ONEMAT(3,3)
      PARAMETER (ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,SIX=6.0D0)
      DATA NEWTON,TOLER/10,1.D-6/
	  
      CALL ONEM(ONEMAT)     
      CALL ZEROM(FTINV)
      CALL ZEROM(AUX1)
      CALL ZEROM(VELGRD)

	  
      CALL M3INV(DFGRD0,FTINV)
      CALL MPROD(DFGRD1,FTINV,AUX1)
      DO 231 I=1,3
        DO 231 J=1,3
          VELGRD(I,J) = (AUX1(I,J)-ONEMAT(I,J))
231   CONTINUE
      
	  NCOUNT=60
       DO I=1,3
       DO J=1,3
	   NCOUNT=NCOUNT+1
	      STATEV(NCOUNT)=VELGRD(I,J)-VELGRD(J,I)
       END DO
      END DO	
c fffffffffffffffffff		 

	  
      DO I=1,3
       DO J=1,3
	      OHM(I,J)=0.5*(VELGRD(I,J)-VELGRD(J,I))
       END DO
      END DO	
c -----------------------------------------
	  
      DO I=1,6
       DStress(I)=0.0
      END DO	
	  
      DO I=1,6
       DO J=1,6
	      NELAS=NELASMAP(I,J)
          DStress(I)=DStress(I)+PROPS(NELAS)*DSTRAN(J)
	      DDSDDE(I,J)=PROPS(NELAS)
       END DO
      END DO	 
	  
c *******************************************
      DO NN=1,6
	     I=TOFULL(1,NN) 
	     J=TOFULL(2,NN)

       DO Q=1,3
	     NELAS=VOIGT(Q,J)
	     DStress(NN)=DStress(NN)+STRESS(NELAS)*OHM(I,Q)
	     NELAS=VOIGT(I,Q)
	     DStress(NN)=DStress(NN)-STRESS(NELAS)*OHM(Q,J)
       END DO

      END DO
c ------------------------------------------------
      NCOUNT=0		
      DO I=1,3
       DO J=1,3
	      NCOUNT=NCOUNT+1
	      STATEV(NCOUNT)=dfgrd0(I,J)
       END DO
      END DO
	  
c ------------------------------------------------	
      DO I=1,3
       DO J=1,3
	      NCOUNT=NCOUNT+1
	      STATEV(NCOUNT)=dfgrd1(I,J)
       END DO
      END DO
c ------------------------------------------------	
      DO I=1,3
       DO J=1,3
	      NCOUNT=NCOUNT+1
	      STATEV(NCOUNT)=drot(I,J)
       END DO
      END DO	  
	  
c ------------------------------------------------	
      DO I=1,6
	      NCOUNT=NCOUNT+1
	      STATEV(NCOUNT)=STRESS(I)
      END DO		  
c ------------------------------------------------	
      DO I=1,6
	      NCOUNT=NCOUNT+1
	      STATEV(NCOUNT)=STRAN(I)
      END DO		  
c ------------------------------------------------	
      DO I=1,6
	      NCOUNT=NCOUNT+1
	      STATEV(NCOUNT)=DSTRAN(I)
      END DO		  
c ------------------------------------------------	
      DO I=1,6
	      NCOUNT=NCOUNT+1
	      STATEV(NCOUNT)=DSTRESS(I)
      END DO	
	  
c ======================================================	
      DO I=1,6
	      NCOUNT=NCOUNT+1
	      STRESS(I)=STRESS(I)+DSTRESS(I)
	      STATEV(NCOUNT)=STRESS(I)
      END DO	  
c ======================================================	
      STATEV(58)=DTIME
      STATEV(59)=TIME(1)
      STATEV(60)=TIME(2)
c ------------------------------------------------		  
      return
      end subroutine UMAT

	  
C----------------------------------------------------------------------


      SUBROUTINE LUDCMP (A, N, NP, INDX, D)

C-----  LU decomposition

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=200, TINY=1.0E-20)
      DIMENSION A(NP,NP), INDX(N), VV(NMAX)

      D=1.
      DO I=1,N
         AAMAX=0.

         DO J=1,N
            IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
         END DO

         IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
         VV(I)=1./AAMAX
      END DO

      DO J=1,N
         DO I=1,J-1
            SUM=A(I,J)

            DO K=1,I-1
               SUM=SUM-A(I,K)*A(K,J)
            END DO

            A(I,J)=SUM
         END DO
         AAMAX=0.

         DO I=J,N
            SUM=A(I,J)

            DO K=1,J-1
               SUM=SUM-A(I,K)*A(K,J)
            END DO

            A(I,J)=SUM
            DUM=VV(I)*ABS(SUM)
            IF (DUM.GE.AAMAX) THEN
               IMAX=I
               AAMAX=DUM
            END IF
         END DO

         IF (J.NE.IMAX) THEN
            DO K=1,N
               DUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=DUM
            END DO

            D=-D
            VV(IMAX)=VV(J)
         END IF

         INDX(J)=IMAX
         IF (A(J,J).EQ.0.) A(J,J)=TINY
         IF (J.NE.N) THEN
            DUM=1./A(J,J)
            DO I=J+1,N
               A(I,J)=A(I,J)*DUM
            END DO
         END IF

      END DO

      RETURN
      END


C----------------------------------------------------------------------


      SUBROUTINE LUBKSB (A, N, NP, INDX, B)

C-----  Linear equation solver based on LU decomposition

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NP,NP), INDX(N), B(N)

      II=0
      DO I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)

         IF (II.NE.0) THEN
            DO J=II,I-1
               SUM=SUM-A(I,J)*B(J)
            END DO
         ELSE IF (SUM.NE.0.) THEN
            II=I
         END IF

         B(I)=SUM
      END DO

      DO I=N,1,-1
         SUM=B(I)

         IF (I.LT.N) THEN
            DO J=I+1,N
               SUM=SUM-A(I,J)*B(J)
            END DO
         END IF

         B(I)=SUM/A(I,I)
      END DO

      RETURN
      END

C     ***************************************************************
C     UTILITY SUBROUTINES TAKEN FROM KALIDINDI CPFEM UMAT-1992
C     Kalidindi, PhD thesis, MIT, 1992, Cambridge
C     ***************************************************************


     	SUBROUTINE ZEROM(A)
C --  THIS SUBROUTINE SETS ALL ENTRIES OF A 3 BY 3 MATRIX TO 0.0

	IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 A(3,3)

      DO 1 I=1,3
        DO 1 J=1,3
             A(I,J) = 0.0
1	CONTINUE
	
	RETURN
	END
C     END OF SUBROUTINE ZEROM	
C     ---------------------------------------------------------------


	SUBROUTINE ONEM(A)
C --	THIS SUBROUTINE STORES THE IDENTITY MATRIX IN 
C --  THE 3 BY 3 MATRIX [A]

	IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 A(3,3)

	DO 1 I = 1,3
	  DO 1 J = 1,3
	    IF (I .EQ. J) THEN      
	      A(I,J)      = 1.0D0
	      ELSE
	        A(I,J)      = 0.0D0
	    END IF
1	CONTINUE

	RETURN
	END
C     END OF SUBROUTINE ONEM	
C     ---------------------------------------------------------------


	SUBROUTINE MPROD(A,B,C)
C --	THIS SUBROUTINE MULTIPLIES TWO 3 BY 3 MATRICES [A] AND [B],
C --	AND PLACE THEIR PRODUCT IN MATRIX [C]. 

	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 A(3,3), B(3,3), C(3,3)

	DO 2 I = 1, 3
	  DO 2 J = 1, 3
		  C(I,J) = 0.0
		  DO 1 K = 1, 3
			C(I,J) = C(I,J) + A(I,K) * B(K,J)                       
1		CONTINUE
2	CONTINUE
	RETURN
	END
C     END OF SUBROUTINE MPROD	
C -----------------------------------------------------------------


	SUBROUTINE MTRANS(A,ATRANS)
C --	THIS SUBROUTINE CALCULATES THE TRANSPOSE OF AN 3 BY 3 
C --	MATRIX [A], AND PLACES THE RESULT IN ATRANS. 

	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8  A(3,3), ATRANS(3,3)

      CALL ONEM(ATRANS)
	DO 1 I = 1, 3
      	DO 1 J = 1, 3
		  ATRANS(I,J) = A(J,I)
1	CONTINUE
	RETURN
	END
C     END OF SUBROUTINE MTRANS	
C -----------------------------------------------------------------


	SUBROUTINE MDET(A,DET)
C --	THIS SUBROUTINE CALCULATES THE DETERMINANT
C --	OF A 3 BY 3 MATRIX [A].

	IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(3,3)

	DET =	  A(1,1)*A(2,2)*A(3,3) 
     +	+ A(1,2)*A(2,3)*A(3,1)
     +	+ A(1,3)*A(2,1)*A(3,2)
     +	- A(3,1)*A(2,2)*A(1,3)
     +	- A(3,2)*A(2,3)*A(1,1)
     +	- A(3,3)*A(2,1)*A(1,2)

	RETURN
	END
C     END OF SUBROUTINE MDET	
C -----------------------------------------------------------------

	SUBROUTINE M3INV(A,AINV)
C --	THIS SUBROUTINE CALCULATES THE THE INVERSE OF A 3 BY 3 MATRIX [A]
C --	AND PLACES THE RESULT IN [AINV]. 
C --	IF DET(A) IS ZERO, THE CALCULATION
C --	IS TERMINATED AND A DIAGNOSTIC STATEMENT IS PRINTED.

	IMPLICIT REAL*8(A-H,O-Z)	
	REAL*8  A(3,3), AINV(3,3), DET, ACOFAC(3,3), AADJ(3,3)
C
C		A(3,3)	-- THE MATRIX WHOSE INVERSE IS DESIRED.
C		DET		-- THE COMPUTED DETERMINANT OF [A].
C		ACOFAC(3,3)	-- THE MATRIX OF COFACTORS OF A(I,J).
C				   THE SIGNED MINOR (-1)**(I+J)*M_IJ
C				   IS CALLED THE COFACTOR OF A(I,J).
C		AADJ(3,3)	-- THE ADJOINT OF [A]. IT IS THE MATRIX
C				   OBTAINED BY REPLACING EACH ELEMENT OF
C				   [A] BY ITS COFACTOR, AND THEN TAKING
C				   TRANSPOSE OF THE RESULTING MATRIX.
C		AINV(3,3)	-- RETURNED AS INVERSE OF [A].
C				   [AINV] = [AADJ]/DET.

	CALL MDET(A,DET)

	IF ( DET .EQ. 0.0 ) THEN
	    WRITE(91,10)
	    CALL XIT
	END IF

	CALL MCOFAC(A,ACOFAC)
	CALL MTRANS(ACOFAC,AADJ)
	DO 1 I = 1,3
	DO 1 J = 1,3
	     AINV(I,J) = AADJ(I,J)/DET
1	CONTINUE

C	FORMAT
 10	FORMAT(5X,'--ERROR IN M3INV--- THE MATRIX IS SINGULAR',/,
     +       10X,'PROGRAM TERMINATED')

	RETURN
	END
C     END OF SUBROUTINE M3INV	
C -----------------------------------------------------------------


	SUBROUTINE MCOFAC(A,ACOFAC)
C --	THIS SUBROUTINE CALCULATES THE COFACTOR OF A 3 BY 3 MATRIX [A],
C --	AND PLACES THE RESULT IN ACOFAC. 

	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8  A(3,3), ACOFAC(3,3)

	ACOFAC(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
	ACOFAC(1,2) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))
	ACOFAC(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
	ACOFAC(2,1) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3))
	ACOFAC(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
	ACOFAC(2,3) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))
	ACOFAC(3,1) = A(1,2)*A(2,3)  - A(2,2)*A(1,3)
	ACOFAC(3,2) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))
	ACOFAC(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)

	RETURN
	END
C     END OF SUBROUTINE MCOFAC	
C -----------------------------------------------------------------