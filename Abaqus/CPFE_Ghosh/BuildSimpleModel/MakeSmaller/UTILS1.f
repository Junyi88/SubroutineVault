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