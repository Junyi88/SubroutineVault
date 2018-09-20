      subroutine InitialiseUMAT(
     1 PROPS, STATEV, 
     1 NPROPS, NSTATV,
     1 KSTEP, KINC, NOEL, NPT,
     1 ROTM
     1 )
	  
      IMPLICIT NONE    
      INTEGER, PARAMETER :: NDIR = 3 
      include 'UserParameters.f'   
      
      INTEGER, INTENT(IN) :: NPROPS, NSTATV
      INTEGER, INTENT(IN) :: KSTEP, KINC, NOEL, NPT
      REAL*8, INTENT(IN) :: PROPS(NPROPS)
      REAL*8, INTENT(INOUT) :: STATEV(NSTATV)

      REAL*8, INTENT(OUT) :: ROTM(NDIR,NDIR)
 
      INTEGER :: I,J,K, NN
      
      REAL*8 :: kgausscoords, kFp, kcurlFp,      
      COMMON/UMPS/kgausscoords(TOTALELEMENTNUM,8,3),
     1 kFp(TOTALELEMENTNUM,8, 9),
     1 kcurlFp(TOTALELEMENTNUM, 8, 9)
     
c ==============================================
      IF ((KSTEP.LE.1).AND.(KINC.LE.1)) THEN
        
        DO I = 1, NSTATV
            STATEV(I) = 0.0
        END DO
c ---
        DO I = 1, NDIR
            DO J = 1, NDIR
              NN = J + (I-1)*NDIR
              STATEV(NN) = PROPS(NN)
            END DO
        END DO
c ---
        STATEV(10) = 1.0
        STATEV(14) = 1.0    
        STATEV(18) = 1.0     
        
      END IF
     
c ==============================================
      IF (KINC.LE.1) THEN
c ---
        DO I = 1, NDIR
            DO J = 1, NDIR
              NN = J + (I-1)*NDIR
              ROTM(I,J) = STATEV(NN)
            END DO
        END DO
c ---

        IF (NPT.EQ.8) THEN
        CALL MUTEXLOCK(1)
            DO I = 1,9
                kFp(NOEL, NPT, I) = STATEV(9+I)
            END DO
        CALL MUTEXUNLOCK(1)
        END IF
c -------      
      END IF

c ==============================================
      IF ((KSTEP.LE.1).AND.(KINC.LE.1)) THEN       
        DO I = 1,9
            STATEV(18+I) = 0.0
        END DO        
      ELSE 
        DO I = 1,9
            STATEV(18+I) = kcurlFp(NOEL, NPT, I)
        END DO  
      END IF


        

      RETURN
      END 