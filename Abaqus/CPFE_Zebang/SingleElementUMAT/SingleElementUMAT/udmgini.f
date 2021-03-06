      SUBROUTINE UDMGINI(FINDEX,NFINDEX,FNORMAL,NDI,NSHR,NTENS,PROPS,
     1 NPROPS,STATEV,NSTATEV,STRESS,STRAIN,STRAINEE,LXFEM,TIME,
     2 DTIME,TEMP,DTEMP,PREDEF,DPRED,NFIELD,COORDS,NOEL,NPT,
     3 KLAYER,KSPT,KSTEP,INC,KDIRCYC,KCYCLELCF,TIMECYC,SSE,SPD,
     4 SCD,SVD,SMD,JMAC,JMATYP,MATLAYO,LACCFLA,CELENT,DROT,ORI)
C
      INCLUDE 'ABA_PARAM.INC'
CC
      DIMENSION FINDEX(NFINDEX),FNORMAL(NDI,NFINDEX),COORDS(*),
     1 STRESS(NTENS),STRAIN(NTENS),STRAINEE(NTENS),PROPS(NPROPS), 
     2 STATEV(NSTATEV),PREDEF(NFIELD),DPRED(NFIELD),TIME(2),
     3 JMAC(*),JMATYP(*),DROR(3,3),ORI(3,3)
      
      DIMENSION PS(3), AN(3,3), WT(6)
      PS(1)=0.0
      PS(2)=0.0
      PS(3)=0.0
C
C ROTATE THE STRESS TO GLOBAL SYSTEM IF THERE IS ORIENTATION
C
      CALL ROTSIG(STRESS,ORI,WT,1,NDI,NSHR)
C
C GET MAXIMUM PRINCIPAL STRESS DIRECTION FOR PROPAGATION
C
      CALL SPRIND(WT,PS,AN,1,NDI,NSHR)
      SIG1 = PS(1)
      KMAX=1
      DO K1 = 2, NDI
         IF(PS(K1).GT.SIG1) THEN
            SIG1 = PS(K1)
            KMAX = K1
         END IF
      END DO
      
      DO K1=1, NDI
      	FNORMAL(K1,1) = AN(KMAX,K1)
      END DO

C THE STORED ENERGY IS CALCULATED IN KMAT, CRITICAL VALUE DEFINED IN INPUT DECK      
      G = STATEV(126)
      FINDEX(1) = G / PROPS(1)
      
      RETURN
      END