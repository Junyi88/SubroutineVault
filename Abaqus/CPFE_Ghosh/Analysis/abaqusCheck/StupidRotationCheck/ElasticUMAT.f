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
	 
      INTEGER I,J, NELAS, NCOUNT
      INTEGER, PARAMETER:: NELASMAP(6,6)=reshape([
     1  1, 7, 8, 9, 10, 11, 
     2  7, 2, 12, 13, 14, 15, 
     3  8, 12, 3, 16, 17, 18, 
     4  9, 13, 16, 4, 19, 20, 
     5  10, 14, 17, 19, 5, 21, 
     6  11, 15, 18, 20, 10, 6 
     1  ], [6,6])
c ------------------------------------------------
      Real*8:: DStress(6)
	  
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

	  




