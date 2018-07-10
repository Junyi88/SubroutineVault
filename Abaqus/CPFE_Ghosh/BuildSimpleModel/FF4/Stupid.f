      SUBROUTINE UMAT(stress,statev,ddsdde,sse,spd,scd,
     1 rpl, ddsddt, drplde, drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
	 
      include 'aba_param.inc'
c
#include <SMAAspUserSubroutines.hdr>
      CHARACTER*8 CMNAME
      EXTERNAL F

      INTEGER:: ISYS, IVAL
      real*8:: MODULUS(6,6), DSTRESS(6)
	  
      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     2 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     3 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
       STATEV(1)=noel
       STATEV(2)=npt
       STATEV(3)=coords(1)
       STATEV(4)=coords(2)
       STATEV(5)=coords(3)	
c---------------------------
      MODULUS(1,1)=0.0
      MODULUS(1,2)=0.0
      MODULUS(1,3)=0.0 
      MODULUS(1,4)=0.0
      MODULUS(1,5)=0.0
      MODULUS(1,6)=0.0

      MODULUS(2,1)=0.0
      MODULUS(2,2)=0.0
      MODULUS(2,3)=0.0 
      MODULUS(2,4)=0.0
      MODULUS(2,5)=0.0
      MODULUS(2,6)=0.0

      MODULUS(3,1)=0.0
      MODULUS(3,2)=0.0
      MODULUS(3,3)=0.0 
      MODULUS(3,4)=0.0
      MODULUS(3,5)=0.0
      MODULUS(3,6)=0.0

      MODULUS(4,1)=0.0
      MODULUS(4,2)=0.0
      MODULUS(4,3)=0.0 
      MODULUS(4,4)=0.0
      MODULUS(4,5)=0.0
      MODULUS(4,6)=0.0

      MODULUS(5,1)=0.0
      MODULUS(5,2)=0.0
      MODULUS(5,3)=0.0 
      MODULUS(5,4)=0.0
      MODULUS(5,5)=0.0
      MODULUS(5,6)=0.0

      MODULUS(6,1)=0.0
      MODULUS(6,2)=0.0
      MODULUS(6,3)=0.0 
      MODULUS(6,4)=0.0
      MODULUS(6,5)=0.0
      MODULUS(6,6)=0.0

	  
      DO ISYS=1,6
        DStress(ISYS)=0.0
        DO IVAL=1,6
           DStress(ISYS)=DStress(ISYS)+
     1     MODULUS(ISYS,IVAL)*
     1     DSTRAN(IVAL)	 
        END DO		
      END DO		 

      DO ISYS=1,6
        Stress(ISYS)=Stress(ISYS)+DStress(ISYS)	
      END DO		
	  
      DO ISYS=1,6
        DO IVAL=1,6
           ddsdde(ISYS,IVAL)=MODULUS(ISYS,IVAL)
        END DO		
      END DO	
	  
      DO ISYS=1,nstatv
	  statev(ISYS)=ISYS
      END DO	  
      return
      end subroutine UMAT