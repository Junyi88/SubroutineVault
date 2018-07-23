       SUBROUTINE CREEP(DECRA,DESWA,STATEV,SERD,EC,ESW,P,QTILD
     1 TEMP,DTEMP,PREDEF,DPRED,TIME,DTIME,CMNAME,LEXIMP,LEND,
     2 COORDS,NSTATV,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
       INCLUDE 'ABA_PARAM.INC'
C
       CHARACTER*80 CMNAME
C
       DIMENSION DECRA(5),DESWA(5),STATEV(4),PREDEF(5),DPRED(5),
     1 TIME(2),EC(2),ESW(2),COORDS(5),SMETAL(13),DYDX(3)
c      SMETAL(13),DYDX(3)
       REAL*8 D1,Tab,RIP
c      REAL*8 RIP
     

**-------INTEGER1---------2---------3---------4---------5---------6---------7**  
**        MATCONS - MATERIAL CONSTANTS
**-------1---------2---------3---------4---------5---------6---------7**
*     MATERIAL CONSTANTS     THE MATERIAL CONSTANTS ARE 31 
*      sk=SMETAL(1)  
*      LK=SMETAL(2)
*      n=SMETAL(3)
*      B=SMETAL(4)
*      gamma=SMETAL(5)
*      A1=SMETAL(6)
*      A2=SMETAL(7)
*      c1=SMETAL(8)
*      eta=SMETAL(9)
*      eta1=SMETAL(10)
*      eta2=SMETAL(11)
*      E=SMETAL(12)
*      Rc=SMETAL(13)
*      


*     It depends on the number of mat constants in your equations
*     STATEV(1)  =PLASTIC STRAIN  ?
*     STATEV(2)  =HARDENING H
*     STATEV(3)  =DISLOCATION DENSITY ? 
*     STATEV(4)  =GLOUBULARIZATION VOLUME FRACTION S
*     
*      
*     DYDX(1)  =PLASTIC STRAIN RATE  ?'
*     DYDX(2)  =DISLOCATION DENSITY RATE ?' 
*     DYDX(3)  =GLOBULARIZATION VOLUME FRACTION RATE S'
*     
*  
     
      SMETAL(1)=5.00D-01            !sk=SMETAL(1)         
      SMETAL(2)=2.346984D03             !LK=SMETAL(2)
      SMETAL(3)=1.3               !n=SMETAL(3)
      SMETAL(4)=2.400067D01             !B=SMETAL(4)
      SMETAL(5)=6.2D-1                !gamma=SMETAL(5)
      SMETAL(6)=4.4915D01             !A1=SMETAL(6)
      SMETAL(7)=2.1745D-02             !A2=SMETAL(7)
      SMETAL(8)=1.45D-02             !c1=SMETAL(8)
      SMETAL(9)=8.4D-02            !eta=SMETAL(9)
      SMETAL(10)=3.43898D-01             !eta1=SMETAL(10)
      SMETAL(11)=1.7581D-02            !eta2=SMETAL(11)
      SMETAL(12)=9.711426D03           !E=SMETAL(12)
      SMETAL(13)=5.0D-01            !Rc=SMETAL(13)
      

  


	                               
!**-------1---------2---------3---------4---------5---------6---------7**
    	IF (STATEV(1) .LE. 0.0) THEN
		STATEV(1) = 0.0
	    END IF
       IF (STATEV(2) .LE. 0.0) THEN
		STATEV(2) = 0.0
1         END IF
  !      IF (STATEV(3) .LE. 0.0) THEN
	!	STATEV(3) = 0.0
        ! END IF
	 ! IF (STATEV(4) .LE. 0.0) THEN
		!STATEV(4) = 0.0
	  !  END IF  

!**	CALCULATE Hardening
      
      STATEV(2)=SMETAL(4)*(STATEV(3)**(SMETAL(5)))
   
      
      
!**	Calculate Plastic Strain Rate DYDX(1)
         D1=(QTILD-STATEV(2)-SMETAL(1))
         IF(D1.LE.1.0D-15) THEN
             DYDX(1)=0.0
	        END IF
         DYDX(1)=((D1/SMETAL(2))**SMETAL(3))*(1/(1-STATEV(4)))
            !END   
             
        
!**	CALCULATE DISLOCATION DENSITY RATE DYDX(2)
     
      Tab=SMETAL(7)*(STATEV(3)**SMETAL(8))
      
      DYDX(2)=SMETAL(6)*(1-STATEV(3))*ABS(DYDX(1))-Tab
      IF (STATEV(3).LE.0.0)   STATEV(3)=1.0D-10
      IF (STATEV(3).GE.0.999) STATEV(3)=0.999
	  
      

!**	CALCULATE Recrsytallization DYDX(3)
       IF (STATEV(4).GE.STATEV(3)) STATEV(4)=STATEV(3)
       IF (STATEV(3).LE.SMETAL(13))  STATEV(3)=5.02D-01
       
       RIP=((STATEV(3)-SMETAL(13))**SMETAL(11))*(STATEV(3)-STATEV(4))
       DYDX(3)=SMETAL(9)*(DYDX(1)**SMETAL(10))*RIP
       
       IF (STATEV(4).LE.0.0)  STATEV(4)=1.0D-10
       IF (STATEV(4).GE.0.999) STATEV(4)=0.999
       
       
      IF (STATEV(4) .LT. 0.0) THEN               !not getting
		DYDX(3) = 0.0
       END IF
      

!**-------1---------2---------3---------4---------5---------6---------7**  
!**  Increment-Euler
!**-------1---------2---------3---------4---------5---------6---------7**
       DECRA(1)=DYDX(1)*DTIME
       IF(LEXIMP.EQ.1) THEN
       T13=((QTILD-STATEV(2)-SMETAL(1))/SMETAL(2))**(SMETAL(3)-1)
        DECRA(5)=(1/SMETAL(2))*SMETAL(3)*(1/(1-STATEV(4)))*T13*DTIME
       END IF
      
      
      STATEV(1)=STATEV(1)+DYDX(1)*DTIME
      STATEV(3)=STATEV(3)+DYDX(2)*DTIME
      STATEV(4)=STATEV(4)+DYDX(3)*DTIME
      
      
      
       RETURN
       END
