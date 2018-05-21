!Slip Rule Re-development
      subroutine kslip6(xNorm,xDir,tau,tauc,burgerv,
     + caratio, dtime,nSys,r,iphase,Lp,tmat,TEMP,gammaDot)
         
      implicit none
      integer,intent(in):: nSys,iphase
      real*8,intent(in) :: TEMP
      real*8,intent(in) :: dtime,r,caratio
      real*8,intent(in) :: xNorm(nSys,3),xDir(nSys,3),tau(nSys),tauc(nSys),
     +    burgerv(nSys)
      real*8,intent(out) :: Lp(3,3),tmat(6,6)
      
      integer :: i
      real*8  :: rhoInitial,xalpha,xbeta,result1,rhoc,
     + xlambdap,xvol,rhossdm,psi,xboltzman,xtemp,xhelmholtz,xfreq,
     + xsnt(3,3),xsnv(6),xnsv(6),xsnnst(6,6),xnst(3,3),result4(6,6),
     + tempNorm(3), tempDir(3),gammaDot(nSys), gammazero
      

C **************************************************************** 
C     CALCULATE THE DERIVATIVE OF PLASTIC STRAIN INCREMENT WITH 
C     RESPECT TO THE STRESS DEFINED AS tmat
C ****************************************************************
                      
      ! Changed 03/03/16 - DJW - Zr Gong Britton, Cuddihy.... Acta 2015
      rhossdm = 5.0                                   
      xhelmholtz = 9.0058345E-20                               
      xfreq = 1E11 ! This is in s^-1 TIME IS IN SECONDS   
      gammazero = 3.53e-4                             
      xboltzman = 1.381E-23 !CONSTANT
      xtemp = 293.0
      rhoInitial = 0.01 !DECOUPLE RHOSSD
      
      
      !Later we can use iphase to prescribe different xhelmholtz, xfreq, etc for the different phase within here!
      !To permanently avoid issues with units in the sinh, let its argument be completely in SI-units!
      !This leads to the required factor of 1.0e-12
C
      tmat=0.; Lp = 0.;result4=0.
	           
      DO I=1,3           
         xlambdap = 1.0/sqrt((rhoInitial)) !overall       
         xvol = xlambdap*burgerv(I)*burgerv(I)*1.0E-12
         xbeta = (xvol*gammazero*0.102)/(xboltzman*xtemp) 
         xalpha = rhossdm*burgerv(I)*burgerv(I)*xfreq*
     +  exp(-xhelmholtz*0.75556/(xboltzman*xtemp))
C
         IF (tau(I) >= tauc(I)) THEN
         
         gammaDot(I)=xalpha*sinh(xbeta*abs(tau(I)-tauc(I)))
         tempNorm = xNorm(I,:); tempDir = xDir(I,:)
         xsnt = spread(tempDir,2,3)*spread(tempNorm,1,3)
         xnst = spread(tempNorm,2,3)*spread(tempDir,1,3)
         CALL KGMATVEC6(xsnt,xsnv)         
         CALL KGMATVEC6(xnst,xnsv) 
         xsnnst = spread(xsnv,2,6)*spread(xnsv,1,6)
         result1 = cosh(xbeta*abs(tau(I)-tauc(I)))
          
         result4 = result4 + xalpha*xbeta*dtime*result1*xsnnst          
         Lp = Lp + gammaDot(I)*xsnt
        
         ELSE
            gammaDot(I)=0.0
         END IF
C
      END DO
	  
	  DO I=4,nSys            
         xlambdap = 1.0/sqrt((rhoInitial)) !overall       
         xvol = xlambdap*burgerv(I)*burgerv(I)*1.0E-12
         xbeta = (xvol*gammazero)/(xboltzman*xtemp) 
         xalpha = rhossdm*burgerv(I)*burgerv(I)*xfreq*
     +  exp(-xhelmholtz/(xboltzman*xtemp))
C
         IF (tau(I) >= tauc(I)) THEN
         
         gammaDot(I)=xalpha*sinh(xbeta*abs(tau(I)-tauc(I)))
         tempNorm = xNorm(I,:); tempDir = xDir(I,:)
         xsnt = spread(tempDir,2,3)*spread(tempNorm,1,3)
         xnst = spread(tempNorm,2,3)*spread(tempDir,1,3)
         CALL KGMATVEC6(xsnt,xsnv)         
         CALL KGMATVEC6(xnst,xnsv) 
         xsnnst = spread(xsnv,2,6)*spread(xnsv,1,6)
         result1 = cosh(xbeta*abs(tau(I)-tauc(I)))
          
         result4 = result4 + xalpha*xbeta*dtime*result1*xsnnst          
         Lp = Lp + gammaDot(I)*xsnt
        
         ELSE
            gammaDot(I)=0.0
         END IF
C
      END DO
      
      tmat = 0.5*(result4+transpose(result4)) !diffenertial of the
											  ! strain tensor increment
             
      RETURN
      END SUBROUTINE kslip6
