!Slip Rule Re-development
      subroutine kslip5(xNorm,xDir,tau,tauc,burgerv,rhossd,gndold,caratio,
     + dtime,L,r,iphase,xlp,tmat)
      implicit none
      integer,intent(in):: L,iphase
      real*8,intent(in) :: rhossd,dtime,r,caratio
      real*8,intent(in) :: xNorm(L,3),xDir(L,3),tau(L),tauc(L),burgerv(L),gndold(L)
      real*8,intent(out) :: xlp(3,3),tmat(6,6)
      
      integer :: i
      real*8  :: rhognd,rhogndold,xalpha,xbeta,result1,rhoc,ags,
     + xlambdap,xvol,rhossdm,psi,xboltzman,xtemp,xhelmholtz,xfreq,
     + xsnt(3,3),xsnv(6),xnsv(6),xsnnst(6,6),xnst(3,3),result4(6,6),
     + tempNorm(3), tempDir(3),gammaDot(L), gamm0
      
C
C  *** CALCULATE THE DERIVATIVE OF PLASTIC STRAIN INCREMENT WITH 
C   RESPECT TO THE STRESS DEFINED AS tmat***
C      
      
      psi = 1.457e-4       !82.8457 if placed above in the slip rule

      
      rhossdm = 5.0e-6
      xhelmholtz = 3.4559E-20 *1E3 ! N.mm = mJ
      xfreq = 1.0E+11
      
      xboltzman = 1.381E-23 *1E3 ! mJ / K
      xtemp = 293.0 !823.0  !for Zr! CRSS was chosen based on T.
      !Later we can use iphase to prescribe different xhelmholtz, xfreq, etc for the different phase within here!
      !L = microns, stress = MPa, F = microN, therefore E = pJ
      gamm0 = 8.33E-6 !8.33E-6 ! this is some reference strain which appears in some papers eg http://dx.doi.org/10.1016/j.ijsolstr.2015.02.023
C
      rhogndold=sum(gndold)
      Do I=1,L
	              
         xlambdap = 1.0/sqrt(psi*(rhogndold+rhossd)) !overall       

         xvol = xlambdap*burgerv(I)*burgerv(I)
         
         xbeta = gamm0*xvol/(xboltzman*xtemp) 

         xalpha = rhossdm*burgerv(I)*burgerv(I)*xfreq*
     +  exp(-xhelmholtz/(xboltzman*xtemp))
C
         if (tau(I) >= tauc(I)) THEN
          gammaDot(I)=xalpha*sinh(xbeta*abs(tau(I)-r-tauc(I)))
         
          tempNorm = xNorm(I,:); tempDir = xDir(I,:)
          xsnt = spread(tempDir,2,3)*spread(tempNorm,1,3)
          xnst = spread(tempNorm,2,3)*spread(tempDir,1,3)
          CALL KGMATVEC6(xsnt,xsnv)         
          CALL KGMATVEC6(xnst,xnsv) 
          xsnnst = spread(xsnv,2,6)*spread(xnsv,1,6)
          result1 = cosh(xbeta*abs(tau(I)-r-tauc(I)))
          
          result4 = result4 + xalpha*xbeta*dtime*result1*xsnnst          
          xlp = xlp + gammaDot(I)*xsnt
         else
            gammaDot(I)=0.0
         end if
C
      END DO
      
      tmat = 0.5*(result4+transpose(result4))
             
      return
      end subroutine kslip5
