!Least Squares Density Minimisation
      subroutine kgndl2(curlfe,xNin,xDin,tau2,ids,burgerv,iphase,L,
     + jelem,kint,time,gndnew,capyramedge,capyramscrew,apyramedge,
     + aprismedge,abasedge,ascrew)

!      implicit real*8(a-h,o-z)
      implicit none

      integer :: i,m,n,iec,isc,ne,ns,nsmax,ialloc,idir,iloc,k,info
      real*8  :: det,tausum,gndesum,gndssum
      real*8,parameter :: zero = 1.0e-6
      
      integer,intent(in):: L,iphase,jelem,kint
      real*8,intent(in) :: time(2)
      real*8,intent(in) :: curlfe(3,3),xNin(L,3),xDin(L,3),tau2(L),
     + burgerv(L)
      integer,intent(in)  :: ids(L)
      
      real*8,dimension(3,3) :: xrot,btdyad,bsdyad
      real*8,dimension(L)   :: gndnew,gnde,gnds,absgnde,absgnds
      real*8,dimension(3)   :: tempn,temps,tempt
      real*8,dimension(9)   :: btvec,bsvec,gv
      
      real*8,dimension(:,:),allocatable :: A,A1,A2,Ainv
      real*8,dimension(:),allocatable :: xvec
      integer,dimension(:,:),allocatable :: iscrews
      
      character (len=*), parameter :: fmt2 = "(24(' ',(I2,1X)))",
     + fmt3="(3(' ',(ES11.3,1X)))",
     + fmt9 = "(9(' ',(ES11.3,1X)))",fmt33 = "(33(' ',(ES11.3,1X)))"
     
      !Ben's experimental data storage groups
      real*8,intent(out):: capyramedge,capyramscrew,apyramedge,
     + aprismedge,abasedge,ascrew
      
     
      EXTERNAL DGELSD

      select case(iphase)
      case (0); nsmax = 9 !HCP       
      case (1); nsmax = 4 !BCC     
      case (2); nsmax = 6 !FCC      
      end select 
      
      gnde = 0.; gnds = 0.; gndnew = 0.0
      
      !Ben's experimental data storage groups
      capyramedge=0.;capyramscrew=0.;apyramedge=0.;aprismedge=0.
      abasedge=0.; ascrew=0.     
      
!      write(6,*) "time,jelem,kint =",time(1),jelem,kint
            
      !Compute the geometric dislocation tensor
      gv = reshape(curlfe,(/9/))!This happens by column.
      
!      write(6,*) "curlfe"; write(6,fmt3) (curlfe(i,:),i=1,ubound(curlfe,1))
           
      !===BEGIN LONG IF
      if(maxval(abs(gv)) <= zero) then 
      !Trying to handle the situation when G = 0.
         gnde = 0.; gnds = 0.
      else         
      
      ne = count(ids /= 0) !Number of edge      
      !Number of screw 
      ns=0      
      select case(iphase)
      case (0)
      !a-slip
      if(ids(1) /= 0 .or. ids(4) /= 0 .or. ids(7) /=0 .or. ids(10) /=0) 
     + ns=ns+1
      if(ids(2) /= 0 .or. ids(5) /= 0 .or. ids(8) /=0 .or. ids(11) /=0)
     +  ns=ns+1
      if(ids(3) /= 0 .or. ids(6) /= 0 .or. ids(9) /=0 .or. ids(12) /=0)
     +  ns=ns+1
      !c+a-slip
      if(ids(13) /= 0 .or. ids(19) /= 0) ns=ns+1
      if(ids(14) /= 0 .or. ids(20) /= 0) ns=ns+1
      if(ids(15) /= 0 .or. ids(21) /= 0) ns=ns+1
      if(ids(16) /= 0 .or. ids(22) /= 0) ns=ns+1
      if(ids(17) /= 0 .or. ids(23) /= 0) ns=ns+1
      if(ids(18) /= 0 .or. ids(24) /= 0) ns=ns+1
      
      case(1)
      if(ids(1)/=0.or.ids(8)/=0.or.ids(11)/=0.or.ids(15)/=0.or.
     + ids(20)/=0.or.ids(21)/=0.or.ids(27)/=0.or.ids(32)/=0.or.
     + ids(36)/=0.or.ids(39)/=0.or.ids(41)/=0.or.ids(45)/=0) ns=ns+1
     
      if(ids(2)/=0.or.ids(4)/=0.or.ids(9)/=0.or.ids(14)/=0.or.
     + ids(17)/=0.or.ids(24)/=0.or.ids(26)/=0.or.ids(29)/=0.or.
     + ids(33)/=0.or.ids(38)/=0.or.ids(44)/=0.or.ids(48)/=0) ns=ns+1
     
      if(ids(3)/=0.or.ids(5)/=0.or.ids(7)/=0.or.ids(16)/=0.or.
     + ids(19)/=0.or.ids(22)/=0.or.ids(28)/=0.or.ids(31)/=0.or.
     + ids(35)/=0.or.ids(40)/=0.or.ids(42)/=0.or.ids(46)/=0) ns=ns+1
     
      if(ids(6)/=0.or.ids(10)/=0.or.ids(12)/=0.or.ids(13)/=0.or.
     + ids(18)/=0.or.ids(23)/=0.or.ids(25)/=0.or.ids(30)/=0.or.
     + ids(34)/=0.or.ids(37)/=0.or.ids(43)/=0.or.ids(47)/=0) ns=ns+1
      
      case(2)
      if(ids(1) /= 0 .or. ids(10) /=0) ns=ns+1
      if(ids(2) /= 0 .or. ids(5) /=0) ns=ns+1
      if(ids(3) /= 0 .or. ids(9) /=0) ns=ns+1
      if(ids(4) /= 0 .or. ids(7) /=0) ns=ns+1
      if(ids(6) /= 0 .or. ids(12) /=0) ns=ns+1
      if(ids(8) /= 0 .or. ids(11) /=0) ns=ns+1
      end select
      
      m = 9; n = ne+ns      
      allocate(A(m,n),Ainv(n,m),xvec(n),iscrews(nsmax,2),STAT=ialloc)
      A = 0.; Ainv = 0.; xvec=0. 
      
      if (m < n) then !right inverse
       allocate(A1(m,m),A2(m,m),STAT=ialloc)
       A1=0.; A2=0.
      else !if (m > n) then !left inverse and the usual inverse
       allocate(A1(n,n),A2(n,n),STAT=ialloc)
       A1=0.; A2=0.
      end if
      
      !Construct the matrix of dyadics
      !Edge
      iec = 0; iscrews = 0
      do i = 1,L
	      if(ids(i) /= 0)then
	      iec=iec+1
	      
	      if(iphase == 0) then
	      !Screw preparation for hcp.
	      select case(i)
	      case(1,4,7,10)
	         if(iscrews(1,1) /= 0) then; iscrews(1,2) = iscrews(1,2)+1
	         else; iscrews(1,1) = i; iscrews(1,2) = 1
	         end if
	      case(2,5,8,11)
	         if(iscrews(2,1) /= 0) then; iscrews(2,2) = iscrews(2,2)+1
	         else; iscrews(2,1) = i; iscrews(2,2) = 1
	         end if
	      case(3,6,9,12)
	         if(iscrews(3,1) /= 0) then; iscrews(3,2) = iscrews(3,2)+1
	         else; iscrews(3,1) = i; iscrews(3,2) = 1
	         end if
	      case(13,19)
	         if(iscrews(4,1) /= 0) then; iscrews(4,2) = iscrews(4,2)+1
	         else; iscrews(4,1) = i; iscrews(4,2) = 1
	         end if
	      case(14,20)
	         if(iscrews(5,1) /= 0) then; iscrews(5,2) = iscrews(5,2)+1
	         else; iscrews(5,1) = i; iscrews(5,2) = 1
	         end if
	      case(15,21)
	         if(iscrews(6,1) /= 0) then; iscrews(6,2) = iscrews(6,2)+1
	         else; iscrews(6,1) = i; iscrews(6,2) = 1
	         end if
	      case(16,22)
	         if(iscrews(7,1) /= 0) then; iscrews(7,2) = iscrews(7,2)+1
	         else; iscrews(7,1) = i; iscrews(7,2) = 1
	         end if
	      case(17,23)
	         if(iscrews(8,1) /= 0) then; iscrews(8,2) = iscrews(8,2)+1
	         else; iscrews(8,1) = i; iscrews(8,2) = 1
	         end if
	      case(18,24)
	         if(iscrews(9,1) /= 0) then; iscrews(9,2) = iscrews(9,2)+1
	         else; iscrews(9,1) = i; iscrews(9,2) = 1
	         end if
	      end select 
	        
	      else if(iphase == 1) then
	      select case(i)
	      case(1,8,11,15,20,21,27,32,36,39,41,45)
	         if(iscrews(1,1) /= 0) then; iscrews(1,2) = iscrews(1,2)+1
	         else; iscrews(1,1) = i; iscrews(1,2) = 1
	         end if
	      case(2,4,9,14,17,24,26,29,33,38,44,48)
	         if(iscrews(2,1) /= 0) then; iscrews(2,2) = iscrews(2,2)+1
	         else; iscrews(2,1) = i; iscrews(2,2) = 1
	         end if
	      case(3,5,7,16,19,22,28,31,35,40,42,46)
	         if(iscrews(3,1) /= 0) then; iscrews(3,2) = iscrews(3,2)+1
	         else; iscrews(3,1) = i; iscrews(3,2) = 1
	         end if
	      case(6,10,12,13,18,23,25,30,34,37,43,47)
	         if(iscrews(4,1) /= 0) then; iscrews(4,2) = iscrews(4,2)+1
	         else; iscrews(4,1) = i; iscrews(4,2) = 1
	         end if
	      end select
	       
	      else !phase 2
	      select case(i)
	      case(1,10)
	         if(iscrews(1,1) /= 0) then; iscrews(1,2) = iscrews(1,2)+1
	         else; iscrews(1,1) = i; iscrews(1,2) = 1
	         end if
	      case(2,5)
	         if(iscrews(2,1) /= 0) then; iscrews(2,2) = iscrews(2,2)+1
	         else; iscrews(2,1) = i; iscrews(2,2) = 1
	         end if
	      case(3,9)
	         if(iscrews(3,1) /= 0) then; iscrews(3,2) = iscrews(3,2)+1
	         else; iscrews(3,1) = i; iscrews(3,2) = 1
	         end if
	      case(4,7)
	         if(iscrews(4,1) /= 0) then; iscrews(4,2) = iscrews(4,2)+1
	         else; iscrews(4,1) = i; iscrews(4,2) = 1
	         end if
	      case(6,12)
	         if(iscrews(5,1) /= 0) then; iscrews(5,2) = iscrews(5,2)+1
	         else; iscrews(5,1) = i; iscrews(5,2) = 1
	         end if
	      case(8,11)
	         if(iscrews(6,1) /= 0) then; iscrews(6,2) = iscrews(6,2)+1
	         else; iscrews(6,1) = i; iscrews(6,2) = 1
	         end if
	      end select
	       
	      end if
         

         tempn = xNin(i,:)
         temps = xDin(i,:)
         CALL KVECPROD(temps,tempn,tempt)

         btdyad = spread(tempt,2,3)*spread(temps,1,3)*burgerv(i)
         A(:,iec) = reshape(btdyad,(/9/))     
         end if         
      end do 
      
      !Screw
      isc = 0
      do i = 1,nsmax
	      if(iscrews(i,1) /= 0)then
	      isc=isc+1; idir = iscrews(i,1) 
         temps = xDin(idir,:)
         bsdyad = spread(temps,2,3)*spread(temps,1,3)*burgerv(idir)
         A(:,ne+isc) = reshape(bsdyad,(/9/))    
         end if         
      end do       
     
!      if ((jelem == 7910 .or. jelem == 8010) .and. kint == 6) then
!         write(6,*) "A, ne",ubound(A,1),ubound(A,2),ne
!         write(6,fmt33) (A(k,:),k=1,ubound(A,1))  !Print rows
!         write(6,*)  
!      end if 

!!      call ksvd(A,m,n,Ainv)
!      
      !three-ifs direct inverses
      if (m < n) then !right inverse
       A1 = matmul(A,transpose(A))
       call lapinverse(A1,m,info,A2)
!       if(info /= 0) write(6,*) "inverse failure: A1 in kgndl2"
       Ainv = matmul(transpose(A),A2)
      else !left inverse
       A1 = matmul(transpose(A),A)
       call lapinverse(A1,n,info,A2)
!       if(info /= 0) write(6,*) "inverse failure: A1 in kgndl2"
       Ainv = matmul(A2,transpose(A))
      end if
      
      !three-ifs solution
      xvec = matmul(Ainv,gv)
      
      do i=1,ubound(xvec,1)
         if (abs(xvec(i)) <= 1.0e-6 .or. xvec(i) /= xvec(i)) 
     +     xvec(i) = 0.0
      end do
      
!      if ((jelem == 7910 .or. jelem == 8010) .and. kint == 6) then
!            write(6,*) "xvec"; write(6,fmt33) xvec
!      end if
      
      !Read back dislocations
      !Edge
      iec = 0
      do i=1,L
       if(ids(i) /= 0) then; iec = iec+1; gnde(i) = xvec(iec); end if !three-ifs
      end do
            
      !Screw
      !For three-ifs solution, replace B(iloc,1) with xvec(iloc)
      isc = 0
      do i = 1,nsmax
	      if(iscrews(i,1) /= 0)then
	      isc = isc+1; iloc = ne+isc
	      if(iphase == 0) then
	      select case(iscrews(i,1))
	      case(1,4,7,10)
	         tausum = tau2(1)+tau2(4)+tau2(7)+tau2(10)
	         gnds(1) = xvec(iloc)*tau2(1)/tausum
	         gnds(4) = xvec(iloc)*tau2(4)/tausum
	         gnds(7) = xvec(iloc)*tau2(7)/tausum
	         gnds(10) = xvec(iloc)*tau2(10)/tausum
	      case(2,5,8,11)
	         tausum = tau2(2)+tau2(5)+tau2(8)+tau2(11)
	         gnds(2) = xvec(iloc)*tau2(2)/tausum
	         gnds(5) = xvec(iloc)*tau2(5)/tausum
	         gnds(8) = xvec(iloc)*tau2(8)/tausum
	         gnds(11) = xvec(iloc)*tau2(11)/tausum
	      case(3,6,9,12)
	         tausum = tau2(3)+tau2(6)+tau2(9)+tau2(12)
	         gnds(3) = xvec(iloc)*tau2(3)/tausum
	         gnds(6) = xvec(iloc)*tau2(6)/tausum
	         gnds(9) = xvec(iloc)*tau2(9)/tausum
	         gnds(12) = xvec(iloc)*tau2(12)/tausum
	      case(13,19)
	         tausum = tau2(13)+tau2(19)
	         gnds(13) = xvec(iloc)*tau2(13)/tausum
	         gnds(19) = xvec(iloc)*tau2(19)/tausum
!	         write(6,*) "tausum 13,19",tausum,gnds(13),gnds(19)
	      case(14,20)
	         tausum = tau2(14)+tau2(20)
	         gnds(14) = xvec(iloc)*tau2(14)/tausum
	         gnds(20) = xvec(iloc)*tau2(20)/tausum
	      case(15,21)
	         tausum = tau2(15)+tau2(21)
	         gnds(15) = xvec(iloc)*tau2(15)/tausum
	         gnds(21) = xvec(iloc)*tau2(21)/tausum
	      case(16,22)
	         tausum = tau2(16)+tau2(22)
	         gnds(16) = xvec(iloc)*tau2(16)/tausum
	         gnds(22) = xvec(iloc)*tau2(22)/tausum
	      case(17,23)
	         tausum = tau2(17)+tau2(23)
	         gnds(17) = xvec(iloc)*tau2(17)/tausum
	         gnds(23) = xvec(iloc)*tau2(23)/tausum
!	         write(6,*) "tausum 17,23",tausum,gnds(17),gnds(23)
	      case(18,24)
	         tausum = tau2(18)+tau2(24)
	         gnds(18) = xvec(iloc)*tau2(18)/tausum
	         gnds(24) = xvec(iloc)*tau2(24)/tausum
	      end select
	      
	      else if(iphase == 1) then
	      select case(i)
	      case(1,8,11,15,20,21,27,32,36,39,41,45)
	         tausum = tau2(1)+tau2(8)+tau2(11)+tau2(15)+tau2(20)+tau2(21)
!     +      +tau2(27)+tau2(32)+tau2(36)+tau2(39)+tau2(41)+tau2(45) !comment for 24
	         gnds(1) = xvec(iloc)*tau2(1)/tausum
	         gnds(8) = xvec(iloc)*tau2(8)/tausum
	         gnds(11) = xvec(iloc)*tau2(11)/tausum
	         gnds(15) = xvec(iloc)*tau2(15)/tausum
	         gnds(20) = xvec(iloc)*tau2(20)/tausum
	         gnds(21) = xvec(iloc)*tau2(21)/tausum
!	         gnds(27) = xvec(iloc)*tau2(27)/tausum !comment from here for 24
!	         gnds(32) = xvec(iloc)*tau2(32)/tausum
!	         gnds(36) = xvec(iloc)*tau2(36)/tausum
!	         gnds(39) = xvec(iloc)*tau2(39)/tausum
!	         gnds(41) = xvec(iloc)*tau2(41)/tausum
!	         gnds(45) = xvec(iloc)*tau2(45)/tausum
	      case(2,4,9,14,17,24,26,29,33,38,44,48)
	         tausum = tau2(2)+tau2(4)+tau2(9)+tau2(14)+tau2(17)+tau2(24)
!     +      +tau2(26)+tau2(29)+tau2(33)+tau2(38)+tau2(44)+tau2(48) !comment for 24
	         gnds(2) = xvec(iloc)*tau2(2)/tausum
	         gnds(4) = xvec(iloc)*tau2(4)/tausum
	         gnds(9) = xvec(iloc)*tau2(9)/tausum
	         gnds(14) = xvec(iloc)*tau2(14)/tausum
	         gnds(17) = xvec(iloc)*tau2(17)/tausum
	         gnds(24) = xvec(iloc)*tau2(24)/tausum
!	         gnds(26) = xvec(iloc)*tau2(26)/tausum !comment from here for 24
!	         gnds(29) = xvec(iloc)*tau2(29)/tausum
!	         gnds(33) = xvec(iloc)*tau2(33)/tausum
!	         gnds(38) = xvec(iloc)*tau2(38)/tausum
!	         gnds(44) = xvec(iloc)*tau2(44)/tausum
!	         gnds(48) = xvec(iloc)*tau2(48)/tausum
	      case(3,5,7,16,19,22,28,31,35,40,42,46)
	         tausum = tau2(3)+tau2(5)+tau2(7)+tau2(16)+tau2(19)+tau2(22)
!     +      +tau2(28)+tau2(31)+tau2(35)+tau2(40)+tau2(42)+tau2(46) !comment for 24
	         gnds(3) = xvec(iloc)*tau2(3)/tausum
	         gnds(5) = xvec(iloc)*tau2(5)/tausum
	         gnds(7) = xvec(iloc)*tau2(7)/tausum
	         gnds(16) = xvec(iloc)*tau2(16)/tausum
	         gnds(19) = xvec(iloc)*tau2(19)/tausum
	         gnds(22) = xvec(iloc)*tau2(22)/tausum
!	         gnds(28) = xvec(iloc)*tau2(28)/tausum !comment from here for 24
!	         gnds(31) = xvec(iloc)*tau2(31)/tausum
!	         gnds(35) = xvec(iloc)*tau2(35)/tausum
!	         gnds(40) = xvec(iloc)*tau2(40)/tausum
!	         gnds(42) = xvec(iloc)*tau2(42)/tausum
!	         gnds(46) = xvec(iloc)*tau2(46)/tausum
	      case(6,10,12,13,18,23,25,30,34,37,43,47)
	         tausum = tau2(6)+tau2(10)+tau2(12)+tau2(13)+tau2(18)+tau2(23)
!     +     +tau2(25)+tau2(30)+tau2(34)+tau2(37)+tau2(43)+tau2(47) !comment for 24
	         gnds(6) = xvec(iloc)*tau2(6)/tausum
	         gnds(10) = xvec(iloc)*tau2(10)/tausum
	         gnds(12) = xvec(iloc)*tau2(12)/tausum
	         gnds(13) = xvec(iloc)*tau2(13)/tausum
	         gnds(18) = xvec(iloc)*tau2(18)/tausum
	         gnds(23) = xvec(iloc)*tau2(23)/tausum
!	         gnds(25) = xvec(iloc)*tau2(25)/tausum !comment from here for 24
!	         gnds(30) = xvec(iloc)*tau2(30)/tausum
!	         gnds(34) = xvec(iloc)*tau2(34)/tausum
!	         gnds(37) = xvec(iloc)*tau2(37)/tausum
!	         gnds(43) = xvec(iloc)*tau2(43)/tausum
!	         gnds(47) = xvec(iloc)*tau2(47)/tausum
	      end select
	      
	      else
	      select case(iscrews(i,1))
	      case(1,10)
	         tausum = tau2(1)+tau2(10)
	         gnds(1) = xvec(iloc)*tau2(1)/tausum
	         gnds(10) = xvec(iloc)*tau2(10)/tausum
	      case(2,5)
	         tausum = tau2(2)+tau2(5)
	         gnds(2) = xvec(iloc)*tau2(2)/tausum
	         gnds(5) = xvec(iloc)*tau2(5)/tausum
	      case(3,9)
	         tausum = tau2(3)+tau2(9)
	         gnds(3) = xvec(iloc)*tau2(3)/tausum
	         gnds(9) = xvec(iloc)*tau2(9)/tausum
	      case(4,7)
	         tausum = tau2(4)+tau2(7)
	         gnds(4) = xvec(iloc)*tau2(4)/tausum
	         gnds(7) = xvec(iloc)*tau2(7)/tausum
	      case(6,12)
	         tausum = tau2(6)+tau2(12)
	         gnds(6) = xvec(iloc)*tau2(6)/tausum
	         gnds(12) = xvec(iloc)*tau2(12)/tausum
	      case(8,11)
	         tausum = tau2(8)+tau2(11)
	         gnds(8) = xvec(iloc)*tau2(8)/tausum
	         gnds(11) = xvec(iloc)*tau2(11)/tausum
	      end select


	      end if
	               
         end if         
      end do
      
      do i=1,L; if(abs(gnds(i)) <= zero .or. gnds(i) /= gnds(i)) 
     + gnds(i) = 0.0; end do
      
      deallocate(A,iscrews,A1,A2,Ainv,xvec)     
      !===END LONG IF
      end if
      
      !Density to output
      gndnew = abs(sqrt(gnde*gnde+gnds*gnds))       
      do i=1,L; if(gndnew(i)/=gndnew(i)) gndnew(i)=0.0; end do  !catching NaN
      where(gndnew > 1.0e7) gndnew = 0.0 !catching infinities. 1.0e7 for microns, 1.0e19 for metres.
      
!      if ((jelem == 7910 .or. jelem == 8010) .and. kint == 6) then
!            write(6,*) "gndnew"; write(6,fmt33) gndnew      
!      end if   
      
      !Ben's experimental data storage groups
      absgnde = abs(gnde); absgnds = abs(gnds)
      if(iphase == 0) then
      !edge
      abasedge = sum(absgnde(1:3))
      aprismedge = sum(absgnde(4:6))
      apyramedge = sum(absgnde(7:12))
      capyramedge = sum(absgnde(13:24))
      !screw
      ascrew = sum(absgnds(1:12))
      capyramscrew = sum(absgnds(13:24))
      end if
       
      return
      end
