      subroutine kdirns(xrot,iphase,L,ddir1,dnor1)
      
      implicit none      
      integer :: i,j,k
      integer,parameter :: m = 3
      integer,intent(in):: iphase,L
      
      real*8,intent(in):: xrot(m,m)
      real*8,intent(out):: ddir1(L,m),dnor1(L,m)
      
      real(kind=8):: ddir(L,m),dnor(L,m)
      real(kind=8):: tdir(m),tnor(m),tdir1(m),tnor1(m)
     
      real(kind=8):: xdirmag,xnormag

         ddir1 = 0.0; dnor1 = 0.0 
      
      if(iphase .eq. 0) then !hcp
         include 'xDir0.f'
         include 'xNorm0.f'
      else if(iphase .eq. 1) then !bcc (24/48)
         include 'xDir1.f'
         include 'xNorm1.f'
      else if(iphase .eq. 2) then !fcc
         include 'xDir2.f'
         include 'xNorm2.f'      
      else !fcc and cubic
         include 'xDir2.f'
         include 'xNorm2.f' 
      end if     

       do k=1,L
        do i=1,m
        tdir(i) = ddir(k,i)
        tnor(i) = dnor(k,i)
        end do
        
        tdir1 = matmul(xrot,tdir)
        tnor1 = matmul(xrot,tnor)
        
        xdirmag = sqrt(dot_product(tdir1,tdir1))
        xnormag = sqrt(dot_product(tnor1,tnor1))
        
        do i=1,m
        ddir1(k,i) = tdir1(i)/xdirmag
        dnor1(k,i) = tnor1(i)/xnormag
        end do
      end do
      
      return
      end 

