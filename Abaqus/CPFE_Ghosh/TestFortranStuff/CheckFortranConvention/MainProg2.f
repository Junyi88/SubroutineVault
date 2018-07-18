      program MainProg2
          implicit none

      Integer:: FourthT(3,3,3,3), X81(81)
      Real*8:: A(3,3),B(3,3),x(3),y(3),C(3,3),B9(9), C2(3,3),F(3,3)
      integer:: I,J,K,L,Z
      A(1,:)=(/1.,2.,3./)
      A(2,:)=(/4.,5.,6./)
      A(3,:)=(/7.,8.,9./)

      B9=(/10.0,13.0,16.0,11.0,14.0,17.0,12.0,15.0,18.0/)
      B=reshape(B9,(/3, 3/))

      x=B(1,:)
      C2=0.
      do I=1,3
      do J=1,3
      do k=1,3
         C2(I,J)=C2(I,J)+A(I,K)*B(K,J)
      end do
      end do
      end do

      y=matmul(A,x)
      C=matmul(A,B)

      do I=1,3
      do J=1,3
         F(I,J)=I*10.0+J
      end do
      end do

      print *,  "F--------------"
      print *,  F(1,1), F(1,2), F(1,3)
      print *,  F(2,1), F(2,2), F(2,3)
      print *,  F(3,1), F(3,2), F(3,3)

      print *,  "A--------------"
      print *,  A(1,1), A(1,2), A(1,3)
      print *,  A(2,1), A(2,2), A(2,3)
      print *,  A(3,1), A(3,2), A(3,3)
      print *,  "  "
      print *,  "B--------------"
      print *,  B(1,1), B(1,2), B(1,3)
      print *,  B(2,1), B(2,2), B(2,3)
      print *,  B(3,1), B(3,2), B(3,3)
      print *,  "  "

      print *,  "y_EXP,y_GET--------------"
      print *,  "68, 167, 266"
      print *,  y(1), y(2), y(3)
      print *,  "  "

      print *,  "C_EXP--------------"
      print *,  "84, 90, 96"
      print *,  "201, 216, 231"
      print *,  "318, 342, 366"
      print *,  "C_Get--------------"
      print *,  C(1,1), C(1,2), C(1,3)
      print *,  C(2,1), C(2,2), C(2,3)
      print *,  C(3,1), C(3,2), C(3,3)
      print *,  "  "
      print *,  "C2--------------"
      print *,  C2(1,1), C2(1,2), C2(1,3)
      print *,  C2(2,1), C2(2,2), C2(2,3)
      print *,  C2(3,1), C2(3,2), C2(3,3)
      print *,  "  "

c ===================================================


      do I=1,81
         X81(I)=I
      end do


      FourthT=reshape(X81,(/3, 3,3 ,3/))
      print *,  "FourthT--------------"
      do I=1,3
      do J=1,3
      do k=1,3
      do l=1,3
         Z = I + (J-1)*3 + (k-1)*9 +(L-1)*27
         print *,  I,J,K,L, FourthT(I,J,K,L),Z
      end do
      end do
      end do
      end do

      end program MainProg2

      include 'utils.f'
      include 'ksvd2.f'
