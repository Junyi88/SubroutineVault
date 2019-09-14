      program MainProg
          implicit none

      Real*8:: My33Mat(3,3), My9Mat(9),
     1 My9Reshape33(3,3),My33Reshape9(9),
     1 A(3,3),AT(3,3),x(3),y(3)

      My33Mat(1,1)=11
      My33Mat(1,2)=12
      My33Mat(1,3)=13
      My33Mat(2,1)=21
      My33Mat(2,2)=22
      My33Mat(2,3)=23
      My33Mat(3,1)=31
      My33Mat(3,2)=32
      My33Mat(3,3)=33

      My9Mat(1)=1.0
      My9Mat(2)=2.0
      My9Mat(3)=3.0
      My9Mat(4)=4.0
      My9Mat(5)=5.0
      My9Mat(6)=6.0
      My9Mat(7)=7.0
      My9Mat(8)=8.0
      My9Mat(9)=9.0



      print *,  "My33Mat---------------"
      print *,  "(1,1),(1,2),(1,3)"
      print *,  "(2,1),(2,2),(2,3)"
      print *,  "(3,1),(3,2),(3,3)"
      print *,  My33Mat(1,1),My33Mat(1,2),My33Mat(1,3)
      print *,  My33Mat(2,1),My33Mat(2,2),My33Mat(2,3)
      print *,  My33Mat(3,1),My33Mat(3,2),My33Mat(3,3)
      print *,  "   "
      print *,  My33Mat
      print *,  "   "

      My9Reshape33= reshape(My9Mat,(/3, 3/))
      My33Reshape9= reshape(My33Mat,(/9/))
      A=My9Reshape33
      x(1)=10.0
      x(2)=11.0
      x(3)=12.0
      AT=transpose(A)
      y=matmul(A,x)



      print *,  "My9Mat---------------"
      print *,  My9Mat(1), My9Mat(2), My9Mat(3)
      print *,  My9Mat(4), My9Mat(5), My9Mat(6)
      print *,  My9Mat(7), My9Mat(8), My9Mat(9)
      print *,  "   "

      print *,  "My33Reshape9---------------"
      print *,  My33Reshape9(1), My33Reshape9(2), My33Reshape9(3)
      print *,  My33Reshape9(4), My33Reshape9(5), My33Reshape9(6)
      print *,  My33Reshape9(7), My33Reshape9(8), My33Reshape9(9)
      print *,  "   "

      print *,  "My9Reshape33---------------"
      print *,  "(1,1),(1,2),(1,3)"
      print *,  "(2,1),(2,2),(2,3)"
      print *,  "(3,1),(3,2),(3,3)"
      print *,  My9Reshape33(1,1),My9Reshape33(1,2),My9Reshape33(1,3)
      print *,  My9Reshape33(2,1),My9Reshape33(2,2),My9Reshape33(2,3)
      print *,  My9Reshape33(3,1),My9Reshape33(3,2),My9Reshape33(3,3)
      print *,  "   "
      print *,  My9Reshape33
      print *,  "   "

      print *,  "y---------------"
      print *,  y(1),y(2),y(3)


      end program MainProg

      include 'utils.f'
      include 'ksvd2.f'
