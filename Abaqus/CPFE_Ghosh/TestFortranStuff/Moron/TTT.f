      program TTT
          implicit none
          real*8  :: xNorm(3),xDir(3)
          real*8  :: xsnt(3,3),xsnv(6),xnsv(6),xsnnst(6,6),xnst(3,3)
          integer, parameter :: Moron=5
          integer, parameter :: Moron2 = Moron*2

          real*8 :: A(3,6), SG(3)

          integer :: I, J

          DO I = 1,3
           DO J = 1,6
              A(I,J)=I*10.0+J
           END DO
          END DO

          SG(1)=1.0
          SG(2)=-1.0
          SG(3)= 1.0

          print *, 'A ='
          print *, A(1,:)
          print *, A(2,:)
          print *, A(3,:)
          print *, 'A2=='
          print *, SG(1)*A(1,:)
          print *, SG(2)*A(2,:)
          print *, SG(3)*A(3,:)


          xNorm =  (/ 1, 2, 3 /)
          xDir =  (/ 4, 5, 6 /)

          print *, 'Spread1 = ', spread(xNorm,1,3)
          print *, 'Spread2 = ', spread(xDir,2,3)
          xsnt = spread(xDir,2,3)*spread(xNorm,1,3)
          print *, xsnt

          print *, 'Spread3 = ', spread(xNorm,2,3)
          print *, 'Spread4 = ', spread(xDir,1,3)
          xnst = spread(xNorm,2,3)*spread(xDir,1,3)
          print *, xnst



      end program TTT
