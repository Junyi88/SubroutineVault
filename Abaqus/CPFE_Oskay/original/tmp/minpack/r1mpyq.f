      subroutine r1mpyq(m,n,a,lda,v,w)
      integer m,n,lda
      double precision a(lda,n),v(n),w(n)
!     **********

!     subroutine r1mpyq

!     given an m by n matrix a, this subroutine computes a*q where
!     q is the product of 2*(n - 1) transformations

!           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)

!     and gv(i), gw(i) are givens rotations in the (i,n) plane which
!     eliminate elements in the i-th and n-th planes, respectively.
!     q itself is not given, rather the information to recover the
!     gv, gw rotations is supplied.

!     the subroutine statement is

!       subroutine r1mpyq(m,n,a,lda,v,w)

!     where

!       m is a positive integer input variable set to the number
!         of rows of a.

!       n is a positive integer input variable set to the number
!         of columns of a.

!       a is an m by n array. on input a must contain the matrix
!         to be postmultiplied by the orthogonal matrix q
!         described above. on output a*q has replaced a.

!       lda is a positive integer input variable not less than m
!         which specifies the leading dimension of the array a.

!       v is an input array of length n. v(i) must contain the
!         information necessary to recover the givens rotation gv(i)
!         described above.

!       w is an input array of length n. w(i) must contain the
!         information necessary to recover the givens rotation gw(i)
!         described above.

!     subroutines called

!       fortran-supplied ... dabs,dsqrt

!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more

!     **********
      integer i,j,nmj,nm1
      double precision cos,one,sin,temp
      data one /1.0d0/

!     apply the first set of givens rotations to a.

      nm1 = n - 1
      if (nm1 .lt. 1) go to 50
      do 20 nmj = 1, nm1
         j = n - nmj
         if (dabs(v(j)) .gt. one) cos = one/v(j)
         if (dabs(v(j)) .gt. one) sin = dsqrt(one-cos**2)
         if (dabs(v(j)) .le. one) sin = v(j)
         if (dabs(v(j)) .le. one) cos = dsqrt(one-sin**2)
         do 10 i = 1, m
            temp = cos*a(i,j) - sin*a(i,n)
            a(i,n) = sin*a(i,j) + cos*a(i,n)
            a(i,j) = temp
   10       continue
   20    continue

!     apply the second set of givens rotations to a.

      do 40 j = 1, nm1
         if (dabs(w(j)) .gt. one) cos = one/w(j)
         if (dabs(w(j)) .gt. one) sin = dsqrt(one-cos**2)
         if (dabs(w(j)) .le. one) sin = w(j)
         if (dabs(w(j)) .le. one) cos = dsqrt(one-sin**2)
         do 30 i = 1, m
            temp = cos*a(i,j) + sin*a(i,n)
            a(i,n) = -sin*a(i,j) + cos*a(i,n)
            a(i,j) = temp
   30       continue
   40    continue
   50 continue
      return

!     last card of subroutine r1mpyq.

      end
