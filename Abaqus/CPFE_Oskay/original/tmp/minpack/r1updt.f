      subroutine r1updt(m,n,s,ls,u,v,w,sing)
      integer m,n,ls
      logical sing
      double precision s(ls),u(m),v(n),w(m)
!     **********

!     subroutine r1updt

!     given an m by n lower trapezoidal matrix s, an m-vector u,
!     and an n-vector v, the problem is to determine an
!     orthogonal matrix q such that

!                   t
!           (s + u*v )*q

!     is again lower trapezoidal.

!     this subroutine determines q as the product of 2*(n - 1)
!     transformations

!           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)

!     where gv(i), gw(i) are givens rotations in the (i,n) plane
!     which eliminate elements in the i-th and n-th planes,
!     respectively. q itself is not accumulated, rather the
!     information to recover the gv, gw rotations is returned.

!     the subroutine statement is

!       subroutine r1updt(m,n,s,ls,u,v,w,sing)

!     where

!       m is a positive integer input variable set to the number
!         of rows of s.

!       n is a positive integer input variable set to the number
!         of columns of s. n must not exceed m.

!       s is an array of length ls. on input s must contain the lower
!         trapezoidal matrix s stored by columns. on output s contains
!         the lower trapezoidal matrix produced as described above.

!       ls is a positive integer input variable not less than
!         (n*(2*m-n+1))/2.

!       u is an input array of length m which must contain the
!         vector u.

!       v is an array of length n. on input v must contain the vector
!         v. on output v(i) contains the information necessary to
!         recover the givens rotation gv(i) described above.

!       w is an output array of length m. w(i) contains information
!         necessary to recover the givens rotation gw(i) described
!         above.

!       sing is a logical output variable. sing is set true if any
!         of the diagonal elements of the output s are zero. otherwise
!         sing is set false.

!     subprograms called

!       minpack-supplied ... dpmpar

!       fortran-supplied ... dabs,dsqrt

!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more,
!     john l. nazareth

!     **********
      integer i,j,jj,l,nmj,nm1
      double precision cos,cotan,giant,one,p5,p25,sin,tan,tau,temp,&
                      zero
      double precision dpmpar
      data one,p5,p25,zero /1.0d0,5.0d-1,2.5d-1,0.0d0/

!     giant is the largest magnitude.

      giant = dpmpar(3)

!     initialize the diagonal element pointer.

      jj = (n*(2*m - n + 1))/2 - (m - n)

!     move the nontrivial part of the last column of s into w.

      l = jj
      do 10 i = n, m
         w(i) = s(l)
         l = l + 1
   10    continue

!     rotate the vector v into a multiple of the n-th unit vector
!     in such a way that a spike is introduced into w.

      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 nmj = 1, nm1
         j = n - nmj
         jj = jj - (m - j + 1)
         w(j) = zero
         if (v(j) .eq. zero) go to 50

!        determine a givens rotation which eliminates the
!        j-th element of v.

         if (dabs(v(n)) .ge. dabs(v(j))) go to 20
            cotan = v(n)/v(j)
            sin = p5/dsqrt(p25+p25*cotan**2)
            cos = sin*cotan
            tau = one
            if (dabs(cos)*giant .gt. one) tau = one/cos
            go to 30
   20    continue
            tan = v(j)/v(n)
            cos = p5/dsqrt(p25+p25*tan**2)
            sin = cos*tan
            tau = sin
   30    continue

!        apply the transformation to v and store the information
!        necessary to recover the givens rotation.

         v(n) = sin*v(j) + cos*v(n)
         v(j) = tau

!        apply the transformation to s and extend the spike in w.

         l = jj
         do 40 i = j, m
            temp = cos*s(l) - sin*w(i)
            w(i) = sin*s(l) + cos*w(i)
            s(l) = temp
            l = l + 1
   40       continue
   50    continue
   60    continue
   70 continue

!     add the spike from the rank 1 update to w.

      do 80 i = 1, m
         w(i) = w(i) + v(n)*u(i)
   80    continue

!     eliminate the spike.

      sing = .false.
      if (nm1 .lt. 1) go to 140
      do 130 j = 1, nm1
         if (w(j) .eq. zero) go to 120

!        determine a givens rotation which eliminates the
!        j-th element of the spike.

         if (dabs(s(jj)) .ge. dabs(w(j))) go to 90
            cotan = s(jj)/w(j)
            sin = p5/dsqrt(p25+p25*cotan**2)
            cos = sin*cotan
            tau = one
            if (dabs(cos)*giant .gt. one) tau = one/cos
            go to 100
   90    continue
            tan = w(j)/s(jj)
            cos = p5/dsqrt(p25+p25*tan**2)
            sin = cos*tan
            tau = sin
  100    continue

!        apply the transformation to s and reduce the spike in w.

         l = jj
         do 110 i = j, m
            temp = cos*s(l) + sin*w(i)
            w(i) = -sin*s(l) + cos*w(i)
            s(l) = temp
            l = l + 1
  110       continue

!        store the information necessary to recover the
!        givens rotation.

         w(j) = tau
  120    continue

!        test for zero diagonal elements in the output s.

         if (s(jj) .eq. zero) sing = .true.
         jj = jj + (m - j + 1)
  130    continue
  140 continue

!     move w back into the last column of the output s.

      l = jj
      do 150 i = n, m
         s(l) = w(i)
         l = l + 1
  150    continue
      if (s(jj) .eq. zero) sing = .true.
      return

!     last card of subroutine r1updt.

      end
