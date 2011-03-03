      subroutine bandv1(nm,n,mbw,a,e21,m,w,z,ierr,nv,rv,rv6)
!
      IMPLICIT NONE
      integer i,j,k,m,n,r,ii,ij,jj,kj,mb,m1,nm,nv,ij1,its,kj1,mbw,m21, &
              ierr,maxj,maxk,group
      double precision a(nm,mbw),w(m),z(nm,m),rv(nv),rv6(n)
      double precision u,v,uk,xu,x0,x1,e21,eps2,eps3,eps4,norm,order, &
             epslon,pythag
      DOUBLE PRECISION ISEED( 4 )
!
!     this subroutine finds those eigenvectors of a real symmetric
!     band matrix corresponding to specified eigenvalues, using inverse
!     iteration.  the subroutine may also be used to solve systems
!     of linear equations with a symmetric or non-symmetric band
!     coefficient matrix.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        mbw is the number of columns of the array a used to store the
!          band matrix.  if the matrix is symmetric, mbw is its (half)
!          band width, denoted mb and defined as the number of adjacent
!          diagonals, including the principal diagonal, required to
!          specify the non-zero portion of the lower triangle of the
!          matrix.  if the subroutine is being used to solve systems
!          of linear equations and the coefficient matrix is not
!          symmetric, it must however have the same number of adjacent
!          diagonals above the main diagonal as below, and in this
!          case, mbw=2*mb-1.
!
!        a contains the lower triangle of the symmetric band input
!          matrix stored as an n by mb array.  its lowest subdiagonal
!          is stored in the last n+1-mb positions of the first column,
!          its next subdiagonal in the last n+2-mb positions of the
!          second column, further subdiagonals similarly, and finally
!          its principal diagonal in the n positions of column mb.
!          if the subroutine is being used to solve systems of linear
!          equations and the coefficient matrix is not symmetric, a is
!          n by 2*mb-1 instead with lower triangle as above and with
!          its first superdiagonal stored in the first n-1 positions of
!          column mb+1, its second superdiagonal in the first n-2
!          positions of column mb+2, further superdiagonals similarly,
!          and finally its highest superdiagonal in the first n+1-mb
!          positions of the last column.
!          contents of storages not part of the matrix are arbitrary.
!
!        e21 specifies the ordering of the eigenvalues and contains
!            0.0d0 if the eigenvalues are in ascending order, or
!            2.0d0 if the eigenvalues are in descending order.
!          if the subroutine is being used to solve systems of linear
!          equations, e21 should be set to 1.0d0 if the coefficient
!          matrix is symmetric and to -1.0d0 if not.
!
!        m is the number of specified eigenvalues or the number of
!          systems of linear equations.
!
!        w contains the m eigenvalues in ascending or descending order.
!          if the subroutine is being used to solve systems of linear
!          equations (a-w(r)*i)*x(r)=b(r), where i is the identity
!          matrix, w(r) should be set accordingly, for r=1,2,...,m.
!
!        z contains the constant matrix columns (b(r),r=1,2,...,m), if
!          the subroutine is used to solve systems of linear equations.
!
!        nv must be set to the dimension of the array parameter rv
!          as declared in the calling program dimension statement.
!
!     on output
!
!        a and w are unaltered.
!
!        z contains the associated set of orthogonal eigenvectors.
!          any vector which fails to converge is set to zero.  if the
!          subroutine is used to solve systems of linear equations,
!          z contains the solution matrix columns (x(r),r=1,2,...,m).
!
!        ierr is set to
!          zero       for normal return,
!          -r         if the eigenvector corresponding to the r-th
!                     eigenvalue fails to converge, or if the r-th
!                     system of linear equations is nearly singular.
!
!        rv and rv6 are temporary storage arrays.  note that rv is
!          of dimension at least n*(2*mb-1).  if the subroutine
!          is being used to solve systems of linear equations, the
!          determinant (up to sign) of a-w(m)*i is available, upon
!          return, as the product of the first n elements of rv.
!
!     calls pythag for  dsqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ISEED(1) = 1
      ISEED(2) = 1
      ISEED(3) = 1
      ISEED(4) = 1
      ierr = 0
      if (m .eq. 0) go to 1001
      mb = mbw
      if (e21 .lt. 0.0d0) mb = (mbw + 1) / 2
      m1 = mb - 1
      m21 = m1 + mb
      order = 1.0d0 - dabs(e21)
!     .......... find vectors by inverse iteration ..........
      do 920 r = 1, m
         its = 1
         x1 = w(r)
         if (r .ne. 1) go to 100
!     .......... compute norm of matrix ..........
         norm = 0.0d0
!
         do 60 j = 1, mb
            jj = mb + 1 - j
            kj = jj + m1
            ij = 1
            v = 0.0d0
!
            do 40 i = jj, n
               v = v + dabs(a(i,j))
               if (e21 .ge. 0.0d0) go to 40
               v = v + dabs(a(ij,kj))
               ij = ij + 1
   40       continue
!
            norm = dmax1(norm,v)
   60    continue
!
         if (e21 .lt. 0.0d0) norm = 0.5d0 * norm
!     .......... eps2 is the criterion for grouping,
!                eps3 replaces zero pivots and equal
!                roots are modified by eps3,
!                eps4 is taken very small to avoid overflow ..........
         if (norm .eq. 0.0d0) norm = 1.0d0
         eps2 = 1.0d-3 * norm * dabs(order)
         eps3 = epslon(norm)
         uk = n
         uk = dsqrt(uk)
         eps4 = uk * eps3
   80    group = 0
         go to 120
!     .......... look for close or coincident roots ..........
  100    if (dabs(x1-x0) .ge. eps2) go to 80
         group = group + 1
         if (order * (x1 - x0) .le. 0.0d0) x1 = x0 + order * eps3
!     .......... expand matrix, subtract eigenvalue,
!                and initialize vector ..........
  120    do 200 i = 1, n
            ij = i + min0(0,i-m1) * n
            kj = ij + mb * n
            ij1 = kj + m1 * n
            if (m1 .eq. 0) go to 180
!
            do 150 j = 1, m1
               if (ij .gt. m1) go to 125
               if (ij .gt. 0) go to 130
               rv(ij1) = 0.0d0
               ij1 = ij1 + n
               go to 130
  125          rv(ij) = a(i,j)
  130          ij = ij + n
               ii = i + j
               if (ii .gt. n) go to 150
               jj = mb - j
               if (e21 .ge. 0.0d0) go to 140
               ii = i
               jj = mb + j
  140          rv(kj) = a(ii,jj)
               kj = kj + n
  150       continue
!
  180       rv(ij) = a(i,mb) - x1
            rv6(i) = eps4
            if (order .eq. 0.0d0) rv6(i) = z(i,r)
  200    continue
!
         if (m1 .eq. 0) go to 600
!     .......... elimination with interchanges ..........
         do 580 i = 1, n
            ii = i + 1
            maxk = min0(i+m1-1,n)
            maxj = min0(n-i,m21-2) * n
!
            do 360 k = i, maxk
               kj1 = k
               j = kj1 + n
               jj = j + maxj
!
               do 340 kj = j, jj, n
                  rv(kj1) = rv(kj)
                  kj1 = kj
  340          continue
!
               rv(kj1) = 0.0d0
  360       continue
!
            if (i .eq. n) go to 580
            u = 0.0d0
            maxk = min0(i+m1,n)
            maxj = min0(n-ii,m21-2) * n
!
            do 450 j = i, maxk
               if (dabs(rv(j)) .lt. dabs(u)) go to 450
               u = rv(j)
               k = j
  450       continue
!
            j = i + n
            jj = j + maxj
            if (k .eq. i) go to 520
            kj = k
!
            do 500 ij = i, jj, n
               v = rv(ij)
               rv(ij) = rv(kj)
               rv(kj) = v
               kj = kj + n
  500       continue
!
            if (order .ne. 0.0d0) go to 520
            v = rv6(i)
            rv6(i) = rv6(k)
            rv6(k) = v
  520       if (u .eq. 0.0d0) go to 580
!
            do 560 k = ii, maxk
               v = rv(k) / u
               kj = k
!
               do 540 ij = j, jj, n
                  kj = kj + n
                  rv(kj) = rv(kj) - v * rv(ij)
  540          continue
!
               if (order .eq. 0.0d0) rv6(k) = rv6(k) - v * rv6(i)
  560       continue
!
  580    continue
!     .......... back substitution
!                for i=n step -1 until 1 do -- ..........
  600    do 630 ii = 1, n
            i = n + 1 - ii
            maxj = min0(ii,m21)
            if (maxj .eq. 1) go to 620
            ij1 = i
            j = ij1 + n
            jj = j + (maxj - 2) * n
!
            do 610 ij = j, jj, n
               ij1 = ij1 + 1
               rv6(i) = rv6(i) - rv(ij) * rv6(ij1)
  610       continue
!
  620       v = rv(i)
            if (dabs(v) .ge. eps3) go to 625
!     .......... set error -- nearly singular linear system ..........
            if (order .eq. 0.0d0) ierr = -r
            v = dsign(eps3,v)
  625       rv6(i) = rv6(i) / v
  630    continue
!
         xu = 1.0d0
         if (order .eq. 0.0d0) go to 870
!     .......... orthogonalize with respect to previous
!                members of group ..........
         if (group .eq. 0) go to 700
!
         do 680 jj = 1, group
            j = r - group - 1 + jj
            xu = 0.0d0
!
            do 640 i = 1, n
  640       xu = xu + rv6(i) * z(i,j)
!
            do 660 i = 1, n
  660       rv6(i) = rv6(i) - xu * z(i,j)
!
  680    continue
!
  700    norm = 0.0d0
!
         do 720 i = 1, n
  720    norm = norm + dabs(rv6(i))
!
         if (norm .ge. 0.1d0) go to 840
!     .......... in-line procedure for choosing
!                a new starting vector ..........
         if (its .ge. n) go to 830
         its = its + 1


         CALL DLARNV( 2, ISEED, N, RV6 )

!         xu = eps4 / (uk + 1.0d0)
!         rv6(1) = eps4
!c
!         do 760 i = 2, n
!  760    rv6(i) = xu
!c
!         rv6(its) = rv6(its) - eps4 * uk

         go to 600
!     .......... set error -- non-converged eigenvector ..........
  830    ierr = -r
         xu = 0.0d0
         go to 870
!     .......... normalize so that sum of squares is
!                1 and expand to full order ..........
  840    u = 0.0d0
!
         do 860 i = 1, n
  860    u = pythag(u,rv6(i))
!
         xu = 1.0d0 / u
!
  870    do 900 i = 1, n
  900    z(i,r) = rv6(i) * xu
!
         x0 = x1
  920 continue
!
 1001 return
      end
