      subroutine dchdc(a,lda,p,work,jpvt,job,info)
      implicit real*8 (a-h,o-z)
      integer lda,p,jpvt(1),job,info
      double precision a(lda,*),work(*)
!
!     Linpack routine
!     dchdc computes the cholesky decomposition of a positive definite
!     matrix.  a pivoting option allows the user to estimate the
!     condition of a positive definite matrix or determine the rank
!     of a positive semidefinite matrix.
!
!     on entry
!
!         a      double precision(lda,p).
!                a contains the matrix whose decomposition is to
!                be computed.  onlt the upper half of a need be stored.
!                the lower part of the array a is not referenced.
!
!         lda    integer.
!                lda is the leading dimension of the array a.
!
!         p      integer.
!                p is the order of the matrix.
!
!         work   double precision.
!                work is a work array.
!
!         jpvt   integer(p).
!                jpvt contains integers that control the selection
!                of the pivot elements, if pivoting has been requested.
!                each diagonal element a(k,k)
!                is placed in one of three classes according to the
!                value of jpvt(k).
!
!                   if jpvt(k) .gt. 0, then x(k) is an initial
!                                      element.
!
!                   if jpvt(k) .eq. 0, then x(k) is a free element.
!
!                   if jpvt(k) .lt. 0, then x(k) is a final element.
!
!                before the decomposition is computed, initial elements
!                are moved by symmetric row and column interchanges to
!                the beginning of the array a and final
!                elements to the end.  both initial and final elements
!                are frozen in place during the computation and only
!                free elements are moved.  at the k-th stage of the
!                reduction, if a(k,k) is occupied by a free element
!                it is interchanged with the largest free element
!                a(l,l) with l .ge. k.  jpvt is not referenced if
!                job .eq. 0.
!
!        job     integer.
!                job is an integer that initiates column pivoting.
!                if job .eq. 0, no pivoting is done.
!                if job .ne. 0, pivoting is done.
!
!     on return
!
!         a      a contains in its upper half the cholesky factor
!                of the matrix a as it has been permuted by pivoting.
!
!         jpvt   jpvt(j) contains the index of the diagonal element
!                of a that was moved into the j-th position,
!                provided pivoting was requested.
!
!         info   contains the index of the last positive diagonal
!                element of the cholesky factor.
!
!     for positive definite matrices info = p is the normal return.
!     for pivoting with positive semidefinite matrices info will
!     in general be less than p.  however, info may be greater than
!     the rank of a, since rounding error can cause an otherwise zero
!     element to be positive. indefinite systems will always cause
!     info to be less than p.
!
!     linpack. this version dated 08/14/78 .
!     j.j. dongarra and g.w. stewart, argonne national laboratory and
!     university of maryland.
!
!
!     blas daxpy,dswap
!     fortran dsqrt
!
!     internal variables
!
      integer pu,pl,plp1,i,j,jp,jt,k,kb,km1,kp1,l,maxl
      double precision temp
      double precision maxdia
      logical swapk,negk
!
      pl = 1
      pu = 0
      info = p
      if (job .eq. 0) go to 160
!
!        pivoting has been requested. rearrange the
!        the elements according to jpvt.
!
         do 70 k = 1, p
            swapk = jpvt(k) .gt. 0
            negk = jpvt(k) .lt. 0
            jpvt(k) = k
            if (negk) jpvt(k) = -jpvt(k)
            if (.not.swapk) go to 60
               if (k .eq. pl) go to 50
                  call dswap(pl-1,a(1,k),1,a(1,pl),1)
                  temp = a(k,k)
                  a(k,k) = a(pl,pl)
                  a(pl,pl) = temp
                  plp1 = pl + 1
                  if (p .lt. plp1) go to 40
                  do 30 j = plp1, p
                     if (j .ge. k) go to 10
                        temp = a(pl,j)
                        a(pl,j) = a(j,k)
                        a(j,k) = temp
                     go to 20
   10                continue
                     if (j .eq. k) go to 20
                        temp = a(k,j)
                        a(k,j) = a(pl,j)
                        a(pl,j) = temp
   20                continue
   30             continue
   40             continue
                  jpvt(k) = jpvt(pl)
                  jpvt(pl) = k
   50          continue
               pl = pl + 1
   60       continue
   70    continue
         pu = p
         if (p .lt. pl) go to 150
         do 140 kb = pl, p
            k = p - kb + pl
            if (jpvt(k) .ge. 0) go to 130
               jpvt(k) = -jpvt(k)
               if (pu .eq. k) go to 120
                  call dswap(k-1,a(1,k),1,a(1,pu),1)
                  temp = a(k,k)
                  a(k,k) = a(pu,pu)
                  a(pu,pu) = temp
                  kp1 = k + 1
                  if (p .lt. kp1) go to 110
                  do 100 j = kp1, p
                     if (j .ge. pu) go to 80
                        temp = a(k,j)
                        a(k,j) = a(j,pu)
                        a(j,pu) = temp
                     go to 90
   80                continue
                     if (j .eq. pu) go to 90
                        temp = a(k,j)
                        a(k,j) = a(pu,j)
                        a(pu,j) = temp
   90                continue
  100             continue
  110             continue
                  jt = jpvt(k)
                  jpvt(k) = jpvt(pu)
                  jpvt(pu) = jt
  120          continue
               pu = pu - 1
  130       continue
  140    continue
  150    continue
  160 continue
      do 270 k = 1, p
!
!        reduction loop.
!
         maxdia = a(k,k)
         kp1 = k + 1
         maxl = k
!
!        determine the pivot element.
!
         if (k .lt. pl .or. k .ge. pu) go to 190
            do 180 l = kp1, pu
               if (a(l,l) .le. maxdia) go to 170
                  maxdia = a(l,l)
                  maxl = l
  170          continue
  180       continue
  190    continue
!
!        quit if the pivot element is not positive.
!
         if (maxdia .gt. 0.0d0) go to 200
            info = k - 1
!     ......exit
            go to 280
  200    continue
         if (k .eq. maxl) go to 210
!
!           start the pivoting and update jpvt.
!
            km1 = k - 1
            call dswap(km1,a(1,k),1,a(1,maxl),1)
            a(maxl,maxl) = a(k,k)
            a(k,k) = maxdia
            jp = jpvt(maxl)
            jpvt(maxl) = jpvt(k)
            jpvt(k) = jp
  210    continue
!
!        reduction step. pivoting is contained across the rows.
!
         work(k) = dsqrt(a(k,k))
         a(k,k) = work(k)
         if (p .lt. kp1) go to 260
         do 250 j = kp1, p
            if (k .eq. maxl) go to 240
               if (j .ge. maxl) go to 220
                  temp = a(k,j)
                  a(k,j) = a(j,maxl)
                  a(j,maxl) = temp
               go to 230
  220          continue
               if (j .eq. maxl) go to 230
                  temp = a(k,j)
                  a(k,j) = a(maxl,j)
                  a(maxl,j) = temp
  230          continue
  240       continue
            a(k,j) = a(k,j)/work(k)
            work(j) = a(k,j)
            temp = -a(k,j)
            call daxpy(j-k,temp,work(kp1),1,a(kp1,j),1)
  250    continue
  260    continue
  270 continue
  280 continue
      return
      end
