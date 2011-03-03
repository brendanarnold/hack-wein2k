      SUBROUTINE DLADIV( A, B, C, D, P, Q )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, P, Q
!     ..
!
!  Purpose
!  =======
!
!  DLADIV performs complex division in  real arithmetic
!
!                        a + i*b
!             p + i*q = ---------
!                        c + i*d
!
!  The algorithm is due to Robert L. Smith and can be found
!  in D. Knuth, The art of Computer Programming, Vol.2, p.195
!
!  Arguments
!  =========
!
!  A       (input) DOUBLE PRECISION
!  B       (input) DOUBLE PRECISION
!  C       (input) DOUBLE PRECISION
!  D       (input) DOUBLE PRECISION
!          The scalars a, b, c, and d in the above expression.
!
!  P       (output) DOUBLE PRECISION
!  Q       (output) DOUBLE PRECISION
!          The scalars p and q in the above expression.
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION   E, F
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
      IF( ABS( D ).LT.ABS( C ) ) THEN
         E = D / C
         F = C + D*E
         P = ( A+B*E ) / F
         Q = ( B-A*E ) / F
      ELSE
         E = C / D
         F = D + C*E
         P = ( B+A*E ) / F
         Q = ( -A+B*E ) / F
      END IF
!
      RETURN
!
!     End of DLADIV
!
      END
      SUBROUTINE DLAE2( A, B, C, RT1, RT2 )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, RT1, RT2
!     ..
!
!  Purpose
!  =======
!
!  DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
!     [  A   B  ]
!     [  B   C  ].
!  On return, RT1 is the eigenvalue of larger absolute value, and RT2
!  is the eigenvalue of smaller absolute value.
!
!  Arguments
!  =========
!
!  A       (input) DOUBLE PRECISION
!          The (1,1) element of the 2-by-2 matrix.
!
!  B       (input) DOUBLE PRECISION
!          The (1,2) and (2,1) elements of the 2-by-2 matrix.
!
!  C       (input) DOUBLE PRECISION
!          The (2,2) element of the 2-by-2 matrix.
!
!  RT1     (output) DOUBLE PRECISION
!          The eigenvalue of larger absolute value.
!
!  RT2     (output) DOUBLE PRECISION
!          The eigenvalue of smaller absolute value.
!
!  Further Details
!  ===============
!
!  RT1 is accurate to a few ulps barring over/underflow.
!
!  RT2 may be inaccurate if there is massive cancellation in the
!  determinant A*C-B*B; higher precision or correctly rounded or
!  correctly truncated arithmetic would be needed to compute RT2
!  accurately in all cases.
!
!  Overflow is possible only if RT1 is within a factor of 5 of overflow.
!  Underflow is harmless if the input data is 0 or exceeds
!     underflow_threshold / macheps.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   AB, ACMN, ACMX, ADF, DF, RT, SM, TB
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
!     ..
!     .. Executable Statements ..
!
!     Compute the eigenvalues
!
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
!
!        Includes case AB=ADF=0
!
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
!
!        Includes case RT1 = RT2 = 0
!
         RT1 = HALF*RT
         RT2 = -HALF*RT
      END IF
      RETURN
!
!     End of DLAE2
!
      END
      SUBROUTINE DLAEBZ( IJOB, NITMAX, N, MMAX, MINP, NBMIN, ABSTOL, &
                         RELTOL, PIVMIN, D, E, E2, NVAL, AB, C, MOUT, &
                         NAB, WORK, IWORK, INFO )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            IJOB, INFO, MINP, MMAX, MOUT, N, NBMIN, NITMAX
      DOUBLE PRECISION   ABSTOL, PIVMIN, RELTOL
!     ..
!     .. Array Arguments ..
      INTEGER            IWORK( * ), NAB( MMAX, * ), NVAL( * )
      DOUBLE PRECISION   AB( MMAX, * ), C( * ), D( * ), E( * ), E2( * ), &
                         WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLAEBZ contains the iteration loops which compute and use the
!  function N(w), which is the count of eigenvalues of a symmetric
!  tridiagonal matrix T less than or equal to its argument  w.  It
!  performs a choice of two types of loops:
!
!  IJOB=1, followed by
!  IJOB=2: It takes as input a list of intervals and returns a list of
!          sufficiently small intervals whose union contains the same
!          eigenvalues as the union of the original intervals.
!          The input intervals are (AB(j,1),AB(j,2)], j=1,...,MINP.
!          The output interval (AB(j,1),AB(j,2)] will contain
!          eigenvalues NAB(j,1)+1,...,NAB(j,2), where 1 <= j <= MOUT.
!
!  IJOB=3: It performs a binary search in each input interval
!          (AB(j,1),AB(j,2)] for a point  w(j)  such that
!          N(w(j))=NVAL(j), and uses  C(j)  as the starting point of
!          the search.  If such a w(j) is found, then on output
!          AB(j,1)=AB(j,2)=w.  If no such w(j) is found, then on output
!          (AB(j,1),AB(j,2)] will be a small interval containing the
!          point where N(w) jumps through NVAL(j), unless that point
!          lies outside the initial interval.
!
!  Note that the intervals are in all cases half-open intervals,
!  i.e., of the form  (a,b] , which includes  b  but not  a .
!
!  To avoid underflow, the matrix should be scaled so that its largest
!  element is no greater than  overflow**(1/2) * underflow**(1/4)
!  in absolute value.  To assure the most accurate computation
!  of small eigenvalues, the matrix should be scaled to be
!  not much smaller than that, either.
!
!  See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
!  Matrix", Report CS41, Computer Science Dept., Stanford
!  University, July 21, 1966
!
!  Note: the arguments are, in general, *not* checked for unreasonable
!  values.
!
!  Arguments
!  =========
!
!  IJOB    (input) INTEGER
!          Specifies what is to be done:
!          = 1:  Compute NAB for the initial intervals.
!          = 2:  Perform bisection iteration to find eigenvalues of T.
!          = 3:  Perform bisection iteration to invert N(w), i.e.,
!                to find a point which has a specified number of
!                eigenvalues of T to its left.
!          Other values will cause DLAEBZ to return with INFO=-1.
!
!  NITMAX  (input) INTEGER
!          The maximum number of "levels" of bisection to be
!          performed, i.e., an interval of width W will not be made
!          smaller than 2^(-NITMAX) * W.  If not all intervals
!          have converged after NITMAX iterations, then INFO is set
!          to the number of non-converged intervals.
!
!  N       (input) INTEGER
!          The dimension n of the tridiagonal matrix T.  It must be at
!          least 1.
!
!  MMAX    (input) INTEGER
!          The maximum number of intervals.  If more than MMAX intervals
!          are generated, then DLAEBZ will quit with INFO=MMAX+1.
!
!  MINP    (input) INTEGER
!          The initial number of intervals.  It may not be greater than
!          MMAX.
!
!  NBMIN   (input) INTEGER
!          The smallest number of intervals that should be processed
!          using a vector loop.  If zero, then only the scalar loop
!          will be used.
!
!  ABSTOL  (input) DOUBLE PRECISION
!          The minimum (absolute) width of an interval.  When an
!          interval is narrower than ABSTOL, or than RELTOL times the
!          larger (in magnitude) endpoint, then it is considered to be
!          sufficiently small, i.e., converged.  This must be at least
!          zero.
!
!  RELTOL  (input) DOUBLE PRECISION
!          The minimum relative width of an interval.  When an interval
!          is narrower than ABSTOL, or than RELTOL times the larger (in
!          magnitude) endpoint, then it is considered to be
!          sufficiently small, i.e., converged.  Note: this should
!          always be at least radix*machine epsilon.
!
!  PIVMIN  (input) DOUBLE PRECISION
!          The minimum absolute value of a "pivot" in the Sturm
!          sequence loop.  This *must* be at least  max |e(j)**2| *
!          safe_min  and at least safe_min, where safe_min is at least
!          the smallest number that can divide one without overflow.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T.
!
!  E       (input) DOUBLE PRECISION array, dimension (N)
!          The offdiagonal elements of the tridiagonal matrix T in
!          positions 1 through N-1.  E(N) is arbitrary.
!
!  E2      (input) DOUBLE PRECISION array, dimension (N)
!          The squares of the offdiagonal elements of the tridiagonal
!          matrix T.  E2(N) is ignored.
!
!  NVAL    (input/output) INTEGER array, dimension (MINP)
!          If IJOB=1 or 2, not referenced.
!          If IJOB=3, the desired values of N(w).  The elements of NVAL
!          will be reordered to correspond with the intervals in AB.
!          Thus, NVAL(j) on output will not, in general be the same as
!          NVAL(j) on input, but it will correspond with the interval
!          (AB(j,1),AB(j,2)] on output.
!
!  AB      (input/output) DOUBLE PRECISION array, dimension (MMAX,2)
!          The endpoints of the intervals.  AB(j,1) is  a(j), the left
!          endpoint of the j-th interval, and AB(j,2) is b(j), the
!          right endpoint of the j-th interval.  The input intervals
!          will, in general, be modified, split, and reordered by the
!          calculation.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (MMAX)
!          If IJOB=1, ignored.
!          If IJOB=2, workspace.
!          If IJOB=3, then on input C(j) should be initialized to the
!          first search point in the binary search.
!
!  MOUT    (output) INTEGER
!          If IJOB=1, the number of eigenvalues in the intervals.
!          If IJOB=2 or 3, the number of intervals output.
!          If IJOB=3, MOUT will equal MINP.
!
!  NAB     (input/output) INTEGER array, dimension (MMAX,2)
!          If IJOB=1, then on output NAB(i,j) will be set to N(AB(i,j)).
!          If IJOB=2, then on input, NAB(i,j) should be set.  It must
!             satisfy the condition:
!             N(AB(i,1)) <= NAB(i,1) <= NAB(i,2) <= N(AB(i,2)),
!             which means that in interval i only eigenvalues
!             NAB(i,1)+1,...,NAB(i,2) will be considered.  Usually,
!             NAB(i,j)=N(AB(i,j)), from a previous call to DLAEBZ with
!             IJOB=1.
!             On output, NAB(i,j) will contain
!             max(na(k),min(nb(k),N(AB(i,j)))), where k is the index of
!             the input interval that the output interval
!             (AB(j,1),AB(j,2)] came from, and na(k) and nb(k) are the
!             the input values of NAB(k,1) and NAB(k,2).
!          If IJOB=3, then on output, NAB(i,j) contains N(AB(i,j)),
!             unless N(w) > NVAL(i) for all search points  w , in which
!             case NAB(i,1) will not be modified, i.e., the output
!             value will be the same as the input value (modulo
!             reorderings -- see NVAL and AB), or unless N(w) < NVAL(i)
!             for all search points  w , in which case NAB(i,2) will
!             not be modified.  Normally, NAB should be set to some
!             distinctive value(s) before DLAEBZ is called.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (MMAX)
!          Workspace.
!
!  IWORK   (workspace) INTEGER array, dimension (MMAX)
!          Workspace.
!
!  INFO    (output) INTEGER
!          = 0:       All intervals converged.
!          = 1--MMAX: The last INFO intervals did not converge.
!          = MMAX+1:  More than MMAX intervals were generated.
!
!  Further Details
!  ===============
!
!      This routine is intended to be called only by other LAPACK
!  routines, thus the interface is less user-friendly.  It is intended
!  for two purposes:
!
!  (a) finding eigenvalues.  In this case, DLAEBZ should have one or
!      more initial intervals set up in AB, and DLAEBZ should be called
!      with IJOB=1.  This sets up NAB, and also counts the eigenvalues.
!      Intervals with no eigenvalues would usually be thrown out at
!      this point.  Also, if not all the eigenvalues in an interval i
!      are desired, NAB(i,1) can be increased or NAB(i,2) decreased.
!      For example, set NAB(i,1)=NAB(i,2)-1 to get the largest
!      eigenvalue.  DLAEBZ is then called with IJOB=2 and MMAX
!      no smaller than the value of MOUT returned by the call with
!      IJOB=1.  After this (IJOB=2) call, eigenvalues NAB(i,1)+1
!      through NAB(i,2) are approximately AB(i,1) (or AB(i,2)) to the
!      tolerance specified by ABSTOL and RELTOL.
!
!  (b) finding an interval (a',b'] containing eigenvalues w(f),...,w(l).
!      In this case, start with a Gershgorin interval  (a,b).  Set up
!      AB to contain 2 search intervals, both initially (a,b).  One
!      NVAL element should contain  f-1  and the other should contain  l
!      , while C should contain a and b, resp.  NAB(i,1) should be -1
!      and NAB(i,2) should be N+1, to flag an error if the desired
!      interval does not lie in (a,b).  DLAEBZ is then called with
!      IJOB=3.  On exit, if w(f-1) < w(f), then one of the intervals --
!      j -- will have AB(j,1)=AB(j,2) and NAB(j,1)=NAB(j,2)=f-1, while
!      if, to the specified tolerance, w(f-k)=...=w(f+r), k > 0 and r
!      >= 0, then the interval will have  N(AB(j,1))=NAB(j,1)=f-k and
!      N(AB(j,2))=NAB(j,2)=f+r.  The cases w(l) < w(l+1) and
!      w(l-r)=...=w(l+k) are handled similarly.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, TWO, HALF
      PARAMETER          ( ZERO = 0.0D0, TWO = 2.0D0, &
                         HALF = 1.0D0 / TWO )
!     ..
!     .. Local Scalars ..
      INTEGER            ITMP1, ITMP2, J, JI, JIT, JP, KF, KFNEW, KL, &
                         KLNEW
      DOUBLE PRECISION   TMP1, TMP2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Check for Errors
!
      INFO = 0
      IF( IJOB.LT.1 .OR. IJOB.GT.3 ) THEN
         INFO = -1
         RETURN
      END IF
!
!     Initialize NAB
!
      IF( IJOB.EQ.1 ) THEN
!
!        Compute the number of eigenvalues in the initial intervals.
!
         MOUT = 0
         DO 30 JI = 1, MINP
            DO 20 JP = 1, 2
               TMP1 = D( 1 ) - AB( JI, JP )
               IF( ABS( TMP1 ).LT.PIVMIN ) &
                  TMP1 = -PIVMIN
               NAB( JI, JP ) = 0
               IF( TMP1.LE.ZERO ) &
                  NAB( JI, JP ) = 1
!
               DO 10 J = 2, N
                  TMP1 = D( J ) - E2( J-1 ) / TMP1 - AB( JI, JP )
                  IF( ABS( TMP1 ).LT.PIVMIN ) &
                     TMP1 = -PIVMIN
                  IF( TMP1.LE.ZERO ) &
                     NAB( JI, JP ) = NAB( JI, JP ) + 1
   10          CONTINUE
   20       CONTINUE
            MOUT = MOUT + NAB( JI, 2 ) - NAB( JI, 1 )
   30    CONTINUE
         RETURN
      END IF
!
!     Initialize for loop
!
!     KF and KL have the following meaning:
!        Intervals 1,...,KF-1 have converged.
!        Intervals KF,...,KL  still need to be refined.
!
      KF = 1
      KL = MINP
!
!     If IJOB=2, initialize C.
!     If IJOB=3, use the user-supplied starting point.
!
      IF( IJOB.EQ.2 ) THEN
         DO 40 JI = 1, MINP
            C( JI ) = HALF*( AB( JI, 1 )+AB( JI, 2 ) )
   40    CONTINUE
      END IF
!
!     Iteration loop
!
      DO 130 JIT = 1, NITMAX
!
!        Loop over intervals
!
         IF( KL-KF+1.GE.NBMIN .AND. NBMIN.GT.0 ) THEN
!
!           Begin of Parallel Version of the loop
!
            DO 60 JI = KF, KL
!
!              Compute N(c), the number of eigenvalues less than c
!
               WORK( JI ) = D( 1 ) - C( JI )
               IWORK( JI ) = 0
               IF( WORK( JI ).LE.PIVMIN ) THEN
                  IWORK( JI ) = 1
                  WORK( JI ) = MIN( WORK( JI ), -PIVMIN )
               END IF
!
               DO 50 J = 2, N
                  WORK( JI ) = D( J ) - E2( J-1 ) / WORK( JI ) - C( JI )
                  IF( WORK( JI ).LE.PIVMIN ) THEN
                     IWORK( JI ) = IWORK( JI ) + 1
                     WORK( JI ) = MIN( WORK( JI ), -PIVMIN )
                  END IF
   50          CONTINUE
   60       CONTINUE
!
            IF( IJOB.LE.2 ) THEN
!
!              IJOB=2: Choose all intervals containing eigenvalues.
!
               KLNEW = KL
               DO 70 JI = KF, KL
!
!                 Insure that N(w) is monotone
!
                  IWORK( JI ) = MIN( NAB( JI, 2 ), &
                                MAX( NAB( JI, 1 ), IWORK( JI ) ) )
!
!                 Update the Queue -- add intervals if both halves
!                 contain eigenvalues.
!
                  IF( IWORK( JI ).EQ.NAB( JI, 2 ) ) THEN
!
!                    No eigenvalue in the upper interval:
!                    just use the lower interval.
!
                     AB( JI, 2 ) = C( JI )
!
                  ELSE IF( IWORK( JI ).EQ.NAB( JI, 1 ) ) THEN
!
!                    No eigenvalue in the lower interval:
!                    just use the upper interval.
!
                     AB( JI, 1 ) = C( JI )
                  ELSE
                     KLNEW = KLNEW + 1
                     IF( KLNEW.LE.MMAX ) THEN
!
!                       Eigenvalue in both intervals -- add upper to
!                       queue.
!
                        AB( KLNEW, 2 ) = AB( JI, 2 )
                        NAB( KLNEW, 2 ) = NAB( JI, 2 )
                        AB( KLNEW, 1 ) = C( JI )
                        NAB( KLNEW, 1 ) = IWORK( JI )
                        AB( JI, 2 ) = C( JI )
                        NAB( JI, 2 ) = IWORK( JI )
                     ELSE
                        INFO = MMAX + 1
                     END IF
                  END IF
   70          CONTINUE
               IF( INFO.NE.0 ) &
                  RETURN
               KL = KLNEW
            ELSE
!
!              IJOB=3: Binary search.  Keep only the interval containing
!                      w   s.t. N(w) = NVAL
!
               DO 80 JI = KF, KL
                  IF( IWORK( JI ).LE.NVAL( JI ) ) THEN
                     AB( JI, 1 ) = C( JI )
                     NAB( JI, 1 ) = IWORK( JI )
                  END IF
                  IF( IWORK( JI ).GE.NVAL( JI ) ) THEN
                     AB( JI, 2 ) = C( JI )
                     NAB( JI, 2 ) = IWORK( JI )
                  END IF
   80          CONTINUE
            END IF
!
         ELSE
!
!           End of Parallel Version of the loop
!
!           Begin of Serial Version of the loop
!
            KLNEW = KL
            DO 100 JI = KF, KL
!
!              Compute N(w), the number of eigenvalues less than w
!
               TMP1 = C( JI )
               TMP2 = D( 1 ) - TMP1
               ITMP1 = 0
               IF( TMP2.LE.PIVMIN ) THEN
                  ITMP1 = 1
                  TMP2 = MIN( TMP2, -PIVMIN )
               END IF
!
!              A series of compiler directives to defeat vectorization
!              for the next loop
!
!$PL$ CMCHAR=' '
!DIR$          NEXTSCALAR
!$DIR          SCALAR
!DIR$          NEXT SCALAR
!VD$L          NOVECTOR
!DEC$          NOVECTOR
!VD$           NOVECTOR
!VDIR          NOVECTOR
!VOCL          LOOP,SCALAR
!IBM           PREFER SCALAR
!$PL$ CMCHAR='*'
!
               DO 90 J = 2, N
                  TMP2 = D( J ) - E2( J-1 ) / TMP2 - TMP1
                  IF( TMP2.LE.PIVMIN ) THEN
                     ITMP1 = ITMP1 + 1
                     TMP2 = MIN( TMP2, -PIVMIN )
                  END IF
   90          CONTINUE
!
               IF( IJOB.LE.2 ) THEN
!
!                 IJOB=2: Choose all intervals containing eigenvalues.
!
!                 Insure that N(w) is monotone
!
                  ITMP1 = MIN( NAB( JI, 2 ), &
                          MAX( NAB( JI, 1 ), ITMP1 ) )
!
!                 Update the Queue -- add intervals if both halves
!                 contain eigenvalues.
!
                  IF( ITMP1.EQ.NAB( JI, 2 ) ) THEN
!
!                    No eigenvalue in the upper interval:
!                    just use the lower interval.
!
                     AB( JI, 2 ) = TMP1
!
                  ELSE IF( ITMP1.EQ.NAB( JI, 1 ) ) THEN
!
!                    No eigenvalue in the lower interval:
!                    just use the upper interval.
!
                     AB( JI, 1 ) = TMP1
                  ELSE IF( KLNEW.LT.MMAX ) THEN
!
!                    Eigenvalue in both intervals -- add upper to queue.
!
                     KLNEW = KLNEW + 1
                     AB( KLNEW, 2 ) = AB( JI, 2 )
                     NAB( KLNEW, 2 ) = NAB( JI, 2 )
                     AB( KLNEW, 1 ) = TMP1
                     NAB( KLNEW, 1 ) = ITMP1
                     AB( JI, 2 ) = TMP1
                     NAB( JI, 2 ) = ITMP1
                  ELSE
                     INFO = MMAX + 1
                     RETURN
                  END IF
               ELSE
!
!                 IJOB=3: Binary search.  Keep only the interval
!                         containing  w  s.t. N(w) = NVAL
!
                  IF( ITMP1.LE.NVAL( JI ) ) THEN
                     AB( JI, 1 ) = TMP1
                     NAB( JI, 1 ) = ITMP1
                  END IF
                  IF( ITMP1.GE.NVAL( JI ) ) THEN
                     AB( JI, 2 ) = TMP1
                     NAB( JI, 2 ) = ITMP1
                  END IF
               END IF
  100       CONTINUE
            KL = KLNEW
!
!           End of Serial Version of the loop
!
         END IF
!
!        Check for convergence
!
         KFNEW = KF
         DO 110 JI = KF, KL
            TMP1 = ABS( AB( JI, 2 )-AB( JI, 1 ) )
            TMP2 = MAX( ABS( AB( JI, 2 ) ), ABS( AB( JI, 1 ) ) )
            IF( TMP1.LT.MAX( ABSTOL, PIVMIN, RELTOL*TMP2 ) .OR. &
                NAB( JI, 1 ).GE.NAB( JI, 2 ) ) THEN
!
!              Converged -- Swap with position KFNEW,
!                           then increment KFNEW
!
               IF( JI.GT.KFNEW ) THEN
                  TMP1 = AB( JI, 1 )
                  TMP2 = AB( JI, 2 )
                  ITMP1 = NAB( JI, 1 )
                  ITMP2 = NAB( JI, 2 )
                  AB( JI, 1 ) = AB( KFNEW, 1 )
                  AB( JI, 2 ) = AB( KFNEW, 2 )
                  NAB( JI, 1 ) = NAB( KFNEW, 1 )
                  NAB( JI, 2 ) = NAB( KFNEW, 2 )
                  AB( KFNEW, 1 ) = TMP1
                  AB( KFNEW, 2 ) = TMP2
                  NAB( KFNEW, 1 ) = ITMP1
                  NAB( KFNEW, 2 ) = ITMP2
                  IF( IJOB.EQ.3 ) THEN
                     ITMP1 = NVAL( JI )
                     NVAL( JI ) = NVAL( KFNEW )
                     NVAL( KFNEW ) = ITMP1
                  END IF
               END IF
               KFNEW = KFNEW + 1
            END IF
  110    CONTINUE
         KF = KFNEW
!
!        Choose Midpoints
!
         DO 120 JI = KF, KL
            C( JI ) = HALF*( AB( JI, 1 )+AB( JI, 2 ) )
  120    CONTINUE
!
!        If no more intervals to refine, quit.
!
         IF( KF.GT.KL ) &
            GO TO 140
  130 CONTINUE
!
!     Converged
!
  140 CONTINUE
      INFO = MAX( KL+1-KF, 0 )
      MOUT = KL
!
      RETURN
!
!     End of DLAEBZ
!
      END
      SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
!     ..
!
!  Purpose
!  =======
!
!  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
!     [  A   B  ]
!     [  B   C  ].
!  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
!  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
!  eigenvector for RT1, giving the decomposition
!
!     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
!     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
!
!  Arguments
!  =========
!
!  A       (input) DOUBLE PRECISION
!          The (1,1) element of the 2-by-2 matrix.
!
!  B       (input) DOUBLE PRECISION
!          The (1,2) element and the conjugate of the (2,1) element of
!          the 2-by-2 matrix.
!
!  C       (input) DOUBLE PRECISION
!          The (2,2) element of the 2-by-2 matrix.
!
!  RT1     (output) DOUBLE PRECISION
!          The eigenvalue of larger absolute value.
!
!  RT2     (output) DOUBLE PRECISION
!          The eigenvalue of smaller absolute value.
!
!  CS1     (output) DOUBLE PRECISION
!  SN1     (output) DOUBLE PRECISION
!          The vector (CS1, SN1) is a unit right eigenvector for RT1.
!
!  Further Details
!  ===============
!
!  RT1 is accurate to a few ulps barring over/underflow.
!
!  RT2 may be inaccurate if there is massive cancellation in the
!  determinant A*C-B*B; higher precision or correctly rounded or
!  correctly truncated arithmetic would be needed to compute RT2
!  accurately in all cases.
!
!  CS1 and SN1 are accurate to a few ulps barring over/underflow.
!
!  Overflow is possible only if RT1 is within a factor of 5 of overflow.
!  Underflow is harmless if the input data is 0 or exceeds
!     underflow_threshold / macheps.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            SGN1, SGN2
      DOUBLE PRECISION   AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM, &
                         TB, TN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
!     ..
!     .. Executable Statements ..
!
!     Compute the eigenvalues
!
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
!
!        Includes case AB=ADF=0
!
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
         SGN1 = -1
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
         SGN1 = 1
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
!
!        Includes case RT1 = RT2 = 0
!
         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF
!
!     Compute the eigenvector
!
      IF( DF.GE.ZERO ) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS( CS )
      IF( ACS.GT.AB ) THEN
         CT = -TB / CS
         SN1 = ONE / SQRT( ONE+CT*CT )
         CS1 = CT*SN1
      ELSE
         IF( AB.EQ.ZERO ) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS / TB
            CS1 = ONE / SQRT( ONE+TN*TN )
            SN1 = TN*CS1
         END IF
      END IF
      IF( SGN1.EQ.SGN2 ) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN
!
!     End of DLAEV2
!
      END
      SUBROUTINE DLAGTF( N, A, LAMBDA, B, C, TOL, D, IN, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, N
      DOUBLE PRECISION   LAMBDA, TOL
!     ..
!     .. Array Arguments ..
      INTEGER            IN( * )
      DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * )
!     ..
!
!  Purpose
!  =======
!
!  DLAGTF factorizes the matrix (T - lambda*I), where T is an n by n
!  tridiagonal matrix and lambda is a scalar, as
!
!     T - lambda*I = PLU,
!
!  where P is a permutation matrix, L is a unit lower tridiagonal matrix
!  with at most one non-zero sub-diagonal elements per column and U is
!  an upper triangular matrix with at most two non-zero super-diagonal
!  elements per column.
!
!  The factorization is obtained by Gaussian elimination with partial
!  pivoting and implicit row scaling.
!
!  The parameter LAMBDA is included in the routine so that DLAGTF may
!  be used, in conjunction with DLAGTS, to obtain eigenvectors of T by
!  inverse iteration.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix T.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, A must contain the diagonal elements of T.
!
!          On exit, A is overwritten by the n diagonal elements of the
!          upper triangular matrix U of the factorization of T.
!
!  LAMBDA  (input) DOUBLE PRECISION
!          On entry, the scalar lambda.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, B must contain the (n-1) super-diagonal elements of
!          T.
!
!          On exit, B is overwritten by the (n-1) super-diagonal
!          elements of the matrix U of the factorization of T.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, C must contain the (n-1) sub-diagonal elements of
!          T.
!
!          On exit, C is overwritten by the (n-1) sub-diagonal elements
!          of the matrix L of the factorization of T.
!
!  TOL     (input) DOUBLE PRECISION
!          On entry, a relative tolerance used to indicate whether or
!          not the matrix (T - lambda*I) is nearly singular. TOL should
!          normally be chose as approximately the largest relative error
!          in the elements of T. For example, if the elements of T are
!          correct to about 4 significant figures, then TOL should be
!          set to about 5*10**(-4). If TOL is supplied as less than eps,
!          where eps is the relative machine precision, then the value
!          eps is used in place of TOL.
!
!  D       (output) DOUBLE PRECISION array, dimension (N-2)
!          On exit, D is overwritten by the (n-2) second super-diagonal
!          elements of the matrix U of the factorization of T.
!
!  IN      (output) INTEGER array, dimension (N)
!          On exit, IN contains details of the permutation matrix P. If
!          an interchange occurred at the kth step of the elimination,
!          then IN(k) = 1, otherwise IN(k) = 0. The element IN(n)
!          returns the smallest positive integer j such that
!
!             abs( u(j,j) ).le. norm( (T - lambda*I)(j) )*TOL,
!
!          where norm( A(j) ) denotes the sum of the absolute values of
!          the jth row of the matrix A. If no such j exists then IN(n)
!          is returned as zero. If IN(n) is returned as positive, then a
!          diagonal element of U is small, indicating that
!          (T - lambda*I) is singular or nearly singular,
!
!  INFO    (output)
!          = 0   : successful exit
!          .lt. 0: if INFO = -k, the kth argument had an illegal value
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            K
      DOUBLE PRECISION   EPS, MULT, PIV1, PIV2, SCALE1, SCALE2, TEMP, TL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DLAGTF', -INFO )
         RETURN
      END IF
!
      IF( N.EQ.0 ) &
         RETURN
!
      A( 1 ) = A( 1 ) - LAMBDA
      IN( N ) = 0
      IF( N.EQ.1 ) THEN
         IF( A( 1 ).EQ.ZERO ) &
            IN( 1 ) = 1
         RETURN
      END IF
!
      EPS = DLAMCH( 'Epsilon' )
!
      TL = MAX( TOL, EPS )
      SCALE1 = ABS( A( 1 ) ) + ABS( B( 1 ) )
      DO 10 K = 1, N - 1
         A( K+1 ) = A( K+1 ) - LAMBDA
         SCALE2 = ABS( C( K ) ) + ABS( A( K+1 ) )
         IF( K.LT.( N-1 ) ) &
            SCALE2 = SCALE2 + ABS( B( K+1 ) )
         IF( A( K ).EQ.ZERO ) THEN
            PIV1 = ZERO
         ELSE
            PIV1 = ABS( A( K ) ) / SCALE1
         END IF
         IF( C( K ).EQ.ZERO ) THEN
            IN( K ) = 0
            PIV2 = ZERO
            SCALE1 = SCALE2
            IF( K.LT.( N-1 ) ) &
               D( K ) = ZERO
         ELSE
            PIV2 = ABS( C( K ) ) / SCALE2
            IF( PIV2.LE.PIV1 ) THEN
               IN( K ) = 0
               SCALE1 = SCALE2
               C( K ) = C( K ) / A( K )
               A( K+1 ) = A( K+1 ) - C( K )*B( K )
               IF( K.LT.( N-1 ) ) &
                  D( K ) = ZERO
            ELSE
               IN( K ) = 1
               MULT = A( K ) / C( K )
               A( K ) = C( K )
               TEMP = A( K+1 )
               A( K+1 ) = B( K ) - MULT*TEMP
               IF( K.LT.( N-1 ) ) THEN
                  D( K ) = B( K+1 )
                  B( K+1 ) = -MULT*D( K )
               END IF
               B( K ) = TEMP
               C( K ) = MULT
            END IF
         END IF
         IF( ( MAX( PIV1, PIV2 ).LE.TL ) .AND. ( IN( N ).EQ.0 ) ) &
            IN( N ) = K
   10 CONTINUE
      IF( ( ABS( A( N ) ).LE.SCALE1*TL ) .AND. ( IN( N ).EQ.0 ) ) &
         IN( N ) = N
!
      RETURN
!
!     End of DLAGTF
!
      END
      SUBROUTINE DLAGTS( JOB, N, A, B, C, D, IN, Y, TOL, INFO )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, JOB, N
      DOUBLE PRECISION   TOL
!     ..
!     .. Array Arguments ..
      INTEGER            IN( * )
      DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  DLAGTS may be used to solve one of the systems of equations
!
!     (T - lambda*I)*x = y   or   (T - lambda*I)'*x = y,
!
!  where T is an n by n tridiagonal matrix, for x, following the
!  factorization of (T - lambda*I) as
!
!     (T - lambda*I) = P*L*U ,
!
!  by routine DLAGTF. The choice of equation to be solved is
!  controlled by the argument JOB, and in each case there is an option
!  to perturb zero or very small diagonal elements of U, this option
!  being intended for use in applications such as inverse iteration.
!
!  Arguments
!  =========
!
!  JOB     (input) INTEGER
!          Specifies the job to be performed by DLAGTS as follows:
!          =  1: The equations  (T - lambda*I)x = y  are to be solved,
!                but diagonal elements of U are not to be perturbed.
!          = -1: The equations  (T - lambda*I)x = y  are to be solved
!                and, if overflow would otherwise occur, the diagonal
!                elements of U are to be perturbed. See argument TOL
!                below.
!          =  2: The equations  (T - lambda*I)'x = y  are to be solved,
!                but diagonal elements of U are not to be perturbed.
!          = -2: The equations  (T - lambda*I)'x = y  are to be solved
!                and, if overflow would otherwise occur, the diagonal
!                elements of U are to be perturbed. See argument TOL
!                below.
!
!  N       (input) INTEGER
!          The order of the matrix T.
!
!  A       (input) DOUBLE PRECISION array, dimension (N)
!          On entry, A must contain the diagonal elements of U as
!          returned from DLAGTF.
!
!  B       (input) DOUBLE PRECISION array, dimension (N-1)
!          On entry, B must contain the first super-diagonal elements of
!          U as returned from DLAGTF.
!
!  C       (input) DOUBLE PRECISION array, dimension (N-1)
!          On entry, C must contain the sub-diagonal elements of L as
!          returned from DLAGTF.
!
!  D       (input) DOUBLE PRECISION array, dimension (N-2)
!          On entry, D must contain the second super-diagonal elements
!          of U as returned from DLAGTF.
!
!  IN      (input) INTEGER array, dimension (N)
!          On entry, IN must contain details of the matrix P as returned
!          from DLAGTF.
!
!  Y       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the right hand side vector y.
!          On exit, Y is overwritten by the solution vector x.
!
!  TOL     (input/output) DOUBLE PRECISION
!          On entry, with  JOB .lt. 0, TOL should be the minimum
!          perturbation to be made to very small diagonal elements of U.
!          TOL should normally be chosen as about eps*norm(U), where eps
!          is the relative machine precision, but if TOL is supplied as
!          non-positive, then it is reset to eps*max( abs( u(i,j) ) ).
!          If  JOB .gt. 0  then TOL is not referenced.
!
!          On exit, TOL is changed as described above, only if TOL is
!          non-positive on entry. Otherwise TOL is unchanged.
!
!  INFO    (output) INTEGER
!          = 0   : successful exit
!          .lt. 0: if INFO = -i, the i-th argument had an illegal value
!          .gt. 0: overflow would occur when computing the INFO(th)
!                  element of the solution vector x. This can only occur
!                  when JOB is supplied as positive and either means
!                  that a diagonal element of U is very small, or that
!                  the elements of the right-hand side vector y are very
!                  large.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            K
      DOUBLE PRECISION   ABSAK, AK, BIGNUM, EPS, PERT, SFMIN, TEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      IF( ( ABS( JOB ).GT.2 ) .OR. ( JOB.EQ.0 ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAGTS', -INFO )
         RETURN
      END IF
!
      IF( N.EQ.0 ) &
         RETURN
!
      EPS = DLAMCH( 'Epsilon' )
      SFMIN = DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SFMIN
!
      IF( JOB.LT.0 ) THEN
         IF( TOL.LE.ZERO ) THEN
            TOL = ABS( A( 1 ) )
            IF( N.GT.1 ) &
               TOL = MAX( TOL, ABS( A( 2 ) ), ABS( B( 1 ) ) )
            DO 10 K = 3, N
               TOL = MAX( TOL, ABS( A( K ) ), ABS( B( K-1 ) ), &
                     ABS( D( K-2 ) ) )
   10       CONTINUE
            TOL = TOL*EPS
            IF( TOL.EQ.ZERO ) &
               TOL = EPS
         END IF
      END IF
!
      IF( ABS( JOB ).EQ.1 ) THEN
         DO 20 K = 2, N
            IF( IN( K-1 ).EQ.0 ) THEN
               Y( K ) = Y( K ) - C( K-1 )*Y( K-1 )
            ELSE
               TEMP = Y( K-1 )
               Y( K-1 ) = Y( K )
               Y( K ) = TEMP - C( K-1 )*Y( K )
            END IF
   20    CONTINUE
         IF( JOB.EQ.1 ) THEN
            DO 30 K = N, 1, -1
               IF( K.LE.N-2 ) THEN
                  TEMP = Y( K ) - B( K )*Y( K+1 ) - D( K )*Y( K+2 )
               ELSE IF( K.EQ.N-1 ) THEN
                  TEMP = Y( K ) - B( K )*Y( K+1 )
               ELSE
                  TEMP = Y( K )
               END IF
               AK = A( K )
               ABSAK = ABS( AK )
               IF( ABSAK.LT.ONE ) THEN
                  IF( ABSAK.LT.SFMIN ) THEN
                     IF( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) &
                          THEN
                        INFO = K
                        RETURN
                     ELSE
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     END IF
                  ELSE IF( ABS( TEMP ).GT.ABSAK*BIGNUM ) THEN
                     INFO = K
                     RETURN
                  END IF
               END IF
               Y( K ) = TEMP / AK
   30       CONTINUE
         ELSE
            DO 50 K = N, 1, -1
               IF( K.LE.N-2 ) THEN
                  TEMP = Y( K ) - B( K )*Y( K+1 ) - D( K )*Y( K+2 )
               ELSE IF( K.EQ.N-1 ) THEN
                  TEMP = Y( K ) - B( K )*Y( K+1 )
               ELSE
                  TEMP = Y( K )
               END IF
               AK = A( K )
               PERT = SIGN( TOL, AK )
   40          CONTINUE
               ABSAK = ABS( AK )
               IF( ABSAK.LT.ONE ) THEN
                  IF( ABSAK.LT.SFMIN ) THEN
                     IF( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) &
                          THEN
                        AK = AK + PERT
                        PERT = 2*PERT
                        GO TO 40
                     ELSE
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     END IF
                  ELSE IF( ABS( TEMP ).GT.ABSAK*BIGNUM ) THEN
                     AK = AK + PERT
                     PERT = 2*PERT
                     GO TO 40
                  END IF
               END IF
               Y( K ) = TEMP / AK
   50       CONTINUE
         END IF
      ELSE
!
!        Come to here if  JOB = 2 or -2
!
         IF( JOB.EQ.2 ) THEN
            DO 60 K = 1, N
               IF( K.GE.3 ) THEN
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 ) - D( K-2 )*Y( K-2 )
               ELSE IF( K.EQ.2 ) THEN
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 )
               ELSE
                  TEMP = Y( K )
               END IF
               AK = A( K )
               ABSAK = ABS( AK )
               IF( ABSAK.LT.ONE ) THEN
                  IF( ABSAK.LT.SFMIN ) THEN
                     IF( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) &
                          THEN
                        INFO = K
                        RETURN
                     ELSE
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     END IF
                  ELSE IF( ABS( TEMP ).GT.ABSAK*BIGNUM ) THEN
                     INFO = K
                     RETURN
                  END IF
               END IF
               Y( K ) = TEMP / AK
   60       CONTINUE
         ELSE
            DO 80 K = 1, N
               IF( K.GE.3 ) THEN
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 ) - D( K-2 )*Y( K-2 )
               ELSE IF( K.EQ.2 ) THEN
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 )
               ELSE
                  TEMP = Y( K )
               END IF
               AK = A( K )
               PERT = SIGN( TOL, AK )
   70          CONTINUE
               ABSAK = ABS( AK )
               IF( ABSAK.LT.ONE ) THEN
                  IF( ABSAK.LT.SFMIN ) THEN
                     IF( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) &
                          THEN
                        AK = AK + PERT
                        PERT = 2*PERT
                        GO TO 70
                     ELSE
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     END IF
                  ELSE IF( ABS( TEMP ).GT.ABSAK*BIGNUM ) THEN
                     AK = AK + PERT
                     PERT = 2*PERT
                     GO TO 70
                  END IF
               END IF
               Y( K ) = TEMP / AK
   80       CONTINUE
         END IF
!
         DO 90 K = N, 2, -1
            IF( IN( K-1 ).EQ.0 ) THEN
               Y( K-1 ) = Y( K-1 ) - C( K-1 )*Y( K )
            ELSE
               TEMP = Y( K-1 )
               Y( K-1 ) = Y( K )
               Y( K ) = TEMP - C( K-1 )*Y( K )
            END IF
   90    CONTINUE
      END IF
!
!     End of DLAGTS
!
      END
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          CMACH
!     ..
!
!  Purpose
!  =======
!
!  DLAMCH determines double precision machine parameters.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER*1
!          Specifies the value to be returned by DLAMCH:
!          = 'E' or 'e',   DLAMCH := eps
!          = 'S' or 's ,   DLAMCH := sfmin
!          = 'B' or 'b',   DLAMCH := base
!          = 'P' or 'p',   DLAMCH := eps*base
!          = 'N' or 'n',   DLAMCH := t
!          = 'R' or 'r',   DLAMCH := rnd
!          = 'M' or 'm',   DLAMCH := emin
!          = 'U' or 'u',   DLAMCH := rmin
!          = 'L' or 'l',   DLAMCH := emax
!          = 'O' or 'o',   DLAMCH := rmax
!
!          where
!
!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps*base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold - base**(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold  - (base**emax)*(1-eps)
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN, &
                         RND, SFMIN, SMALL, T
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAMC2
!     ..
!     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN, &
                         EMAX, RMAX, PREC
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
!
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
!
      DLAMCH = RMACH
      RETURN
!
!     End of DLAMCH
!
      END
!
!***********************************************************************
!
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
!     ..
!
!  Purpose
!  =======
!
!  DLAMC1 determines the machine parameters given by BETA, T, RND, and
!  IEEE1.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  IEEE1   (output) LOGICAL
!          Specifies whether rounding appears to be done in the IEEE
!          'round to nearest' style.
!
!  Further Details
!  ===============
!
!  The routine is based on the routine  ENVRON  by Malcolm and
!  incorporates suggestions by Gentleman and Marovich. See
!
!     Malcolm M. A. (1972) Algorithms to reveal properties of
!        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
!
!     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
!        that reveal properties of floating point arithmetic units.
!        Comms. of the ACM, 17, 276-277.
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
!     ..
!     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
!
!        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
!        IEEE1, T and RND.
!
!        Throughout this routine  we use the function  DLAMC3  to ensure
!        that relevant values are  stored and not held in registers,  or
!        are not affected by optimizers.
!
!        Compute  a = 2.0**m  with the  smallest positive integer m such
!        that
!
!           fl( a + 1.0 ) = a.
!
         A = 1
         C = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
!+       END WHILE
!
!        Now compute  b = 2.0**m  with the smallest positive integer m
!        such that
!
!           fl( a + b ) .gt. a.
!
         B = 1
         C = DLAMC3( A, B )
!
!+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF
!+       END WHILE
!
!        Now compute the base.  a and c  are neighbouring floating point
!        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
!        their difference is beta. Adding 0.25 to c is to ensure that it
!        is truncated to beta and not ( beta - 1 ).
!
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = C + QTR
!
!        Now determine whether rounding or chopping occurs,  by adding a
!        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
!
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) ) &
            LRND = .FALSE.
!
!        Try and decide whether rounding is done in the  IEEE  'round to
!        nearest' style. B/2 is half a unit in the last place of the two
!        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
!        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
!        A, but adding B/2 to SAVEC should change SAVEC.
!
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
!
!        Now find  the  mantissa, t.  It should  be the  integer part of
!        log to the base beta of a,  however it is safer to determine  t
!        by powering.  So we find t as the smallest positive integer for
!        which
!
!           fl( beta**t + 1.0 ) = 1.0.
!
         LT = 0
         A = 1
         C = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
!+       END WHILE
!
      END IF
!
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
!
!     End of DLAMC1
!
      END
!
!***********************************************************************
!
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
!     ..
!
!  Purpose
!  =======
!
!  DLAMC2 determines the machine parameters specified in its argument
!  list.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  EPS     (output) DOUBLE PRECISION
!          The smallest positive number such that
!
!             fl( 1.0 - EPS ) .LT. 1.0,
!
!          where fl denotes the computed value.
!
!  EMIN    (output) INTEGER
!          The minimum exponent before (gradual) underflow occurs.
!
!  RMIN    (output) DOUBLE PRECISION
!          The smallest normalized number for the machine, given by
!          BASE**( EMIN - 1 ), where  BASE  is the floating point value
!          of BETA.
!
!  EMAX    (output) INTEGER
!          The maximum exponent before overflow occurs.
!
!  RMAX    (output) DOUBLE PRECISION
!          The largest positive number for the machine, given by
!          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
!          value of BETA.
!
!  Further Details
!  ===============
!
!  The computation of  EPS  is based on a routine PARANOIA by
!  W. Kahan of the University of California at Berkeley.
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT, &
                         NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE, &
                         SIXTH, SMALL, THIRD, TWO, ZERO
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX, &
                         LRMIN, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
!
!        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
!        BETA, T, RND, EPS, EMIN and RMIN.
!
!        Throughout this routine  we use the function  DLAMC3  to ensure
!        that relevant values are stored  and not held in registers,  or
!        are not affected by optimizers.
!
!        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
!
!        Start to find EPS.
!
         B = LBETA
         A = B**( -LT )
         LEPS = A
!
!        Try some tricks to see whether or not this is the correct  EPS.
!
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS ) &
            B = LEPS
!
         LEPS = 1
!
!+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF
!+       END WHILE
!
         IF( A.LT.LEPS ) &
            LEPS = A
!
!        Computation of EPS complete.
!
!        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
!        Keep dividing  A by BETA until (gradual) underflow occurs. This
!        is detected when we cannot recover the previous A.
!
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
!
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
!            ( Non twos-complement machines, no gradual underflow;
!              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
!            ( Non twos-complement machines, with gradual underflow;
!              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
!            ( Twos-complement machines, no gradual underflow;
!              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND. &
                  ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
!            ( Twos-complement machines with gradual underflow;
!              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
!         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
!**
! Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
!**
!
!        Assume IEEE arithmetic if we found denormalised  numbers above,
!        or if arithmetic seems to round in the  IEEE style,  determined
!        in routine DLAMC1. A true IEEE machine should have both  things
!        true; however, faulty machines may have one or the other.
!
         IEEE = IEEE .OR. LIEEE1
!
!        Compute  RMIN by successive division by  BETA. We could compute
!        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
!        this computation.
!
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
!
!        Finally, call DLAMC5 to compute EMAX and RMAX.
!
         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
!
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
!
      RETURN
!
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-', &
            '  EMIN = ', I8, / &
            ' If, after inspection, the value EMIN looks', &
            ' acceptable please comment out ', &
            / ' the IF block as marked within the code of routine', &
            ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
!
!     End of DLAMC2
!
      END
!
!***********************************************************************
!
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
!     ..
!
!  Purpose
!  =======
!
!  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
!  the addition of  A  and  B ,  for use in situations where optimizers
!  might hold one of these in a register.
!
!  Arguments
!  =========
!
!  A, B    (input) DOUBLE PRECISION
!          The values A and B.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      DLAMC3 = A + B
!
      RETURN
!
!     End of DLAMC3
!
      END
!
!***********************************************************************
!
      SUBROUTINE DLAMC4( EMIN, START, BASE )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
!     ..
!
!  Purpose
!  =======
!
!  DLAMC4 is a service routine for DLAMC2.
!
!  Arguments
!  =========
!
!  EMIN    (output) EMIN
!          The minimum exponent before (gradual) underflow, computed by
!          setting A = START and dividing by BASE until the previous A
!          can not be recovered.
!
!  START   (input) DOUBLE PRECISION
!          The starting point for determining EMIN.
!
!  BASE    (input) INTEGER
!          The base of the machine.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
!     ..
!     .. Executable Statements ..
!
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
!+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
!    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND. &
          ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
!+    END WHILE
!
      RETURN
!
!     End of DLAMC4
!
      END
!
!***********************************************************************
!
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
!     ..
!
!  Purpose
!  =======
!
!  DLAMC5 attempts to compute RMAX, the largest machine floating-point
!  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
!  approximately to a power of 2.  It will fail on machines where this
!  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
!  EMAX = 28718).  It will also fail if the value supplied for EMIN is
!  too large (i.e. too close to zero), probably with overflow.
!
!  Arguments
!  =========
!
!  BETA    (input) INTEGER
!          The base of floating-point arithmetic.
!
!  P       (input) INTEGER
!          The number of base BETA digits in the mantissa of a
!          floating-point value.
!
!  EMIN    (input) INTEGER
!          The minimum exponent before (gradual) underflow.
!
!  IEEE    (input) LOGICAL
!          A logical flag specifying whether or not the arithmetic
!          system is thought to comply with the IEEE standard.
!
!  EMAX    (output) INTEGER
!          The largest exponent before overflow
!
!  RMAX    (output) DOUBLE PRECISION
!          The largest machine floating-point number.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
!
!     First compute LEXP and UEXP, two powers of 2 that bound
!     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
!     approximately to the bound that is closest to abs(EMIN).
!     (EMAX is the exponent of the required number RMAX).
!
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
!
!     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
!     than or equal to EMIN. EXBITS is the number of bits needed to
!     store the exponent.
!
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
!
!     EXPSUM is the exponent range, approximately equal to
!     EMAX - EMIN + 1 .
!
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
!
!     NBITS is the total number of bits needed to store a
!     floating-point number.
!
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
!
!        Either there are an odd number of bits used to store a
!        floating-point number, which is unlikely, or some bits are
!        not used in the representation of numbers, which is possible,
!        (e.g. Cray machines) or the mantissa has an implicit bit,
!        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
!        most likely. We have to assume the last alternative.
!        If this is true, then we need to reduce EMAX by one because
!        there must be some way of representing zero in an implicit-bit
!        system. On machines like Cray, we are reducing EMAX by one
!        unnecessarily.
!
         EMAX = EMAX - 1
      END IF
!
      IF( IEEE ) THEN
!
!        Assume we are on an IEEE machine which reserves one exponent
!        for infinity and NaN.
!
         EMAX = EMAX - 1
      END IF
!
!     Now create RMAX, the largest machine number, which should
!     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
!
!     First compute 1.0 - BETA**(-P), being careful that the
!     result is less than 1.0 .
!
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE ) &
            OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE ) &
         Y = OLDY
!
!     Now multiply by BETA**EMAX to get RMAX.
!
      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
!
      RMAX = Y
      RETURN
!
!     End of DLAMC5
!
      END
      DOUBLE PRECISION FUNCTION DLANSP( NORM, UPLO, N, AP, WORK )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AP( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLANSP  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real symmetric matrix A,  supplied in packed form.
!
!  Description
!  ===========
!
!  DLANSP returns the value
!
!     DLANSP = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in DLANSP as described
!          above.
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is supplied.
!          = 'U':  Upper triangular part of A is supplied
!          = 'L':  Lower triangular part of A is supplied
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, DLANSP is
!          set to zero.
!
!  AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          The upper or lower triangle of the symmetric matrix A, packed
!          columnwise in a linear array.  The j-th column of A is stored
!          in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
!          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
!          WORK is not referenced.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, K
      DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASSQ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            K = 1
            DO 20 J = 1, N
               DO 10 I = K, K + J - 1
                  VALUE = MAX( VALUE, ABS( AP( I ) ) )
   10          CONTINUE
               K = K + J
   20       CONTINUE
         ELSE
            K = 1
            DO 40 J = 1, N
               DO 30 I = K, K + N - J
                  VALUE = MAX( VALUE, ABS( AP( I ) ) )
   30          CONTINUE
               K = K + N - J + 1
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. &
               ( NORM.EQ.'1' ) ) THEN
!
!        Find normI(A) ( = norm1(A), since A is symmetric).
!
         VALUE = ZERO
         K = 1
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( AP( K ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
                  K = K + 1
   50          CONTINUE
               WORK( J ) = SUM + ABS( AP( K ) )
               K = K + 1
   60       CONTINUE
            DO 70 I = 1, N
               VALUE = MAX( VALUE, WORK( I ) )
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( AP( K ) )
               K = K + 1
               DO 90 I = J + 1, N
                  ABSA = ABS( AP( K ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
                  K = K + 1
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         K = 2
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 2, N
               CALL DLASSQ( J-1, AP( K ), 1, SCALE, SUM )
               K = K + J
  110       CONTINUE
         ELSE
            DO 120 J = 1, N - 1
               CALL DLASSQ( N-J, AP( K ), 1, SCALE, SUM )
               K = K + N - J + 1
  120       CONTINUE
         END IF
         SUM = 2*SUM
         K = 1
         DO 130 I = 1, N
            IF( AP( K ).NE.ZERO ) THEN
               ABSA = ABS( AP( K ) )
               IF( SCALE.LT.ABSA ) THEN
                  SUM = ONE + SUM*( SCALE / ABSA )**2
                  SCALE = ABSA
               ELSE
                  SUM = SUM + ( ABSA / SCALE )**2
               END IF
            END IF
            IF( LSAME( UPLO, 'U' ) ) THEN
               K = K + I + 1
            ELSE
               K = K + N - I + 1
            END IF
  130    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
!
      DLANSP = VALUE
      RETURN
!
!     End of DLANSP
!
      END
      DOUBLE PRECISION FUNCTION DLANST( NORM, N, D, E )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
!     ..
!
!  Purpose
!  =======
!
!  DLANST  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real symmetric tridiagonal matrix A.
!
!  Description
!  ===========
!
!  DLANST returns the value
!
!     DLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in DLANST as described
!          above.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, DLANST is
!          set to zero.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of A.
!
!  E       (input) DOUBLE PRECISION array, dimension (N-1)
!          The (n-1) sub-diagonal or super-diagonal elements of A.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   ANORM, SCALE, SUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASSQ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( N.LE.0 ) THEN
         ANORM = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         ANORM = ABS( D( N ) )
         DO 10 I = 1, N - 1
            ANORM = MAX( ANORM, ABS( D( I ) ) )
            ANORM = MAX( ANORM, ABS( E( I ) ) )
   10    CONTINUE
      ELSE IF( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' .OR. &
               LSAME( NORM, 'I' ) ) THEN
!
!        Find norm1(A).
!
         IF( N.EQ.1 ) THEN
            ANORM = ABS( D( 1 ) )
         ELSE
            ANORM = MAX( ABS( D( 1 ) )+ABS( E( 1 ) ), &
                    ABS( E( N-1 ) )+ABS( D( N ) ) )
            DO 20 I = 2, N - 1
               ANORM = MAX( ANORM, ABS( D( I ) )+ABS( E( I ) )+ &
                       ABS( E( I-1 ) ) )
   20       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         IF( N.GT.1 ) THEN
            CALL DLASSQ( N-1, E, 1, SCALE, SUM )
            SUM = 2*SUM
         END IF
         CALL DLASSQ( N, D, 1, SCALE, SUM )
         ANORM = SCALE*SQRT( SUM )
      END IF
!
      DLANST = ANORM
      RETURN
!
!     End of DLANST
!
      END
      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
!     ..
!
!  Purpose
!  =======
!
!  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!  overflow.
!
!  Arguments
!  =========
!
!  X       (input) DOUBLE PRECISION
!  Y       (input) DOUBLE PRECISION
!          X and Y specify the values x and y.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, Z
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         DLAPY2 = W
      ELSE
         DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
!
!     End of DLAPY2
!
      END
      DOUBLE PRECISION FUNCTION DLAPY3( X, Y, Z )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y, Z
!     ..
!
!  Purpose
!  =======
!
!  DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
!  unnecessary overflow.
!
!  Arguments
!  =========
!
!  X       (input) DOUBLE PRECISION
!  Y       (input) DOUBLE PRECISION
!  Z       (input) DOUBLE PRECISION
!          X, Y and Z specify the values x, y and z.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, ZABS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      XABS = ABS( X )
      YABS = ABS( Y )
      ZABS = ABS( Z )
      W = MAX( XABS, YABS, ZABS )
      IF( W.EQ.ZERO ) THEN
         DLAPY3 = ZERO
      ELSE
         DLAPY3 = W*SQRT( ( XABS / W )**2+( YABS / W )**2+ &
                  ( ZABS / W )**2 )
      END IF
      RETURN
!
!     End of DLAPY3
!
      END
      SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLARF applies a real elementary reflector H to a real m by n matrix
!  C, from either the left or the right. H is represented in the form
!
!        H = I - tau * v * v'
!
!  where tau is a real scalar and v is a real vector.
!
!  If tau = 0, then H is taken to be the unit matrix.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': form  H * C
!          = 'R': form  C * H
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  V       (input) DOUBLE PRECISION array, dimension
!                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!          The vector v in the representation of H. V is not used if
!          TAU = 0.
!
!  INCV    (input) INTEGER
!          The increment between elements of v. INCV <> 0.
!
!  TAU     (input) DOUBLE PRECISION
!          The value tau in the representation of H.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!          or C * H if SIDE = 'R'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                         (N) if SIDE = 'L'
!                      or (M) if SIDE = 'R'
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  H * C
!
         IF( TAU.NE.ZERO ) THEN
!
!           w := C' * v
!
            CALL DGEMV( 'Transpose', M, N, ONE, C, LDC, V, INCV, ZERO, &
                        WORK, 1 )
!
!           C := C - v * w'
!
            CALL DGER( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
!
!        Form  C * H
!
         IF( TAU.NE.ZERO ) THEN
!
!           w := C * v
!
            CALL DGEMV( 'No transpose', M, N, ONE, C, LDC, V, INCV, &
                        ZERO, WORK, 1 )
!
!           C := C - w * v'
!
            CALL DGER( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
!
!     End of DLARF
!
      END
      SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ALPHA, TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
!     ..
!
!  Purpose
!  =======
!
!  DLARFG generates a real elementary reflector H of order n, such
!  that
!
!        H * ( alpha ) = ( beta ),   H' * H = I.
!            (   x   )   (   0  )
!
!  where alpha and beta are scalars, and x is an (n-1)-element real
!  vector. H is represented in the form
!
!        H = I - tau * ( 1 ) * ( 1 v' ) ,
!                      ( v )
!
!  where tau is a real scalar and v is a real (n-1)-element
!  vector.
!
!  If the elements of x are all zero, then tau = 0 and H is taken to be
!  the unit matrix.
!
!  Otherwise  1 <= tau <= 2.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the elementary reflector.
!
!  ALPHA   (input/output) DOUBLE PRECISION
!          On entry, the value alpha.
!          On exit, it is overwritten with the value beta.
!
!  X       (input/output) DOUBLE PRECISION array, dimension
!                         (1+(N-2)*abs(INCX))
!          On entry, the vector x.
!          On exit, it is overwritten with the vector v.
!
!  INCX    (input) INTEGER
!          The increment between elements of X. INCX > 0.
!
!  TAU     (output) DOUBLE PRECISION
!          The value tau.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
      EXTERNAL           DLAMCH, DLAPY2, DNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL
!     ..
!     .. Executable Statements ..
!
      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF
!
      XNORM = DNRM2( N-1, X, INCX )
!
      IF( XNORM.EQ.ZERO ) THEN
!
!        H  =  I
!
         TAU = ZERO
      ELSE
!
!        general case
!
         BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         IF( ABS( BETA ).LT.SAFMIN ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
            RSAFMN = ONE / SAFMIN
            KNT = 0
   10       CONTINUE
            KNT = KNT + 1
            CALL DSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN ) &
               GO TO 10
!
!           New BETA is at most 1, at least SAFMIN
!
            XNORM = DNRM2( N-1, X, INCX )
            BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
            TAU = ( BETA-ALPHA ) / BETA
            CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
!
!           If ALPHA is subnormal, it may lose relative accuracy
!
            ALPHA = BETA
            DO 20 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   20       CONTINUE
         ELSE
            TAU = ( BETA-ALPHA ) / BETA
            CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
            ALPHA = BETA
         END IF
      END IF
!
      RETURN
!
!     End of DLARFG
!
      END
      SUBROUTINE DLARNV( IDIST, ISEED, N, X )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994 
!
!     .. Scalar Arguments ..
      INTEGER            IDIST, N
!     ..
!     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   X( * )
!     ..
!
!  Purpose
!  =======
!
!  DLARNV returns a vector of n random real numbers from a uniform or
!  normal distribution.
!
!  Arguments
!  =========
!
!  IDIST   (input) INTEGER
!          Specifies the distribution of the random numbers:
!          = 1:  uniform (0,1)
!          = 2:  uniform (-1,1)
!          = 3:  normal (0,1)
!
!  ISEED   (input/output) INTEGER array, dimension (4)
!          On entry, the seed of the random number generator; the array
!          elements must be between 0 and 4095, and ISEED(4) must be
!          odd.
!          On exit, the seed is updated.
!
!  N       (input) INTEGER
!          The number of random numbers to be generated.
!
!  X       (output) DOUBLE PRECISION array, dimension (N)
!          The generated random numbers.
!
!  Further Details
!  ===============
!
!  This routine calls the auxiliary routine DLARUV to generate random
!  real numbers from a uniform (0,1) distribution, in batches of up to
!  128 using vectorisable code. The Box-Muller method is used to
!  transform numbers from a uniform to a normal distribution.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, TWO
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0 )
      INTEGER            LV
      PARAMETER          ( LV = 128 )
      DOUBLE PRECISION   TWOPI
      PARAMETER          ( TWOPI = 6.2831853071795864769252867663D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IL, IL2, IV
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   U( LV )
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          COS, LOG, MIN, SQRT
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARUV
!     ..
!     .. Executable Statements ..
!
      DO 40 IV = 1, N, LV / 2
         IL = MIN( LV / 2, N-IV+1 )
         IF( IDIST.EQ.3 ) THEN
            IL2 = 2*IL
         ELSE
            IL2 = IL
         END IF
!
!        Call DLARUV to generate IL2 numbers from a uniform (0,1)
!        distribution (IL2 <= LV)
!
         CALL DLARUV( ISEED, IL2, U )
!
         IF( IDIST.EQ.1 ) THEN
!
!           Copy generated numbers
!
            DO 10 I = 1, IL
               X( IV+I-1 ) = U( I )
   10       CONTINUE
         ELSE IF( IDIST.EQ.2 ) THEN
!
!           Convert generated numbers to uniform (-1,1) distribution
!
            DO 20 I = 1, IL
               X( IV+I-1 ) = TWO*U( I ) - ONE
   20       CONTINUE
         ELSE IF( IDIST.EQ.3 ) THEN
!
!           Convert generated numbers to normal (0,1) distribution
!
            DO 30 I = 1, IL
               X( IV+I-1 ) = SQRT( -TWO*LOG( U( 2*I-1 ) ) )* &
                             COS( TWOPI*U( 2*I ) )
   30       CONTINUE
         END IF
   40 CONTINUE
      RETURN
!
!     End of DLARNV
!
      END
      SUBROUTINE DLARTG( F, G, CS, SN, R )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   CS, F, G, R, SN
!     ..
!
!  Purpose
!  =======
!
!  DLARTG generate a plane rotation so that
!
!     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
!     [ -SN  CS  ]     [ G ]     [ 0 ]
!
!  This is a slower, more accurate version of the BLAS1 routine DROTG,
!  with the following other differences:
!     F and G are unchanged on return.
!     If G=0, then CS=1 and SN=0.
!     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
!        floating point operations (saves work in DBDSQR when
!        there are zeros on the diagonal).
!
!  If F exceeds G in magnitude, CS will be positive.
!
!  Arguments
!  =========
!
!  F       (input) DOUBLE PRECISION
!          The first component of vector to be rotated.
!
!  G       (input) DOUBLE PRECISION
!          The second component of vector to be rotated.
!
!  CS      (output) DOUBLE PRECISION
!          The cosine of the rotation.
!
!  SN      (output) DOUBLE PRECISION
!          The sine of the rotation.
!
!  R       (output) DOUBLE PRECISION
!          The nonzero component of the rotated vector.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            FIRST
      INTEGER            COUNT, I
      DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, SQRT
!     ..
!     .. Save statement ..
      SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         SAFMIN = DLAMCH( 'S' )
         EPS = DLAMCH( 'E' )
         SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) / &
                  LOG( DLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
      END IF
      IF( G.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         F1 = F
         G1 = G
         SCALE = MAX( ABS( F1 ), ABS( G1 ) )
         IF( SCALE.GE.SAFMX2 ) THEN
            COUNT = 0
   10       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.GE.SAFMX2 ) &
               GO TO 10
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 20 I = 1, COUNT
               R = R*SAFMX2
   20       CONTINUE
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
            COUNT = 0
   30       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.LE.SAFMN2 ) &
               GO TO 30
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 40 I = 1, COUNT
               R = R*SAFMN2
   40       CONTINUE
         ELSE
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         END IF
         IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
            CS = -CS
            SN = -SN
            R = -R
         END IF
      END IF
      RETURN
!
!     End of DLARTG
!
      END
      SUBROUTINE DLARUV( ISEED, N, X )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            N
!     ..
!     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   X( N )
!     ..
!
!  Purpose
!  =======
!
!  DLARUV returns a vector of n random real numbers from a uniform (0,1)
!  distribution (n <= 128).
!
!  This is an auxiliary routine called by DLARNV and ZLARNV.
!
!  Arguments
!  =========
!
!  ISEED   (input/output) INTEGER array, dimension (4)
!          On entry, the seed of the random number generator; the array
!          elements must be between 0 and 4095, and ISEED(4) must be
!          odd.
!          On exit, the seed is updated.
!
!  N       (input) INTEGER
!          The number of random numbers to be generated. N <= 128.
!
!  X       (output) DOUBLE PRECISION array, dimension (N)
!          The generated random numbers.
!
!  Further Details
!  ===============
!
!  This routine uses a multiplicative congruential method with modulus
!  2**48 and multiplier 33952834046453 (see G.S.Fishman,
!  'Multiplicative congruential random number generators with modulus
!  2**b: an exhaustive analysis for b = 32 and a partial analysis for
!  b = 48', Math. Comp. 189, pp 331-344, 1990).
!
!  48-bit integers are stored in 4 integer array elements with 12 bits
!  per element. Hence the routine is portable across machines with
!  integers of 32 bits or more.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      INTEGER            LV, IPW2
      DOUBLE PRECISION   R
      PARAMETER          ( LV = 128, IPW2 = 4096, R = ONE / IPW2 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I1, I2, I3, I4, IT1, IT2, IT3, IT4, J
!     ..
!     .. Local Arrays ..
      INTEGER            MM( LV, 4 )
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MIN, MOD
!     ..
!     .. Data statements ..
      DATA               ( MM( 1, J ), J = 1, 4 ) / 494, 322, 2508, &
                         2549 /
      DATA               ( MM( 2, J ), J = 1, 4 ) / 2637, 789, 3754, &
                         1145 /
      DATA               ( MM( 3, J ), J = 1, 4 ) / 255, 1440, 1766, &
                         2253 /
      DATA               ( MM( 4, J ), J = 1, 4 ) / 2008, 752, 3572, &
                         305 /
      DATA               ( MM( 5, J ), J = 1, 4 ) / 1253, 2859, 2893, &
                         3301 /
      DATA               ( MM( 6, J ), J = 1, 4 ) / 3344, 123, 307, &
                         1065 /
      DATA               ( MM( 7, J ), J = 1, 4 ) / 4084, 1848, 1297, &
                         3133 /
      DATA               ( MM( 8, J ), J = 1, 4 ) / 1739, 643, 3966, &
                         2913 /
      DATA               ( MM( 9, J ), J = 1, 4 ) / 3143, 2405, 758, &
                         3285 /
      DATA               ( MM( 10, J ), J = 1, 4 ) / 3468, 2638, 2598, &
                         1241 /
      DATA               ( MM( 11, J ), J = 1, 4 ) / 688, 2344, 3406, &
                         1197 /
      DATA               ( MM( 12, J ), J = 1, 4 ) / 1657, 46, 2922, &
                         3729 /
      DATA               ( MM( 13, J ), J = 1, 4 ) / 1238, 3814, 1038, &
                         2501 /
      DATA               ( MM( 14, J ), J = 1, 4 ) / 3166, 913, 2934, &
                         1673 /
      DATA               ( MM( 15, J ), J = 1, 4 ) / 1292, 3649, 2091, &
                         541 /
      DATA               ( MM( 16, J ), J = 1, 4 ) / 3422, 339, 2451, &
                         2753 /
      DATA               ( MM( 17, J ), J = 1, 4 ) / 1270, 3808, 1580, &
                         949 /
      DATA               ( MM( 18, J ), J = 1, 4 ) / 2016, 822, 1958, &
                         2361 /
      DATA               ( MM( 19, J ), J = 1, 4 ) / 154, 2832, 2055, &
                         1165 /
      DATA               ( MM( 20, J ), J = 1, 4 ) / 2862, 3078, 1507, &
                         4081 /
      DATA               ( MM( 21, J ), J = 1, 4 ) / 697, 3633, 1078, &
                         2725 /
      DATA               ( MM( 22, J ), J = 1, 4 ) / 1706, 2970, 3273, &
                         3305 /
      DATA               ( MM( 23, J ), J = 1, 4 ) / 491, 637, 17, &
                         3069 /
      DATA               ( MM( 24, J ), J = 1, 4 ) / 931, 2249, 854, &
                         3617 /
      DATA               ( MM( 25, J ), J = 1, 4 ) / 1444, 2081, 2916, &
                         3733 /
      DATA               ( MM( 26, J ), J = 1, 4 ) / 444, 4019, 3971, &
                         409 /
      DATA               ( MM( 27, J ), J = 1, 4 ) / 3577, 1478, 2889, &
                         2157 /
      DATA               ( MM( 28, J ), J = 1, 4 ) / 3944, 242, 3831, &
                         1361 /
      DATA               ( MM( 29, J ), J = 1, 4 ) / 2184, 481, 2621, &
                         3973 /
      DATA               ( MM( 30, J ), J = 1, 4 ) / 1661, 2075, 1541, &
                         1865 /
      DATA               ( MM( 31, J ), J = 1, 4 ) / 3482, 4058, 893, &
                         2525 /
      DATA               ( MM( 32, J ), J = 1, 4 ) / 657, 622, 736, &
                         1409 /
      DATA               ( MM( 33, J ), J = 1, 4 ) / 3023, 3376, 3992, &
                         3445 /
      DATA               ( MM( 34, J ), J = 1, 4 ) / 3618, 812, 787, &
                         3577 /
      DATA               ( MM( 35, J ), J = 1, 4 ) / 1267, 234, 2125, &
                         77 /
      DATA               ( MM( 36, J ), J = 1, 4 ) / 1828, 641, 2364, &
                         3761 /
      DATA               ( MM( 37, J ), J = 1, 4 ) / 164, 4005, 2460, &
                         2149 /
      DATA               ( MM( 38, J ), J = 1, 4 ) / 3798, 1122, 257, &
                         1449 /
      DATA               ( MM( 39, J ), J = 1, 4 ) / 3087, 3135, 1574, &
                         3005 /
      DATA               ( MM( 40, J ), J = 1, 4 ) / 2400, 2640, 3912, &
                         225 /
      DATA               ( MM( 41, J ), J = 1, 4 ) / 2870, 2302, 1216, &
                         85 /
      DATA               ( MM( 42, J ), J = 1, 4 ) / 3876, 40, 3248, &
                         3673 /
      DATA               ( MM( 43, J ), J = 1, 4 ) / 1905, 1832, 3401, &
                         3117 /
      DATA               ( MM( 44, J ), J = 1, 4 ) / 1593, 2247, 2124, &
                         3089 /
      DATA               ( MM( 45, J ), J = 1, 4 ) / 1797, 2034, 2762, &
                         1349 /
      DATA               ( MM( 46, J ), J = 1, 4 ) / 1234, 2637, 149, &
                         2057 /
      DATA               ( MM( 47, J ), J = 1, 4 ) / 3460, 1287, 2245, &
                         413 /
      DATA               ( MM( 48, J ), J = 1, 4 ) / 328, 1691, 166, &
                         65 /
      DATA               ( MM( 49, J ), J = 1, 4 ) / 2861, 496, 466, &
                         1845 /
      DATA               ( MM( 50, J ), J = 1, 4 ) / 1950, 1597, 4018, &
                         697 /
      DATA               ( MM( 51, J ), J = 1, 4 ) / 617, 2394, 1399, &
                         3085 /
      DATA               ( MM( 52, J ), J = 1, 4 ) / 2070, 2584, 190, &
                         3441 /
      DATA               ( MM( 53, J ), J = 1, 4 ) / 3331, 1843, 2879, &
                         1573 /
      DATA               ( MM( 54, J ), J = 1, 4 ) / 769, 336, 153, &
                         3689 /
      DATA               ( MM( 55, J ), J = 1, 4 ) / 1558, 1472, 2320, &
                         2941 /
      DATA               ( MM( 56, J ), J = 1, 4 ) / 2412, 2407, 18, &
                         929 /
      DATA               ( MM( 57, J ), J = 1, 4 ) / 2800, 433, 712, &
                         533 /
      DATA               ( MM( 58, J ), J = 1, 4 ) / 189, 2096, 2159, &
                         2841 /
      DATA               ( MM( 59, J ), J = 1, 4 ) / 287, 1761, 2318, &
                         4077 /
      DATA               ( MM( 60, J ), J = 1, 4 ) / 2045, 2810, 2091, &
                         721 /
      DATA               ( MM( 61, J ), J = 1, 4 ) / 1227, 566, 3443, &
                         2821 /
      DATA               ( MM( 62, J ), J = 1, 4 ) / 2838, 442, 1510, &
                         2249 /
      DATA               ( MM( 63, J ), J = 1, 4 ) / 209, 41, 449, &
                         2397 /
      DATA               ( MM( 64, J ), J = 1, 4 ) / 2770, 1238, 1956, &
                         2817 /
      DATA               ( MM( 65, J ), J = 1, 4 ) / 3654, 1086, 2201, &
                         245 /
      DATA               ( MM( 66, J ), J = 1, 4 ) / 3993, 603, 3137, &
                         1913 /
      DATA               ( MM( 67, J ), J = 1, 4 ) / 192, 840, 3399, &
                         1997 /
      DATA               ( MM( 68, J ), J = 1, 4 ) / 2253, 3168, 1321, &
                         3121 /
      DATA               ( MM( 69, J ), J = 1, 4 ) / 3491, 1499, 2271, &
                         997 /
      DATA               ( MM( 70, J ), J = 1, 4 ) / 2889, 1084, 3667, &
                         1833 /
      DATA               ( MM( 71, J ), J = 1, 4 ) / 2857, 3438, 2703, &
                         2877 /
      DATA               ( MM( 72, J ), J = 1, 4 ) / 2094, 2408, 629, &
                         1633 /
      DATA               ( MM( 73, J ), J = 1, 4 ) / 1818, 1589, 2365, &
                         981 /
      DATA               ( MM( 74, J ), J = 1, 4 ) / 688, 2391, 2431, &
                         2009 /
      DATA               ( MM( 75, J ), J = 1, 4 ) / 1407, 288, 1113, &
                         941 /
      DATA               ( MM( 76, J ), J = 1, 4 ) / 634, 26, 3922, &
                         2449 /
      DATA               ( MM( 77, J ), J = 1, 4 ) / 3231, 512, 2554, &
                         197 /
      DATA               ( MM( 78, J ), J = 1, 4 ) / 815, 1456, 184, &
                         2441 /
      DATA               ( MM( 79, J ), J = 1, 4 ) / 3524, 171, 2099, &
                         285 /
      DATA               ( MM( 80, J ), J = 1, 4 ) / 1914, 1677, 3228, &
                         1473 /
      DATA               ( MM( 81, J ), J = 1, 4 ) / 516, 2657, 4012, &
                         2741 /
      DATA               ( MM( 82, J ), J = 1, 4 ) / 164, 2270, 1921, &
                         3129 /
      DATA               ( MM( 83, J ), J = 1, 4 ) / 303, 2587, 3452, &
                         909 /
      DATA               ( MM( 84, J ), J = 1, 4 ) / 2144, 2961, 3901, &
                         2801 /
      DATA               ( MM( 85, J ), J = 1, 4 ) / 3480, 1970, 572, &
                         421 /
      DATA               ( MM( 86, J ), J = 1, 4 ) / 119, 1817, 3309, &
                         4073 /
      DATA               ( MM( 87, J ), J = 1, 4 ) / 3357, 676, 3171, &
                         2813 /
      DATA               ( MM( 88, J ), J = 1, 4 ) / 837, 1410, 817, &
                         2337 /
      DATA               ( MM( 89, J ), J = 1, 4 ) / 2826, 3723, 3039, &
                         1429 /
      DATA               ( MM( 90, J ), J = 1, 4 ) / 2332, 2803, 1696, &
                         1177 /
      DATA               ( MM( 91, J ), J = 1, 4 ) / 2089, 3185, 1256, &
                         1901 /
      DATA               ( MM( 92, J ), J = 1, 4 ) / 3780, 184, 3715, &
                         81 /
      DATA               ( MM( 93, J ), J = 1, 4 ) / 1700, 663, 2077, &
                         1669 /
      DATA               ( MM( 94, J ), J = 1, 4 ) / 3712, 499, 3019, &
                         2633 /
      DATA               ( MM( 95, J ), J = 1, 4 ) / 150, 3784, 1497, &
                         2269 /
      DATA               ( MM( 96, J ), J = 1, 4 ) / 2000, 1631, 1101, &
                         129 /
      DATA               ( MM( 97, J ), J = 1, 4 ) / 3375, 1925, 717, &
                         1141 /
      DATA               ( MM( 98, J ), J = 1, 4 ) / 1621, 3912, 51, &
                         249 /
      DATA               ( MM( 99, J ), J = 1, 4 ) / 3090, 1398, 981, &
                         3917 /
      DATA               ( MM( 100, J ), J = 1, 4 ) / 3765, 1349, 1978, &
                         2481 /
      DATA               ( MM( 101, J ), J = 1, 4 ) / 1149, 1441, 1813, &
                         3941 /
      DATA               ( MM( 102, J ), J = 1, 4 ) / 3146, 2224, 3881, &
                         2217 /
      DATA               ( MM( 103, J ), J = 1, 4 ) / 33, 2411, 76, &
                         2749 /
      DATA               ( MM( 104, J ), J = 1, 4 ) / 3082, 1907, 3846, &
                         3041 /
      DATA               ( MM( 105, J ), J = 1, 4 ) / 2741, 3192, 3694, &
                         1877 /
      DATA               ( MM( 106, J ), J = 1, 4 ) / 359, 2786, 1682, &
                         345 /
      DATA               ( MM( 107, J ), J = 1, 4 ) / 3316, 382, 124, &
                         2861 /
      DATA               ( MM( 108, J ), J = 1, 4 ) / 1749, 37, 1660, &
                         1809 /
      DATA               ( MM( 109, J ), J = 1, 4 ) / 185, 759, 3997, &
                         3141 /
      DATA               ( MM( 110, J ), J = 1, 4 ) / 2784, 2948, 479, &
                         2825 /
      DATA               ( MM( 111, J ), J = 1, 4 ) / 2202, 1862, 1141, &
                         157 /
      DATA               ( MM( 112, J ), J = 1, 4 ) / 2199, 3802, 886, &
                         2881 /
      DATA               ( MM( 113, J ), J = 1, 4 ) / 1364, 2423, 3514, &
                         3637 /
      DATA               ( MM( 114, J ), J = 1, 4 ) / 1244, 2051, 1301, &
                         1465 /
      DATA               ( MM( 115, J ), J = 1, 4 ) / 2020, 2295, 3604, &
                         2829 /
      DATA               ( MM( 116, J ), J = 1, 4 ) / 3160, 1332, 1888, &
                         2161 /
      DATA               ( MM( 117, J ), J = 1, 4 ) / 2785, 1832, 1836, &
                         3365 /
      DATA               ( MM( 118, J ), J = 1, 4 ) / 2772, 2405, 1990, &
                         361 /
      DATA               ( MM( 119, J ), J = 1, 4 ) / 1217, 3638, 2058, &
                         2685 /
      DATA               ( MM( 120, J ), J = 1, 4 ) / 1822, 3661, 692, &
                         3745 /
      DATA               ( MM( 121, J ), J = 1, 4 ) / 1245, 327, 1194, &
                         2325 /
      DATA               ( MM( 122, J ), J = 1, 4 ) / 2252, 3660, 20, &
                         3609 /
      DATA               ( MM( 123, J ), J = 1, 4 ) / 3904, 716, 3285, &
                         3821 /
      DATA               ( MM( 124, J ), J = 1, 4 ) / 2774, 1842, 2046, &
                         3537 /
      DATA               ( MM( 125, J ), J = 1, 4 ) / 997, 3987, 2107, &
                         517 /
      DATA               ( MM( 126, J ), J = 1, 4 ) / 2573, 1368, 3508, &
                         3017 /
      DATA               ( MM( 127, J ), J = 1, 4 ) / 1148, 1848, 3525, &
                         2141 /
      DATA               ( MM( 128, J ), J = 1, 4 ) / 545, 2366, 3801, &
                         1537 /
!     ..
!     .. Executable Statements ..
!
      I1 = ISEED( 1 )
      I2 = ISEED( 2 )
      I3 = ISEED( 3 )
      I4 = ISEED( 4 )
!
      DO 10 I = 1, MIN( N, LV )
!
!        Multiply the seed by i-th power of the multiplier modulo 2**48
!
         IT4 = I4*MM( I, 4 )
         IT3 = IT4 / IPW2
         IT4 = IT4 - IPW2*IT3
         IT3 = IT3 + I3*MM( I, 4 ) + I4*MM( I, 3 )
         IT2 = IT3 / IPW2
         IT3 = IT3 - IPW2*IT2
         IT2 = IT2 + I2*MM( I, 4 ) + I3*MM( I, 3 ) + I4*MM( I, 2 )
         IT1 = IT2 / IPW2
         IT2 = IT2 - IPW2*IT1
         IT1 = IT1 + I1*MM( I, 4 ) + I2*MM( I, 3 ) + I3*MM( I, 2 ) + &
               I4*MM( I, 1 )
         IT1 = MOD( IT1, IPW2 )
!
!        Convert 48-bit integer to a real number in the interval (0,1)
!
         X( I ) = R*( DBLE( IT1 )+R*( DBLE( IT2 )+R*( DBLE( IT3 )+R* &
                  DBLE( IT4 ) ) ) )
   10 CONTINUE
!
!     Return final value of seed
!
      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4
      RETURN
!
!     End of DLARUV
!
      END
      SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASCL multiplies the M by N real matrix A by the real scalar
!  CTO/CFROM.  This is done without over/underflow as long as the final
!  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
!  A may be full, upper triangular, lower triangular, upper Hessenberg,
!  or banded.
!
!  Arguments
!  =========
!
!  TYPE    (input) CHARACTER*1
!          TYPE indices the storage type of the input matrix.
!          = 'G':  A is a full matrix.
!          = 'L':  A is a lower triangular matrix.
!          = 'U':  A is an upper triangular matrix.
!          = 'H':  A is an upper Hessenberg matrix.
!          = 'B':  A is a symmetric band matrix with lower bandwidth KL
!                  and upper bandwidth KU and with the only the lower
!                  half stored.
!          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
!                  and upper bandwidth KU and with the only the upper
!                  half stored.
!          = 'Z':  A is a band matrix with lower bandwidth KL and upper
!                  bandwidth KU.
!
!  KL      (input) INTEGER
!          The lower bandwidth of A.  Referenced only if TYPE = 'B',
!          'Q' or 'Z'.
!
!  KU      (input) INTEGER
!          The upper bandwidth of A.  Referenced only if TYPE = 'B',
!          'Q' or 'Z'.
!
!  CFROM   (input) DOUBLE PRECISION
!  CTO     (input) DOUBLE PRECISION
!          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
!          without over/underflow if the final result CTO*A(I,J)/CFROM
!          can be represented without over/underflow.  CFROM must be
!          nonzero.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,M)
!          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
!          storage type.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  INFO    (output) INTEGER
!          0  - successful exit
!          <0 - if INFO = -i, the i-th argument had an illegal value.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
!
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
!
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR. &
               ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR. &
                  ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) ) &
                   THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR. &
                  ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR. &
                  ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASCL', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. M.EQ.0 ) &
         RETURN
!
!     Get machine parameters
!
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
!
      CFROMC = CFROM
      CTOC = CTO
!
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      CTO1 = CTOC / BIGNUM
      IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
         MUL = SMLNUM
         DONE = .FALSE.
         CFROMC = CFROM1
      ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
         MUL = BIGNUM
         DONE = .FALSE.
         CTOC = CTO1
      ELSE
         MUL = CTOC / CFROMC
         DONE = .TRUE.
      END IF
!
      IF( ITYPE.EQ.0 ) THEN
!
!        Full matrix
!
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
!
      ELSE IF( ITYPE.EQ.1 ) THEN
!
!        Lower triangular matrix
!
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
!
      ELSE IF( ITYPE.EQ.2 ) THEN
!
!        Upper triangular matrix
!
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
!
      ELSE IF( ITYPE.EQ.3 ) THEN
!
!        Upper Hessenberg matrix
!
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
!
      ELSE IF( ITYPE.EQ.4 ) THEN
!
!        Lower half of a symmetric band matrix
!
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
!
      ELSE IF( ITYPE.EQ.5 ) THEN
!
!        Upper half of a symmetric band matrix
!
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
!
      ELSE IF( ITYPE.EQ.6 ) THEN
!
!        Band matrix
!
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
!
      END IF
!
      IF( .NOT.DONE ) &
         GO TO 10
!
      RETURN
!
!     End of DLASCL
!
      END
      SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      DOUBLE PRECISION   ALPHA, BETA
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASET initializes an m-by-n matrix A to BETA on the diagonal and
!  ALPHA on the offdiagonals.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be set.
!          = 'U':      Upper triangular part is set; the strictly lower
!                      triangular part of A is not changed.
!          = 'L':      Lower triangular part is set; the strictly upper
!                      triangular part of A is not changed.
!          Otherwise:  All of the matrix A is set.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  ALPHA   (input) DOUBLE PRECISION
!          The constant to which the offdiagonal elements are to be set.
!
!  BETA    (input) DOUBLE PRECISION
!          The constant to which the diagonal elements are to be set.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On exit, the leading m-by-n submatrix of A is set as follows:
!
!          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
!          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
!          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
!
!          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Set the strictly upper triangular or trapezoidal part of the
!        array to ALPHA.
!
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
!
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
!
!        Set the strictly lower triangular or trapezoidal part of the
!        array to ALPHA.
!
         DO 40 J = 1, MIN( M, N )
            DO 30 I = J + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
   40    CONTINUE
!
      ELSE
!
!        Set the leading m-by-n submatrix to ALPHA.
!
         DO 60 J = 1, N
            DO 50 I = 1, M
               A( I, J ) = ALPHA
   50       CONTINUE
   60    CONTINUE
      END IF
!
!     Set the first min(M,N) diagonal elements to BETA.
!
      DO 70 I = 1, MIN( M, N )
         A( I, I ) = BETA
   70 CONTINUE
!
      RETURN
!
!     End of DLASET
!
      END
      SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASR   performs the transformation
!
!     A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )
!
!     A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )
!
!  where A is an m by n real matrix and P is an orthogonal matrix,
!  consisting of a sequence of plane rotations determined by the
!  parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'
!  and z = n when SIDE = 'R' or 'r' ):
!
!  When  DIRECT = 'F' or 'f'  ( Forward sequence ) then
!
!     P = P( z - 1 )*...*P( 2 )*P( 1 ),
!
!  and when DIRECT = 'B' or 'b'  ( Backward sequence ) then
!
!     P = P( 1 )*P( 2 )*...*P( z - 1 ),
!
!  where  P( k ) is a plane rotation matrix for the following planes:
!
!     when  PIVOT = 'V' or 'v'  ( Variable pivot ),
!        the plane ( k, k + 1 )
!
!     when  PIVOT = 'T' or 't'  ( Top pivot ),
!        the plane ( 1, k + 1 )
!
!     when  PIVOT = 'B' or 'b'  ( Bottom pivot ),
!        the plane ( k, z )
!
!  c( k ) and s( k )  must contain the  cosine and sine that define the
!  matrix  P( k ).  The two by two plane rotation part of the matrix
!  P( k ), R( k ), is assumed to be of the form
!
!     R( k ) = (  c( k )  s( k ) ).
!              ( -s( k )  c( k ) )
!
!  This version vectorises across rows of the array A when SIDE = 'L'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          Specifies whether the plane rotation matrix P is applied to
!          A on the left or the right.
!          = 'L':  Left, compute A := P*A
!          = 'R':  Right, compute A:= A*P'
!
!  DIRECT  (input) CHARACTER*1
!          Specifies whether P is a forward or backward sequence of
!          plane rotations.
!          = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )
!          = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )
!
!  PIVOT   (input) CHARACTER*1
!          Specifies the plane for which P(k) is a plane rotation
!          matrix.
!          = 'V':  Variable pivot, the plane (k,k+1)
!          = 'T':  Top pivot, the plane (1,k+1)
!          = 'B':  Bottom pivot, the plane (k,z)
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  If m <= 1, an immediate
!          return is effected.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  If n <= 1, an
!          immediate return is effected.
!
!  C, S    (input) DOUBLE PRECISION arrays, dimension
!                  (M-1) if SIDE = 'L'
!                  (N-1) if SIDE = 'R'
!          c(k) and s(k) contain the cosine and sine that define the
!          matrix P(k).  The two by two plane rotation part of the
!          matrix P(k), R(k), is assumed to be of the form
!          R( k ) = (  c( k )  s( k ) ).
!                   ( -s( k )  c( k ) )
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          The m by n matrix A.  On exit, A is overwritten by P*A if
!          SIDE = 'R' or by A*P' if SIDE = 'L'.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   CTEMP, STEMP, TEMP
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT, &
               'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) ) &
                THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASR ', INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) ) &
         RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  P * A
!
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!        Form A * P'
!
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of DLASR
!
      END
      SUBROUTINE DLASRT( ID, N, D, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * )
!     ..
!
!  Purpose
!  =======
!
!  Sort the numbers in D in increasing order (if ID = 'I') or
!  in decreasing order (if ID = 'D' ).
!
!  Use Quick Sort, reverting to Insertion sort on arrays of
!  size <= 20. Dimension of STACK limits N to about 2**32.
!
!  Arguments
!  =========
!
!  ID      (input) CHARACTER*1
!          = 'I': sort D in increasing order;
!          = 'D': sort D in decreasing order.
!
!  N       (input) INTEGER
!          The length of the array D.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the array to be sorted.
!          On exit, D has been sorted into increasing order
!          (D(1) <= ... <= D(N) ) or into decreasing order
!          (D(1) >= ... >= D(N) ), depending on ID.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
!     ..
!     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, J, START, STKPNT
      DOUBLE PRECISION   D1, D2, D3, DMNMX, TMP
!     ..
!     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input paramters.
!
      INFO = 0
      DIR = -1
      IF( LSAME( ID, 'D' ) ) THEN
         DIR = 0
      ELSE IF( LSAME( ID, 'I' ) ) THEN
         DIR = 1
      END IF
      IF( DIR.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASRT', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.1 ) &
         RETURN
!
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
!
!        Do Insertion sort on D( START:ENDD )
!
         IF( DIR.EQ.0 ) THEN
!
!           Sort into decreasing order
!
            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  IF( D( J ).GT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 30
                  END IF
   20          CONTINUE
   30       CONTINUE
!
         ELSE
!
!           Sort into increasing order
!
            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  IF( D( J ).LT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 50
                  END IF
   40          CONTINUE
   50       CONTINUE
!
         END IF
!
      ELSE IF( ENDD-START.GT.SELECT ) THEN
!
!        Partition D( START:ENDD ) and stack parts, largest one first
!
!        Choose partition entry as median of 3
!
         D1 = D( START )
         D2 = D( ENDD )
         I = ( START+ENDD ) / 2
         D3 = D( I )
         IF( D1.LT.D2 ) THEN
            IF( D3.LT.D1 ) THEN
               DMNMX = D1
            ELSE IF( D3.LT.D2 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            END IF
         ELSE
            IF( D3.LT.D2 ) THEN
               DMNMX = D2
            ELSE IF( D3.LT.D1 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            END IF
         END IF
!
         IF( DIR.EQ.0 ) THEN
!
!           Sort into decreasing order
!
            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF( D( J ).LT.DMNMX ) &
               GO TO 70
   80       CONTINUE
            I = I + 1
            IF( D( I ).GT.DMNMX ) &
               GO TO 80
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 60
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         ELSE
!
!           Sort into increasing order
!
            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF( D( J ).GT.DMNMX ) &
               GO TO 100
  110       CONTINUE
            I = I + 1
            IF( D( I ).LT.DMNMX ) &
               GO TO 110
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 90
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         END IF
      END IF
      IF( STKPNT.GT.0 ) &
         GO TO 10
      RETURN
!
!     End of DLASRT
!
      END
      SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASSQ  returns the values  scl  and  smsq  such that
!
!     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
!
!  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
!  assumed to be non-negative and  scl  returns the value
!
!     scl = max( scale, abs( x( i ) ) ).
!
!  scale and sumsq must be supplied in SCALE and SUMSQ and
!  scl and smsq are overwritten on SCALE and SUMSQ respectively.
!
!  The routine makes only one pass through the vector x.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements to be used from the vector X.
!
!  X       (input) DOUBLE PRECISION
!          The vector for which a scaled sum of squares is computed.
!             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
!
!  INCX    (input) INTEGER
!          The increment between successive values of the vector X.
!          INCX > 0.
!
!  SCALE   (input/output) DOUBLE PRECISION
!          On entry, the value  scale  in the equation above.
!          On exit, SCALE is overwritten with  scl , the scaling factor
!          for the sum of squares.
!
!  SUMSQ   (input/output) DOUBLE PRECISION
!          On entry, the value  sumsq  in the equation above.
!          On exit, SUMSQ is overwritten with  smsq , the basic sum of
!          squares from which  scl  has been factored out.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   ABSXI
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( X( IX ).NE.ZERO ) THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
!
!     End of DLASSQ
!
      END
      SUBROUTINE DOPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDQ, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AP( * ), Q( LDQ, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DOPGTR generates a real orthogonal matrix Q which is defined as the
!  product of n-1 elementary reflectors H(i) of order n, as returned by
!  DSPTRD using packed storage:
!
!  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
!
!  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U': Upper triangular packed storage used in previous
!                 call to DSPTRD;
!          = 'L': Lower triangular packed storage used in previous
!                 call to DSPTRD.
!
!  N       (input) INTEGER
!          The order of the matrix Q. N >= 0.
!
!  AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          The vectors which define the elementary reflectors, as
!          returned by DSPTRD.
!
!  TAU     (input) DOUBLE PRECISION array, dimension (N-1)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DSPTRD.
!
!  Q       (output) DOUBLE PRECISION array, dimension (LDQ,N)
!          The N-by-N orthogonal matrix Q.
!
!  LDQ     (input) INTEGER
!          The leading dimension of the array Q. LDQ >= max(1,N).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N-1)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IINFO, IJ, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DORG2L, DORG2R, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DOPGTR', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF( UPPER ) THEN
!
!        Q was determined by a call to DSPTRD with UPLO = 'U'
!
!        Unpack the vectors which define the elementary reflectors and
!        set the last row and column of Q equal to those of the unit
!        matrix
!
         IJ = 2
         DO 20 J = 1, N - 1
            DO 10 I = 1, J - 1
               Q( I, J ) = AP( IJ )
               IJ = IJ + 1
   10       CONTINUE
            IJ = IJ + 2
            Q( N, J ) = ZERO
   20    CONTINUE
         DO 30 I = 1, N - 1
            Q( I, N ) = ZERO
   30    CONTINUE
         Q( N, N ) = ONE
!
!        Generate Q(1:n-1,1:n-1)
!
         CALL DORG2L( N-1, N-1, N-1, Q, LDQ, TAU, WORK, IINFO )
!
      ELSE
!
!        Q was determined by a call to DSPTRD with UPLO = 'L'.
!
!        Unpack the vectors which define the elementary reflectors and
!        set the first row and column of Q equal to those of the unit
!        matrix
!
         Q( 1, 1 ) = ONE
         DO 40 I = 2, N
            Q( I, 1 ) = ZERO
   40    CONTINUE
         IJ = 3
         DO 60 J = 2, N
            Q( 1, J ) = ZERO
            DO 50 I = J + 1, N
               Q( I, J ) = AP( IJ )
               IJ = IJ + 1
   50       CONTINUE
            IJ = IJ + 2
   60    CONTINUE
         IF( N.GT.1 ) THEN
!
!           Generate Q(2:n,2:n)
!
            CALL DORG2R( N-1, N-1, N-1, Q( 2, 2 ), LDQ, TAU, WORK, &
                         IINFO )
         END IF
      END IF
      RETURN
!
!     End of DOPGTR
!
      END
      SUBROUTINE DOPMTR( SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK, &
                         INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS, UPLO
      INTEGER            INFO, LDC, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AP( * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DOPMTR overwrites the general real M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'T':      Q**T * C       C * Q**T
!
!  where Q is a real orthogonal matrix of order nq, with nq = m if
!  SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
!  nq-1 elementary reflectors, as returned by DSPTRD using packed
!  storage:
!
!  if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);
!
!  if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**T from the Left;
!          = 'R': apply Q or Q**T from the Right.
!
!  UPLO    (input) CHARACTER*1
!          = 'U': Upper triangular packed storage used in previous
!                 call to DSPTRD;
!          = 'L': Lower triangular packed storage used in previous
!                 call to DSPTRD.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'T':  Transpose, apply Q**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  AP      (input) DOUBLE PRECISION array, dimension
!                               (M*(M+1)/2) if SIDE = 'L'
!                               (N*(N+1)/2) if SIDE = 'R'
!          The vectors which define the elementary reflectors, as
!          returned by DSPTRD.  AP is modified by the routine but
!          restored on exit.
!
!  TAU     (input) DOUBLE PRECISION array, dimension (M-1) if SIDE = 'L'
!                                     or (N-1) if SIDE = 'R'
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DSPTRD.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                                   (N) if SIDE = 'L'
!                                   (M) if SIDE = 'R'
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            FORWRD, LEFT, NOTRAN, UPPER
      INTEGER            I, I1, I2, I3, IC, II, JC, MI, NI, NQ
      DOUBLE PRECISION   AII
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARF, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      UPPER = LSAME( UPLO, 'U' )
!
!     NQ is the order of Q
!
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DOPMTR', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) &
         RETURN
!
      IF( UPPER ) THEN
!
!        Q was determined by a call to DSPTRD with UPLO = 'U'
!
         FORWRD = ( LEFT .AND. NOTRAN ) .OR. &
                  ( .NOT.LEFT .AND. .NOT.NOTRAN )
!
         IF( FORWRD ) THEN
            I1 = 1
            I2 = NQ - 1
            I3 = 1
            II = 2
         ELSE
            I1 = NQ - 1
            I2 = 1
            I3 = -1
            II = NQ*( NQ+1 ) / 2 - 1
         END IF
!
         IF( LEFT ) THEN
            NI = N
         ELSE
            MI = M
         END IF
!
         DO 10 I = I1, I2, I3
            IF( LEFT ) THEN
!
!              H(i) is applied to C(1:i,1:n)
!
               MI = I
            ELSE
!
!              H(i) is applied to C(1:m,1:i)
!
               NI = I
            END IF
!
!           Apply H(i)
!
            AII = AP( II )
            AP( II ) = ONE
            CALL DLARF( SIDE, MI, NI, AP( II-I+1 ), 1, TAU( I ), C, LDC, &
                        WORK )
            AP( II ) = AII
!
            IF( FORWRD ) THEN
               II = II + I + 2
            ELSE
               II = II - I - 1
            END IF
   10    CONTINUE
      ELSE
!
!        Q was determined by a call to DSPTRD with UPLO = 'L'.
!
         FORWRD = ( LEFT .AND. .NOT.NOTRAN ) .OR. &
                  ( .NOT.LEFT .AND. NOTRAN )
!
         IF( FORWRD ) THEN
            I1 = 1
            I2 = NQ - 1
            I3 = 1
            II = 2
         ELSE
            I1 = NQ - 1
            I2 = 1
            I3 = -1
            II = NQ*( NQ+1 ) / 2 - 1
         END IF
!
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
!
         DO 20 I = I1, I2, I3
            AII = AP( II )
            AP( II ) = ONE
            IF( LEFT ) THEN
!
!              H(i) is applied to C(i+1:m,1:n)
!
               MI = M - I
               IC = I + 1
            ELSE
!
!              H(i) is applied to C(1:m,i+1:n)
!
               NI = N - I
               JC = I + 1
            END IF
!
!           Apply H(i)
!
            CALL DLARF( SIDE, MI, NI, AP( II ), 1, TAU( I ), &
                        C( IC, JC ), LDC, WORK )
            AP( II ) = AII
!
            IF( FORWRD ) THEN
               II = II + NQ - I + 1
            ELSE
               II = II - NQ + I - 2
            END IF
   20    CONTINUE
      END IF
      RETURN
!
!     End of DOPMTR
!
      END
      SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORG2L generates an m by n real matrix Q with orthonormal columns,
!  which is defined as the last n columns of a product of k elementary
!  reflectors of order m
!
!        Q  =  H(k) . . . H(2) H(1)
!
!  as returned by DGEQLF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the (n-k+i)-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQLF in the last k columns of its array
!          argument A.
!          On exit, the m by n matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQLF.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, II, J, L
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARF, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORG2L', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) &
         RETURN
!
!     Initialise columns 1:n-k to columns of the unit matrix
!
      DO 20 J = 1, N - K
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( M-N+J, J ) = ONE
   20 CONTINUE
!
      DO 40 I = 1, K
         II = N - K + I
!
!        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
!
         A( M-N+II, II ) = ONE
         CALL DLARF( 'Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A, &
                     LDA, WORK )
         CALL DSCAL( M-N+II-1, -TAU( I ), A( 1, II ), 1 )
         A( M-N+II, II ) = ONE - TAU( I )
!
!        Set A(m-k+i+1:m,n-k+i) to zero
!
         DO 30 L = M - N + II + 1, M
            A( L, II ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
!
!     End of DORG2L
!
      END
      SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORG2R generates an m by n real matrix Q with orthonormal columns,
!  which is defined as the first n columns of a product of k elementary
!  reflectors of order m
!
!        Q  =  H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQRF in the first k columns of its array
!          argument A.
!          On exit, the m-by-n matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, L
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLARF, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORG2R', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) &
         RETURN
!
!     Initialise columns k+1:n to columns of the unit matrix
!
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
!
      DO 40 I = K, 1, -1
!
!        Apply H(i) to A(i:m,i:n) from the left
!
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), &
                        A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M ) &
            CALL DSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
!
!        Set A(1:i-1,i) to zero
!
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
!
!     End of DORG2R
!
      END
      SUBROUTINE DPOTF2( UPLO, N, A, LDA, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DPOTF2 computes the Cholesky factorization of a real symmetric
!  positive definite matrix A.
!
!  The factorization has the form
!     A = U' * U ,  if UPLO = 'U', or
!     A = L  * L',  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n by n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n by n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the factor U or L from the Cholesky
!          factorization A = U'*U  or A = L*L'.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, the leading minor of order k is not
!               positive definite, and the factorization could not be
!               completed.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMV, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF( UPPER ) THEN
!
!        Compute the Cholesky factorization A = U'*U.
!
         DO 10 J = 1, N
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
            AJJ = A( J, J ) - DDOT( J-1, A( 1, J ), 1, A( 1, J ), 1 )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
!
!           Compute elements J+1:N of row J.
!
            IF( J.LT.N ) THEN
               CALL DGEMV( 'Transpose', J-1, N-J, -ONE, A( 1, J+1 ), &
                           LDA, A( 1, J ), 1, ONE, A( J, J+1 ), LDA )
               CALL DSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            END IF
   10    CONTINUE
      ELSE
!
!        Compute the Cholesky factorization A = L*L'.
!
         DO 20 J = 1, N
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
            AJJ = A( J, J ) - DDOT( J-1, A( J, 1 ), LDA, A( J, 1 ), &
                  LDA )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
!
!           Compute elements J+1:N of column J.
!
            IF( J.LT.N ) THEN
               CALL DGEMV( 'No transpose', N-J, J-1, -ONE, A( J+1, 1 ), &
                           LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1 )
               CALL DSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
      GO TO 40
!
   30 CONTINUE
      INFO = J
!
   40 CONTINUE
      RETURN
!
!     End of DPOTF2
!
      END
      SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DPOTRF computes the Cholesky factorization of a real symmetric
!  positive definite matrix A.
!
!  The factorization has the form
!     A = U**T * U,  if UPLO = 'U', or
!     A = L  * L**T,  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  This is the block version of the algorithm, calling Level 3 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the factor U or L from the Cholesky
!          factorization A = U**T*U or A = L*L**T.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i is not
!                positive definite, and the factorization could not be
!                completed.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JB, NB
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMM, DPOTF2, DSYRK, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'DPOTRF', UPLO, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
!
!        Use unblocked code.
!
         CALL DPOTF2( UPLO, N, A, LDA, INFO )
      ELSE
!
!        Use blocked code.
!
         IF( UPPER ) THEN
!
!           Compute the Cholesky factorization A = U'*U.
!
            DO 10 J = 1, N, NB
!
!              Update and factorize the current diagonal block and test
!              for non-positive-definiteness.
!
               JB = MIN( NB, N-J+1 )
               CALL DSYRK( 'Upper', 'Transpose', JB, J-1, -ONE, &
                           A( 1, J ), LDA, ONE, A( J, J ), LDA )
               CALL DPOTF2( 'Upper', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 ) &
                  GO TO 30
               IF( J+JB.LE.N ) THEN
!
!                 Compute the current block row.
!
                  CALL DGEMM( 'Transpose', 'No transpose', JB, N-J-JB+1, &
                              J-1, -ONE, A( 1, J ), LDA, A( 1, J+JB ), &
                              LDA, ONE, A( J, J+JB ), LDA )
                  CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', &
                              JB, N-J-JB+1, ONE, A( J, J ), LDA, &
                              A( J, J+JB ), LDA )
               END IF
   10       CONTINUE
!
         ELSE
!
!           Compute the Cholesky factorization A = L*L'.
!
            DO 20 J = 1, N, NB
!
!              Update and factorize the current diagonal block and test
!              for non-positive-definiteness.
!
               JB = MIN( NB, N-J+1 )
               CALL DSYRK( 'Lower', 'No transpose', JB, J-1, -ONE, &
                           A( J, 1 ), LDA, ONE, A( J, J ), LDA )
               CALL DPOTF2( 'Lower', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 ) &
                  GO TO 30
               IF( J+JB.LE.N ) THEN
!
!                 Compute the current block column.
!
                  CALL DGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB, &
                              J-1, -ONE, A( J+JB, 1 ), LDA, A( J, 1 ), &
                              LDA, ONE, A( J+JB, J ), LDA )
                  CALL DTRSM( 'Right', 'Lower', 'Transpose', 'Non-unit', &
                              N-J-JB+1, JB, ONE, A( J, J ), LDA, &
                              A( J+JB, J ), LDA )
               END IF
   20       CONTINUE
         END IF
      END IF
      GO TO 40
!
   30 CONTINUE
      INFO = INFO + J - 1
!
   40 CONTINUE
      RETURN
!
!     End of DPOTRF
!
      END
      SUBROUTINE DPPTRF( UPLO, N, AP, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AP( * )
!     ..
!
!  Purpose
!  =======
!
!  DPPTRF computes the Cholesky factorization of a real symmetric
!  positive definite matrix A stored in packed format.
!
!  The factorization has the form
!     A = U**T * U,  if UPLO = 'U', or
!     A = L  * L**T,  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          On entry, the upper or lower triangle of the symmetric matrix
!          A, packed columnwise in a linear array.  The j-th column of A
!          is stored in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!          See below for further details.
!
!          On exit, if INFO = 0, the triangular factor U or L from the
!          Cholesky factorization A = U**T*U or A = L*L**T, in the same
!          storage format as A.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i is not
!                positive definite, and the factorization could not be
!                completed.
!
!  Further Details
!  ======= =======
!
!  The packed storage scheme is illustrated by the following example
!  when N = 4, UPLO = 'U':
!
!  Two-dimensional storage of the symmetric matrix A:
!
!     a11 a12 a13 a14
!         a22 a23 a24
!             a33 a34     (aij = aji)
!                 a44
!
!  Packed storage of the upper triangle of A:
!
!  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JC, JJ
      DOUBLE PRECISION   AJJ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL, DSPR, DTPSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPPTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF( UPPER ) THEN
!
!        Compute the Cholesky factorization A = U'*U.
!
         JJ = 0
         DO 10 J = 1, N
            JC = JJ + 1
            JJ = JJ + J
!
!           Compute elements 1:J-1 of column J.
!
            IF( J.GT.1 ) &
               CALL DTPSV( 'Upper', 'Transpose', 'Non-unit', J-1, AP, &
                           AP( JC ), 1 )
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
            AJJ = AP( JJ ) - DDOT( J-1, AP( JC ), 1, AP( JC ), 1 )
            IF( AJJ.LE.ZERO ) THEN
               AP( JJ ) = AJJ
               GO TO 30
            END IF
            AP( JJ ) = SQRT( AJJ )
   10    CONTINUE
      ELSE
!
!        Compute the Cholesky factorization A = L*L'.
!
         JJ = 1
         DO 20 J = 1, N
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
            AJJ = AP( JJ )
            IF( AJJ.LE.ZERO ) THEN
               AP( JJ ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            AP( JJ ) = AJJ
!
!           Compute elements J+1:N of column J and update the trailing
!           submatrix.
!
            IF( J.LT.N ) THEN
               CALL DSCAL( N-J, ONE / AJJ, AP( JJ+1 ), 1 )
               CALL DSPR( 'Lower', N-J, -ONE, AP( JJ+1 ), 1, &
                          AP( JJ+N-J+1 ) )
               JJ = JJ + N - J + 1
            END IF
   20    CONTINUE
      END IF
      GO TO 40
!
   30 CONTINUE
      INFO = J
!
   40 CONTINUE
      RETURN
!
!     End of DPPTRF
!
      END
      SUBROUTINE DSPGST( ITYPE, UPLO, N, AP, BP, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AP( * ), BP( * )
!     ..
!
!  Purpose
!  =======
!
!  DSPGST reduces a real symmetric-definite generalized eigenproblem
!  to standard form, using packed storage.
!
!  If ITYPE = 1, the problem is A*x = lambda*B*x,
!  and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
!
!  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
!  B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
!
!  B must have been previously factorized as U**T*U or L*L**T by DPPTRF.
!
!  Arguments
!  =========
!
!  ITYPE   (input) INTEGER
!          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);
!          = 2 or 3: compute U*A*U**T or L**T*A*L.
!
!  UPLO    (input) CHARACTER
!          = 'U':  Upper triangle of A is stored and B is factored as
!                  U**T*U;
!          = 'L':  Lower triangle of A is stored and B is factored as
!                  L*L**T.
!
!  N       (input) INTEGER
!          The order of the matrices A and B.  N >= 0.
!
!  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          On entry, the upper or lower triangle of the symmetric matrix
!          A, packed columnwise in a linear array.  The j-th column of A
!          is stored in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!
!          On exit, if INFO = 0, the transformed matrix, stored in the
!          same format as A.
!
!  BP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          The triangular factor from the Cholesky factorization of B,
!          stored in the same format as A, as returned by DPPTRF.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, HALF
      PARAMETER          ( ONE = 1.0D0, HALF = 0.5D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, J1, J1J1, JJ, K, K1, K1K1, KK
      DOUBLE PRECISION   AJJ, AKK, BJJ, BKK, CT
!     ..
!     .. External Subroutines ..
      EXTERNAL           DAXPY, DSCAL, DSPMV, DSPR2, DTPMV, DTPSV, &
                         XERBLA
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSPGST', -INFO )
         RETURN
      END IF
!
      IF( ITYPE.EQ.1 ) THEN
         IF( UPPER ) THEN
!
!           Compute inv(U')*A*inv(U)
!
!           J1 and JJ are the indices of A(1,j) and A(j,j)
!
            JJ = 0
            DO 10 J = 1, N
               J1 = JJ + 1
               JJ = JJ + J
!
!              Compute the j-th column of the upper triangle of A
!
               BJJ = BP( JJ )
               CALL DTPSV( UPLO, 'Transpose', 'Nonunit', J, BP, &
                           AP( J1 ), 1 )
               CALL DSPMV( UPLO, J-1, -ONE, AP, BP( J1 ), 1, ONE, &
                           AP( J1 ), 1 )
               CALL DSCAL( J-1, ONE / BJJ, AP( J1 ), 1 )
               AP( JJ ) = ( AP( JJ )-DDOT( J-1, AP( J1 ), 1, BP( J1 ), &
                          1 ) ) / BJJ
   10       CONTINUE
         ELSE
!
!           Compute inv(L)*A*inv(L')
!
!           KK and K1K1 are the indices of A(k,k) and A(k+1,k+1)
!
            KK = 1
            DO 20 K = 1, N
               K1K1 = KK + N - K + 1
!
!              Update the lower triangle of A(k:n,k:n)
!
               AKK = AP( KK )
               BKK = BP( KK )
               AKK = AKK / BKK**2
               AP( KK ) = AKK
               IF( K.LT.N ) THEN
                  CALL DSCAL( N-K, ONE / BKK, AP( KK+1 ), 1 )
                  CT = -HALF*AKK
                  CALL DAXPY( N-K, CT, BP( KK+1 ), 1, AP( KK+1 ), 1 )
                  CALL DSPR2( UPLO, N-K, -ONE, AP( KK+1 ), 1, &
                              BP( KK+1 ), 1, AP( K1K1 ) )
                  CALL DAXPY( N-K, CT, BP( KK+1 ), 1, AP( KK+1 ), 1 )
                  CALL DTPSV( UPLO, 'No transpose', 'Non-unit', N-K, &
                              BP( K1K1 ), AP( KK+1 ), 1 )
               END IF
               KK = K1K1
   20       CONTINUE
         END IF
      ELSE
         IF( UPPER ) THEN
!
!           Compute U*A*U'
!
!           K1 and KK are the indices of A(1,k) and A(k,k)
!
            KK = 0
            DO 30 K = 1, N
               K1 = KK + 1
               KK = KK + K
!
!              Update the upper triangle of A(1:k,1:k)
!
               AKK = AP( KK )
               BKK = BP( KK )
               CALL DTPMV( UPLO, 'No transpose', 'Non-unit', K-1, BP, &
                           AP( K1 ), 1 )
               CT = HALF*AKK
               CALL DAXPY( K-1, CT, BP( K1 ), 1, AP( K1 ), 1 )
               CALL DSPR2( UPLO, K-1, ONE, AP( K1 ), 1, BP( K1 ), 1, &
                           AP )
               CALL DAXPY( K-1, CT, BP( K1 ), 1, AP( K1 ), 1 )
               CALL DSCAL( K-1, BKK, AP( K1 ), 1 )
               AP( KK ) = AKK*BKK**2
   30       CONTINUE
         ELSE
!
!           Compute L'*A*L
!
!           JJ and J1J1 are the indices of A(j,j) and A(j+1,j+1)
!
            JJ = 1
            DO 40 J = 1, N
               J1J1 = JJ + N - J + 1
!
!              Compute the j-th column of the lower triangle of A
!
               AJJ = AP( JJ )
               BJJ = BP( JJ )
               AP( JJ ) = AJJ*BJJ + DDOT( N-J, AP( JJ+1 ), 1, &
                          BP( JJ+1 ), 1 )
               CALL DSCAL( N-J, BJJ, AP( JJ+1 ), 1 )
               CALL DSPMV( UPLO, N-J, ONE, AP( J1J1 ), BP( JJ+1 ), 1, &
                           ONE, AP( JJ+1 ), 1 )
               CALL DTPMV( UPLO, 'Transpose', 'Non-unit', N-J+1, &
                           BP( JJ ), AP( JJ ), 1 )
               JJ = J1J1
   40       CONTINUE
         END IF
      END IF
      RETURN
!
!     End of DSPGST
!
      END
      SUBROUTINE DSPTRD( UPLO, N, AP, D, E, TAU, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   AP( * ), D( * ), E( * ), TAU( * )
!     ..
!
!  Purpose
!  =======
!
!  DSPTRD reduces a real symmetric matrix A stored in packed form to
!  symmetric tridiagonal form T by an orthogonal similarity
!  transformation: Q**T * A * Q = T.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          On entry, the upper or lower triangle of the symmetric matrix
!          A, packed columnwise in a linear array.  The j-th column of A
!          is stored in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!          of A are overwritten by the corresponding elements of the
!          tridiagonal matrix T, and the elements above the first
!          superdiagonal, with the array TAU, represent the orthogonal
!          matrix Q as a product of elementary reflectors; if UPLO
!          = 'L', the diagonal and first subdiagonal of A are over-
!          written by the corresponding elements of the tridiagonal
!          matrix T, and the elements below the first subdiagonal, with
!          the array TAU, represent the orthogonal matrix Q as a product
!          of elementary reflectors. See Further Details.
!
!  D       (output) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = A(i,i).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!
!  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n-1) . . . H(2) H(1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in AP,
!  overwriting A(1:i-1,i+1), and tau is stored in TAU(i).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in AP,
!  overwriting A(i+2:n,i), and tau is stored in TAU(i).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO, HALF
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0, &
                         HALF = 1.0D0 / 2.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, I1, I1I1, II
      DOUBLE PRECISION   ALPHA, TAUI
!     ..
!     .. External Subroutines ..
      EXTERNAL           DAXPY, DLARFG, DSPMV, DSPR2, XERBLA
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSPTRD', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) &
         RETURN
!
      IF( UPPER ) THEN
!
!        Reduce the upper triangle of A.
!        I1 is the index in AP of A(1,I+1).
!
         I1 = N*( N-1 ) / 2 + 1
         DO 10 I = N - 1, 1, -1
!
!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(1:i-1,i+1)
!
            CALL DLARFG( I, AP( I1+I-1 ), AP( I1 ), 1, TAUI )
            E( I ) = AP( I1+I-1 )
!
            IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(1:i,1:i)
!
               AP( I1+I-1 ) = ONE
!
!              Compute  y := tau * A * v  storing y in TAU(1:i)
!
               CALL DSPMV( UPLO, I, TAUI, AP, AP( I1 ), 1, ZERO, TAU, &
                           1 )
!
!              Compute  w := y - 1/2 * tau * (y'*v) * v
!
               ALPHA = -HALF*TAUI*DDOT( I, TAU, 1, AP( I1 ), 1 )
               CALL DAXPY( I, ALPHA, AP( I1 ), 1, TAU, 1 )
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'
!
               CALL DSPR2( UPLO, I, -ONE, AP( I1 ), 1, TAU, 1, AP )
!
               AP( I1+I-1 ) = E( I )
            END IF
            D( I+1 ) = AP( I1+I )
            TAU( I ) = TAUI
            I1 = I1 - I
   10    CONTINUE
         D( 1 ) = AP( 1 )
      ELSE
!
!        Reduce the lower triangle of A. II is the index in AP of
!        A(i,i) and I1I1 is the index of A(i+1,i+1).
!
         II = 1
         DO 20 I = 1, N - 1
            I1I1 = II + N - I + 1
!
!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(i+2:n,i)
!
            CALL DLARFG( N-I, AP( II+1 ), AP( II+2 ), 1, TAUI )
            E( I ) = AP( II+1 )
!
            IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(i+1:n,i+1:n)
!
               AP( II+1 ) = ONE
!
!              Compute  y := tau * A * v  storing y in TAU(i:n-1)
!
               CALL DSPMV( UPLO, N-I, TAUI, AP( I1I1 ), AP( II+1 ), 1, &
                           ZERO, TAU( I ), 1 )
!
!              Compute  w := y - 1/2 * tau * (y'*v) * v
!
               ALPHA = -HALF*TAUI*DDOT( N-I, TAU( I ), 1, AP( II+1 ), &
                       1 )
               CALL DAXPY( N-I, ALPHA, AP( II+1 ), 1, TAU( I ), 1 )
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'
!
               CALL DSPR2( UPLO, N-I, -ONE, AP( II+1 ), 1, TAU( I ), 1, &
                           AP( I1I1 ) )
!
               AP( II+1 ) = E( I )
            END IF
            D( I ) = AP( II )
            TAU( I ) = TAUI
            II = I1I1
   20    CONTINUE
         D( N ) = AP( II )
      END IF
!
      RETURN
!
!     End of DSPTRD
!
      END
      SUBROUTINE DSTEBZ( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E, &
                         M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK, &
                         INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          ORDER, RANGE
      INTEGER            IL, INFO, IU, M, N, NSPLIT
      DOUBLE PRECISION   ABSTOL, VL, VU
!     ..
!     .. Array Arguments ..
      INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DSTEBZ computes the eigenvalues of a symmetric tridiagonal
!  matrix T.  The user may ask for all eigenvalues, all eigenvalues
!  in the half-open interval (VL, VU], or the IL-th through IU-th
!  eigenvalues.
!
!  To avoid overflow, the matrix must be scaled so that its
!  largest element is no greater than overflow**(1/2) *
!  underflow**(1/4) in absolute value, and for greatest
!  accuracy, it should not be much smaller than that.
!
!  See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
!  Matrix", Report CS41, Computer Science Dept., Stanford
!  University, July 21, 1966.
!
!  Arguments
!  =========
!
!  RANGE   (input) CHARACTER
!          = 'A': ("All")   all eigenvalues will be found.
!          = 'V': ("Value") all eigenvalues in the half-open interval
!                           (VL, VU] will be found.
!          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the
!                           entire matrix) will be found.
!
!  ORDER   (input) CHARACTER
!          = 'B': ("By Block") the eigenvalues will be grouped by
!                              split-off block (see IBLOCK, ISPLIT) and
!                              ordered from smallest to largest within
!                              the block.
!          = 'E': ("Entire matrix")
!                              the eigenvalues for the entire matrix
!                              will be ordered from smallest to
!                              largest.
!
!  N       (input) INTEGER
!          The order of the tridiagonal matrix T.  N >= 0.
!
!  VL      (input) DOUBLE PRECISION
!  VU      (input) DOUBLE PRECISION
!          If RANGE='V', the lower and upper bounds of the interval to
!          be searched for eigenvalues.  Eigenvalues less than or equal
!          to VL, or greater than VU, will not be returned.  VL < VU.
!          Not referenced if RANGE = 'A' or 'I'.
!
!  IL      (input) INTEGER
!  IU      (input) INTEGER
!          If RANGE='I', the indices (in ascending order) of the
!          smallest and largest eigenvalues to be returned.
!          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!          Not referenced if RANGE = 'A' or 'V'.
!
!  ABSTOL  (input) DOUBLE PRECISION
!          The absolute tolerance for the eigenvalues.  An eigenvalue
!          (or cluster) is considered to be located if it has been
!          determined to lie in an interval whose width is ABSTOL or
!          less.  If ABSTOL is less than or equal to zero, then ULP*|T|
!          will be used, where |T| means the 1-norm of T.
!
!          Eigenvalues will be computed most accurately when ABSTOL is
!          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The n diagonal elements of the tridiagonal matrix T.
!
!  E       (input) DOUBLE PRECISION array, dimension (N-1)
!          The (n-1) off-diagonal elements of the tridiagonal matrix T.
!
!  M       (output) INTEGER
!          The actual number of eigenvalues found. 0 <= M <= N.
!          (See also the description of INFO=2,3.)
!
!  NSPLIT  (output) INTEGER
!          The number of diagonal blocks in the matrix T.
!          1 <= NSPLIT <= N.
!
!  W       (output) DOUBLE PRECISION array, dimension (N)
!          On exit, the first M elements of W will contain the
!          eigenvalues.  (DSTEBZ may use the remaining N-M elements as
!          workspace.)
!
!  IBLOCK  (output) INTEGER array, dimension (N)
!          At each row/column j where E(j) is zero or small, the
!          matrix T is considered to split into a block diagonal
!          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which
!          block (from 1 to the number of blocks) the eigenvalue W(i)
!          belongs.  (DSTEBZ may use the remaining N-M elements as
!          workspace.)
!
!  ISPLIT  (output) INTEGER array, dimension (N)
!          The splitting points, at which T breaks up into submatrices.
!          The first submatrix consists of rows/columns 1 to ISPLIT(1),
!          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
!          etc., and the NSPLIT-th consists of rows/columns
!          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
!          (Only the first NSPLIT elements will actually be used, but
!          since the user cannot know a priori what value NSPLIT will
!          have, N words must be reserved for ISPLIT.)
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)
!
!  IWORK   (workspace) INTEGER array, dimension (3*N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  some or all of the eigenvalues failed to converge or
!                were not computed:
!                =1 or 3: Bisection failed to converge for some
!                        eigenvalues; these eigenvalues are flagged by a
!                        negative block number.  The effect is that the
!                        eigenvalues may not be as accurate as the
!                        absolute and relative tolerances.  This is
!                        generally caused by unexpectedly inaccurate
!                        arithmetic.
!                =2 or 3: RANGE='I' only: Not all of the eigenvalues
!                        IL:IU were found.
!                        Effect: M < IU+1-IL
!                        Cause:  non-monotonic arithmetic, causing the
!                                Sturm sequence to be non-monotonic.
!                        Cure:   recalculate, using RANGE='A', and pick
!                                out eigenvalues IL:IU.  In some cases,
!                                increasing the PARAMETER "FUDGE" may
!                                make things work.
!                = 4:    RANGE='I', and the Gershgorin interval
!                        initially used was too small.  No eigenvalues
!                        were computed.
!                        Probable cause: your machine has sloppy
!                                        floating-point arithmetic.
!                        Cure: Increase the PARAMETER "FUDGE",
!                              recompile, and try again.
!
!  Internal Parameters
!  ===================
!
!  RELFAC  DOUBLE PRECISION, default = 2.0e0
!          The relative tolerance.  An interval (a,b] lies within
!          "relative tolerance" if  b-a < RELFAC*ulp*max(|a|,|b|),
!          where "ulp" is the machine precision (distance from 1 to
!          the next larger floating point number.)
!
!  FUDGE   DOUBLE PRECISION, default = 2
!          A "fudge factor" to widen the Gershgorin intervals.  Ideally,
!          a value of 1 should work, but on machines with sloppy
!          arithmetic, this needs to be larger.  The default for
!          publicly released versions should be large enough to handle
!          the worst machine around.  Note that this has no effect
!          on accuracy of the solution.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, HALF
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, &
                         HALF = 1.0D0 / TWO )
      DOUBLE PRECISION   FUDGE, RELFAC
      PARAMETER          ( FUDGE = 2.0D0, RELFAC = 2.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NCNVRG, TOOFEW
      INTEGER            IB, IBEGIN, IDISCL, IDISCU, IE, IEND, IINFO, &
                         IM, IN, IOFF, IORDER, IOUT, IRANGE, ITMAX, &
                         ITMP1, IW, IWOFF, J, JB, JDISC, JE, NB, NWL, &
                         NWU
      DOUBLE PRECISION   ATOLI, BNORM, GL, GU, PIVMIN, RTOLI, SAFEMN, &
                         TMP1, TMP2, TNORM, ULP, WKILL, WL, WLU, WU, WUL
!     ..
!     .. Local Arrays ..
      INTEGER            IDUMMA( 1 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, ILAENV, DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAEBZ, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
!     Decode RANGE
!
      IF( LSAME( RANGE, 'A' ) ) THEN
         IRANGE = 1
      ELSE IF( LSAME( RANGE, 'V' ) ) THEN
         IRANGE = 2
      ELSE IF( LSAME( RANGE, 'I' ) ) THEN
         IRANGE = 3
      ELSE
         IRANGE = 0
      END IF
!
!     Decode ORDER
!
      IF( LSAME( ORDER, 'B' ) ) THEN
         IORDER = 2
      ELSE IF( LSAME( ORDER, 'E' ) ) THEN
         IORDER = 1
      ELSE
         IORDER = 0
      END IF
!
!     Check for Errors
!
      IF( IRANGE.LE.0 ) THEN
         INFO = -1
      ELSE IF( IORDER.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( IRANGE.EQ.2 .AND. VL.GE.VU ) THEN
         INFO = -5
      ELSE IF( IRANGE.EQ.3 .AND. ( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) ) &
                THEN
         INFO = -6
      ELSE IF( IRANGE.EQ.3 .AND. ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) ) &
                THEN
         INFO = -7
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEBZ', -INFO )
         RETURN
      END IF
!
!     Initialize error flags
!
      INFO = 0
      NCNVRG = .FALSE.
      TOOFEW = .FALSE.
!
!     Quick return if possible
!
      M = 0
      IF( N.EQ.0 ) &
         RETURN
!
!     Simplifications:
!
      IF( IRANGE.EQ.3 .AND. IL.EQ.1 .AND. IU.EQ.N ) &
         IRANGE = 1
!
!     Get machine constants
!     NB is the minimum vector length for vector bisection, or 0
!     if only scalar is to be done.
!
      SAFEMN = DLAMCH( 'S' )
      ULP = DLAMCH( 'P' )
      RTOLI = ULP*RELFAC
      NB = ILAENV( 1, 'DSTEBZ', ' ', N, -1, -1, -1 )
      IF( NB.LE.1 ) &
         NB = 0
!
!     Special Case when N=1
!
      IF( N.EQ.1 ) THEN
         NSPLIT = 1
         ISPLIT( 1 ) = 1
         IF( IRANGE.EQ.2 .AND. ( VL.GE.D( 1 ) .OR. VU.LT.D( 1 ) ) ) THEN
            M = 0
         ELSE
            W( 1 ) = D( 1 )
            IBLOCK( 1 ) = 1
            M = 1
         END IF
         RETURN
      END IF
!
!     Compute Splitting Points
!
      NSPLIT = 1
      WORK( N ) = ZERO
      PIVMIN = ONE
!
      DO 10 J = 2, N
         TMP1 = E( J-1 )**2
         IF( ABS( D( J )*D( J-1 ) )*ULP**2+SAFEMN.GT.TMP1 ) THEN
            ISPLIT( NSPLIT ) = J - 1
            NSPLIT = NSPLIT + 1
            WORK( J-1 ) = ZERO
         ELSE
            WORK( J-1 ) = TMP1
            PIVMIN = MAX( PIVMIN, TMP1 )
         END IF
   10 CONTINUE
      ISPLIT( NSPLIT ) = N
      PIVMIN = PIVMIN*SAFEMN
!
!     Compute Interval and ATOLI
!
      IF( IRANGE.EQ.3 ) THEN
!
!        RANGE='I': Compute the interval containing eigenvalues
!                   IL through IU.
!
!        Compute Gershgorin interval for entire (split) matrix
!        and use it as the initial interval
!
         GU = D( 1 )
         GL = D( 1 )
         TMP1 = ZERO
!
         DO 20 J = 1, N - 1
            TMP2 = SQRT( WORK( J ) )
            GU = MAX( GU, D( J )+TMP1+TMP2 )
            GL = MIN( GL, D( J )-TMP1-TMP2 )
            TMP1 = TMP2
   20    CONTINUE
!
         GU = MAX( GU, D( N )+TMP1 )
         GL = MIN( GL, D( N )-TMP1 )
         TNORM = MAX( ABS( GL ), ABS( GU ) )
         GL = GL - FUDGE*TNORM*ULP*N - FUDGE*TWO*PIVMIN
         GU = GU + FUDGE*TNORM*ULP*N + FUDGE*PIVMIN
!
!        Compute Iteration parameters
!
         ITMAX = INT( ( LOG( TNORM+PIVMIN )-LOG( PIVMIN ) ) / &
                 LOG( TWO ) ) + 2
         IF( ABSTOL.LE.ZERO ) THEN
            ATOLI = ULP*TNORM
         ELSE
            ATOLI = ABSTOL
         END IF
!
         WORK( N+1 ) = GL
         WORK( N+2 ) = GL
         WORK( N+3 ) = GU
         WORK( N+4 ) = GU
         WORK( N+5 ) = GL
         WORK( N+6 ) = GU
         IWORK( 1 ) = -1
         IWORK( 2 ) = -1
         IWORK( 3 ) = N + 1
         IWORK( 4 ) = N + 1
         IWORK( 5 ) = IL - 1
         IWORK( 6 ) = IU
!
         CALL DLAEBZ( 3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN, D, E, &
                      WORK, IWORK( 5 ), WORK( N+1 ), WORK( N+5 ), IOUT, &
                      IWORK, W, IBLOCK, IINFO )
!
         IF( IWORK( 6 ).EQ.IU ) THEN
            WL = WORK( N+1 )
            WLU = WORK( N+3 )
            NWL = IWORK( 1 )
            WU = WORK( N+4 )
            WUL = WORK( N+2 )
            NWU = IWORK( 4 )
         ELSE
            WL = WORK( N+2 )
            WLU = WORK( N+4 )
            NWL = IWORK( 2 )
            WU = WORK( N+3 )
            WUL = WORK( N+1 )
            NWU = IWORK( 3 )
         END IF
!
         IF( NWL.LT.0 .OR. NWL.GE.N .OR. NWU.LT.1 .OR. NWU.GT.N ) THEN
            INFO = 4
            RETURN
         END IF
      ELSE
!
!        RANGE='A' or 'V' -- Set ATOLI
!
         TNORM = MAX( ABS( D( 1 ) )+ABS( E( 1 ) ), &
                 ABS( D( N ) )+ABS( E( N-1 ) ) )
!
         DO 30 J = 2, N - 1
            TNORM = MAX( TNORM, ABS( D( J ) )+ABS( E( J-1 ) )+ &
                    ABS( E( J ) ) )
   30    CONTINUE
!
         IF( ABSTOL.LE.ZERO ) THEN
            ATOLI = ULP*TNORM
         ELSE
            ATOLI = ABSTOL
         END IF
!
         IF( IRANGE.EQ.2 ) THEN
            WL = VL
            WU = VU
         END IF
      END IF
!
!     Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU.
!     NWL accumulates the number of eigenvalues .le. WL,
!     NWU accumulates the number of eigenvalues .le. WU
!
      M = 0
      IEND = 0
      INFO = 0
      NWL = 0
      NWU = 0
!
      DO 70 JB = 1, NSPLIT
         IOFF = IEND
         IBEGIN = IOFF + 1
         IEND = ISPLIT( JB )
         IN = IEND - IOFF
!
         IF( IN.EQ.1 ) THEN
!
!           Special Case -- IN=1
!
            IF( IRANGE.EQ.1 .OR. WL.GE.D( IBEGIN )-PIVMIN ) &
               NWL = NWL + 1
            IF( IRANGE.EQ.1 .OR. WU.GE.D( IBEGIN )-PIVMIN ) &
               NWU = NWU + 1
            IF( IRANGE.EQ.1 .OR. ( WL.LT.D( IBEGIN )-PIVMIN .AND. WU.GE. &
                D( IBEGIN )-PIVMIN ) ) THEN
               M = M + 1
               W( M ) = D( IBEGIN )
               IBLOCK( M ) = JB
            END IF
         ELSE
!
!           General Case -- IN > 1
!
!           Compute Gershgorin Interval
!           and use it as the initial interval
!
            GU = D( IBEGIN )
            GL = D( IBEGIN )
            TMP1 = ZERO
!
            DO 40 J = IBEGIN, IEND - 1
               TMP2 = ABS( E( J ) )
               GU = MAX( GU, D( J )+TMP1+TMP2 )
               GL = MIN( GL, D( J )-TMP1-TMP2 )
               TMP1 = TMP2
   40       CONTINUE
!
            GU = MAX( GU, D( IEND )+TMP1 )
            GL = MIN( GL, D( IEND )-TMP1 )
            BNORM = MAX( ABS( GL ), ABS( GU ) )
            GL = GL - FUDGE*BNORM*ULP*IN - FUDGE*PIVMIN
            GU = GU + FUDGE*BNORM*ULP*IN + FUDGE*PIVMIN
!
!           Compute ATOLI for the current submatrix
!
            IF( ABSTOL.LE.ZERO ) THEN
               ATOLI = ULP*MAX( ABS( GL ), ABS( GU ) )
            ELSE
               ATOLI = ABSTOL
            END IF
!
            IF( IRANGE.GT.1 ) THEN
               IF( GU.LT.WL ) THEN
                  NWL = NWL + IN
                  NWU = NWU + IN
                  GO TO 70
               END IF
               GL = MAX( GL, WL )
               GU = MIN( GU, WU )
               IF( GL.GE.GU ) &
                  GO TO 70
            END IF
!
!           Set Up Initial Interval
!
            WORK( N+1 ) = GL
            WORK( N+IN+1 ) = GU
            CALL DLAEBZ( 1, 0, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN, &
                         D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ), &
                         IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IM, &
                         IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
!
            NWL = NWL + IWORK( 1 )
            NWU = NWU + IWORK( IN+1 )
            IWOFF = M - IWORK( 1 )
!
!           Compute Eigenvalues
!
            ITMAX = INT( ( LOG( GU-GL+PIVMIN )-LOG( PIVMIN ) ) / &
                    LOG( TWO ) ) + 2
            CALL DLAEBZ( 2, ITMAX, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN, &
                         D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ), &
                         IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IOUT, &
                         IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
!
!           Copy Eigenvalues Into W and IBLOCK
!           Use -JB for block number for unconverged eigenvalues.
!
            DO 60 J = 1, IOUT
               TMP1 = HALF*( WORK( J+N )+WORK( J+IN+N ) )
!
!              Flag non-convergence.
!
               IF( J.GT.IOUT-IINFO ) THEN
                  NCNVRG = .TRUE.
                  IB = -JB
               ELSE
                  IB = JB
               END IF
               DO 50 JE = IWORK( J ) + 1 + IWOFF, &
                       IWORK( J+IN ) + IWOFF
                  W( JE ) = TMP1
                  IBLOCK( JE ) = IB
   50          CONTINUE
   60       CONTINUE
!
            M = M + IM
         END IF
   70 CONTINUE
!
!     If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
!     If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
!
      IF( IRANGE.EQ.3 ) THEN
         IM = 0
         IDISCL = IL - 1 - NWL
         IDISCU = NWU - IU
!
         IF( IDISCL.GT.0 .OR. IDISCU.GT.0 ) THEN
            DO 80 JE = 1, M
               IF( W( JE ).LE.WLU .AND. IDISCL.GT.0 ) THEN
                  IDISCL = IDISCL - 1
               ELSE IF( W( JE ).GE.WUL .AND. IDISCU.GT.0 ) THEN
                  IDISCU = IDISCU - 1
               ELSE
                  IM = IM + 1
                  W( IM ) = W( JE )
                  IBLOCK( IM ) = IBLOCK( JE )
               END IF
   80       CONTINUE
            M = IM
         END IF
         IF( IDISCL.GT.0 .OR. IDISCU.GT.0 ) THEN
!
!           Code to deal with effects of bad arithmetic:
!           Some low eigenvalues to be discarded are not in (WL,WLU],
!           or high eigenvalues to be discarded are not in (WUL,WU]
!           so just kill off the smallest IDISCL/largest IDISCU
!           eigenvalues, by simply finding the smallest/largest
!           eigenvalue(s).
!
!           (If N(w) is monotone non-decreasing, this should never
!               happen.)
!
            IF( IDISCL.GT.0 ) THEN
               WKILL = WU
               DO 100 JDISC = 1, IDISCL
                  IW = 0
                  DO 90 JE = 1, M
                     IF( IBLOCK( JE ).NE.0 .AND. &
                         ( W( JE ).LT.WKILL .OR. IW.EQ.0 ) ) THEN
                        IW = JE
                        WKILL = W( JE )
                     END IF
   90             CONTINUE
                  IBLOCK( IW ) = 0
  100          CONTINUE
            END IF
            IF( IDISCU.GT.0 ) THEN
!
               WKILL = WL
               DO 120 JDISC = 1, IDISCU
                  IW = 0
                  DO 110 JE = 1, M
                     IF( IBLOCK( JE ).NE.0 .AND. &
                         ( W( JE ).GT.WKILL .OR. IW.EQ.0 ) ) THEN
                        IW = JE
                        WKILL = W( JE )
                     END IF
  110             CONTINUE
                  IBLOCK( IW ) = 0
  120          CONTINUE
            END IF
            IM = 0
            DO 130 JE = 1, M
               IF( IBLOCK( JE ).NE.0 ) THEN
                  IM = IM + 1
                  W( IM ) = W( JE )
                  IBLOCK( IM ) = IBLOCK( JE )
               END IF
  130       CONTINUE
            M = IM
         END IF
         IF( IDISCL.LT.0 .OR. IDISCU.LT.0 ) THEN
            TOOFEW = .TRUE.
         END IF
      END IF
!
!     If ORDER='B', do nothing -- the eigenvalues are already sorted
!        by block.
!     If ORDER='E', sort the eigenvalues from smallest to largest
!
      IF( IORDER.EQ.1 .AND. NSPLIT.GT.1 ) THEN
         DO 150 JE = 1, M - 1
            IE = 0
            TMP1 = W( JE )
            DO 140 J = JE + 1, M
               IF( W( J ).LT.TMP1 ) THEN
                  IE = J
                  TMP1 = W( J )
               END IF
  140       CONTINUE
!
            IF( IE.NE.0 ) THEN
               ITMP1 = IBLOCK( IE )
               W( IE ) = W( JE )
               IBLOCK( IE ) = IBLOCK( JE )
               W( JE ) = TMP1
               IBLOCK( JE ) = ITMP1
            END IF
  150    CONTINUE
      END IF
!
      INFO = 0
      IF( NCNVRG ) &
         INFO = INFO + 1
      IF( TOOFEW ) &
         INFO = INFO + 2
      RETURN
!
!     End of DSTEBZ
!
      END
      SUBROUTINE DSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK, &
                         IWORK, IFAIL, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDZ, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IBLOCK( * ), IFAIL( * ), ISPLIT( * ), &
                         IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  DSTEIN computes the eigenvectors of a real symmetric tridiagonal
!  matrix T corresponding to specified eigenvalues, using inverse
!  iteration.
!
!  The maximum number of iterations allowed for each eigenvector is
!  specified by an internal parameter MAXITS (currently set to 5).
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The n diagonal elements of the tridiagonal matrix T.
!
!  E       (input) DOUBLE PRECISION array, dimension (N)
!          The (n-1) subdiagonal elements of the tridiagonal matrix
!          T, in elements 1 to N-1.  E(N) need not be set.
!
!  M       (input) INTEGER
!          The number of eigenvectors to be found.  0 <= M <= N.
!
!  W       (input) DOUBLE PRECISION array, dimension (N)
!          The first M elements of W contain the eigenvalues for
!          which eigenvectors are to be computed.  The eigenvalues
!          should be grouped by split-off block and ordered from
!          smallest to largest within the block.  ( The output array
!          W from DSTEBZ with ORDER = 'B' is expected here. )
!
!  IBLOCK  (input) INTEGER array, dimension (N)
!          The submatrix indices associated with the corresponding
!          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to
!          the first submatrix from the top, =2 if W(i) belongs to
!          the second submatrix, etc.  ( The output array IBLOCK
!          from DSTEBZ is expected here. )
!
!  ISPLIT  (input) INTEGER array, dimension (N)
!          The splitting points, at which T breaks up into submatrices.
!          The first submatrix consists of rows/columns 1 to
!          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1
!          through ISPLIT( 2 ), etc.
!          ( The output array ISPLIT from DSTEBZ is expected here. )
!
!  Z       (output) DOUBLE PRECISION array, dimension (LDZ, M)
!          The computed eigenvectors.  The eigenvector associated
!          with the eigenvalue W(i) is stored in the i-th column of
!          Z.  Any vector which fails to converge is set to its current
!          iterate after MAXITS iterations.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= max(1,N).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (5*N)
!
!  IWORK   (workspace) INTEGER array, dimension (N)
!
!  IFAIL   (output) INTEGER array, dimension (M)
!          On normal exit, all elements of IFAIL are zero.
!          If one or more eigenvectors fail to converge after
!          MAXITS iterations, then their indices are stored in
!          array IFAIL.
!
!  INFO    (output) INTEGER
!          = 0: successful exit.
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, then i eigenvectors failed to converge
!               in MAXITS iterations.  Their indices are stored in
!               array IFAIL.
!
!  Internal Parameters
!  ===================
!
!  MAXITS  INTEGER, default = 5
!          The maximum number of iterations performed.
!
!  EXTRA   INTEGER, default = 2
!          The number of iterations performed after norm growth
!          criterion is satisfied, should be at least 1.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TEN, ODM3, ODM1
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TEN = 1.0D+1, &
                         ODM3 = 1.0D-3, ODM1 = 1.0D-1 )
      INTEGER            MAXITS, EXTRA
      PARAMETER          ( MAXITS = 5, EXTRA = 2 )
!     ..
!     .. Local Scalars ..
      INTEGER            B1, BLKSIZ, BN, GPIND, I, IINFO, INDRV1, &
                         INDRV2, INDRV3, INDRV4, INDRV5, ITS, J, J1, &
                         JBLK, JMAX, NBLK, NRMCHK
      DOUBLE PRECISION   DTPCRT, EPS, EPS1, NRM, ONENRM, ORTOL, PERTOL, &
                         SCL, SEP, TOL, XJ, XJM, ZTR
!     ..
!     .. Local Arrays ..
      INTEGER            ISEED( 4 )
!     ..
!     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DASUM, DDOT, DLAMCH, DNRM2
      EXTERNAL           IDAMAX, DASUM, DDOT, DLAMCH, DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DLAGTF, DLAGTS, DLARNV, DSCAL, &
                         XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      DO 10 I = 1, M
         IFAIL( I ) = 0
   10 CONTINUE
!
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 .OR. M.GT.N ) THEN
         INFO = -4
      ELSE IF( LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE
         DO 20 J = 2, M
            IF( IBLOCK( J ).LT.IBLOCK( J-1 ) ) THEN
               INFO = -6
               GO TO 30
            END IF
            IF( IBLOCK( J ).EQ.IBLOCK( J-1 ) .AND. W( J ).LT.W( J-1 ) ) &
                 THEN
               INFO = -5
               GO TO 30
            END IF
   20    CONTINUE
   30    CONTINUE
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEIN', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. M.EQ.0 ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
         Z( 1, 1 ) = ONE
         RETURN
      END IF
!
!     Get machine constants.
!
      EPS = DLAMCH( 'Precision' )
!
!     Initialize seed for random number generator DLARNV.
!
      DO 40 I = 1, 4
         ISEED( I ) = 1
   40 CONTINUE
!
!     Initialize pointers.
!
      INDRV1 = 0
      INDRV2 = INDRV1 + N
      INDRV3 = INDRV2 + N
      INDRV4 = INDRV3 + N
      INDRV5 = INDRV4 + N
!
!     Compute eigenvectors of matrix blocks.
!
      J1 = 1
      DO 160 NBLK = 1, IBLOCK( M )
!
!        Find starting and ending indices of block nblk.
!
         IF( NBLK.EQ.1 ) THEN
            B1 = 1
         ELSE
            B1 = ISPLIT( NBLK-1 ) + 1
         END IF
         BN = ISPLIT( NBLK )
         BLKSIZ = BN - B1 + 1
         IF( BLKSIZ.EQ.1 ) &
            GO TO 60
         GPIND = B1
!
!        Compute reorthogonalization criterion and stopping criterion.
!
         ONENRM = ABS( D( B1 ) ) + ABS( E( B1 ) )
         ONENRM = MAX( ONENRM, ABS( D( BN ) )+ABS( E( BN-1 ) ) )
         DO 50 I = B1 + 1, BN - 1
            ONENRM = MAX( ONENRM, ABS( D( I ) )+ABS( E( I-1 ) )+ &
                     ABS( E( I ) ) )
   50    CONTINUE
         ORTOL = ODM3*ONENRM
!
         DTPCRT = SQRT( ODM1 / BLKSIZ )
!
!        Loop through eigenvalues of block nblk.
!
   60    CONTINUE
         JBLK = 0
         DO 150 J = J1, M
            IF( IBLOCK( J ).NE.NBLK ) THEN
               J1 = J
               GO TO 160
            END IF
            JBLK = JBLK + 1
            XJ = W( J )
!
!           Skip all the work if the block size is one.
!
            IF( BLKSIZ.EQ.1 ) THEN
               WORK( INDRV1+1 ) = ONE
               GO TO 120
            END IF
!
!           If eigenvalues j and j-1 are too close, add a relatively
!           small perturbation.
!
            IF( JBLK.GT.1 ) THEN
               EPS1 = ABS( EPS*XJ )
               PERTOL = TEN*EPS1
               SEP = XJ - XJM
               IF( SEP.LT.PERTOL ) &
                  XJ = XJM + PERTOL
            END IF
!
            ITS = 0
            NRMCHK = 0
!
!           Get random starting vector.
!
            CALL DLARNV( 2, ISEED, BLKSIZ, WORK( INDRV1+1 ) )
!
!           Copy the matrix T so it won't be destroyed in factorization.
!
            CALL DCOPY( BLKSIZ, D( B1 ), 1, WORK( INDRV4+1 ), 1 )
            CALL DCOPY( BLKSIZ-1, E( B1 ), 1, WORK( INDRV2+2 ), 1 )
            CALL DCOPY( BLKSIZ-1, E( B1 ), 1, WORK( INDRV3+1 ), 1 )
!
!           Compute LU factors with partial pivoting  ( PT = LU )
!
            TOL = ZERO
            CALL DLAGTF( BLKSIZ, WORK( INDRV4+1 ), XJ, WORK( INDRV2+2 ), &
                         WORK( INDRV3+1 ), TOL, WORK( INDRV5+1 ), IWORK, &
                         IINFO )
!
!           Update iteration count.
!
   70       CONTINUE
            ITS = ITS + 1
            IF( ITS.GT.MAXITS ) &
               GO TO 100
!
!           Normalize and scale the righthand side vector Pb.
!
            SCL = BLKSIZ*ONENRM*MAX( EPS, &
                  ABS( WORK( INDRV4+BLKSIZ ) ) ) / &
                  DASUM( BLKSIZ, WORK( INDRV1+1 ), 1 )
            CALL DSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
!
!           Solve the system LU = Pb.
!
            CALL DLAGTS( -1, BLKSIZ, WORK( INDRV4+1 ), WORK( INDRV2+2 ), &
                         WORK( INDRV3+1 ), WORK( INDRV5+1 ), IWORK, &
                         WORK( INDRV1+1 ), TOL, IINFO )
!
!           Reorthogonalize by modified Gram-Schmidt if eigenvalues are
!           close enough.
!
            IF( JBLK.EQ.1 ) &
               GO TO 90
            IF( ABS( XJ-XJM ).GT.ORTOL ) &
               GPIND = J
            IF( GPIND.NE.J ) THEN
               DO 80 I = GPIND, J - 1
                  ZTR = -DDOT( BLKSIZ, WORK( INDRV1+1 ), 1, Z( B1, I ), &
                        1 )
                  CALL DAXPY( BLKSIZ, ZTR, Z( B1, I ), 1, &
                              WORK( INDRV1+1 ), 1 )
   80          CONTINUE
            END IF
!
!           Check the infinity norm of the iterate.
!
   90       CONTINUE
            JMAX = IDAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
            NRM = ABS( WORK( INDRV1+JMAX ) )
!
!           Continue for additional iterations after norm reaches
!           stopping criterion.
!
            IF( NRM.LT.DTPCRT ) &
               GO TO 70
            NRMCHK = NRMCHK + 1
            IF( NRMCHK.LT.EXTRA+1 ) &
               GO TO 70
!
            GO TO 110
!
!           If stopping criterion was not satisfied, update info and
!           store eigenvector number in array ifail.
!
  100       CONTINUE
            INFO = INFO + 1
            IFAIL( INFO ) = J
!
!           Accept iterate as jth eigenvector.
!
  110       CONTINUE
            SCL = ONE / DNRM2( BLKSIZ, WORK( INDRV1+1 ), 1 )
            JMAX = IDAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
            IF( WORK( INDRV1+JMAX ).LT.ZERO ) &
               SCL = -SCL
            CALL DSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
  120       CONTINUE
            DO 130 I = 1, N
               Z( I, J ) = ZERO
  130       CONTINUE
            DO 140 I = 1, BLKSIZ
               Z( B1+I-1, J ) = WORK( INDRV1+I )
  140       CONTINUE
!
!           Save the shift to check eigenvalue spacing at next
!           iteration.
!
            XJM = XJ
!
  150    CONTINUE
  160 CONTINUE
!
      RETURN
!
!     End of DSTEIN
!
      END
      SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  DSTEQR computes all eigenvalues and, optionally, eigenvectors of a
!  symmetric tridiagonal matrix using the implicit QL or QR method.
!  The eigenvectors of a full or band symmetric matrix can also be found
!  if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to
!  tridiagonal form.
!
!  Arguments
!  =========
!
!  COMPZ   (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only.
!          = 'V':  Compute eigenvalues and eigenvectors of the original
!                  symmetric matrix.  On entry, Z must contain the
!                  orthogonal matrix used to reduce the original matrix
!                  to tridiagonal form.
!          = 'I':  Compute eigenvalues and eigenvectors of the
!                  tridiagonal matrix.  Z is initialized to the identity
!                  matrix.
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the diagonal elements of the tridiagonal matrix.
!          On exit, if INFO = 0, the eigenvalues in ascending order.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix.
!          On exit, E has been destroyed.
!
!  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
!          On entry, if  COMPZ = 'V', then Z contains the orthogonal
!          matrix used in the reduction to tridiagonal form.
!          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
!          orthonormal eigenvectors of the original symmetric matrix,
!          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!          of the symmetric tridiagonal matrix.
!          If COMPZ = 'N', then Z is not referenced.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          eigenvectors are desired, then  LDZ >= max(1,N).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
!          If COMPZ = 'N', then WORK is not referenced.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  the algorithm has failed to find all the eigenvalues in
!                a total of 30*N iterations; if INFO = i, then i
!                elements of E have not converged to zero; on exit, D
!                and E contain the elements of a symmetric tridiagonal
!                matrix which is orthogonally similar to the original
!                matrix.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, &
                         THREE = 3.0D0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND, &
                         LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1, &
                         NM1, NMAXIT
      DOUBLE PRECISION   ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2, &
                         S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASET, DLASR, &
                         DLASRT, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1, &
               N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEQR', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.EQ.2 ) &
            Z( 1, 1 ) = ONE
         RETURN
      END IF
!
!     Determine the unit roundoff and over/underflow thresholds.
!
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
!
!     Compute the eigenvalues and eigenvectors of the tridiagonal
!     matrix.
!
      IF( ICOMPZ.EQ.2 ) &
         CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
!
      NMAXIT = N*MAXIT
      JTOT = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      L1 = 1
      NM1 = N - 1
!
   10 CONTINUE
      IF( L1.GT.N ) &
         GO TO 160
      IF( L1.GT.1 ) &
         E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO ) &
               GO TO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+ &
                1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
!
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L ) &
         GO TO 10
!
!     Scale submatrix in rows and columns L to LEND
!
      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO ) &
         GO TO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N, &
                      INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N, &
                      INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N, &
                      INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N, &
                      INFO )
      END IF
!
!     Choose between QL and QR iteration
!
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
!
      IF( LEND.GT.L ) THEN
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+ &
                   SAFMIN )GO TO 60
   50       CONTINUE
         END IF
!
         M = LEND
!
   60    CONTINUE
         IF( M.LT.LEND ) &
            E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L ) &
            GO TO 80
!
!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.
!
         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL DLASR( 'R', 'V', 'B', N, 2, WORK( L ), &
                           WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND ) &
               GO TO 40
            GO TO 140
         END IF
!
         IF( JTOT.EQ.NMAXIT ) &
            GO TO 140
         JTOT = JTOT + 1
!
!        Form shift.
!
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
!
         S = ONE
         C = ONE
         P = ZERO
!
!        Inner loop
!
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M-1 ) &
               E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF
!
   70    CONTINUE
!
!        If eigenvectors are desired, then apply saved rotations.
!
         IF( ICOMPZ.GT.0 ) THEN
            MM = M - L + 1
            CALL DLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ), &
                        Z( 1, L ), LDZ )
         END IF
!
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
!
!        Eigenvalue found.
!
   80    CONTINUE
         D( L ) = P
!
         L = L + 1
         IF( L.LE.LEND ) &
            GO TO 40
         GO TO 140
!
      ELSE
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
   90    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+ &
                   SAFMIN )GO TO 110
  100       CONTINUE
         END IF
!
         M = LEND
!
  110    CONTINUE
         IF( M.GT.LEND ) &
            E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L ) &
            GO TO 130
!
!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.
!
         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL DLASR( 'R', 'V', 'F', N, 2, WORK( M ), &
                           WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND ) &
               GO TO 90
            GO TO 140
         END IF
!
         IF( JTOT.EQ.NMAXIT ) &
            GO TO 140
         JTOT = JTOT + 1
!
!        Form shift.
!
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
!
         S = ONE
         C = ONE
         P = ZERO
!
!        Inner loop
!
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M ) &
               E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF
!
  120    CONTINUE
!
!        If eigenvectors are desired, then apply saved rotations.
!
         IF( ICOMPZ.GT.0 ) THEN
            MM = L - M + 1
            CALL DLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ), &
                        Z( 1, M ), LDZ )
         END IF
!
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
!
!        Eigenvalue found.
!
  130    CONTINUE
         D( L ) = P
!
         L = L - 1
         IF( L.GE.LEND ) &
            GO TO 90
         GO TO 140
!
      END IF
!
!     Undo scaling if necessary
!
  140 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1, &
                      D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ), &
                      N, INFO )
      ELSE IF( ISCALE.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1, &
                      D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ), &
                      N, INFO )
      END IF
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
      IF( JTOT.LT.NMAXIT ) &
         GO TO 10
      DO 150 I = 1, N - 1
         IF( E( I ).NE.ZERO ) &
            INFO = INFO + 1
  150 CONTINUE
      GO TO 190
!
!     Order eigenvalues and eigenvectors.
!
  160 CONTINUE
      IF( ICOMPZ.EQ.0 ) THEN
!
!        Use Quick Sort
!
         CALL DLASRT( 'I', N, D, INFO )
!
      ELSE
!
!        Use Selection Sort to minimize swaps of eigenvectors
!
         DO 180 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 170 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
  170       CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL DSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
  180    CONTINUE
      END IF
!
  190 CONTINUE
      RETURN
!
!     End of DSTEQR
!
      END
      SUBROUTINE DSTERF( N, D, E, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
!     ..
!
!  Purpose
!  =======
!
!  DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
!  using the Pal-Walker-Kahan variant of the QL or QR algorithm.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the n diagonal elements of the tridiagonal matrix.
!          On exit, if INFO = 0, the eigenvalues in ascending order.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix.
!          On exit, E has been destroyed.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  the algorithm failed to find all of the eigenvalues in
!                a total of 30*N iterations; if INFO = i, then i
!                elements of E have not converged to zero.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, &
                         THREE = 3.0D0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, ISCALE, JTOT, L, L1, LEND, LENDM1, LENDP1, &
                         LENDSV, LM1, LSV, M, MM1, NM1, NMAXIT
      DOUBLE PRECISION   ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC, &
                         OLDGAM, P, R, RT1, RT2, RTE, S, SAFMAX, SAFMIN, &
                         SIGMA, SSFMAX, SSFMIN, TST
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           DLAMCH, DLANST, DLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAE2, DLASCL, DLASRT, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
!     Quick return if possible
!
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DSTERF', -INFO )
         RETURN
      END IF
      IF( N.LE.1 ) &
         RETURN
!
!     Determine the unit roundoff for this environment.
!
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
!
!     Compute the eigenvalues of the tridiagonal matrix.
!
      NMAXIT = N*MAXIT
      SIGMA = ZERO
      JTOT = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      L1 = 1
      NM1 = N - 1
!
   10 CONTINUE
      IF( L1.GT.N ) &
         GO TO 170
      IF( L1.GT.1 ) &
         E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO ) &
               GO TO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+ &
                1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
!
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L ) &
         GO TO 10
!
!     Scale submatrix in rows and columns L to LEND
!
      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N, &
                      INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N, &
                      INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N, &
                      INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N, &
                      INFO )
      END IF
!
      DO 40 I = L, LEND - 1
         E( I ) = E( I )**2
   40 CONTINUE
!
!     Choose between QL and QR iteration
!
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
!
      IF( LEND.GE.L ) THEN
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
   50    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 60 M = L, LENDM1
               TST = ABS( E( M ) )
               IF( TST.LE.EPS2*ABS( D( M )*D( M+1 ) ) ) &
                  GO TO 70
   60       CONTINUE
         END IF
!
         M = LEND
!
   70    CONTINUE
         IF( M.LT.LEND ) &
            E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L ) &
            GO TO 90
!
!        If remaining matrix is 2 by 2, use DLAE2 to compute its
!        eigenvalues.
!
         IF( M.EQ.L+1 ) THEN
            RTE = SQRT( E( L ) )
            CALL DLAE2( D( L ), RTE, D( L+1 ), RT1, RT2 )
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND ) &
               GO TO 50
            GO TO 150
         END IF
!
         IF( JTOT.EQ.NMAXIT ) &
            GO TO 150
         JTOT = JTOT + 1
!
!        Form shift.
!
         RTE = SQRT( E( L ) )
         SIGMA = ( D( L+1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
!
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
!
!        Inner loop
!
         MM1 = M - 1
         DO 80 I = MM1, L, -1
            BB = E( I )
            R = P + BB
            IF( I.NE.M-1 ) &
               E( I+1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I+1 ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
   80    CONTINUE
!
         E( L ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 50
!
!        Eigenvalue found.
!
   90    CONTINUE
         D( L ) = P
!
         L = L + 1
         IF( L.LE.LEND ) &
            GO TO 50
         GO TO 150
!
      ELSE
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
  100    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 110 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )
               IF( TST.LE.EPS2*ABS( D( M )*D( M-1 ) ) ) &
                  GO TO 120
  110       CONTINUE
         END IF
!
         M = LEND
!
  120    CONTINUE
         IF( M.GT.LEND ) &
            E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L ) &
            GO TO 140
!
!        If remaining matrix is 2 by 2, use DLAE2 to compute its
!        eigenvalues.
!
         IF( M.EQ.L-1 ) THEN
            RTE = SQRT( E( L-1 ) )
            CALL DLAE2( D( L ), RTE, D( L-1 ), RT1, RT2 )
            D( L ) = RT1
            D( L-1 ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND ) &
               GO TO 100
            GO TO 150
         END IF
!
         IF( JTOT.EQ.NMAXIT ) &
            GO TO 150
         JTOT = JTOT + 1
!
!        Form shift.
!
         RTE = SQRT( E( L-1 ) )
         SIGMA = ( D( L-1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
!
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
!
!        Inner loop
!
         LM1 = L - 1
         DO 130 I = M, LM1
            BB = E( I )
            R = P + BB
            IF( I.NE.M ) &
               E( I-1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I+1 )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
  130    CONTINUE
!
         E( LM1 ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 100
!
!        Eigenvalue found.
!
  140    CONTINUE
         D( L ) = P
!
         L = L - 1
         IF( L.GE.LEND ) &
            GO TO 100
         GO TO 150
!
      END IF
!
!     Undo scaling if necessary
!
  150 CONTINUE
      IF( ISCALE.EQ.1 ) &
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1, &
                      D( LSV ), N, INFO )
      IF( ISCALE.EQ.2 ) &
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1, &
                      D( LSV ), N, INFO )
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
      IF( JTOT.EQ.NMAXIT ) THEN
         DO 160 I = 1, N - 1
            IF( E( I ).NE.ZERO ) &
               INFO = INFO + 1
  160    CONTINUE
         RETURN
      END IF
      GO TO 10
!
!     Sort eigenvalues in increasing order.
!
  170 CONTINUE
      CALL DLASRT( 'I', N, D, INFO )
!
      RETURN
!
!     End of DSTERF
!
      END
      SUBROUTINE DSYGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DSYGS2 reduces a real symmetric-definite generalized eigenproblem
!  to standard form.
!
!  If ITYPE = 1, the problem is A*x = lambda*B*x,
!  and A is overwritten by inv(U')*A*inv(U) or inv(L)*A*inv(L')
!
!  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
!  B*A*x = lambda*x, and A is overwritten by U*A*U` or L'*A*L.
!
!  B must have been previously factorized as U'*U or L*L' by DPOTRF.
!
!  Arguments
!  =========
!
!  ITYPE   (input) INTEGER
!          = 1: compute inv(U')*A*inv(U) or inv(L)*A*inv(L');
!          = 2 or 3: compute U*A*U' or L'*A*L.
!
!  UPLO    (input) CHARACTER
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored, and how B has been factorized.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrices A and B.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n by n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n by n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the transformed matrix, stored in the
!          same format as A.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  B       (input) DOUBLE PRECISION array, dimension (LDB,N)
!          The triangular factor from the Cholesky factorization of B,
!          as returned by DPOTRF.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, HALF
      PARAMETER          ( ONE = 1.0D0, HALF = 0.5D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K
      DOUBLE PRECISION   AKK, BKK, CT
!     ..
!     .. External Subroutines ..
      EXTERNAL           DAXPY, DSCAL, DSYR2, DTRMV, DTRSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYGS2', -INFO )
         RETURN
      END IF
!
      IF( ITYPE.EQ.1 ) THEN
         IF( UPPER ) THEN
!
!           Compute inv(U')*A*inv(U)
!
            DO 10 K = 1, N
!
!              Update the upper triangle of A(k:n,k:n)
!
               AKK = A( K, K )
               BKK = B( K, K )
               AKK = AKK / BKK**2
               A( K, K ) = AKK
               IF( K.LT.N ) THEN
                  CALL DSCAL( N-K, ONE / BKK, A( K, K+1 ), LDA )
                  CT = -HALF*AKK
                  CALL DAXPY( N-K, CT, B( K, K+1 ), LDB, A( K, K+1 ), &
                              LDA )
                  CALL DSYR2( UPLO, N-K, -ONE, A( K, K+1 ), LDA, &
                              B( K, K+1 ), LDB, A( K+1, K+1 ), LDA )
                  CALL DAXPY( N-K, CT, B( K, K+1 ), LDB, A( K, K+1 ), &
                              LDA )
                  CALL DTRSV( UPLO, 'Transpose', 'Non-unit', N-K, &
                              B( K+1, K+1 ), LDB, A( K, K+1 ), LDA )
               END IF
   10       CONTINUE
         ELSE
!
!           Compute inv(L)*A*inv(L')
!
            DO 20 K = 1, N
!
!              Update the lower triangle of A(k:n,k:n)
!
               AKK = A( K, K )
               BKK = B( K, K )
               AKK = AKK / BKK**2
               A( K, K ) = AKK
               IF( K.LT.N ) THEN
                  CALL DSCAL( N-K, ONE / BKK, A( K+1, K ), 1 )
                  CT = -HALF*AKK
                  CALL DAXPY( N-K, CT, B( K+1, K ), 1, A( K+1, K ), 1 )
                  CALL DSYR2( UPLO, N-K, -ONE, A( K+1, K ), 1, &
                              B( K+1, K ), 1, A( K+1, K+1 ), LDA )
                  CALL DAXPY( N-K, CT, B( K+1, K ), 1, A( K+1, K ), 1 )
                  CALL DTRSV( UPLO, 'No transpose', 'Non-unit', N-K, &
                              B( K+1, K+1 ), LDB, A( K+1, K ), 1 )
               END IF
   20       CONTINUE
         END IF
      ELSE
         IF( UPPER ) THEN
!
!           Compute U*A*U'
!
            DO 30 K = 1, N
!
!              Update the upper triangle of A(1:k,1:k)
!
               AKK = A( K, K )
               BKK = B( K, K )
               CALL DTRMV( UPLO, 'No transpose', 'Non-unit', K-1, B, &
                           LDB, A( 1, K ), 1 )
               CT = HALF*AKK
               CALL DAXPY( K-1, CT, B( 1, K ), 1, A( 1, K ), 1 )
               CALL DSYR2( UPLO, K-1, ONE, A( 1, K ), 1, B( 1, K ), 1, &
                           A, LDA )
               CALL DAXPY( K-1, CT, B( 1, K ), 1, A( 1, K ), 1 )
               CALL DSCAL( K-1, BKK, A( 1, K ), 1 )
               A( K, K ) = AKK*BKK**2
   30       CONTINUE
         ELSE
!
!           Compute L'*A*L
!
            DO 40 K = 1, N
!
!              Update the lower triangle of A(1:k,1:k)
!
               AKK = A( K, K )
               BKK = B( K, K )
               CALL DTRMV( UPLO, 'Transpose', 'Non-unit', K-1, B, LDB, &
                           A( K, 1 ), LDA )
               CT = HALF*AKK
               CALL DAXPY( K-1, CT, B( K, 1 ), LDB, A( K, 1 ), LDA )
               CALL DSYR2( UPLO, K-1, ONE, A( K, 1 ), LDA, B( K, 1 ), &
                           LDB, A, LDA )
               CALL DAXPY( K-1, CT, B( K, 1 ), LDB, A( K, 1 ), LDA )
               CALL DSCAL( K-1, BKK, A( K, 1 ), LDA )
               A( K, K ) = AKK*BKK**2
   40       CONTINUE
         END IF
      END IF
      RETURN
!
!     End of DSYGS2
!
      END
      SUBROUTINE DSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DSYGST reduces a real symmetric-definite generalized eigenproblem
!  to standard form.
!
!  If ITYPE = 1, the problem is A*x = lambda*B*x,
!  and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
!
!  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
!  B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
!
!  B must have been previously factorized as U**T*U or L*L**T by DPOTRF.
!
!  Arguments
!  =========
!
!  ITYPE   (input) INTEGER
!          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);
!          = 2 or 3: compute U*A*U**T or L**T*A*L.
!
!  UPLO    (input) CHARACTER
!          = 'U':  Upper triangle of A is stored and B is factored as
!                  U**T*U;
!          = 'L':  Lower triangle of A is stored and B is factored as
!                  L*L**T.
!
!  N       (input) INTEGER
!          The order of the matrices A and B.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the transformed matrix, stored in the
!          same format as A.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  B       (input) DOUBLE PRECISION array, dimension (LDB,N)
!          The triangular factor from the Cholesky factorization of B,
!          as returned by DPOTRF.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, HALF
      PARAMETER          ( ONE = 1.0D0, HALF = 0.5D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K, KB, NB
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSYGS2, DSYMM, DSYR2K, DTRMM, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYGST', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'DSYGST', UPLO, N, -1, -1, -1 )
!
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
!
!        Use unblocked code
!
         CALL DSYGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
      ELSE
!
!        Use blocked code
!
         IF( ITYPE.EQ.1 ) THEN
            IF( UPPER ) THEN
!
!              Compute inv(U')*A*inv(U)
!
               DO 10 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
!
!                 Update the upper triangle of A(k:n,k:n)
!
                  CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA, &
                               B( K, K ), LDB, INFO )
                  IF( K+KB.LE.N ) THEN
                     CALL DTRSM( 'Left', UPLO, 'Transpose', 'Non-unit', &
                                 KB, N-K-KB+1, ONE, B( K, K ), LDB, &
                                 A( K, K+KB ), LDA )
                     CALL DSYMM( 'Left', UPLO, KB, N-K-KB+1, -HALF, &
                                 A( K, K ), LDA, B( K, K+KB ), LDB, ONE, &
                                 A( K, K+KB ), LDA )
                     CALL DSYR2K( UPLO, 'Transpose', N-K-KB+1, KB, -ONE, &
                                  A( K, K+KB ), LDA, B( K, K+KB ), LDB, &
                                  ONE, A( K+KB, K+KB ), LDA )
                     CALL DSYMM( 'Left', UPLO, KB, N-K-KB+1, -HALF, &
                                 A( K, K ), LDA, B( K, K+KB ), LDB, ONE, &
                                 A( K, K+KB ), LDA )
                     CALL DTRSM( 'Right', UPLO, 'No transpose', &
                                 'Non-unit', KB, N-K-KB+1, ONE, &
                                 B( K+KB, K+KB ), LDB, A( K, K+KB ), &
                                 LDA )
                  END IF
   10          CONTINUE
            ELSE
!
!              Compute inv(L)*A*inv(L')
!
               DO 20 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
!
!                 Update the lower triangle of A(k:n,k:n)
!
                  CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA, &
                               B( K, K ), LDB, INFO )
                  IF( K+KB.LE.N ) THEN
                     CALL DTRSM( 'Right', UPLO, 'Transpose', 'Non-unit', &
                                 N-K-KB+1, KB, ONE, B( K, K ), LDB, &
                                 A( K+KB, K ), LDA )
                     CALL DSYMM( 'Right', UPLO, N-K-KB+1, KB, -HALF, &
                                 A( K, K ), LDA, B( K+KB, K ), LDB, ONE, &
                                 A( K+KB, K ), LDA )
                     CALL DSYR2K( UPLO, 'No transpose', N-K-KB+1, KB, &
                                  -ONE, A( K+KB, K ), LDA, B( K+KB, K ), &
                                  LDB, ONE, A( K+KB, K+KB ), LDA )
                     CALL DSYMM( 'Right', UPLO, N-K-KB+1, KB, -HALF, &
                                 A( K, K ), LDA, B( K+KB, K ), LDB, ONE, &
                                 A( K+KB, K ), LDA )
                     CALL DTRSM( 'Left', UPLO, 'No transpose', &
                                 'Non-unit', N-K-KB+1, KB, ONE, &
                                 B( K+KB, K+KB ), LDB, A( K+KB, K ), &
                                 LDA )
                  END IF
   20          CONTINUE
            END IF
         ELSE
            IF( UPPER ) THEN
!
!              Compute U*A*U'
!
               DO 30 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
!
!                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1)
!
                  CALL DTRMM( 'Left', UPLO, 'No transpose', 'Non-unit', &
                              K-1, KB, ONE, B, LDB, A( 1, K ), LDA )
                  CALL DSYMM( 'Right', UPLO, K-1, KB, HALF, A( K, K ), &
                              LDA, B( 1, K ), LDB, ONE, A( 1, K ), LDA )
                  CALL DSYR2K( UPLO, 'No transpose', K-1, KB, ONE, &
                               A( 1, K ), LDA, B( 1, K ), LDB, ONE, A, &
                               LDA )
                  CALL DSYMM( 'Right', UPLO, K-1, KB, HALF, A( K, K ), &
                              LDA, B( 1, K ), LDB, ONE, A( 1, K ), LDA )
                  CALL DTRMM( 'Right', UPLO, 'Transpose', 'Non-unit', &
                              K-1, KB, ONE, B( K, K ), LDB, A( 1, K ), &
                              LDA )
                  CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA, &
                               B( K, K ), LDB, INFO )
   30          CONTINUE
            ELSE
!
!              Compute L'*A*L
!
               DO 40 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
!
!                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1)
!
                  CALL DTRMM( 'Right', UPLO, 'No transpose', 'Non-unit', &
                              KB, K-1, ONE, B, LDB, A( K, 1 ), LDA )
                  CALL DSYMM( 'Left', UPLO, KB, K-1, HALF, A( K, K ), &
                              LDA, B( K, 1 ), LDB, ONE, A( K, 1 ), LDA )
                  CALL DSYR2K( UPLO, 'Transpose', K-1, KB, ONE, &
                               A( K, 1 ), LDA, B( K, 1 ), LDB, ONE, A, &
                               LDA )
                  CALL DSYMM( 'Left', UPLO, KB, K-1, HALF, A( K, K ), &
                              LDA, B( K, 1 ), LDB, ONE, A( K, 1 ), LDA )
                  CALL DTRMM( 'Left', UPLO, 'Transpose', 'Non-unit', KB, &
                              K-1, ONE, B( K, K ), LDB, A( K, 1 ), LDA )
                  CALL DSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA, &
                               B( K, K ), LDB, INFO )
   40          CONTINUE
            END IF
         END IF
      END IF
      RETURN
!
!     End of DSYGST
!
      END
      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, &
                       N4 )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. Executable Statements ..
!
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
!
!     Invalid value for ISPEC
!
      ILAENV = -1
      RETURN
!
  100 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!        ASCII character set
!
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 ) &
                  SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!        EBCDIC character set
!
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
             ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
             ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
                   ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
                   ( IC.GE.162 .AND. IC.LE.169 ) ) &
                  SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 ) &
                  SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
!
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) ) &
         RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
!
      GO TO ( 110, 200, 300 ) ISPEC
!
  110 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
                  C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
!
  200 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
             C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
!
  300 CONTINUE
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
             C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
!
  400 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
      ILAENV = 6
      RETURN
!
  500 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ILAENV = 2
      RETURN
!
  600 CONTINUE 
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
!
  700 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
  800 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN
!
!     End of ILAENV
!
      END
      SUBROUTINE ZHEGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZHEGS2 reduces a complex Hermitian-definite generalized
!  eigenproblem to standard form.
!
!  If ITYPE = 1, the problem is A*x = lambda*B*x,
!  and A is overwritten by inv(U')*A*inv(U) or inv(L)*A*inv(L')
!
!  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
!  B*A*x = lambda*x, and A is overwritten by U*A*U` or L'*A*L.
!
!  B must have been previously factorized as U'*U or L*L' by ZPOTRF.
!
!  Arguments
!  =========
!
!  ITYPE   (input) INTEGER
!          = 1: compute inv(U')*A*inv(U) or inv(L)*A*inv(L');
!          = 2 or 3: compute U*A*U' or L'*A*L.
!
!  UPLO    (input) CHARACTER
!          Specifies whether the upper or lower triangular part of the
!          Hermitian matrix A is stored, and how B has been factorized.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrices A and B.  N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!          n by n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n by n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the transformed matrix, stored in the
!          same format as A.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  B       (input) COMPLEX*16 array, dimension (LDB,N)
!          The triangular factor from the Cholesky factorization of B,
!          as returned by ZPOTRF.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, HALF
      PARAMETER          ( ONE = 1.0D+0, HALF = 0.5D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K
      DOUBLE PRECISION   AKK, BKK
      COMPLEX*16         CT
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZAXPY, ZDSCAL, ZHER2, ZLACGV, ZTRMV, &
                         ZTRSV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHEGS2', -INFO )
         RETURN
      END IF
!
      IF( ITYPE.EQ.1 ) THEN
         IF( UPPER ) THEN
!
!           Compute inv(U')*A*inv(U)
!
            DO 10 K = 1, N
!
!              Update the upper triangle of A(k:n,k:n)
!
               AKK = A( K, K )
               BKK = B( K, K )
               AKK = AKK / BKK**2
               A( K, K ) = AKK
               IF( K.LT.N ) THEN
                  CALL ZDSCAL( N-K, ONE / BKK, A( K, K+1 ), LDA )
                  CT = -HALF*AKK
                  CALL ZLACGV( N-K, A( K, K+1 ), LDA )
                  CALL ZLACGV( N-K, B( K, K+1 ), LDB )
                  CALL ZAXPY( N-K, CT, B( K, K+1 ), LDB, A( K, K+1 ), &
                              LDA )
                  CALL ZHER2( UPLO, N-K, -CONE, A( K, K+1 ), LDA, &
                              B( K, K+1 ), LDB, A( K+1, K+1 ), LDA )
                  CALL ZAXPY( N-K, CT, B( K, K+1 ), LDB, A( K, K+1 ), &
                              LDA )
                  CALL ZLACGV( N-K, B( K, K+1 ), LDB )
                  CALL ZTRSV( UPLO, 'Conjugate transpose', 'Non-unit', &
                              N-K, B( K+1, K+1 ), LDB, A( K, K+1 ), &
                              LDA )
                  CALL ZLACGV( N-K, A( K, K+1 ), LDA )
               END IF
   10       CONTINUE
         ELSE
!
!           Compute inv(L)*A*inv(L')
!
            DO 20 K = 1, N
!
!              Update the lower triangle of A(k:n,k:n)
!
               AKK = A( K, K )
               BKK = B( K, K )
               AKK = AKK / BKK**2
               A( K, K ) = AKK
               IF( K.LT.N ) THEN
                  CALL ZDSCAL( N-K, ONE / BKK, A( K+1, K ), 1 )
                  CT = -HALF*AKK
                  CALL ZAXPY( N-K, CT, B( K+1, K ), 1, A( K+1, K ), 1 )
                  CALL ZHER2( UPLO, N-K, -CONE, A( K+1, K ), 1, &
                              B( K+1, K ), 1, A( K+1, K+1 ), LDA )
                  CALL ZAXPY( N-K, CT, B( K+1, K ), 1, A( K+1, K ), 1 )
                  CALL ZTRSV( UPLO, 'No transpose', 'Non-unit', N-K, &
                              B( K+1, K+1 ), LDB, A( K+1, K ), 1 )
               END IF
   20       CONTINUE
         END IF
      ELSE
         IF( UPPER ) THEN
!
!           Compute U*A*U'
!
            DO 30 K = 1, N
!
!              Update the upper triangle of A(1:k,1:k)
!
               AKK = A( K, K )
               BKK = B( K, K )
               CALL ZTRMV( UPLO, 'No transpose', 'Non-unit', K-1, B, &
                           LDB, A( 1, K ), 1 )
               CT = HALF*AKK
               CALL ZAXPY( K-1, CT, B( 1, K ), 1, A( 1, K ), 1 )
               CALL ZHER2( UPLO, K-1, CONE, A( 1, K ), 1, B( 1, K ), 1, &
                           A, LDA )
               CALL ZAXPY( K-1, CT, B( 1, K ), 1, A( 1, K ), 1 )
               CALL ZDSCAL( K-1, BKK, A( 1, K ), 1 )
               A( K, K ) = AKK*BKK**2
   30       CONTINUE
         ELSE
!
!           Compute L'*A*L
!
            DO 40 K = 1, N
!
!              Update the lower triangle of A(1:k,1:k)
!
               AKK = A( K, K )
               BKK = B( K, K )
               CALL ZLACGV( K-1, A( K, 1 ), LDA )
               CALL ZTRMV( UPLO, 'Conjugate transpose', 'Non-unit', K-1, &
                           B, LDB, A( K, 1 ), LDA )
               CT = HALF*AKK
               CALL ZLACGV( K-1, B( K, 1 ), LDB )
               CALL ZAXPY( K-1, CT, B( K, 1 ), LDB, A( K, 1 ), LDA )
               CALL ZHER2( UPLO, K-1, CONE, A( K, 1 ), LDA, B( K, 1 ), &
                           LDB, A, LDA )
               CALL ZAXPY( K-1, CT, B( K, 1 ), LDB, A( K, 1 ), LDA )
               CALL ZLACGV( K-1, B( K, 1 ), LDB )
               CALL ZDSCAL( K-1, BKK, A( K, 1 ), LDA )
               CALL ZLACGV( K-1, A( K, 1 ), LDA )
               A( K, K ) = AKK*BKK**2
   40       CONTINUE
         END IF
      END IF
      RETURN
!
!     End of ZHEGS2
!
      END
      SUBROUTINE ZHEGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZHEGST reduces a complex Hermitian-definite generalized
!  eigenproblem to standard form.
!
!  If ITYPE = 1, the problem is A*x = lambda*B*x,
!  and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
!
!  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
!  B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.
!
!  B must have been previously factorized as U**H*U or L*L**H by ZPOTRF.
!
!  Arguments
!  =========
!
!  ITYPE   (input) INTEGER
!          = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);
!          = 2 or 3: compute U*A*U**H or L**H*A*L.
!
!  UPLO    (input) CHARACTER
!          = 'U':  Upper triangle of A is stored and B is factored as
!                  U**H*U;
!          = 'L':  Lower triangle of A is stored and B is factored as
!                  L*L**H.
!
!  N       (input) INTEGER
!          The order of the matrices A and B.  N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the transformed matrix, stored in the
!          same format as A.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  B       (input) COMPLEX*16 array, dimension (LDB,N)
!          The triangular factor from the Cholesky factorization of B,
!          as returned by ZPOTRF.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      COMPLEX*16         CONE, HALF
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ), &
                         HALF = ( 0.5D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K, KB, NB
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZHEGS2, ZHEMM, ZHER2K, ZTRMM, ZTRSM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHEGST', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'ZHEGST', UPLO, N, -1, -1, -1 )
!
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
!
!        Use unblocked code
!
         CALL ZHEGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
      ELSE
!
!        Use blocked code
!
         IF( ITYPE.EQ.1 ) THEN
            IF( UPPER ) THEN
!
!              Compute inv(U')*A*inv(U)
!
               DO 10 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
!
!                 Update the upper triangle of A(k:n,k:n)
!
                  CALL ZHEGS2( ITYPE, UPLO, KB, A( K, K ), LDA, &
                               B( K, K ), LDB, INFO )
                  IF( K+KB.LE.N ) THEN
                     CALL ZTRSM( 'Left', UPLO, 'Conjugate transpose', &
                                 'Non-unit', KB, N-K-KB+1, CONE, &
                                 B( K, K ), LDB, A( K, K+KB ), LDA )
                     CALL ZHEMM( 'Left', UPLO, KB, N-K-KB+1, -HALF, &
                                 A( K, K ), LDA, B( K, K+KB ), LDB, &
                                 CONE, A( K, K+KB ), LDA )
                     CALL ZHER2K( UPLO, 'Conjugate transpose', N-K-KB+1, &
                                  KB, -CONE, A( K, K+KB ), LDA, &
                                  B( K, K+KB ), LDB, ONE, &
                                  A( K+KB, K+KB ), LDA )
                     CALL ZHEMM( 'Left', UPLO, KB, N-K-KB+1, -HALF, &
                                 A( K, K ), LDA, B( K, K+KB ), LDB, &
                                 CONE, A( K, K+KB ), LDA )
                     CALL ZTRSM( 'Right', UPLO, 'No transpose', &
                                 'Non-unit', KB, N-K-KB+1, CONE, &
                                 B( K+KB, K+KB ), LDB, A( K, K+KB ), &
                                 LDA )
                  END IF
   10          CONTINUE
            ELSE
!
!              Compute inv(L)*A*inv(L')
!
               DO 20 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
!
!                 Update the lower triangle of A(k:n,k:n)
!
                  CALL ZHEGS2( ITYPE, UPLO, KB, A( K, K ), LDA, &
                               B( K, K ), LDB, INFO )
                  IF( K+KB.LE.N ) THEN
                     CALL ZTRSM( 'Right', UPLO, 'Conjugate transpose', &
                                 'Non-unit', N-K-KB+1, KB, CONE, &
                                 B( K, K ), LDB, A( K+KB, K ), LDA )
                     CALL ZHEMM( 'Right', UPLO, N-K-KB+1, KB, -HALF, &
                                 A( K, K ), LDA, B( K+KB, K ), LDB, &
                                 CONE, A( K+KB, K ), LDA )
                     CALL ZHER2K( UPLO, 'No transpose', N-K-KB+1, KB, &
                                  -CONE, A( K+KB, K ), LDA, &
                                  B( K+KB, K ), LDB, ONE, &
                                  A( K+KB, K+KB ), LDA )
                     CALL ZHEMM( 'Right', UPLO, N-K-KB+1, KB, -HALF, &
                                 A( K, K ), LDA, B( K+KB, K ), LDB, &
                                 CONE, A( K+KB, K ), LDA )
                     CALL ZTRSM( 'Left', UPLO, 'No transpose', &
                                 'Non-unit', N-K-KB+1, KB, CONE, &
                                 B( K+KB, K+KB ), LDB, A( K+KB, K ), &
                                 LDA )
                  END IF
   20          CONTINUE
            END IF
         ELSE
            IF( UPPER ) THEN
!
!              Compute U*A*U'
!
               DO 30 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
!
!                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1)
!
                  CALL ZTRMM( 'Left', UPLO, 'No transpose', 'Non-unit', &
                              K-1, KB, CONE, B, LDB, A( 1, K ), LDA )
                  CALL ZHEMM( 'Right', UPLO, K-1, KB, HALF, A( K, K ), &
                              LDA, B( 1, K ), LDB, CONE, A( 1, K ), &
                              LDA )
                  CALL ZHER2K( UPLO, 'No transpose', K-1, KB, CONE, &
                               A( 1, K ), LDA, B( 1, K ), LDB, ONE, A, &
                               LDA )
                  CALL ZHEMM( 'Right', UPLO, K-1, KB, HALF, A( K, K ), &
                              LDA, B( 1, K ), LDB, CONE, A( 1, K ), &
                              LDA )
                  CALL ZTRMM( 'Right', UPLO, 'Conjugate transpose', &
                              'Non-unit', K-1, KB, CONE, B( K, K ), LDB, &
                              A( 1, K ), LDA )
                  CALL ZHEGS2( ITYPE, UPLO, KB, A( K, K ), LDA, &
                               B( K, K ), LDB, INFO )
   30          CONTINUE
            ELSE
!
!              Compute L'*A*L
!
               DO 40 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
!
!                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1)
!
                  CALL ZTRMM( 'Right', UPLO, 'No transpose', 'Non-unit', &
                              KB, K-1, CONE, B, LDB, A( K, 1 ), LDA )
                  CALL ZHEMM( 'Left', UPLO, KB, K-1, HALF, A( K, K ), &
                              LDA, B( K, 1 ), LDB, CONE, A( K, 1 ), &
                              LDA )
                  CALL ZHER2K( UPLO, 'Conjugate transpose', K-1, KB, &
                               CONE, A( K, 1 ), LDA, B( K, 1 ), LDB, &
                               ONE, A, LDA )
                  CALL ZHEMM( 'Left', UPLO, KB, K-1, HALF, A( K, K ), &
                              LDA, B( K, 1 ), LDB, CONE, A( K, 1 ), &
                              LDA )
                  CALL ZTRMM( 'Left', UPLO, 'Conjugate transpose', &
                              'Non-unit', KB, K-1, CONE, B( K, K ), LDB, &
                              A( K, 1 ), LDA )
                  CALL ZHEGS2( ITYPE, UPLO, KB, A( K, K ), LDA, &
                               B( K, K ), LDB, INFO )
   40          CONTINUE
            END IF
         END IF
      END IF
      RETURN
!
!     End of ZHEGST
!
      END
      SUBROUTINE ZHPGST( ITYPE, UPLO, N, AP, BP, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         AP( * ), BP( * )
!     ..
!
!  Purpose
!  =======
!
!  ZHPGST reduces a complex Hermitian-definite generalized
!  eigenproblem to standard form, using packed storage.
!
!  If ITYPE = 1, the problem is A*x = lambda*B*x,
!  and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
!
!  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
!  B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.
!
!  B must have been previously factorized as U**H*U or L*L**H by ZPPTRF.
!
!  Arguments
!  =========
!
!  ITYPE   (input) INTEGER
!          = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);
!          = 2 or 3: compute U*A*U**H or L**H*A*L.
!
!  UPLO    (input) CHARACTER
!          = 'U':  Upper triangle of A is stored and B is factored as
!                  U**H*U;
!          = 'L':  Lower triangle of A is stored and B is factored as
!                  L*L**H.
!
!  N       (input) INTEGER
!          The order of the matrices A and B.  N >= 0.
!
!  AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)
!          On entry, the upper or lower triangle of the Hermitian matrix
!          A, packed columnwise in a linear array.  The j-th column of A
!          is stored in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!
!          On exit, if INFO = 0, the transformed matrix, stored in the
!          same format as A.
!
!  BP      (input) COMPLEX*16 array, dimension (N*(N+1)/2)
!          The triangular factor from the Cholesky factorization of B,
!          stored in the same format as A, as returned by ZPPTRF.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, HALF
      PARAMETER          ( ONE = 1.0D+0, HALF = 0.5D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, J1, J1J1, JJ, K, K1, K1K1, KK
      DOUBLE PRECISION   AJJ, AKK, BJJ, BKK
      COMPLEX*16         CT
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZAXPY, ZDSCAL, ZHPMV, ZHPR2, ZTPMV, &
                         ZTPSV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, ZDOTC
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHPGST', -INFO )
         RETURN
      END IF
!
      IF( ITYPE.EQ.1 ) THEN
         IF( UPPER ) THEN
!
!           Compute inv(U')*A*inv(U)
!
!           J1 and JJ are the indices of A(1,j) and A(j,j)
!
            JJ = 0
            DO 10 J = 1, N
               J1 = JJ + 1
               JJ = JJ + J
!
!              Compute the j-th column of the upper triangle of A
!
               AP( JJ ) = DBLE( AP( JJ ) )
               BJJ = BP( JJ )
               CALL ZTPSV( UPLO, 'Conjugate transpose', 'Non-unit', J, &
                           BP, AP( J1 ), 1 )
               CALL ZHPMV( UPLO, J-1, -CONE, AP, BP( J1 ), 1, CONE, &
                           AP( J1 ), 1 )
               CALL ZDSCAL( J-1, ONE / BJJ, AP( J1 ), 1 )
               AP( JJ ) = ( AP( JJ )-ZDOTC( J-1, AP( J1 ), 1, BP( J1 ), &
                          1 ) ) / BJJ
   10       CONTINUE
         ELSE
!
!           Compute inv(L)*A*inv(L')
!
!           KK and K1K1 are the indices of A(k,k) and A(k+1,k+1)
!
            KK = 1
            DO 20 K = 1, N
               K1K1 = KK + N - K + 1
!
!              Update the lower triangle of A(k:n,k:n)
!
               AKK = AP( KK )
               BKK = BP( KK )
               AKK = AKK / BKK**2
               AP( KK ) = AKK
               IF( K.LT.N ) THEN
                  CALL ZDSCAL( N-K, ONE / BKK, AP( KK+1 ), 1 )
                  CT = -HALF*AKK
                  CALL ZAXPY( N-K, CT, BP( KK+1 ), 1, AP( KK+1 ), 1 )
                  CALL ZHPR2( UPLO, N-K, -CONE, AP( KK+1 ), 1, &
                              BP( KK+1 ), 1, AP( K1K1 ) )
                  CALL ZAXPY( N-K, CT, BP( KK+1 ), 1, AP( KK+1 ), 1 )
                  CALL ZTPSV( UPLO, 'No transpose', 'Non-unit', N-K, &
                              BP( K1K1 ), AP( KK+1 ), 1 )
               END IF
               KK = K1K1
   20       CONTINUE
         END IF
      ELSE
         IF( UPPER ) THEN
!
!           Compute U*A*U'
!
!           K1 and KK are the indices of A(1,k) and A(k,k)
!
            KK = 0
            DO 30 K = 1, N
               K1 = KK + 1
               KK = KK + K
!
!              Update the upper triangle of A(1:k,1:k)
!
               AKK = AP( KK )
               BKK = BP( KK )
               CALL ZTPMV( UPLO, 'No transpose', 'Non-unit', K-1, BP, &
                           AP( K1 ), 1 )
               CT = HALF*AKK
               CALL ZAXPY( K-1, CT, BP( K1 ), 1, AP( K1 ), 1 )
               CALL ZHPR2( UPLO, K-1, CONE, AP( K1 ), 1, BP( K1 ), 1, &
                           AP )
               CALL ZAXPY( K-1, CT, BP( K1 ), 1, AP( K1 ), 1 )
               CALL ZDSCAL( K-1, BKK, AP( K1 ), 1 )
               AP( KK ) = AKK*BKK**2
   30       CONTINUE
         ELSE
!
!           Compute L'*A*L
!
!           JJ and J1J1 are the indices of A(j,j) and A(j+1,j+1)
!
            JJ = 1
            DO 40 J = 1, N
               J1J1 = JJ + N - J + 1
!
!              Compute the j-th column of the lower triangle of A
!
               AJJ = AP( JJ )
               BJJ = BP( JJ )
               AP( JJ ) = AJJ*BJJ + ZDOTC( N-J, AP( JJ+1 ), 1, &
                          BP( JJ+1 ), 1 )
               CALL ZDSCAL( N-J, BJJ, AP( JJ+1 ), 1 )
               CALL ZHPMV( UPLO, N-J, CONE, AP( J1J1 ), BP( JJ+1 ), 1, &
                           CONE, AP( JJ+1 ), 1 )
               CALL ZTPMV( UPLO, 'Conjugate transpose', 'Non-unit', &
                           N-J+1, BP( JJ ), AP( JJ ), 1 )
               JJ = J1J1
   40       CONTINUE
         END IF
      END IF
      RETURN
!
!     End of ZHPGST
!
      END
      SUBROUTINE ZHPTRD( UPLO, N, AP, D, E, TAU, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         AP( * ), TAU( * )
!     ..
!
!  Purpose
!  =======
!
!  ZHPTRD reduces a complex Hermitian matrix A stored in packed form to
!  real symmetric tridiagonal form T by a unitary similarity
!  transformation: Q**H * A * Q = T.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)
!          On entry, the upper or lower triangle of the Hermitian matrix
!          A, packed columnwise in a linear array.  The j-th column of A
!          is stored in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!          of A are overwritten by the corresponding elements of the
!          tridiagonal matrix T, and the elements above the first
!          superdiagonal, with the array TAU, represent the unitary
!          matrix Q as a product of elementary reflectors; if UPLO
!          = 'L', the diagonal and first subdiagonal of A are over-
!          written by the corresponding elements of the tridiagonal
!          matrix T, and the elements below the first subdiagonal, with
!          the array TAU, represent the unitary matrix Q as a product
!          of elementary reflectors. See Further Details.
!
!  D       (output) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = A(i,i).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!
!  TAU     (output) COMPLEX*16 array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n-1) . . . H(2) H(1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a complex scalar, and v is a complex vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in AP,
!  overwriting A(1:i-1,i+1), and tau is stored in TAU(i).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a complex scalar, and v is a complex vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in AP,
!  overwriting A(i+2:n,i), and tau is stored in TAU(i).
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO, HALF
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), &
                         ZERO = ( 0.0D+0, 0.0D+0 ), &
                         HALF = ( 0.5D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, I1, I1I1, II
      COMPLEX*16         ALPHA, TAUI
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZAXPY, ZHPMV, ZHPR2, ZLARFG
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, ZDOTC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHPTRD', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) &
         RETURN
!
      IF( UPPER ) THEN
!
!        Reduce the upper triangle of A.
!        I1 is the index in AP of A(1,I+1).
!
         I1 = N*( N-1 ) / 2 + 1
         AP( I1+N-1 ) = DBLE( AP( I1+N-1 ) )
         DO 10 I = N - 1, 1, -1
!
!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(1:i-1,i+1)
!
            ALPHA = AP( I1+I-1 )
            CALL ZLARFG( I, ALPHA, AP( I1 ), 1, TAUI )
            E( I ) = ALPHA
!
            IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(1:i,1:i)
!
               AP( I1+I-1 ) = ONE
!
!              Compute  y := tau * A * v  storing y in TAU(1:i)
!
               CALL ZHPMV( UPLO, I, TAUI, AP, AP( I1 ), 1, ZERO, TAU, &
                           1 )
!
!              Compute  w := y - 1/2 * tau * (y'*v) * v
!
               ALPHA = -HALF*TAUI*ZDOTC( I, TAU, 1, AP( I1 ), 1 )
               CALL ZAXPY( I, ALPHA, AP( I1 ), 1, TAU, 1 )
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'
!
               CALL ZHPR2( UPLO, I, -ONE, AP( I1 ), 1, TAU, 1, AP )
!
            END IF
            AP( I1+I-1 ) = E( I )
            D( I+1 ) = AP( I1+I )
            TAU( I ) = TAUI
            I1 = I1 - I
   10    CONTINUE
         D( 1 ) = AP( 1 )
      ELSE
!
!        Reduce the lower triangle of A. II is the index in AP of
!        A(i,i) and I1I1 is the index of A(i+1,i+1).
!
         II = 1
         AP( 1 ) = DBLE( AP( 1 ) )
         DO 20 I = 1, N - 1
            I1I1 = II + N - I + 1
!
!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(i+2:n,i)
!
            ALPHA = AP( II+1 )
            CALL ZLARFG( N-I, ALPHA, AP( II+2 ), 1, TAUI )
            E( I ) = ALPHA
!
            IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(i+1:n,i+1:n)
!
               AP( II+1 ) = ONE
!
!              Compute  y := tau * A * v  storing y in TAU(i:n-1)
!
               CALL ZHPMV( UPLO, N-I, TAUI, AP( I1I1 ), AP( II+1 ), 1, &
                           ZERO, TAU( I ), 1 )
!
!              Compute  w := y - 1/2 * tau * (y'*v) * v
!
               ALPHA = -HALF*TAUI*ZDOTC( N-I, TAU( I ), 1, AP( II+1 ), &
                       1 )
               CALL ZAXPY( N-I, ALPHA, AP( II+1 ), 1, TAU( I ), 1 )
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'
!
               CALL ZHPR2( UPLO, N-I, -ONE, AP( II+1 ), 1, TAU( I ), 1, &
                           AP( I1I1 ) )
!
            END IF
            AP( II+1 ) = E( I )
            D( I ) = AP( II )
            TAU( I ) = TAUI
            II = I1I1
   20    CONTINUE
         D( N ) = AP( II )
      END IF
!
      RETURN
!
!     End of ZHPTRD
!
      END
      SUBROUTINE ZLACGV( N, X, INCX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         X( * )
!     ..
!
!  Purpose
!  =======
!
!  ZLACGV conjugates a complex vector of length N.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The length of the vector X.  N >= 0.
!
!  X       (input/output) COMPLEX*16 array, dimension
!                         (1+(N-1)*abs(INCX))
!          On entry, the vector of length N to be conjugated.
!          On exit, X is overwritten with conjg(X).
!
!  INCX    (input) INTEGER
!          The spacing between successive elements of X.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IOFF
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
!     ..
!     .. Executable Statements ..
!
      IF( INCX.EQ.1 ) THEN
         DO 10 I = 1, N
            X( I ) = DCONJG( X( I ) )
   10    CONTINUE
      ELSE
         IOFF = 1
         IF( INCX.LT.0 ) &
            IOFF = 1 - ( N-1 )*INCX
         DO 20 I = 1, N
            X( IOFF ) = DCONJG( X( IOFF ) )
            IOFF = IOFF + INCX
   20    CONTINUE
      END IF
      RETURN
!
!     End of ZLACGV
!
      END
      DOUBLE COMPLEX   FUNCTION ZLADIV( X, Y )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      COMPLEX*16         X, Y
!     ..
!
!  Purpose
!  =======
!
!  ZLADIV := X / Y, where X and Y are complex.  The computation of X / Y
!  will not overflow on an intermediary step unless the results
!  overflows.
!
!  Arguments
!  =========
!
!  X       (input) COMPLEX*16
!  Y       (input) COMPLEX*16
!          The complex scalars X and Y.
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION   ZI, ZR
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLADIV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, DIMAG
!     ..
!     .. Executable Statements ..
!
      CALL DLADIV( DBLE( X ), DIMAG( X ), DBLE( Y ), DIMAG( Y ), ZR, &
                   ZI )
      ZLADIV = DCMPLX( ZR, ZI )
!
      RETURN
!
!     End of ZLADIV
!
      END
      DOUBLE PRECISION FUNCTION ZLANHP( NORM, UPLO, N, AP, WORK )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   WORK( * )
      COMPLEX*16         AP( * )
!     ..
!
!  Purpose
!  =======
!
!  ZLANHP  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  complex hermitian matrix A,  supplied in packed form.
!
!  Description
!  ===========
!
!  ZLANHP returns the value
!
!     ZLANHP = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in ZLANHP as described
!          above.
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          hermitian matrix A is supplied.
!          = 'U':  Upper triangular part of A is supplied
!          = 'L':  Lower triangular part of A is supplied
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, ZLANHP is
!          set to zero.
!
!  AP      (input) COMPLEX*16 array, dimension (N*(N+1)/2)
!          The upper or lower triangle of the hermitian matrix A, packed
!          columnwise in a linear array.  The j-th column of A is stored
!          in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!          Note that the  imaginary parts of the diagonal elements need
!          not be set and are assumed to be zero.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
!          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
!          WORK is not referenced.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, K
      DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZLASSQ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            K = 0
            DO 20 J = 1, N
               DO 10 I = K + 1, K + J - 1
                  VALUE = MAX( VALUE, ABS( AP( I ) ) )
   10          CONTINUE
               K = K + J
               VALUE = MAX( VALUE, ABS( DBLE( AP( K ) ) ) )
   20       CONTINUE
         ELSE
            K = 1
            DO 40 J = 1, N
               VALUE = MAX( VALUE, ABS( DBLE( AP( K ) ) ) )
               DO 30 I = K + 1, K + N - J
                  VALUE = MAX( VALUE, ABS( AP( I ) ) )
   30          CONTINUE
               K = K + N - J + 1
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. &
               ( NORM.EQ.'1' ) ) THEN
!
!        Find normI(A) ( = norm1(A), since A is hermitian).
!
         VALUE = ZERO
         K = 1
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( AP( K ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
                  K = K + 1
   50          CONTINUE
               WORK( J ) = SUM + ABS( DBLE( AP( K ) ) )
               K = K + 1
   60       CONTINUE
            DO 70 I = 1, N
               VALUE = MAX( VALUE, WORK( I ) )
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( DBLE( AP( K ) ) )
               K = K + 1
               DO 90 I = J + 1, N
                  ABSA = ABS( AP( K ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
                  K = K + 1
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         K = 2
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 2, N
               CALL ZLASSQ( J-1, AP( K ), 1, SCALE, SUM )
               K = K + J
  110       CONTINUE
         ELSE
            DO 120 J = 1, N - 1
               CALL ZLASSQ( N-J, AP( K ), 1, SCALE, SUM )
               K = K + N - J + 1
  120       CONTINUE
         END IF
         SUM = 2*SUM
         K = 1
         DO 130 I = 1, N
            IF( DBLE( AP( K ) ).NE.ZERO ) THEN
               ABSA = ABS( DBLE( AP( K ) ) )
               IF( SCALE.LT.ABSA ) THEN
                  SUM = ONE + SUM*( SCALE / ABSA )**2
                  SCALE = ABSA
               ELSE
                  SUM = SUM + ( ABSA / SCALE )**2
               END IF
            END IF
            IF( LSAME( UPLO, 'U' ) ) THEN
               K = K + I + 1
            ELSE
               K = K + N - I + 1
            END IF
  130    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
!
      ZLANHP = VALUE
      RETURN
!
!     End of ZLANHP
!
      END
      SUBROUTINE ZLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      COMPLEX*16         TAU
!     ..
!     .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZLARF applies a complex elementary reflector H to a complex M-by-N
!  matrix C, from either the left or the right. H is represented in the
!  form
!
!        H = I - tau * v * v'
!
!  where tau is a complex scalar and v is a complex vector.
!
!  If tau = 0, then H is taken to be the unit matrix.
!
!  To apply H' (the conjugate transpose of H), supply conjg(tau) instead
!  tau.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': form  H * C
!          = 'R': form  C * H
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  V       (input) COMPLEX*16 array, dimension
!                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!          The vector v in the representation of H. V is not used if
!          TAU = 0.
!
!  INCV    (input) INTEGER
!          The increment between elements of v. INCV <> 0.
!
!  TAU     (input) COMPLEX*16
!          The value tau in the representation of H.
!
!  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!          or C * H if SIDE = 'R'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) COMPLEX*16 array, dimension
!                         (N) if SIDE = 'L'
!                      or (M) if SIDE = 'R'
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), &
                         ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZGEMV, ZGERC
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  H * C
!
         IF( TAU.NE.ZERO ) THEN
!
!           w := C' * v
!
            CALL ZGEMV( 'Conjugate transpose', M, N, ONE, C, LDC, V, &
                        INCV, ZERO, WORK, 1 )
!
!           C := C - v * w'
!
            CALL ZGERC( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
!
!        Form  C * H
!
         IF( TAU.NE.ZERO ) THEN
!
!           w := C * v
!
            CALL ZGEMV( 'No transpose', M, N, ONE, C, LDC, V, INCV, &
                        ZERO, WORK, 1 )
!
!           C := C - w * v'
!
            CALL ZGERC( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
!
!     End of ZLARF
!
      END
      SUBROUTINE ZLARFG( N, ALPHA, X, INCX, TAU )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX*16         ALPHA, TAU
!     ..
!     .. Array Arguments ..
      COMPLEX*16         X( * )
!     ..
!
!  Purpose
!  =======
!
!  ZLARFG generates a complex elementary reflector H of order n, such
!  that
!
!        H' * ( alpha ) = ( beta ),   H' * H = I.
!             (   x   )   (   0  )
!
!  where alpha and beta are scalars, with beta real, and x is an
!  (n-1)-element complex vector. H is represented in the form
!
!        H = I - tau * ( 1 ) * ( 1 v' ) ,
!                      ( v )
!
!  where tau is a complex scalar and v is a complex (n-1)-element
!  vector. Note that H is not hermitian.
!
!  If the elements of x are all zero and alpha is real, then tau = 0
!  and H is taken to be the unit matrix.
!
!  Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the elementary reflector.
!
!  ALPHA   (input/output) COMPLEX*16
!          On entry, the value alpha.
!          On exit, it is overwritten with the value beta.
!
!  X       (input/output) COMPLEX*16 array, dimension
!                         (1+(N-2)*abs(INCX))
!          On entry, the vector x.
!          On exit, it is overwritten with the vector v.
!
!  INCX    (input) INTEGER
!          The increment between elements of X. INCX > 0.
!
!  TAU     (output) COMPLEX*16
!          The value tau.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY3, DZNRM2
      COMPLEX*16         ZLADIV
      EXTERNAL           DLAMCH, DLAPY3, DZNRM2, ZLADIV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DIMAG, SIGN
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZDSCAL, ZSCAL
!     ..
!     .. Executable Statements ..
!
      IF( N.LE.0 ) THEN
         TAU = ZERO
         RETURN
      END IF
!
      XNORM = DZNRM2( N-1, X, INCX )
      ALPHR = DBLE( ALPHA )
      ALPHI = DIMAG( ALPHA )
!
      IF( XNORM.EQ.ZERO .AND. ALPHI.EQ.ZERO ) THEN
!
!        H  =  I
!
         TAU = ZERO
      ELSE
!
!        general case
!
         BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         RSAFMN = ONE / SAFMIN
!
         IF( ABS( BETA ).LT.SAFMIN ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
            KNT = 0
   10       CONTINUE
            KNT = KNT + 1
            CALL ZDSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHI = ALPHI*RSAFMN
            ALPHR = ALPHR*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN ) &
               GO TO 10
!
!           New BETA is at most 1, at least SAFMIN
!
            XNORM = DZNRM2( N-1, X, INCX )
            ALPHA = DCMPLX( ALPHR, ALPHI )
            BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
            TAU = DCMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
            ALPHA = ZLADIV( DCMPLX( ONE ), ALPHA-BETA )
            CALL ZSCAL( N-1, ALPHA, X, INCX )
!
!           If ALPHA is subnormal, it may lose relative accuracy
!
            ALPHA = BETA
            DO 20 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   20       CONTINUE
         ELSE
            TAU = DCMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
            ALPHA = ZLADIV( DCMPLX( ONE ), ALPHA-BETA )
            CALL ZSCAL( N-1, ALPHA, X, INCX )
            ALPHA = BETA
         END IF
      END IF
!
      RETURN
!
!     End of ZLARFG
!
      END
      SUBROUTINE ZLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      COMPLEX*16         ALPHA, BETA
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZLASET initializes a 2-D array A to BETA on the diagonal and
!  ALPHA on the offdiagonals.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be set.
!          = 'U':      Upper triangular part is set. The lower triangle
!                      is unchanged.
!          = 'L':      Lower triangular part is set. The upper triangle
!                      is unchanged.
!          Otherwise:  All of the matrix A is set.
!
!  M       (input) INTEGER
!          On entry, M specifies the number of rows of A.
!
!  N       (input) INTEGER
!          On entry, N specifies the number of columns of A.
!
!  ALPHA   (input) COMPLEX*16
!          All the offdiagonal array elements are set to ALPHA.
!
!  BETA    (input) COMPLEX*16
!          All the diagonal array elements are set to BETA.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the m by n matrix A.
!          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
!                   A(i,i) = BETA , 1 <= i <= min(m,n)
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Set the diagonal to BETA and the strictly upper triangular
!        part of the array to ALPHA.
!
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
         DO 30 I = 1, MIN( N, M )
            A( I, I ) = BETA
   30    CONTINUE
!
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
!
!        Set the diagonal to BETA and the strictly lower triangular
!        part of the array to ALPHA.
!
         DO 50 J = 1, MIN( M, N )
            DO 40 I = J + 1, M
               A( I, J ) = ALPHA
   40       CONTINUE
   50    CONTINUE
         DO 60 I = 1, MIN( N, M )
            A( I, I ) = BETA
   60    CONTINUE
!
      ELSE
!
!        Set the array to BETA on the diagonal and ALPHA on the
!        offdiagonal.
!
         DO 80 J = 1, N
            DO 70 I = 1, M
               A( I, J ) = ALPHA
   70       CONTINUE
   80    CONTINUE
         DO 90 I = 1, MIN( M, N )
            A( I, I ) = BETA
   90    CONTINUE
      END IF
!
      RETURN
!
!     End of ZLASET
!
      END
      SUBROUTINE ZLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( * ), S( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZLASR   performs the transformation
!
!     A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )
!
!     A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )
!
!  where A is an m by n complex matrix and P is an orthogonal matrix,
!  consisting of a sequence of plane rotations determined by the
!  parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'
!  and z = n when SIDE = 'R' or 'r' ):
!
!  When  DIRECT = 'F' or 'f'  ( Forward sequence ) then
!
!     P = P( z - 1 )*...*P( 2 )*P( 1 ),
!
!  and when DIRECT = 'B' or 'b'  ( Backward sequence ) then
!
!     P = P( 1 )*P( 2 )*...*P( z - 1 ),
!
!  where  P( k ) is a plane rotation matrix for the following planes:
!
!     when  PIVOT = 'V' or 'v'  ( Variable pivot ),
!        the plane ( k, k + 1 )
!
!     when  PIVOT = 'T' or 't'  ( Top pivot ),
!        the plane ( 1, k + 1 )
!
!     when  PIVOT = 'B' or 'b'  ( Bottom pivot ),
!        the plane ( k, z )
!
!  c( k ) and s( k )  must contain the  cosine and sine that define the
!  matrix  P( k ).  The two by two plane rotation part of the matrix
!  P( k ), R( k ), is assumed to be of the form
!
!     R( k ) = (  c( k )  s( k ) ).
!              ( -s( k )  c( k ) )
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          Specifies whether the plane rotation matrix P is applied to
!          A on the left or the right.
!          = 'L':  Left, compute A := P*A
!          = 'R':  Right, compute A:= A*P'
!
!  DIRECT  (input) CHARACTER*1
!          Specifies whether P is a forward or backward sequence of
!          plane rotations.
!          = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )
!          = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )
!
!  PIVOT   (input) CHARACTER*1
!          Specifies the plane for which P(k) is a plane rotation
!          matrix.
!          = 'V':  Variable pivot, the plane (k,k+1)
!          = 'T':  Top pivot, the plane (1,k+1)
!          = 'B':  Bottom pivot, the plane (k,z)
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  If m <= 1, an immediate
!          return is effected.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  If n <= 1, an
!          immediate return is effected.
!
!  C, S    (input) DOUBLE PRECISION arrays, dimension
!                  (M-1) if SIDE = 'L'
!                  (N-1) if SIDE = 'R'
!          c(k) and s(k) contain the cosine and sine that define the
!          matrix P(k).  The two by two plane rotation part of the
!          matrix P(k), R(k), is assumed to be of the form
!          R( k ) = (  c( k )  s( k ) ).
!                   ( -s( k )  c( k ) )
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          The m by n matrix A.  On exit, A is overwritten by P*A if
!          SIDE = 'R' or by A*P' if SIDE = 'L'.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   CTEMP, STEMP
      COMPLEX*16         TEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT, &
               'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) ) &
                THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLASR ', INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) ) &
         RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  P * A
!
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!        Form A * P'
!
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of ZLASR
!
      END
      SUBROUTINE ZLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
!     ..
!     .. Array Arguments ..
      COMPLEX*16         X( * )
!     ..
!
!  Purpose
!  =======
!
!  ZLASSQ returns the values scl and ssq such that
!
!     ( scl**2 )*ssq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
!
!  where x( i ) = abs( X( 1 + ( i - 1 )*INCX ) ). The value of sumsq is
!  assumed to be at least unity and the value of ssq will then satisfy
!
!     1.0 .le. ssq .le. ( sumsq + 2*n ).
!
!  scale is assumed to be non-negative and scl returns the value
!
!     scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),
!            i
!
!  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
!  SCALE and SUMSQ are overwritten by scl and ssq respectively.
!
!  The routine makes only one pass through the vector X.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements to be used from the vector X.
!
!  X       (input) DOUBLE PRECISION
!          The vector x as described above.
!             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
!
!  INCX    (input) INTEGER
!          The increment between successive values of the vector X.
!          INCX > 0.
!
!  SCALE   (input/output) DOUBLE PRECISION
!          On entry, the value  scale  in the equation above.
!          On exit, SCALE is overwritten with the value  scl .
!
!  SUMSQ   (input/output) DOUBLE PRECISION
!          On entry, the value  sumsq  in the equation above.
!          On exit, SUMSQ is overwritten with the value  ssq .
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   TEMP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG
!     ..
!     .. Executable Statements ..
!
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( DBLE( X( IX ) ).NE.ZERO ) THEN
               TEMP1 = ABS( DBLE( X( IX ) ) )
               IF( SCALE.LT.TEMP1 ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
               END IF
            END IF
            IF( DIMAG( X( IX ) ).NE.ZERO ) THEN
               TEMP1 = ABS( DIMAG( X( IX ) ) )
               IF( SCALE.LT.TEMP1 ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
!
      RETURN
!
!     End of ZLASSQ
!
      END
      SUBROUTINE ZPOTF2( UPLO, N, A, LDA, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZPOTF2 computes the Cholesky factorization of a complex Hermitian
!  positive definite matrix A.
!
!  The factorization has the form
!     A = U' * U ,  if UPLO = 'U', or
!     A = L  * L',  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          Hermitian matrix A is stored.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!          n by n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n by n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the factor U or L from the Cholesky
!          factorization A = U'*U  or A = L*L'.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, the leading minor of order k is not
!               positive definite, and the factorization could not be
!               completed.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, ZDOTC
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZDSCAL, ZGEMV, ZLACGV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPOTF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF( UPPER ) THEN
!
!        Compute the Cholesky factorization A = U'*U.
!
         DO 10 J = 1, N
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
            AJJ = DBLE( A( J, J ) ) - ZDOTC( J-1, A( 1, J ), 1, &
                  A( 1, J ), 1 )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
!
!           Compute elements J+1:N of row J.
!
            IF( J.LT.N ) THEN
               CALL ZLACGV( J-1, A( 1, J ), 1 )
               CALL ZGEMV( 'Transpose', J-1, N-J, -CONE, A( 1, J+1 ), &
                           LDA, A( 1, J ), 1, CONE, A( J, J+1 ), LDA )
               CALL ZLACGV( J-1, A( 1, J ), 1 )
               CALL ZDSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            END IF
   10    CONTINUE
      ELSE
!
!        Compute the Cholesky factorization A = L*L'.
!
         DO 20 J = 1, N
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
            AJJ = DBLE( A( J, J ) ) - ZDOTC( J-1, A( J, 1 ), LDA, &
                  A( J, 1 ), LDA )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
!
!           Compute elements J+1:N of column J.
!
            IF( J.LT.N ) THEN
               CALL ZLACGV( J-1, A( J, 1 ), LDA )
               CALL ZGEMV( 'No transpose', N-J, J-1, -CONE, A( J+1, 1 ), &
                           LDA, A( J, 1 ), LDA, CONE, A( J+1, J ), 1 )
               CALL ZLACGV( J-1, A( J, 1 ), LDA )
               CALL ZDSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
      GO TO 40
!
   30 CONTINUE
      INFO = J
!
   40 CONTINUE
      RETURN
!
!     End of ZPOTF2
!
      END
      SUBROUTINE ZPOTRF( UPLO, N, A, LDA, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZPOTRF computes the Cholesky factorization of a complex Hermitian
!  positive definite matrix A.
!
!  The factorization has the form
!     A = U**H * U,  if UPLO = 'U', or
!     A = L  * L**H,  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  This is the block version of the algorithm, calling Level 3 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the factor U or L from the Cholesky
!          factorization A = U**H*U or A = L*L**H.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i is not
!                positive definite, and the factorization could not be
!                completed.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      COMPLEX*16         CONE
      PARAMETER          ( ONE = 1.0D+0, CONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JB, NB
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZHERK, ZPOTF2, ZTRSM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPOTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'ZPOTRF', UPLO, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
!
!        Use unblocked code.
!
         CALL ZPOTF2( UPLO, N, A, LDA, INFO )
      ELSE
!
!        Use blocked code.
!
         IF( UPPER ) THEN
!
!           Compute the Cholesky factorization A = U'*U.
!
            DO 10 J = 1, N, NB
!
!              Update and factorize the current diagonal block and test
!              for non-positive-definiteness.
!
               JB = MIN( NB, N-J+1 )
               CALL ZHERK( 'Upper', 'Conjugate transpose', JB, J-1, &
                           -ONE, A( 1, J ), LDA, ONE, A( J, J ), LDA )
               CALL ZPOTF2( 'Upper', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 ) &
                  GO TO 30
               IF( J+JB.LE.N ) THEN
!
!                 Compute the current block row.
!
                  CALL ZGEMM( 'Conjugate transpose', 'No transpose', JB, &
                              N-J-JB+1, J-1, -CONE, A( 1, J ), LDA, &
                              A( 1, J+JB ), LDA, CONE, A( J, J+JB ), &
                              LDA )
                  CALL ZTRSM( 'Left', 'Upper', 'Conjugate transpose', &
                              'Non-unit', JB, N-J-JB+1, CONE, A( J, J ), &
                              LDA, A( J, J+JB ), LDA )
               END IF
   10       CONTINUE
!
         ELSE
!
!           Compute the Cholesky factorization A = L*L'.
!
            DO 20 J = 1, N, NB
!
!              Update and factorize the current diagonal block and test
!              for non-positive-definiteness.
!
               JB = MIN( NB, N-J+1 )
               CALL ZHERK( 'Lower', 'No transpose', JB, J-1, -ONE, &
                           A( J, 1 ), LDA, ONE, A( J, J ), LDA )
               CALL ZPOTF2( 'Lower', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 ) &
                  GO TO 30
               IF( J+JB.LE.N ) THEN
!
!                 Compute the current block column.
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', &
                              N-J-JB+1, JB, J-1, -CONE, A( J+JB, 1 ), &
                              LDA, A( J, 1 ), LDA, CONE, A( J+JB, J ), &
                              LDA )
                  CALL ZTRSM( 'Right', 'Lower', 'Conjugate transpose', &
                              'Non-unit', N-J-JB+1, JB, CONE, A( J, J ), &
                              LDA, A( J+JB, J ), LDA )
               END IF
   20       CONTINUE
         END IF
      END IF
      GO TO 40
!
   30 CONTINUE
      INFO = INFO + J - 1
!
   40 CONTINUE
      RETURN
!
!     End of ZPOTRF
!
      END
      SUBROUTINE ZPPTRF( UPLO, N, AP, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         AP( * )
!     ..
!
!  Purpose
!  =======
!
!  ZPPTRF computes the Cholesky factorization of a complex Hermitian
!  positive definite matrix A stored in packed format.
!
!  The factorization has the form
!     A = U**H * U,  if UPLO = 'U', or
!     A = L  * L**H,  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)
!          On entry, the upper or lower triangle of the Hermitian matrix
!          A, packed columnwise in a linear array.  The j-th column of A
!          is stored in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!          See below for further details.
!
!          On exit, if INFO = 0, the triangular factor U or L from the
!          Cholesky factorization A = U**H*U or A = L*L**H, in the same
!          storage format as A.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i is not
!                positive definite, and the factorization could not be
!                completed.
!
!  Further Details
!  ===============
!
!  The packed storage scheme is illustrated by the following example
!  when N = 4, UPLO = 'U':
!
!  Two-dimensional storage of the Hermitian matrix A:
!
!     a11 a12 a13 a14
!         a22 a23 a24
!             a33 a34     (aij = conjg(aji))
!                 a44
!
!  Packed storage of the upper triangle of A:
!
!  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JC, JJ
      DOUBLE PRECISION   AJJ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, ZDOTC
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZDSCAL, ZHPR, ZTPSV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPPTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF( UPPER ) THEN
!
!        Compute the Cholesky factorization A = U'*U.
!
         JJ = 0
         DO 10 J = 1, N
            JC = JJ + 1
            JJ = JJ + J
!
!           Compute elements 1:J-1 of column J.
!
            IF( J.GT.1 ) &
               CALL ZTPSV( 'Upper', 'Conjugate transpose', 'Non-unit', &
                           J-1, AP, AP( JC ), 1 )
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
            AJJ = DBLE( AP( JJ ) ) - ZDOTC( J-1, AP( JC ), 1, AP( JC ), &
                  1 )
            IF( AJJ.LE.ZERO ) THEN
               AP( JJ ) = AJJ
               GO TO 30
            END IF
            AP( JJ ) = SQRT( AJJ )
   10    CONTINUE
      ELSE
!
!        Compute the Cholesky factorization A = L*L'.
!
         JJ = 1
         DO 20 J = 1, N
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
            AJJ = DBLE( AP( JJ ) )
            IF( AJJ.LE.ZERO ) THEN
               AP( JJ ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            AP( JJ ) = AJJ
!
!           Compute elements J+1:N of column J and update the trailing
!           submatrix.
!
            IF( J.LT.N ) THEN
               CALL ZDSCAL( N-J, ONE / AJJ, AP( JJ+1 ), 1 )
               CALL ZHPR( 'Lower', N-J, -ONE, AP( JJ+1 ), 1, &
                          AP( JJ+N-J+1 ) )
               JJ = JJ + N - J + 1
            END IF
   20    CONTINUE
      END IF
      GO TO 40
!
   30 CONTINUE
      INFO = J
!
   40 CONTINUE
      RETURN
!
!     End of ZPPTRF
!
      END
      SUBROUTINE ZSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK, &
                         IWORK, IFAIL, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDZ, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IBLOCK( * ), IFAIL( * ), ISPLIT( * ), &
                         IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )
      COMPLEX*16         Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  ZSTEIN computes the eigenvectors of a real symmetric tridiagonal
!  matrix T corresponding to specified eigenvalues, using inverse
!  iteration.
!
!  The maximum number of iterations allowed for each eigenvector is
!  specified by an internal parameter MAXITS (currently set to 5).
!
!  Although the eigenvectors are real, they are stored in a complex
!  array, which may be passed to ZUNMTR or ZUPMTR for back
!  transformation to the eigenvectors of a complex Hermitian matrix
!  which was reduced to tridiagonal form.
!
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The n diagonal elements of the tridiagonal matrix T.
!
!  E       (input) DOUBLE PRECISION array, dimension (N)
!          The (n-1) subdiagonal elements of the tridiagonal matrix
!          T, stored in elements 1 to N-1; E(N) need not be set.
!
!  M       (input) INTEGER
!          The number of eigenvectors to be found.  0 <= M <= N.
!
!  W       (input) DOUBLE PRECISION array, dimension (N)
!          The first M elements of W contain the eigenvalues for
!          which eigenvectors are to be computed.  The eigenvalues
!          should be grouped by split-off block and ordered from
!          smallest to largest within the block.  ( The output array
!          W from DSTEBZ with ORDER = 'B' is expected here. )
!
!  IBLOCK  (input) INTEGER array, dimension (N)
!          The submatrix indices associated with the corresponding
!          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to
!          the first submatrix from the top, =2 if W(i) belongs to
!          the second submatrix, etc.  ( The output array IBLOCK
!          from DSTEBZ is expected here. )
!
!  ISPLIT  (input) INTEGER array, dimension (N)
!          The splitting points, at which T breaks up into submatrices.
!          The first submatrix consists of rows/columns 1 to
!          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1
!          through ISPLIT( 2 ), etc.
!          ( The output array ISPLIT from DSTEBZ is expected here. )
!
!  Z       (output) COMPLEX*16 array, dimension (LDZ, M)
!          The computed eigenvectors.  The eigenvector associated
!          with the eigenvalue W(i) is stored in the i-th column of
!          Z.  Any vector which fails to converge is set to its current
!          iterate after MAXITS iterations.
!          The imaginary parts of the eigenvectors are set to zero.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= max(1,N).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (5*N)
!
!  IWORK   (workspace) INTEGER array, dimension (N)
!
!  IFAIL   (output) INTEGER array, dimension (M)
!          On normal exit, all elements of IFAIL are zero.
!          If one or more eigenvectors fail to converge after
!          MAXITS iterations, then their indices are stored in
!          array IFAIL.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, then i eigenvectors failed to converge
!               in MAXITS iterations.  Their indices are stored in
!               array IFAIL.
!
!  Internal Parameters
!  ===================
!
!  MAXITS  INTEGER, default = 5
!          The maximum number of iterations performed.
!
!  EXTRA   INTEGER, default = 2
!          The number of iterations performed after norm growth
!          criterion is satisfied, should be at least 1.
!
! =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), &
                         CONE = ( 1.0D+0, 0.0D+0 ) )
      DOUBLE PRECISION   ZERO, ONE, TEN, ODM3, ODM1
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TEN = 1.0D+1, &
                         ODM3 = 1.0D-3, ODM1 = 1.0D-1 )
      INTEGER            MAXITS, EXTRA
      PARAMETER          ( MAXITS = 5, EXTRA = 2 )
!     ..
!     .. Local Scalars ..
      INTEGER            B1, BLKSIZ, BN, GPIND, I, IINFO, INDRV1, &
                         INDRV2, INDRV3, INDRV4, INDRV5, ITS, J, J1, &
                         JBLK, JMAX, JR, NBLK, NRMCHK
      DOUBLE PRECISION   DTPCRT, EPS, EPS1, NRM, ONENRM, ORTOL, PERTOL, &
                         SCL, SEP, TOL, XJ, XJM, ZTR
!     ..
!     .. Local Arrays ..
      INTEGER            ISEED( 4 )
!     ..
!     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DASUM, DLAMCH, DNRM2
      EXTERNAL           IDAMAX, DASUM, DLAMCH, DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAGTF, DLAGTS, DLARNV, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      DO 10 I = 1, M
         IFAIL( I ) = 0
   10 CONTINUE
!
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 .OR. M.GT.N ) THEN
         INFO = -4
      ELSE IF( LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE
         DO 20 J = 2, M
            IF( IBLOCK( J ).LT.IBLOCK( J-1 ) ) THEN
               INFO = -6
               GO TO 30
            END IF
            IF( IBLOCK( J ).EQ.IBLOCK( J-1 ) .AND. W( J ).LT.W( J-1 ) ) &
                 THEN
               INFO = -5
               GO TO 30
            END IF
   20    CONTINUE
   30    CONTINUE
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSTEIN', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. M.EQ.0 ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
         Z( 1, 1 ) = CONE
         RETURN
      END IF
!
!     Get machine constants.
!
      EPS = DLAMCH( 'Precision' )
!
!     Initialize seed for random number generator DLARNV.
!
      DO 40 I = 1, 4
         ISEED( I ) = 1
   40 CONTINUE
!
!     Initialize pointers.
!
      INDRV1 = 0
      INDRV2 = INDRV1 + N
      INDRV3 = INDRV2 + N
      INDRV4 = INDRV3 + N
      INDRV5 = INDRV4 + N
!
!     Compute eigenvectors of matrix blocks.
!
      J1 = 1
      DO 180 NBLK = 1, IBLOCK( M )
!
!        Find starting and ending indices of block nblk.
!
         IF( NBLK.EQ.1 ) THEN
            B1 = 1
         ELSE
            B1 = ISPLIT( NBLK-1 ) + 1
         END IF
         BN = ISPLIT( NBLK )
         BLKSIZ = BN - B1 + 1
         IF( BLKSIZ.EQ.1 ) &
            GO TO 60
         GPIND = B1
!
!        Compute reorthogonalization criterion and stopping criterion.
!
         ONENRM = ABS( D( B1 ) ) + ABS( E( B1 ) )
         ONENRM = MAX( ONENRM, ABS( D( BN ) )+ABS( E( BN-1 ) ) )
         DO 50 I = B1 + 1, BN - 1
            ONENRM = MAX( ONENRM, ABS( D( I ) )+ABS( E( I-1 ) )+ &
                     ABS( E( I ) ) )
   50    CONTINUE
         ORTOL = ODM3*ONENRM
!
         DTPCRT = SQRT( ODM1 / BLKSIZ )
!
!        Loop through eigenvalues of block nblk.
!
   60    CONTINUE
         JBLK = 0
         DO 170 J = J1, M
            IF( IBLOCK( J ).NE.NBLK ) THEN
               J1 = J
               GO TO 180
            END IF
            JBLK = JBLK + 1
            XJ = W( J )
!
!           Skip all the work if the block size is one.
!
            IF( BLKSIZ.EQ.1 ) THEN
               WORK( INDRV1+1 ) = ONE
               GO TO 140
            END IF
!
!           If eigenvalues j and j-1 are too close, add a relatively
!           small perturbation.
!
            IF( JBLK.GT.1 ) THEN
               EPS1 = ABS( EPS*XJ )
               PERTOL = TEN*EPS1
               SEP = XJ - XJM
               IF( SEP.LT.PERTOL ) &
                  XJ = XJM + PERTOL
            END IF
!
            ITS = 0
            NRMCHK = 0
!
!           Get random starting vector.
!
            CALL DLARNV( 2, ISEED, BLKSIZ, WORK( INDRV1+1 ) )
!
!           Copy the matrix T so it won't be destroyed in factorization.
!
            CALL DCOPY( BLKSIZ, D( B1 ), 1, WORK( INDRV4+1 ), 1 )
            CALL DCOPY( BLKSIZ-1, E( B1 ), 1, WORK( INDRV2+2 ), 1 )
            CALL DCOPY( BLKSIZ-1, E( B1 ), 1, WORK( INDRV3+1 ), 1 )
!
!           Compute LU factors with partial pivoting  ( PT = LU )
!
            TOL = ZERO
            CALL DLAGTF( BLKSIZ, WORK( INDRV4+1 ), XJ, WORK( INDRV2+2 ), &
                         WORK( INDRV3+1 ), TOL, WORK( INDRV5+1 ), IWORK, &
                         IINFO )
!
!           Update iteration count.
!
   70       CONTINUE
            ITS = ITS + 1
            IF( ITS.GT.MAXITS ) &
               GO TO 120
!
!           Normalize and scale the righthand side vector Pb.
!
            SCL = BLKSIZ*ONENRM*MAX( EPS, &
                  ABS( WORK( INDRV4+BLKSIZ ) ) ) / &
                  DASUM( BLKSIZ, WORK( INDRV1+1 ), 1 )
            CALL DSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
!
!           Solve the system LU = Pb.
!
            CALL DLAGTS( -1, BLKSIZ, WORK( INDRV4+1 ), WORK( INDRV2+2 ), &
                         WORK( INDRV3+1 ), WORK( INDRV5+1 ), IWORK, &
                         WORK( INDRV1+1 ), TOL, IINFO )
!
!           Reorthogonalize by modified Gram-Schmidt if eigenvalues are
!           close enough.
!
            IF( JBLK.EQ.1 ) &
               GO TO 110
            IF( ABS( XJ-XJM ).GT.ORTOL ) &
               GPIND = J
            IF( GPIND.NE.J ) THEN
               DO 100 I = GPIND, J - 1
                  ZTR = ZERO
                  DO 80 JR = 1, BLKSIZ
                     ZTR = ZTR + WORK( INDRV1+JR )* &
                           DBLE( Z( B1-1+JR, I ) )
   80             CONTINUE
                  DO 90 JR = 1, BLKSIZ
                     WORK( INDRV1+JR ) = WORK( INDRV1+JR ) - &
                                         ZTR*DBLE( Z( B1-1+JR, I ) )
   90             CONTINUE
  100          CONTINUE
            END IF
!
!           Check the infinity norm of the iterate.
!
  110       CONTINUE
            JMAX = IDAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
            NRM = ABS( WORK( INDRV1+JMAX ) )
!
!           Continue for additional iterations after norm reaches
!           stopping criterion.
!
            IF( NRM.LT.DTPCRT ) &
               GO TO 70
            NRMCHK = NRMCHK + 1
            IF( NRMCHK.LT.EXTRA+1 ) &
               GO TO 70
!
            GO TO 130
!
!           If stopping criterion was not satisfied, update info and
!           store eigenvector number in array ifail.
!
  120       CONTINUE
            INFO = INFO + 1
            IFAIL( INFO ) = J
!
!           Accept iterate as jth eigenvector.
!
  130       CONTINUE
            SCL = ONE / DNRM2( BLKSIZ, WORK( INDRV1+1 ), 1 )
            JMAX = IDAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
            IF( WORK( INDRV1+JMAX ).LT.ZERO ) &
               SCL = -SCL
            CALL DSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
  140       CONTINUE
            DO 150 I = 1, N
               Z( I, J ) = CZERO
  150       CONTINUE
            DO 160 I = 1, BLKSIZ
               Z( B1+I-1, J ) = DCMPLX( WORK( INDRV1+I ), ZERO )
  160       CONTINUE
!
!           Save the shift to check eigenvalue spacing at next
!           iteration.
!
            XJM = XJ
!
  170    CONTINUE
  180 CONTINUE
!
      RETURN
!
!     End of ZSTEIN
!
      END
      SUBROUTINE ZSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), WORK( * )
      COMPLEX*16         Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  ZSTEQR computes all eigenvalues and, optionally, eigenvectors of a
!  symmetric tridiagonal matrix using the implicit QL or QR method.
!  The eigenvectors of a full or band complex Hermitian matrix can also
!  be found if ZHETRD or ZHPTRD or ZHBTRD has been used to reduce this
!  matrix to tridiagonal form.
!
!  Arguments
!  =========
!
!  COMPZ   (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only.
!          = 'V':  Compute eigenvalues and eigenvectors of the original
!                  Hermitian matrix.  On entry, Z must contain the
!                  unitary matrix used to reduce the original matrix
!                  to tridiagonal form.
!          = 'I':  Compute eigenvalues and eigenvectors of the
!                  tridiagonal matrix.  Z is initialized to the identity
!                  matrix.
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the diagonal elements of the tridiagonal matrix.
!          On exit, if INFO = 0, the eigenvalues in ascending order.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix.
!          On exit, E has been destroyed.
!
!  Z       (input/output) COMPLEX*16 array, dimension (LDZ, N)
!          On entry, if  COMPZ = 'V', then Z contains the unitary
!          matrix used in the reduction to tridiagonal form.
!          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
!          orthonormal eigenvectors of the original Hermitian matrix,
!          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!          of the symmetric tridiagonal matrix.
!          If COMPZ = 'N', then Z is not referenced.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          eigenvectors are desired, then  LDZ >= max(1,N).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
!          If COMPZ = 'N', then WORK is not referenced.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  the algorithm has failed to find all the eigenvalues in
!                a total of 30*N iterations; if INFO = i, then i
!                elements of E have not converged to zero; on exit, D
!                and E contain the elements of a symmetric tridiagonal
!                matrix which is unitarily similar to the original
!                matrix.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, &
                         THREE = 3.0D0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ), &
                         CONE = ( 1.0D0, 0.0D0 ) )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND, &
                         LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1, &
                         NM1, NMAXIT
      DOUBLE PRECISION   ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2, &
                         S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASRT, XERBLA, &
                         ZLASET, ZLASR, ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1, &
               N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSTEQR', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.EQ.2 ) &
            Z( 1, 1 ) = CONE
         RETURN
      END IF
!
!     Determine the unit roundoff and over/underflow thresholds.
!
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
!
!     Compute the eigenvalues and eigenvectors of the tridiagonal
!     matrix.
!
      IF( ICOMPZ.EQ.2 ) &
         CALL ZLASET( 'Full', N, N, CZERO, CONE, Z, LDZ )
!
      NMAXIT = N*MAXIT
      JTOT = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      L1 = 1
      NM1 = N - 1
!
   10 CONTINUE
      IF( L1.GT.N ) &
         GO TO 160
      IF( L1.GT.1 ) &
         E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO ) &
               GO TO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+ &
                1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
!
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L ) &
         GO TO 10
!
!     Scale submatrix in rows and columns L to LEND
!
      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO ) &
         GO TO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N, &
                      INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N, &
                      INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N, &
                      INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N, &
                      INFO )
      END IF
!
!     Choose between QL and QR iteration
!
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
!
      IF( LEND.GT.L ) THEN
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+ &
                   SAFMIN )GO TO 60
   50       CONTINUE
         END IF
!
         M = LEND
!
   60    CONTINUE
         IF( M.LT.LEND ) &
            E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L ) &
            GO TO 80
!
!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.
!
         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL ZLASR( 'R', 'V', 'B', N, 2, WORK( L ), &
                           WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND ) &
               GO TO 40
            GO TO 140
         END IF
!
         IF( JTOT.EQ.NMAXIT ) &
            GO TO 140
         JTOT = JTOT + 1
!
!        Form shift.
!
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
!
         S = ONE
         C = ONE
         P = ZERO
!
!        Inner loop
!
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M-1 ) &
               E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF
!
   70    CONTINUE
!
!        If eigenvectors are desired, then apply saved rotations.
!
         IF( ICOMPZ.GT.0 ) THEN
            MM = M - L + 1
            CALL ZLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ), &
                        Z( 1, L ), LDZ )
         END IF
!
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
!
!        Eigenvalue found.
!
   80    CONTINUE
         D( L ) = P
!
         L = L + 1
         IF( L.LE.LEND ) &
            GO TO 40
         GO TO 140
!
      ELSE
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
   90    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+ &
                   SAFMIN )GO TO 110
  100       CONTINUE
         END IF
!
         M = LEND
!
  110    CONTINUE
         IF( M.GT.LEND ) &
            E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L ) &
            GO TO 130
!
!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.
!
         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL ZLASR( 'R', 'V', 'F', N, 2, WORK( M ), &
                           WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND ) &
               GO TO 90
            GO TO 140
         END IF
!
         IF( JTOT.EQ.NMAXIT ) &
            GO TO 140
         JTOT = JTOT + 1
!
!        Form shift.
!
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
!
         S = ONE
         C = ONE
         P = ZERO
!
!        Inner loop
!
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M ) &
               E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF
!
  120    CONTINUE
!
!        If eigenvectors are desired, then apply saved rotations.
!
         IF( ICOMPZ.GT.0 ) THEN
            MM = L - M + 1
            CALL ZLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ), &
                        Z( 1, M ), LDZ )
         END IF
!
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
!
!        Eigenvalue found.
!
  130    CONTINUE
         D( L ) = P
!
         L = L - 1
         IF( L.GE.LEND ) &
            GO TO 90
         GO TO 140
!
      END IF
!
!     Undo scaling if necessary
!
  140 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1, &
                      D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ), &
                      N, INFO )
      ELSE IF( ISCALE.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1, &
                      D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ), &
                      N, INFO )
      END IF
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
      IF( JTOT.EQ.NMAXIT ) THEN
         DO 150 I = 1, N - 1
            IF( E( I ).NE.ZERO ) &
               INFO = INFO + 1
  150    CONTINUE
         RETURN
      END IF
      GO TO 10
!
!     Order eigenvalues and eigenvectors.
!
  160 CONTINUE
      IF( ICOMPZ.EQ.0 ) THEN
!
!        Use Quick Sort
!
         CALL DLASRT( 'I', N, D, INFO )
!
      ELSE
!
!        Use Selection Sort to minimize swaps of eigenvectors
!
         DO 180 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 170 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
  170       CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL ZSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
  180    CONTINUE
      END IF
      RETURN
!
!     End of ZSTEQR
!
      END
      SUBROUTINE ZUNG2L( M, N, K, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZUNG2L generates an m by n complex matrix Q with orthonormal columns,
!  which is defined as the last n columns of a product of k elementary
!  reflectors of order m
!
!        Q  =  H(k) . . . H(2) H(1)
!
!  as returned by ZGEQLF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the (n-k+i)-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by ZGEQLF in the last k columns of its array
!          argument A.
!          On exit, the m-by-n matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) COMPLEX*16 array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by ZGEQLF.
!
!  WORK    (workspace) COMPLEX*16 array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), &
                         ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, II, J, L
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF, ZSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNG2L', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) &
         RETURN
!
!     Initialise columns 1:n-k to columns of the unit matrix
!
      DO 20 J = 1, N - K
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( M-N+J, J ) = ONE
   20 CONTINUE
!
      DO 40 I = 1, K
         II = N - K + I
!
!        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
!
         A( M-N+II, II ) = ONE
         CALL ZLARF( 'Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A, &
                     LDA, WORK )
         CALL ZSCAL( M-N+II-1, -TAU( I ), A( 1, II ), 1 )
         A( M-N+II, II ) = ONE - TAU( I )
!
!        Set A(m-k+i+1:m,n-k+i) to zero
!
         DO 30 L = M - N + II + 1, M
            A( L, II ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
!
!     End of ZUNG2L
!
      END
      SUBROUTINE ZUNG2R( M, N, K, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZUNG2R generates an m by n complex matrix Q with orthonormal columns,
!  which is defined as the first n columns of a product of k elementary
!  reflectors of order m
!
!        Q  =  H(1) H(2) . . . H(k)
!
!  as returned by ZGEQRF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by ZGEQRF in the first k columns of its array
!          argument A.
!          On exit, the m by n matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) COMPLEX*16 array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by ZGEQRF.
!
!  WORK    (workspace) COMPLEX*16 array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), &
                         ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, L
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF, ZSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNG2R', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) &
         RETURN
!
!     Initialise columns k+1:n to columns of the unit matrix
!
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
!
      DO 40 I = K, 1, -1
!
!        Apply H(i) to A(i:m,i:n) from the left
!
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL ZLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), &
                        A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M ) &
            CALL ZSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
!
!        Set A(1:i-1,i) to zero
!
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
!
!     End of ZUNG2R
!
      END
      SUBROUTINE ZUPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDQ, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         AP( * ), Q( LDQ, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZUPGTR generates a complex unitary matrix Q which is defined as the
!  product of n-1 elementary reflectors H(i) of order n, as returned by
!  ZHPTRD using packed storage:
!
!  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
!
!  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U': Upper triangular packed storage used in previous
!                 call to ZHPTRD;
!          = 'L': Lower triangular packed storage used in previous
!                 call to ZHPTRD.
!
!  N       (input) INTEGER
!          The order of the matrix Q. N >= 0.
!
!  AP      (input) COMPLEX*16 array, dimension (N*(N+1)/2)
!          The vectors which define the elementary reflectors, as
!          returned by ZHPTRD.
!
!  TAU     (input) COMPLEX*16 array, dimension (N-1)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by ZHPTRD.
!
!  Q       (output) COMPLEX*16 array, dimension (LDQ,N)
!          The N-by-N unitary matrix Q.
!
!  LDQ     (input) INTEGER
!          The leading dimension of the array Q. LDQ >= max(1,N).
!
!  WORK    (workspace) COMPLEX*16 array, dimension (N-1)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), &
                         CONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IINFO, IJ, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZUNG2L, ZUNG2R
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUPGTR', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF( UPPER ) THEN
!
!        Q was determined by a call to ZHPTRD with UPLO = 'U'
!
!        Unpack the vectors which define the elementary reflectors and
!        set the last row and column of Q equal to those of the unit
!        matrix
!
         IJ = 2
         DO 20 J = 1, N - 1
            DO 10 I = 1, J - 1
               Q( I, J ) = AP( IJ )
               IJ = IJ + 1
   10       CONTINUE
            IJ = IJ + 2
            Q( N, J ) = CZERO
   20    CONTINUE
         DO 30 I = 1, N - 1
            Q( I, N ) = CZERO
   30    CONTINUE
         Q( N, N ) = CONE
!
!        Generate Q(1:n-1,1:n-1)
!
         CALL ZUNG2L( N-1, N-1, N-1, Q, LDQ, TAU, WORK, IINFO )
!
      ELSE
!
!        Q was determined by a call to ZHPTRD with UPLO = 'L'.
!
!        Unpack the vectors which define the elementary reflectors and
!        set the first row and column of Q equal to those of the unit
!        matrix
!
         Q( 1, 1 ) = CONE
         DO 40 I = 2, N
            Q( I, 1 ) = CZERO
   40    CONTINUE
         IJ = 3
         DO 60 J = 2, N
            Q( 1, J ) = CZERO
            DO 50 I = J + 1, N
               Q( I, J ) = AP( IJ )
               IJ = IJ + 1
   50       CONTINUE
            IJ = IJ + 2
   60    CONTINUE
         IF( N.GT.1 ) THEN
!
!           Generate Q(2:n,2:n)
!
            CALL ZUNG2R( N-1, N-1, N-1, Q( 2, 2 ), LDQ, TAU, WORK, &
                         IINFO )
         END IF
      END IF
      RETURN
!
!     End of ZUPGTR
!
      END
      SUBROUTINE ZUPMTR( SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK, &
                         INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS, UPLO
      INTEGER            INFO, LDC, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         AP( * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  ZUPMTR overwrites the general complex M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'C':      Q**H * C       C * Q**H
!
!  where Q is a complex unitary matrix of order nq, with nq = m if
!  SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
!  nq-1 elementary reflectors, as returned by ZHPTRD using packed
!  storage:
!
!  if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);
!
!  if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**H from the Left;
!          = 'R': apply Q or Q**H from the Right.
!
!  UPLO    (input) CHARACTER*1
!          = 'U': Upper triangular packed storage used in previous
!                 call to ZHPTRD;
!          = 'L': Lower triangular packed storage used in previous
!                 call to ZHPTRD.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'C':  Conjugate transpose, apply Q**H.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  AP      (input) COMPLEX*16 array, dimension
!                               (M*(M+1)/2) if SIDE = 'L'
!                               (N*(N+1)/2) if SIDE = 'R'
!          The vectors which define the elementary reflectors, as
!          returned by ZHPTRD.  AP is modified by the routine but
!          restored on exit.
!
!  TAU     (input) COMPLEX*16 array, dimension (M-1) if SIDE = 'L'
!                                     or (N-1) if SIDE = 'R'
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by ZHPTRD.
!
!  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) COMPLEX*16 array, dimension
!                                   (N) if SIDE = 'L'
!                                   (M) if SIDE = 'R'
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            FORWRD, LEFT, NOTRAN, UPPER
      INTEGER            I, I1, I2, I3, IC, II, JC, MI, NI, NQ
      COMPLEX*16         AII, TAUI
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      UPPER = LSAME( UPLO, 'U' )
!
!     NQ is the order of Q
!
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUPMTR', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) &
         RETURN
!
      IF( UPPER ) THEN
!
!        Q was determined by a call to ZHPTRD with UPLO = 'U'
!
         FORWRD = ( LEFT .AND. NOTRAN ) .OR. &
                  ( .NOT.LEFT .AND. .NOT.NOTRAN )
!
         IF( FORWRD ) THEN
            I1 = 1
            I2 = NQ - 1
            I3 = 1
            II = 2
         ELSE
            I1 = NQ - 1
            I2 = 1
            I3 = -1
            II = NQ*( NQ+1 ) / 2 - 1
         END IF
!
         IF( LEFT ) THEN
            NI = N
         ELSE
            MI = M
         END IF
!
         DO 10 I = I1, I2, I3
            IF( LEFT ) THEN
!
!              H(i) or H(i)' is applied to C(1:i,1:n)
!
               MI = I
            ELSE
!
!              H(i) or H(i)' is applied to C(1:m,1:i)
!
               NI = I
            END IF
!
!           Apply H(i) or H(i)'
!
            IF( NOTRAN ) THEN
               TAUI = TAU( I )
            ELSE
               TAUI = DCONJG( TAU( I ) )
            END IF
            AII = AP( II )
            AP( II ) = ONE
            CALL ZLARF( SIDE, MI, NI, AP( II-I+1 ), 1, TAUI, C, LDC, &
                        WORK )
            AP( II ) = AII
!
            IF( FORWRD ) THEN
               II = II + I + 2
            ELSE
               II = II - I - 1
            END IF
   10    CONTINUE
      ELSE
!
!        Q was determined by a call to ZHPTRD with UPLO = 'L'.
!
         FORWRD = ( LEFT .AND. .NOT.NOTRAN ) .OR. &
                  ( .NOT.LEFT .AND. NOTRAN )
!
         IF( FORWRD ) THEN
            I1 = 1
            I2 = NQ - 1
            I3 = 1
            II = 2
         ELSE
            I1 = NQ - 1
            I2 = 1
            I3 = -1
            II = NQ*( NQ+1 ) / 2 - 1
         END IF
!
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
!
         DO 20 I = I1, I2, I3
            AII = AP( II )
            AP( II ) = ONE
            IF( LEFT ) THEN
!
!              H(i) or H(i)' is applied to C(i+1:m,1:n)
!
               MI = M - I
               IC = I + 1
            ELSE
!
!              H(i) or H(i)' is applied to C(1:m,i+1:n)
!
               NI = N - I
               JC = I + 1
            END IF
!
!           Apply H(i) or H(i)'
!
            IF( NOTRAN ) THEN
               TAUI = TAU( I )
            ELSE
               TAUI = DCONJG( TAU( I ) )
            END IF
            CALL ZLARF( SIDE, MI, NI, AP( II ), 1, TAUI, C( IC, JC ), &
                        LDC, WORK )
            AP( II ) = AII
!
            IF( FORWRD ) THEN
               II = II + NQ - I + 1
            ELSE
               II = II - NQ + I - 2
            END IF
   20    CONTINUE
      END IF
      RETURN
!
!     End of ZUPMTR
!
      END

      SUBROUTINE XERBLA( SRNAME, INFO )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      WRITE( *, FMT = 9999 )SRNAME, INFO
!
      STOP
!
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ', &
            'an illegal value' )
!
!     End of XERBLA
!
      END

      LOGICAL          FUNCTION LSAME( CA, CB )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          CA, CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = CA.EQ.CB
      IF( LSAME ) &
         RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR. &
             INTA.GE.145 .AND. INTA.LE.153 .OR. &
             INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR. &
             INTB.GE.145 .AND. INTB.LE.153 .OR. &
             INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
!
!     RETURN
!
!     End of LSAME
!
      END