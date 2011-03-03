      SUBROUTINE DB2INK(X,NX,Y,NY,FCN,LDF,KX,KY,TX,TY,BCOEF,WORK,IFLAG)
!***BEGIN PROLOGUE  DB2INK
!***DATE WRITTEN   25 MAY 1982
!***REVISION DATE  25 MAY 1982
!***REVISION HISTORY  (YYMMDD)
!   000330  Modified array declarations.  (JEC)
!
!***CATEGORY NO.  E1A
!***KEYWORDS  INTERPOLATION, TWO-DIMENSIONS, GRIDDED DATA, SPLINES,
!             PIECEWISE POLYNOMIALS
!***AUTHOR  BOISVERT, RONALD, NBS
!             SCIENTIFIC COMPUTING DIVISION
!             NATIONAL BUREAU OF STANDARDS
!             WASHINGTON, DC 20234
!***PURPOSE  DOUBLE PRECISION VERSION OF B2INK.
!            DB2INK DETERMINES A PIECEWISE POLYNOMIAL FUNCTION THAT
!            INTERPOLATES TWO-DIMENSIONAL GRIDDED DATA. USERS SPECIFY
!            THE POLYNOMIAL ORDER (DEGREE+1) OF THE INTERPOLANT AND
!            (OPTIONALLY) THE KNOT SEQUENCE.
!***DESCRIPTION
!
!   DB2INK determines the parameters of a  function  that  interpolates
!   the two-dimensional gridded data (X(i),Y(j),FCN(i,j)) for i=1,..,NX
!   and j=1,..,NY. The interpolating function and its  derivatives  may
!   subsequently be evaluated by the function DB2VAL.
!
!   The interpolating  function  is  a  piecewise  polynomial  function
!   represented as a tensor product of one-dimensional  B-splines.  The
!   form of this function is
!
!                          NX   NY
!              S(x,y)  =  SUM  SUM  a   U (x) V (y)
!                         i=1  j=1   ij  i     j
!
!   where the functions U(i)  and  V(j)  are  one-dimensional  B-spline
!   basis functions. The coefficients a(i,j) are chosen so that
!
!         S(X(i),Y(j)) = FCN(i,j)   for i=1,..,NX and j=1,..,NY
!
!   Note that  for  each  fixed  value  of  y  S(x,y)  is  a  piecewise
!   polynomial function of x alone, and for each fixed value of x  S(x,
!   y) is a piecewise polynomial function of y alone. In one  dimension
!   a piecewise polynomial may  be  created  by  partitioning  a  given
!   interval into subintervals and defining a distinct polynomial piece
!   on each one. The points where adjacent subintervals meet are called
!   knots. Each of the functions U(i) and V(j)  above  is  a  piecewise
!   polynomial.
!
!   Users of DB2INK choose  the  order  (degree+1)  of  the  polynomial
!   pieces used to define the piecewise polynomial in each of the x and
!   y directions (KX and KY). Users also  may  define  their  own  knot
!   sequence in x and y separately (TX and TY).  If  IFLAG=0,  however,
!   DB2INK will choose sequences of knots that result  in  a  piecewise
!   polynomial interpolant with KX-2 continuous partial derivatives  in
!   x and KY-2 continuous partial derivatives in y. (KX knots are taken
!   near each endpoint in the x direction,  not-a-knot  end  conditions
!   are used, and the remaining knots are placed at data points  if  KX
!   is even or at midpoints between data points if KX  is  odd.  The  y
!   direction is treated similarly.)
!
!   After a call to DB2INK, all information  necessary  to  define  the
!   interpolating function are contained in the parameters NX, NY,  KX,
!   KY, TX, TY, and BCOEF. These quantities should not be altered until
!   after the last call of the evaluation routine DB2VAL.
!
!
!   I N P U T
!   ---------
!
!   X       Double precision 1D array (size NX)
!           Array of x abcissae. Must be strictly increasing.
!
!   NX      Integer scalar (.GE. 3)
!           Number of x abcissae.
!
!   Y       Double precision 1D array (size NY)
!           Array of y abcissae. Must be strictly increasing.
!
!   NY      Integer scalar (.GE. 3)
!           Number of y abcissae.
!
!   FCN     Double precision 2D array (size LDF by NY)
!           Array of function values to interpolate. FCN(I,J) should
!           contain the function value at the point (X(I),Y(J))
!
!   LDF     Integer scalar (.GE. NX)
!           The actual leading dimension of FCN used in the calling
!           calling program.
!
!   KX      Integer scalar (.GE. 2, .LT. NX)
!           The order of spline pieces in x.
!           (Order = polynomial degree + 1)
!
!   KY      Integer scalar (.GE. 2, .LT. NY)
!           The order of spline pieces in y.
!           (Order = polynomial degree + 1)
!
!
!   I N P U T   O R   O U T P U T
!   -----------------------------
!
!   TX      Double precision 1D array (size NX+KX)
!           The knots in the x direction for the spline interpolant.
!           If IFLAG=0 these are chosen by DB2INK.
!           If IFLAG=1 these are specified by the user.
!                      (Must be non-decreasing.)
!
!   TY      Double precision 1D array (size NY+KY)
!           The knots in the y direction for the spline interpolant.
!           If IFLAG=0 these are chosen by DB2INK.
!           If IFLAG=1 these are specified by the user.
!                      (Must be non-decreasing.)
!
!
!   O U T P U T
!   -----------
!
!   BCOEF   Double precision 2D array (size NX by NY)
!           Array of coefficients of the B-spline interpolant.
!           This may be the same array as FCN.
!
!
!   M I S C E L L A N E O U S
!   -------------------------
!
!   WORK    Double precision 1D array (size NX*NY + max( 2*KX*(NX+1),
!                                             2*KY*(NY+1) ))
!           Array of working storage.
!
!   IFLAG   Integer scalar.
!           On input:  0 == knot sequence chosen by DB2INK
!                      1 == knot sequence chosen by user.
!           On output: 1 == successful execution
!                      2 == IFLAG out of range
!                      3 == NX out of range
!                      4 == KX out of range
!                      5 == X not strictly increasing
!                      6 == TX not non-decreasing
!                      7 == NY out of range
!                      8 == KY out of range
!                      9 == Y not strictly increasing
!                     10 == TY not non-decreasing
!
!***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES,
!                 SPRINGER-VERLAG, NEW YORK, 1978.
!               CARL DE BOOR, EFFICIENT COMPUTER MANIPULATION OF TENSOR
!                 PRODUCTS, ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!                 VOL. 5 (1979), PP. 173-182.
!***ROUTINES CALLED  DBTPCF,DBKNOT
!***END PROLOGUE  DB2INK
!
!  ------------
!  DECLARATIONS
!  ------------
!
!  PARAMETERS
!
      INTEGER NX, NY, LDF, KX, KY, IFLAG
      DOUBLE PRECISION &
          X(NX), Y(NY), FCN(LDF,NY), TX(*), TY(*), BCOEF(NX,NY),  &
          WORK(*)
!
!  LOCAL VARIABLES
!
      INTEGER I, IW, NPK
!
!  -----------------------
!  CHECK VALIDITY OF INPUT
!  -----------------------
!
!***FIRST EXECUTABLE STATEMENT
      IF ((IFLAG .LT. 0) .OR. (IFLAG .GT. 1))  GO TO 920
      IF (NX .LT. 3)  GO TO 930
      IF (NY .LT. 3)  GO TO 970
      IF ((KX .LT. 2) .OR. (KX .GE. NX))  GO TO 940
      IF ((KY .LT. 2) .OR. (KY .GE. NY))  GO TO 980
      DO 10 I=2,NX
         IF (X(I) .LE. X(I-1))  GO TO 950
   10 CONTINUE
      DO 20 I=2,NY
         IF (Y(I) .LE. Y(I-1))  GO TO 990
   20 CONTINUE
      IF (IFLAG .EQ. 0)  GO TO 50
         NPK = NX + KX
         DO 30 I=2,NPK
            IF (TX(I) .LT. TX(I-1))  GO TO 960
   30    CONTINUE
         NPK = NY + KY
         DO 40 I=2,NPK
            IF (TY(I) .LT. TY(I-1))  GO TO 1000
   40    CONTINUE
   50 CONTINUE
!
!  ------------
!  CHOOSE KNOTS
!  ------------
!
      IF (IFLAG .NE. 0)  GO TO 100
         CALL DBKNOT(X,NX,KX,TX)
         CALL DBKNOT(Y,NY,KY,TY)
  100 CONTINUE
!
!  -------------------------------
!  CONSTRUCT B-SPLINE COEFFICIENTS
!  -------------------------------
!
      IFLAG = 1
      IW = NX*NY + 1
      CALL DBTPCF(X,NX,FCN,LDF,NY,TX,KX,WORK,WORK(IW))
      CALL DBTPCF(Y,NY,WORK,NY,NX,TY,KY,BCOEF,WORK(IW))
      GO TO 9999
!
!  -----
!  EXITS
!  -----
!
  920 CONTINUE
      CALL XERRWV('DB2INK -  IFLAG=I1 IS OUT OF RANGE.',  &
                 36,2,1,1,IFLAG,I2,0,R1,R2)
      IFLAG = 2
      GO TO 9999
!
  930 CONTINUE
      IFLAG = 3
      CALL XERRWV('DB2INK -  NX=I1 IS OUT OF RANGE.',  &
                 32,IFLAG,1,1,NX,I2,0,R1,R2)
      GO TO 9999
!
  940 CONTINUE
      IFLAG = 4
      CALL XERRWV('DB2INK -  KX=I1 IS OUT OF RANGE.',   &
                 32,IFLAG,1,1,KX,I2,0,R1,R2)
      GO TO 9999
!
  950 CONTINUE
      IFLAG = 5
      CALL XERRWV('DB2INK -  X ARRAY MUST BE STRICTLY INCREASING.',  &
                 46,IFLAG,1,0,I1,I2,0,R1,R2)
      GO TO 9999
!
  960 CONTINUE
      IFLAG = 6
      CALL XERRWV('DB2INK -  TX ARRAY MUST BE NON-DECREASING.',  &
                 42,IFLAG,1,0,I1,I2,0,R1,R2)
      GO TO 9999
!
  970 CONTINUE
      IFLAG = 7
      CALL XERRWV('DB2INK -  NY=I1 IS OUT OF RANGE.', &
                 32,IFLAG,1,1,NY,I2,0,R1,R2)
      GO TO 9999
!
  980 CONTINUE
      IFLAG = 8
      CALL XERRWV('DB2INK -  KY=I1 IS OUT OF RANGE.',  &
                 32,IFLAG,1,1,KY,I2,0,R1,R2)
      GO TO 9999
!
  990 CONTINUE
      IFLAG = 9
      CALL XERRWV('DB2INK -  Y ARRAY MUST BE STRICTLY INCREASING.',  &
                 46,IFLAG,1,0,I1,I2,0,R1,R2)
      GO TO 9999
!
 1000 CONTINUE
      IFLAG = 10
      CALL XERRWV('DB2INK -  TY ARRAY MUST BE NON-DECREASING.', &
                 42,IFLAG,1,0,I1,I2,0,R1,R2)
      GO TO 9999
!
 9999 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DB2VAL(XVAL,YVAL,IDX,IDY,TX,TY,NX,NY, &
       KX,KY,BCOEF,WORK)
!***BEGIN PROLOGUE  DB2VAL
!***DATE WRITTEN   25 MAY 1982
!***REVISION DATE  25 MAY 1982
!***REVISION HISTORY  (YYMMDD)
!   000330  Modified array declarations.  (JEC)
!
!***CATEGORY NO.  E1A
!***KEYWORDS  INTERPOLATION, TWO-DIMENSIONS, GRIDDED DATA, SPLINES,
!             PIECEWISE POLYNOMIALS
!***AUTHOR  BOISVERT, RONALD, NBS
!             SCIENTIFIC COMPUTING DIVISION
!             NATIONAL BUREAU OF STANDARDS
!             WASHINGTON, DC 20234
!***PURPOSE  DB2VAL EVALUATES THE PIECEWISE POLYNOMIAL INTERPOLATING
!            FUNCTION CONSTRUCTED BY THE ROUTINE DB2INK OR ONE OF ITS
!            PARTIAL DERIVATIVES.
!            DOUBLE PRECISION VERSION OF B2VAL.
!***DESCRIPTION
!
!   DB2VAL  evaluates   the   tensor   product   piecewise   polynomial
!   interpolant constructed  by  the  routine  DB2INK  or  one  of  its
!   derivatives at the point (XVAL,YVAL). To evaluate  the  interpolant
!   itself, set IDX=IDY=0, to evaluate the first partial  with  respect
!   to x, set IDX=1,IDY=0, and so on.
!
!   DB2VAL returns 0.0E0 if (XVAL,YVAL) is out of range. That is, if
!            XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR.
!            YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY)
!   If the knots TX  and  TY  were  chosen  by  DB2INK,  then  this  is
!   equivalent to
!            XVAL.LT.X(1) .OR. XVAL.GT.X(NX)+EPSX .OR.
!            YVAL.LT.Y(1) .OR. YVAL.GT.Y(NY)+EPSY
!   where EPSX = 0.1*(X(NX)-X(NX-1)) and EPSY = 0.1*(Y(NY)-Y(NY-1)).
!
!   The input quantities TX, TY, NX, NY, KX, KY, and  BCOEF  should  be
!   unchanged since the last call of DB2INK.
!
!
!   I N P U T
!   ---------
!
!   XVAL    Double precision scalar
!           X coordinate of evaluation point.
!
!   YVAL    Double precision scalar
!           Y coordinate of evaluation point.
!
!   IDX     Integer scalar
!           X derivative of piecewise polynomial to evaluate.
!
!   IDY     Integer scalar
!           Y derivative of piecewise polynomial to evaluate.
!
!   TX      Double precision 1D array (size NX+KX)
!           Sequence of knots defining the piecewise polynomial in
!           the x direction.  (Same as in last call to DB2INK.)
!
!   TY      Double precision 1D array (size NY+KY)
!           Sequence of knots defining the piecewise polynomial in
!           the y direction.  (Same as in last call to DB2INK.)
!
!   NX      Integer scalar
!           The number of interpolation points in x.
!           (Same as in last call to DB2INK.)
!
!   NY      Integer scalar
!           The number of interpolation points in y.
!           (Same as in last call to DB2INK.)
!
!   KX      Integer scalar
!           Order of polynomial pieces in x.
!           (Same as in last call to DB2INK.)
!
!   KY      Integer scalar
!           Order of polynomial pieces in y.
!           (Same as in last call to DB2INK.)
!
!   BCOEF   Double precision 2D array (size NX by NY)
!           The B-spline coefficients computed by DB2INK.
!
!   WORK    Double precision 1D array (size 3*max(KX,KY) + KY)
!           A working storage array.
!
!***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES,
!                 SPRINGER-VERLAG, NEW YORK, 1978.
!***ROUTINES CALLED  DINTRV,DBVALU
!***END PROLOGUE  DB2VAL
!
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
!   MODIFICATION
!   ------------
!
!   ADDED CHECK TO SEE IF X OR Y IS OUT OF RANGE, IF SO, RETURN 0.0
!
!   R.F. BOISVERT, NIST
!   22 FEB 00
!
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!  ------------
!  DECLARATIONS
!  ------------
!
!  PARAMETERS
!
      INTEGER IDX, IDY, NX, NY, KX, KY
      DOUBLE PRECISION  XVAL, YVAL, TX(*), TY(*), BCOEF(NX,NY), WORK(*)
!
!  LOCAL VARIABLES
!
      INTEGER  ILOY, INBVX, INBV, K, LEFTY, MFLAG, KCOL, IW
      DOUBLE PRECISION  DBVALU
!
!      DATA ILOY /1/,  INBVX /1/
      ILOY=1
      INBVX=1
!      SAVE ILOY    ,  INBVX
!
!
!***FIRST EXECUTABLE STATEMENT
      DB2VAL = 0.0D0
!  NEXT STATEMENT - RFB MOD
      IF (XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR. &
         YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY))      GO TO 100
      CALL DINTRV(TY,NY+KY,YVAL,ILOY,LEFTY,MFLAG)
      IF (MFLAG .NE. 0)  GO TO 100
         IW = KY + 1
         KCOL = LEFTY - KY
         DO 50 K=1,KY
            KCOL = KCOL + 1
            WORK(K) = DBVALU(TX,BCOEF(1,KCOL),NX,KX,IDX,XVAL,INBVX, &
                           WORK(IW))
   50    CONTINUE
         INBV = 1
         KCOL = LEFTY - KY + 1
         DB2VAL = DBVALU(TY(KCOL),WORK,KY,KY,IDY,YVAL,INBV,WORK(IW))
  100 CONTINUE
      RETURN
      END
      SUBROUTINE DBINTK(X,Y,T,N,K,BCOEF,Q,WORK)
!***BEGIN PROLOGUE  DBINTK
!***DATE WRITTEN   800901   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***REVISION HISTORY  (YYMMDD)
!   000330  Modified array declarations.  (JEC)
!
!***CATEGORY NO.  E1A
!***KEYWORDS  B-SPLINE,DATA FITTING,DOUBLE PRECISION,INTERPOLATION,
!             SPLINE
!***AUTHOR  AMOS, D. E., (SNLA)
!***PURPOSE  Produces the B-spline coefficients, BCOEF, of the
!            B-spline of order K with knots T(I), I=1,...,N+K, which
!            takes on the value Y(I) at X(I), I=1,...,N.
!***DESCRIPTION
!
!     Written by Carl de Boor and modified by D. E. Amos
!
!     References
!
!         A Practical Guide to Splines by C. de Boor, Applied
!         Mathematics Series 27, Springer, 1979.
!
!     Abstract    **** a double precision routine ****
!
!         DBINTK is the SPLINT routine of the reference.
!
!         DBINTK produces the B-spline coefficients, BCOEF, of the
!         B-spline of order K with knots T(I), I=1,...,N+K, which
!         takes on the value Y(I) at X(I), I=1,...,N.  The spline or
!         any of its derivatives can be evaluated by calls to DBVALU.
!
!         The I-th equation of the linear system A*BCOEF = B for the
!         coefficients of the interpolant enforces interpolation at
!         X(I), I=1,...,N.  Hence, B(I) = Y(I), for all I, and A is
!         a band matrix with 2K-1 bands if A is invertible.  The matrix
!         A is generated row by row and stored, diagonal by diagonal,
!         in the rows of Q, with the main diagonal going into row K.
!         The banded system is then solved by a call to DBNFAC (which
!         constructs the triangular factorization for A and stores it
!         again in Q), followed by a call to DBNSLV (which then
!         obtains the solution BCOEF by substitution).  DBNFAC does no
!         pivoting, since the total positivity of the matrix A makes
!         this unnecessary.  The linear system to be solved is
!         (theoretically) invertible if and only if
!                 T(I) .LT. X(I) .LT. T(I+K),        for all I.
!         Equality is permitted on the left for I=1 and on the right
!         for I=N when K knots are used at X(1) or X(N).  Otherwise,
!         violation of this condition is certain to lead to an error.
!
!         DBINTK calls DBSPVN, DBNFAC, DBNSLV, XERROR
!
!     Description of Arguments
!
!         Input       X,Y,T are double precision
!           X       - vector of length N containing data point abscissa
!                     in strictly increasing order.
!           Y       - corresponding vector of length N containing data
!                     point ordinates.
!           T       - knot vector of length N+K
!                     Since T(1),..,T(K) .LE. X(1) and T(N+1),..,T(N+K)
!                     .GE. X(N), this leaves only N-K knots (not nec-
!                     essarily X(I) values) interior to (X(1),X(N))
!           N       - number of data points, N .GE. K
!           K       - order of the spline, K .GE. 1
!
!         Output      BCOEF,Q,WORK are double precision
!           BCOEF   - a vector of length N containing the B-spline
!                     coefficients
!           Q       - a work vector of length (2*K-1)*N, containing
!                     the triangular factorization of the coefficient
!                     matrix of the linear system being solved.  The
!                     coefficients for the interpolant of an
!                     additional data set (X(I),yY(I)), I=1,...,N
!                     with the same abscissa can be obtained by loading
!                     YY into BCOEF and then executing
!                         CALL DBNSLV(Q,2K-1,N,K-1,K-1,BCOEF)
!           WORK    - work vector of length 2*K
!
!     Error Conditions
!         Improper input is a fatal error
!         Singular system of equations is a fatal error
!***REFERENCES  C. DE BOOR, *A PRACTICAL GUIDE TO SPLINES*, APPLIED
!                 MATHEMATICS SERIES 27, SPRINGER, 1979.
!               D.E. AMOS, *COMPUTATION WITH SPLINES AND B-SPLINES*,
!                 SAND78-1968,SANDIA LABORATORIES,MARCH,1979.
!***ROUTINES CALLED  DBNFAC,DBNSLV,DBSPVN,XERROR
!***END PROLOGUE  DBINTK
!
!
      INTEGER IFLAG, IWORK, K, N, I, ILP1MX, J, JJ, KM1, KPKM2, LEFT, &
       LENQ, NP1
      DOUBLE PRECISION BCOEF(N), Y(N), Q(*), T(*), X(N), XI, WORK(*)
!     DIMENSION Q(2*K-1,N), T(N+K)
!***FIRST EXECUTABLE STATEMENT  DBINTK
      IF(K.LT.1) GO TO 100
      IF(N.LT.K) GO TO 105
      JJ = N - 1
      IF(JJ.EQ.0) GO TO 6
      DO 5 I=1,JJ
      IF(X(I).GE.X(I+1)) GO TO 110
    5 CONTINUE
    6 CONTINUE
      NP1 = N + 1
      KM1 = K - 1
      KPKM2 = 2*KM1
      LEFT = K
!                ZERO OUT ALL ENTRIES OF Q
      LENQ = N*(K+KM1)
      DO 10 I=1,LENQ
        Q(I) = 0.0D0
   10 CONTINUE
!
!  ***   LOOP OVER I TO CONSTRUCT THE  N  INTERPOLATION EQUATIONS
      DO 50 I=1,N
        XI = X(I)
        ILP1MX = MIN0(I+K,NP1)
!        *** FIND  LEFT  IN THE CLOSED INTERVAL (I,I+K-1) SUCH THAT
!                T(LEFT) .LE. X(I) .LT. T(LEFT+1)
!        MATRIX IS SINGULAR IF THIS IS NOT POSSIBLE
        LEFT = MAX0(LEFT,I)
        IF (XI.LT.T(LEFT)) GO TO 80
   20   IF (XI.LT.T(LEFT+1)) GO TO 30
        LEFT = LEFT + 1
        IF (LEFT.LT.ILP1MX) GO TO 20
        LEFT = LEFT - 1
        IF (XI.GT.T(LEFT+1)) GO TO 80
!        *** THE I-TH EQUATION ENFORCES INTERPOLATION AT XI, HENCE
!        A(I,J) = B(J,K,T)(XI), ALL J. ONLY THE  K  ENTRIES WITH  J =
!        LEFT-K+1,...,LEFT ACTUALLY MIGHT BE NONZERO. THESE  K  NUMBERS
!        ARE RETURNED, IN  BCOEF (USED FOR TEMP.STORAGE HERE), BY THE
!        FOLLOWING
   30   CALL DBSPVN(T, K, K, 1, XI, LEFT, BCOEF, WORK, IWORK)
!        WE THEREFORE WANT  BCOEF(J) = B(LEFT-K+J)(XI) TO GO INTO
!        A(I,LEFT-K+J), I.E., INTO  Q(I-(LEFT+J)+2*K,(LEFT+J)-K) SINCE
!        A(I+J,J)  IS TO GO INTO  Q(I+K,J), ALL I,J,  IF WE CONSIDER  Q
!        AS A TWO-DIM. ARRAY , WITH  2*K-1  ROWS (SEE COMMENTS IN
!        DBNFAC). IN THE PRESENT PROGRAM, WE TREAT  Q  AS AN EQUIVALENT
!        ONE-DIMENSIONAL ARRAY (BECAUSE OF FORTRAN RESTRICTIONS ON
!        DIMENSION STATEMENTS) . WE THEREFORE WANT  BCOEF(J) TO GO INTO
!        ENTRY
!            I -(LEFT+J) + 2*K + ((LEFT+J) - K-1)*(2*K-1)
!                   =  I-LEFT+1 + (LEFT -K)*(2*K-1) + (2*K-2)*J
!        OF  Q .
        JJ = I - LEFT + 1 + (LEFT-K)*(K+KM1)
        DO 40 J=1,K
          JJ = JJ + KPKM2
          Q(JJ) = BCOEF(J)
   40   CONTINUE
   50 CONTINUE
!
!     ***OBTAIN FACTORIZATION OF  A  , STORED AGAIN IN  Q.
      CALL DBNFAC(Q, K+KM1, N, KM1, KM1, IFLAG)
      GO TO (60, 90), IFLAG
!     *** SOLVE  A*BCOEF = Y  BY BACKSUBSTITUTION
   60 DO 70 I=1,N
        BCOEF(I) = Y(I)
   70 CONTINUE
      CALL DBNSLV(Q, K+KM1, N, KM1, KM1, BCOEF)
      RETURN
!
!
   80 CONTINUE
      CALL XERROR( &
      ' DBINTK,  SOME ABSCISSA WAS NOT IN THE SUPPORT OF THE CORRESPONDING BASIS FUNCTION AND THE SYSTEM IS SINGULAR.' &
      ,109,2, 1)
      RETURN
   90 CONTINUE
      CALL XERROR(  &
    ' DBINTK,  THE SYSTEM OF SOLVER DETECTS A SINGULAR SYSTEM ALTHOUGH THE THEORETICAL CONDITIONS FOR A SOLUTION WERE SATISFIED.',&
      123,8,1)
      RETURN
  100 CONTINUE
      CALL XERROR( ' DBINTK,  K DOES NOT SATISFY K.GE.1', 35, 2, 1)
      RETURN
  105 CONTINUE
      CALL XERROR( ' DBINTK,  N DOES NOT SATISFY N.GE.K', 35, 2, 1)
      RETURN
  110 CONTINUE
      CALL XERROR( ' DBINTK,  X(I) DOES NOT SATISFY X(I).LT.X(I+1) FOR SOME I', &
      57, 2, 1)
      RETURN
      END
      SUBROUTINE DBKNOT(X,N,K,T)
!***BEGIN PROLOGUE  DBKNOT
!***REFER TO  DB2INK,DB3INK
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   000330  Modified array declarations.  (JEC)
!
!***END PROLOGUE  DBKNOT
!
!  --------------------------------------------------------------------
!  DBKNOT CHOOSES A KNOT SEQUENCE FOR INTERPOLATION OF ORDER K AT THE
!  DATA POINTS X(I), I=1,..,N.  THE N+K KNOTS ARE PLACED IN THE ARRAY
!  T.  K KNOTS ARE PLACED AT EACH ENDPOINT AND NOT-A-KNOT END
!  CONDITIONS ARE USED.  THE REMAINING KNOTS ARE PLACED AT DATA POINTS
!  IF N IS EVEN AND BETWEEN DATA POINTS IF N IS ODD.  THE RIGHTMOST
!  KNOT IS SHIFTED SLIGHTLY TO THE RIGHT TO INSURE PROPER INTERPOLATION
!  AT X(N) (SEE PAGE 350 OF THE REFERENCE).
!  DOUBLE PRECISION VERSION OF BKNOT.
!  --------------------------------------------------------------------
!
!  ------------
!  DECLARATIONS
!  ------------
!
!  PARAMETERS
!
      INTEGER N, K
      DOUBLE PRECISION  X(N), T(*)
!
!  LOCAL VARIABLES
!
      INTEGER   I, J, IPJ, NPJ, IP1
      DOUBLE PRECISION  RNOT
!
!
!  ----------------------------
!  PUT K KNOTS AT EACH ENDPOINT
!  ----------------------------
!
!     (SHIFT RIGHT ENPOINTS SLIGHTLY -- SEE PG 350 OF REFERENCE)
      RNOT = X(N) + 0.10D0*( X(N)-X(N-1) )
      DO 110 J=1,K
         T(J) = X(1)
         NPJ = N + J
         T(NPJ) = RNOT
  110 CONTINUE
!
!  --------------------------
!  DISTRIBUTE REMAINING KNOTS
!  --------------------------
!
      IF (MOD(K,2) .EQ. 1)  GO TO 150
!
!     CASE OF EVEN K --  KNOTS AT DATA POINTS
!
      I = (K/2) - K
      JSTRT = K+1
      DO 120 J=JSTRT,N
         IPJ = I + J
         T(J) = X(IPJ)
  120 CONTINUE
      GO TO 200
!
!     CASE OF ODD K --  KNOTS BETWEEN DATA POINTS
!
  150 CONTINUE
      I = (K-1)/2 - K
      IP1 = I + 1
      JSTRT = K + 1
      DO 160 J=JSTRT,N
         IPJ = I + J
         T(J) = 0.50D0*( X(IPJ) + X(IPJ+1) )
  160 CONTINUE
  200 CONTINUE
!
      RETURN
      END
      SUBROUTINE DBNFAC(W,NROWW,NROW,NBANDL,NBANDU,IFLAG)
!***BEGIN PROLOGUE  DBNFAC
!***REFER TO  DBINT4,DBINTK
!
!  DBNFAC is the BANFAC routine from
!        * A Practical Guide to Splines *  by C. de Boor
!
!  DBNFAC is a double precision routine
!
!  Returns in  W  the LU-factorization (without pivoting) of the banded
!  matrix  A  of order  NROW  with  (NBANDL + 1 + NBANDU) bands or diag-
!  onals in the work array  W .
!
! *****  I N P U T  ****** W is double precision
!  W.....Work array of size  (NROWW,NROW)  containing the interesting
!        part of a banded matrix  A , with the diagonals or bands of  A
!        stored in the rows of  W , while columns of  A  correspond to
!        columns of  W . This is the storage mode used in  LINPACK  and
!        results in efficient innermost loops.
!           Explicitly,  A  has  NBANDL  bands below the diagonal
!                            +     1     (main) diagonal
!                            +   NBANDU  bands above the diagonal
!        and thus, with    MIDDLE = NBANDU + 1,
!          A(I+J,J)  is in  W(I+MIDDLE,J)  for I=-NBANDU,...,NBANDL
!                                              J=1,...,NROW .
!        For example, the interesting entries of A (1,2)-banded matrix
!        of order  9  would appear in the first  1+1+2 = 4  rows of  W
!        as follows.
!                          13 24 35 46 57 68 79
!                       12 23 34 45 56 67 78 89
!                    11 22 33 44 55 66 77 88 99
!                    21 32 43 54 65 76 87 98
!
!        All other entries of  W  not identified in this way with an en-
!        try of  A  are never referenced .
!  NROWW.....Row dimension of the work array  W .
!        must be  .GE.  NBANDL + 1 + NBANDU  .
!  NBANDL.....Number of bands of  A  below the main diagonal
!  NBANDU.....Number of bands of  A  above the main diagonal .
!
! *****  O U T P U T  ****** W is double precision
!  IFLAG.....Integer indicating success( = 1) or failure ( = 2) .
!     If  IFLAG = 1, then
!  W.....contains the LU-factorization of  A  into a unit lower triangu-
!        lar matrix  L  and an upper triangular matrix  U (both banded)
!        and stored in customary fashion over the corresponding entries
!        of  A . This makes it possible to solve any particular linear
!        system  A*X = B  for  X  by a
!              CALL DBNSLV ( W, NROWW, NROW, NBANDL, NBANDU, B )
!        with the solution X  contained in  B  on return .
!     If  IFLAG = 2, then
!        one of  NROW-1, NBANDL,NBANDU failed to be nonnegative, or else
!        one of the potential pivots was found to be zero indicating
!        that  A  does not have an LU-factorization. This implies that
!        A  is singular in case it is totally positive .
!
! *****  M E T H O D  ******
!     Gauss elimination  W I T H O U T  pivoting is used. The routine is
!  intended for use with matrices  A  which do not require row inter-
!  changes during factorization, especially for the  T O T A L L Y
!  P O S I T I V E  matrices which occur in spline calculations.
!     The routine should NOT be used for an arbitrary banded matrix.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DBNFAC
!
      INTEGER IFLAG, NBANDL, NBANDU, NROW, NROWW, I, IPK, J, JMAX, K, &
       KMAX, MIDDLE, MIDMK, NROWM1
      DOUBLE PRECISION W(NROWW,NROW), FACTOR, PIVOT
!
!***FIRST EXECUTABLE STATEMENT  DBNFAC
      IFLAG = 1
      MIDDLE = NBANDU + 1
!                         W(MIDDLE,.) CONTAINS THE MAIN DIAGONAL OF  A .
      NROWM1 = NROW - 1
      IF (NROWM1) 120, 110, 10
   10 IF (NBANDL.GT.0) GO TO 30
!                A IS UPPER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO .
      DO 20 I=1,NROWM1
        IF (W(MIDDLE,I).EQ.0.0D0) GO TO 120
   20 CONTINUE
      GO TO 110
   30 IF (NBANDU.GT.0) GO TO 60
!              A IS LOWER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO AND
!                 DIVIDE EACH COLUMN BY ITS DIAGONAL .
      DO 50 I=1,NROWM1
        PIVOT = W(MIDDLE,I)
        IF (PIVOT.EQ.0.0D0) GO TO 120
        JMAX = MIN0(NBANDL,NROW-I)
        DO 40 J=1,JMAX
          W(MIDDLE+J,I) = W(MIDDLE+J,I)/PIVOT
   40   CONTINUE
   50 CONTINUE
      RETURN
!
!        A  IS NOT JUST A TRIANGULAR MATRIX. CONSTRUCT LU FACTORIZATION
   60 DO 100 I=1,NROWM1
!                                  W(MIDDLE,I)  IS PIVOT FOR I-TH STEP .
        PIVOT = W(MIDDLE,I)
        IF (PIVOT.EQ.0.0D0) GO TO 120
!                 JMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN COLUMN  I
!                     BELOW THE DIAGONAL .
        JMAX = MIN0(NBANDL,NROW-I)
!              DIVIDE EACH ENTRY IN COLUMN  I  BELOW DIAGONAL BY PIVOT .
        DO 70 J=1,JMAX
          W(MIDDLE+J,I) = W(MIDDLE+J,I)/PIVOT
   70   CONTINUE
!                 KMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN ROW  I  TO
!                     THE RIGHT OF THE DIAGONAL .
        KMAX = MIN0(NBANDU,NROW-I)
!                  SUBTRACT  A(I,I+K)*(I-TH COLUMN) FROM (I+K)-TH COLUMN
!                  (BELOW ROW  I ) .
        DO 90 K=1,KMAX
          IPK = I + K
          MIDMK = MIDDLE - K
          FACTOR = W(MIDMK,IPK)
          DO 80 J=1,JMAX
            W(MIDMK+J,IPK) = W(MIDMK+J,IPK) - W(MIDDLE+J,I)*FACTOR
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE
!                                       CHECK THE LAST DIAGONAL ENTRY .
  110 IF (W(MIDDLE,NROW).NE.0.0D0) RETURN
  120 IFLAG = 2
      RETURN
      END
      SUBROUTINE DBNSLV(W,NROWW,NROW,NBANDL,NBANDU,B)
!***BEGIN PROLOGUE  DBNSLV
!***REFER TO  DBINT4,DBINTK
!
!  DBNSLV is the BANSLV routine from
!        * A Practical Guide to Splines *  by C. de Boor
!
!  DBNSLV is a double precision routine
!
!  Companion routine to  DBNFAC . It returns the solution  X  of the
!  linear system  A*X = B  in place of  B , given the LU-factorization
!  for  A  in the work array  W from DBNFAC.
!
! *****  I N P U T  ****** W,B are DOUBLE PRECISION
!  W, NROWW,NROW,NBANDL,NBANDU.....Describe the LU-factorization of a
!        banded matrix  A  of order  NROW  as constructed in  DBNFAC .
!        For details, see  DBNFAC .
!  B.....Right side of the system to be solved .
!
! *****  O U T P U T  ****** B is DOUBLE PRECISION
!  B.....Contains the solution  X , of order  NROW .
!
! *****  M E T H O D  ******
!     (With  A = L*U, as stored in  W,) the unit lower triangular system
!  L(U*X) = B  is solved for  Y = U*X, and  Y  stored in  B . Then the
!  upper triangular system  U*X = Y  is solved for  X  . The calcul-
!  ations are so arranged that the innermost loops stay within columns.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DBNSLV
!
      INTEGER NBANDL, NBANDU, NROW, NROWW, I, J, JMAX, MIDDLE, NROWM1
      DOUBLE PRECISION W(NROWW,NROW), B(NROW)
!***FIRST EXECUTABLE STATEMENT  DBNSLV
      MIDDLE = NBANDU + 1
      IF (NROW.EQ.1) GO TO 80
      NROWM1 = NROW - 1
      IF (NBANDL.EQ.0) GO TO 30
!                                 FORWARD PASS
!            FOR I=1,2,...,NROW-1, SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN
!            OF  L )  FROM RIGHT SIDE  (BELOW I-TH ROW) .
      DO 20 I=1,NROWM1
        JMAX = MIN0(NBANDL,NROW-I)
        DO 10 J=1,JMAX
          B(I+J) = B(I+J) - B(I)*W(MIDDLE+J,I)
   10   CONTINUE
   20 CONTINUE
!                                 BACKWARD PASS
!            FOR I=NROW,NROW-1,...,1, DIVIDE RIGHT SIDE(I) BY I-TH DIAG-
!            ONAL ENTRY OF  U, THEN SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN
!            OF  U)  FROM RIGHT SIDE  (ABOVE I-TH ROW).
   30 IF (NBANDU.GT.0) GO TO 50
!                                A  IS LOWER TRIANGULAR .
      DO 40 I=1,NROW
        B(I) = B(I)/W(1,I)
   40 CONTINUE
      RETURN
   50 I = NROW
   60 B(I) = B(I)/W(MIDDLE,I)
      JMAX = MIN0(NBANDU,I-1)
      DO 70 J=1,JMAX
        B(I-J) = B(I-J) - B(I)*W(MIDDLE-J,I)
   70 CONTINUE
      I = I - 1
      IF (I.GT.1) GO TO 60
   80 B(1) = B(1)/W(MIDDLE,1)
      RETURN
      END
      SUBROUTINE DBSPVN(T,JHIGH,K,INDEX,X,ILEFT,VNIKX,WORK,IWORK)
!***BEGIN PROLOGUE  DBSPVN
!***DATE WRITTEN   800901   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***REVISION HISTORY  (YYMMDD)
!   000330  Modified array declarations.  (JEC)
!
!***CATEGORY NO.  E3,K6
!***KEYWORDS  B-SPLINE,DATA FITTING,DOUBLE PRECISION,INTERPOLATION,
!             SPLINE
!***AUTHOR  AMOS, D. E., (SNLA)
!***PURPOSE  Calculates the value of all (possibly) nonzero basis
!            functions at X.
!***DESCRIPTION
!
!     Written by Carl de Boor and modified by D. E. Amos
!
!     Reference
!         SIAM J. Numerical Analysis, 14, No. 3, June, 1977, pp.441-472.
!
!     Abstract    **** a double precision routine ****
!         DBSPVN is the BSPLVN routine of the reference.
!
!         DBSPVN calculates the value of all (possibly) nonzero basis
!         functions at X of order MAX(JHIGH,(J+1)*(INDEX-1)), where T(K)
!         .LE. X .LE. T(N+1) and J=IWORK is set inside the routine on
!         the first call when INDEX=1.  ILEFT is such that T(ILEFT) .LE.
!         X .LT. T(ILEFT+1).  A call to DINTRV(T,N+1,X,ILO,ILEFT,MFLAG)
!         produces the proper ILEFT.  DBSPVN calculates using the basic
!         algorithm needed in DBSPVD.  If only basis functions are
!         desired, setting JHIGH=K and INDEX=1 can be faster than
!         calling DBSPVD, but extra coding is required for derivatives
!         (INDEX=2) and DBSPVD is set up for this purpose.
!
!         Left limiting values are set up as described in DBSPVD.
!
!     Description of Arguments
!
!         Input      T,X are double precision
!          T       - knot vector of length N+K, where
!                    N = number of B-spline basis functions
!                    N = sum of knot multiplicities-K
!          JHIGH   - order of B-spline, 1 .LE. JHIGH .LE. K
!          K       - highest possible order
!          INDEX   - INDEX = 1 gives basis functions of order JHIGH
!                          = 2 denotes previous entry with work, IWORK
!                              values saved for subsequent calls to
!                              DBSPVN.
!          X       - argument of basis functions,
!                    T(K) .LE. X .LE. T(N+1)
!          ILEFT   - largest integer such that
!                    T(ILEFT) .LE. X .LT.  T(ILEFT+1)
!
!         Output     VNIKX, WORK are double precision
!          VNIKX   - vector of length K for spline values.
!          WORK    - a work vector of length 2*K
!          IWORK   - a work parameter.  Both WORK and IWORK contain
!                    information necessary to continue for INDEX = 2.
!                    When INDEX = 1 exclusively, these are scratch
!                    variables and can be used for other purposes.
!
!     Error Conditions
!         Improper input is a fatal error.
!***REFERENCES  C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,
!                 SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3,
!                 JUNE 1977, PP. 441-472.
!***ROUTINES CALLED  XERROR
!***END PROLOGUE  DBSPVN
!
!
      INTEGER ILEFT, IMJP1, INDEX, IPJ, IWORK, JHIGH, JP1, JP1ML, K, L
      DOUBLE PRECISION T, VM, VMPREV, VNIKX, WORK, X
!     DIMENSION T(ILEFT+JHIGH)
      DIMENSION T(*), VNIKX(K), WORK(*)
!     CONTENT OF J, DELTAM, DELTAP IS EXPECTED UNCHANGED BETWEEN CALLS.
!     WORK(I) = DELTAP(I), WORK(K+I) = DELTAM(I), I = 1,K
!***FIRST EXECUTABLE STATEMENT  DBSPVN
      IF(K.LT.1) GO TO 90
      IF(JHIGH.GT.K .OR. JHIGH.LT.1) GO TO 100
      IF(INDEX.LT.1 .OR. INDEX.GT.2) GO TO 105
      IF(X.LT.T(ILEFT) .OR. X.GT.T(ILEFT+1)) GO TO 110
      GO TO (10, 20), INDEX
   10 IWORK = 1
      VNIKX(1) = 1.0D0
      IF (IWORK.GE.JHIGH) GO TO 40
!
   20 IPJ = ILEFT + IWORK
      WORK(IWORK) = T(IPJ) - X
      IMJP1 = ILEFT - IWORK + 1
      WORK(K+IWORK) = X - T(IMJP1)
      VMPREV = 0.0D0
      JP1 = IWORK + 1
      DO 30 L=1,IWORK
        JP1ML = JP1 - L
        VM = VNIKX(L)/(WORK(L)+WORK(K+JP1ML))
        VNIKX(L) = VM*WORK(L) + VMPREV
        VMPREV = VM*WORK(K+JP1ML)
   30 CONTINUE
      VNIKX(JP1) = VMPREV
      IWORK = JP1
      IF (IWORK.LT.JHIGH) GO TO 20
!
   40 RETURN
!
!
   90 CONTINUE
      CALL XERROR( ' DBSPVN,  K DOES NOT SATISFY K.GE.1', 35, 2, 1)
      RETURN
  100 CONTINUE
      CALL XERROR( ' DBSPVN,  JHIGH DOES NOT SATISFY 1.LE.JHIGH.LE.K', &
      48, 2, 1)
      RETURN
  105 CONTINUE
      CALL XERROR( ' DBSPVN,  INDEX IS NOT 1 OR 2',29,2,1)
      RETURN
  110 CONTINUE
      CALL XERROR( ' DBSPVN,  X DOES NOT SATISFY T(ILEFT).LE.X.LE.T(ILEFT+1)', &
      56, 2, 1)
      RETURN
      END
      SUBROUTINE DBTPCF(X,N,FCN,LDF,NF,T,K,BCOEF,WORK)
!***BEGIN PROLOGUE  DBTPCF
!***REFER TO  DB2INK,DB3INK
!***ROUTINES CALLED  DBINTK,DBNSLV
!***REVISION HISTORY  (YYMMDD)
!   000330  Modified array declarations.  (JEC)
!
!***END PROLOGUE  DBTPCF
!
!  -----------------------------------------------------------------
!  DBTPCF COMPUTES B-SPLINE INTERPOLATION COEFFICIENTS FOR NF SETS
!  OF DATA STORED IN THE COLUMNS OF THE ARRAY FCN. THE B-SPLINE
!  COEFFICIENTS ARE STORED IN THE ROWS OF BCOEF HOWEVER.
!  EACH INTERPOLATION IS BASED ON THE N ABCISSA STORED IN THE
!  ARRAY X, AND THE N+K KNOTS STORED IN THE ARRAY T. THE ORDER
!  OF EACH INTERPOLATION IS K. THE WORK ARRAY MUST BE OF LENGTH
!  AT LEAST 2*K*(N+1).
!  DOUBLE PRECISION VERSION OF BTPCF.
!  -----------------------------------------------------------------
!
!  ------------
!  DECLARATIONS
!  ------------
!
!  PARAMETERS
!
      INTEGER  N, LDF, K
      DOUBLE PRECISION  X(N), FCN(LDF,NF), T(*), BCOEF(NF,N), WORK(*)
!
!  LOCAL VARIABLES
!
      INTEGER  I, J, K1, K2, IQ, IW
!
!  ---------------------------------------------
!  CHECK FOR NULL INPUT AND PARTITION WORK ARRAY
!  ---------------------------------------------
!
!***FIRST EXECUTABLE STATEMENT
      IF (NF .LE. 0)  GO TO 500
      K1 = K - 1
      K2 = K1 + K
      IQ = 1 + N
      IW = IQ + K2*N+1
!
!  -----------------------------
!  COMPUTE B-SPLINE COEFFICIENTS
!  -----------------------------
!
!
!   FIRST DATA SET
!
      CALL DBINTK(X,FCN,T,N,K,WORK,WORK(IQ),WORK(IW))
      DO 20 I=1,N
         BCOEF(1,I) = WORK(I)
   20 CONTINUE
!
!  ALL REMAINING DATA SETS BY BACK-SUBSTITUTION
!
      IF (NF .EQ. 1)  GO TO 500
      DO 100 J=2,NF
         DO 50 I=1,N
            WORK(I) = FCN(I,J)
   50    CONTINUE
         CALL DBNSLV(WORK(IQ),K2,N,K1,K1,WORK)
         DO 60 I=1,N
            BCOEF(J,I) = WORK(I)
   60    CONTINUE
  100 CONTINUE
!
!  ----
!  EXIT
!  ----
!
  500 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DBVALU(T,A,N,K,IDERIV,X,INBV,WORK)
!***BEGIN PROLOGUE  DBVALU
!***DATE WRITTEN   800901   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***REVISION HISTORY  (YYMMDD)
!   000330  Modified array declarations.  (JEC)
!
!***CATEGORY NO.  E3,K6
!***KEYWORDS  B-SPLINE,DATA FITTING,DOUBLE PRECISION,INTERPOLATION,
!             SPLINE
!***AUTHOR  AMOS, D. E., (SNLA)
!***PURPOSE  Evaluates the B-representation of a B-spline at X for the
!            function value or any of its derivatives.
!***DESCRIPTION
!
!     Written by Carl de Boor and modified by D. E. Amos
!
!     Reference
!         SIAM J. Numerical Analysis, 14, No. 3, June, 1977, pp.441-472.
!
!     Abstract   **** a double precision routine ****
!         DBVALU is the BVALUE function of the reference.
!
!         DBVALU evaluates the B-representation (T,A,N,K) of a B-spline
!         at X for the function value on IDERIV=0 or any of its
!         derivatives on IDERIV=1,2,...,K-1.  Right limiting values
!         (right derivatives) are returned except at the right end
!         point X=T(N+1) where left limiting values are computed.  The
!         spline is defined on T(K) .LE. X .LE. T(N+1).  DBVALU returns
!         a fatal error message when X is outside of this interval.
!
!         To compute left derivatives or left limiting values at a
!         knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1.
!
!         DBVALU calls DINTRV
!
!     Description of Arguments
!
!         Input      T,A,X are double precision
!          T       - knot vector of length N+K
!          A       - B-spline coefficient vector of length N
!          N       - number of B-spline coefficients
!                    N = sum of knot multiplicities-K
!          K       - order of the B-spline, K .GE. 1
!          IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1
!                    IDERIV = 0 returns the B-spline value
!          X       - argument, T(K) .LE. X .LE. T(N+1)
!          INBV    - an initialization parameter which must be set
!                    to 1 the first time DBVALU is called.
!
!         Output     WORK,DBVALU are double precision
!          INBV    - INBV contains information for efficient process-
!                    ing after the initial call and INBV must not
!                    be changed by the user.  Distinct splines require
!                    distinct INBV parameters.
!          WORK    - work vector of length 3*K.
!          DBVALU  - value of the IDERIV-th derivative at X
!
!     Error Conditions
!         An improper input is a fatal error
!***REFERENCES  C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,
!                 SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3,
!                 JUNE 1977, PP. 441-472.
!***ROUTINES CALLED  DINTRV,XERROR
!***END PROLOGUE  DBVALU
!
!
      INTEGER I,IDERIV,IDERP1,IHI,IHMKMJ,ILO,IMK,IMKPJ, INBV, IPJ, &
      IP1, IP1MJ, J, JJ, J1, J2, K, KMIDER, KMJ, KM1, KPK, MFLAG, N
      DOUBLE PRECISION A, FKMJ, T, WORK, X
      DIMENSION T(*), A(N), WORK(*)
!***FIRST EXECUTABLE STATEMENT  DBVALU
      DBVALU = 0.0D0
      IF(K.LT.1) GO TO 102
      IF(N.LT.K) GO TO 101
      IF(IDERIV.LT.0 .OR. IDERIV.GE.K) GO TO 110
      KMIDER = K - IDERIV
!
! *** FIND *I* IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1)
!     (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)).
      KM1 = K - 1
      CALL DINTRV(T, N+1, X, INBV, I, MFLAG)
      IF (X.LT.T(K)) GO TO 120
      IF (MFLAG.EQ.0) GO TO 20
      IF (X.GT.T(I)) GO TO 130
   10 IF (I.EQ.K) GO TO 140
      I = I - 1
      IF (X.EQ.T(I)) GO TO 10
!
! *** DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES
!     WORK(I) = AJ(I), WORK(K+I) = DP(I), WORK(K+K+I) = DM(I), I=1.K
!
   20 IMK = I - K
      DO 30 J=1,K
        IMKPJ = IMK + J
        WORK(J) = A(IMKPJ)
   30 CONTINUE
      IF (IDERIV.EQ.0) GO TO 60
      DO 50 J=1,IDERIV
        KMJ = K - J
        FKMJ = DBLE(FLOAT(KMJ))
        DO 40 JJ=1,KMJ
          IHI = I + JJ
          IHMKMJ = IHI - KMJ
          WORK(JJ) = (WORK(JJ+1)-WORK(JJ))/(T(IHI)-T(IHMKMJ))*FKMJ
   40   CONTINUE
   50 CONTINUE
!
! *** COMPUTE VALUE AT *X* IN (T(I),(T(I+1)) OF IDERIV-TH DERIVATIVE,
!     GIVEN ITS RELEVANT B-SPLINE COEFF. IN AJ(1),...,AJ(K-IDERIV).
   60 IF (IDERIV.EQ.KM1) GO TO 100
      IP1 = I + 1
      KPK = K + K
      J1 = K + 1
      J2 = KPK + 1
      DO 70 J=1,KMIDER
        IPJ = I + J
        WORK(J1) = T(IPJ) - X
        IP1MJ = IP1 - J
        WORK(J2) = X - T(IP1MJ)
        J1 = J1 + 1
        J2 = J2 + 1
   70 CONTINUE
      IDERP1 = IDERIV + 1
      DO 90 J=IDERP1,KM1
        KMJ = K - J
        ILO = KMJ
        DO 80 JJ=1,KMJ
          WORK(JJ) = (WORK(JJ+1)*WORK(KPK+ILO)+WORK(JJ)  &
                   *WORK(K+JJ))/(WORK(KPK+ILO)+WORK(K+JJ))
          ILO = ILO - 1
   80   CONTINUE
   90 CONTINUE
  100 DBVALU = WORK(1)
      RETURN
!
!
  101 CONTINUE
      CALL XERROR( ' DBVALU,  N DOES NOT SATISFY N.GE.K',35,2,1)
      RETURN
  102 CONTINUE
      CALL XERROR( ' DBVALU,  K DOES NOT SATISFY K.GE.1',35,2,1)
      RETURN
  110 CONTINUE
      CALL XERROR( ' DBVALU,  IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K', &
      50, 2, 1)
      RETURN
  120 CONTINUE
      CALL XERROR( ' DBVALU,  X IS N0T GREATER THAN OR EQUAL TO T(K)', &
      48, 2, 1)
      RETURN
  130 CONTINUE
      CALL XERROR( ' DBVALU,  X IS NOT LESS THAN OR EQUAL TO T(N+1)', &
      47, 2, 1)
      RETURN
  140 CONTINUE
      CALL XERROR( ' DBVALU,  A LEFT LIMITING VALUE CANN0T BE OBTAINED AT T(K)',&
          58, 2, 1)
      RETURN
      END
      SUBROUTINE DCHFEV(X1,X2,F1,F2,D1,D2,NE,XE,FE,NEXT,IERR)
!***BEGIN PROLOGUE  DCHFEV
!***DATE WRITTEN   811019   (YYMMDD)
!***REVISION DATE  870707   (YYMMDD)
!***CATEGORY NO.  E3,H1
!***KEYWORDS  LIBRARY=SLATEC(PCHIP),
!            TYPE=DOUBLE PRECISION(CHFEV-S DCHFEV-D),
!            CUBIC HERMITE EVALUATION,CUBIC POLYNOMIAL EVALUATION
!***AUTHOR  FRITSCH, F. N., (LLNL)
!            MATHEMATICS AND STATISTICS DIVISION
!            LAWRENCE LIVERMORE NATIONAL LABORATORY
!            P.O. BOX 808  (L-316)
!            LIVERMORE, CA  94550
!            FTS 532-4275, (415) 422-4275
!***PURPOSE  Evaluate a cubic polynomial given in Hermite form at an
!           array of points.  While designed for use by DPCHFE, it may
!           be useful directly as an evaluator for a piecewise cubic
!           Hermite function in applications, such as graphing, where
!           the interval is known in advance.
!***DESCRIPTION
!
!      **** Double Precision version of CHFEV ****
!
!         DCHFEV:  Cubic Hermite Function EValuator
!
!    Evaluates the cubic polynomial determined by function values
!    F1,F2 and derivatives D1,D2 on interval (X1,X2) at the points
!    XE(J), J=1(1)NE.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!       INTEGER  NE, NEXT(2), IERR
!       DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE)
!
!       CALL  DCHFEV (X1,X2, F1,F2, D1,D2, NE, XE, FE, NEXT, IERR)
!
!   Parameters:
!
!    X1,X2 -- (input) endpoints of interval of definition of cubic.
!          (Error return if  X1.EQ.X2 .)
!
!    F1,F2 -- (input) values of function at X1 and X2, respectively.
!
!    D1,D2 -- (input) values of derivative at X1 and X2, respectively.
!
!    NE -- (input) number of evaluation points.  (Error return if
!          NE.LT.1 .)
!
!    XE -- (input) real*8 array of points at which the function is to
!          be evaluated.  If any of the XE are outside the interval
!          [X1,X2], a warning error is returned in NEXT.
!
!    FE -- (output) real*8 array of values of the cubic function
!          defined by  X1,X2, F1,F2, D1,D2  at the points  XE.
!
!    NEXT -- (output) integer array indicating number of extrapolation
!          points:
!           NEXT(1) = number of evaluation points to left of interval.
!           NEXT(2) = number of evaluation points to right of interval.
!
!    IERR -- (output) error flag.
!          Normal return:
!             IERR = 0  (no errors).
!          "Recoverable" errors:
!             IERR = -1  if NE.LT.1 .
!             IERR = -2  if X1.EQ.X2 .
!               (The FE-array has not been changed in either case.)
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERROR
!***END PROLOGUE  DCHFEV
!
! ----------------------------------------------------------------------
!
!  Change record:
!    82-08-03   Minor cosmetic changes for release 1.
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!    To produce a single precision version, simply:
!       a. Change DCHFEV to CHFEV wherever it occurs,
!       b. Change the double precision declaration to real,
!       c. Change the constant ZERO to single precision, and
!       d. Change the names of the Fortran functions:  AMAX1, AMIN1.
!
!  DECLARE ARGUMENTS.
!
      INTEGER NE, NEXT(2), IERR
      DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE)
!
!  DECLARE LOCAL VARIABLES.
!
      INTEGER I
      DOUBLE PRECISION  C2, C3, DEL1, DEL2, DELTA, H, X, XMI, XMA, &
       ZERO
!      DATA  ZERO /0.D0/
       parameter (ZERO=0.D0)
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  DCHFEV
      IF (NE .LT. 1)  GO TO 5001
      H = X2 - X1
      IF (H .EQ. ZERO)  GO TO 5002
!
!  INITIALIZE.
!
      IERR = 0
      NEXT(1) = 0
      NEXT(2) = 0
      XMI = DMIN1(ZERO, H)
      XMA = DMAX1(ZERO, H)
!
!  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
!
      DELTA = (F2 - F1)/H
      DEL1 = (D1 - DELTA)/H
      DEL2 = (D2 - DELTA)/H
!                                          (DELTA IS NO LONGER NEEDED.)
      C2 = -(DEL1+DEL1 + DEL2)
      C3 = (DEL1 + DEL2)/H
!                              (H, DEL1 AND DEL2 ARE NO LONGER NEEDED.)
!
!  EVALUATION LOOP.
!
      DO 500  I = 1, NE
         X = XE(I) - X1
         FE(I) = F1 + X*(D1 + X*(C2 + X*C3))
!         COUNT EXTRAPOLATION POINTS.
         IF ( X.LT.XMI )  NEXT(1) = NEXT(1) + 1
         IF ( X.GT.XMA )  NEXT(2) = NEXT(2) + 1
!       (NOTE REDUNDANCY--IF EITHER CONDITION IS TRUE, OTHER IS FALSE.)
  500 CONTINUE
!
!  NORMAL RETURN.
!
      RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!    NE.LT.1 RETURN.
      IERR = -1
      CALL XERROR ('DCHFEV -- NUMBER OF EVALUATION POINTS LESS THAN ONE'  &
                , 51, IERR, 1)
      RETURN
!
 5002 CONTINUE
!    X1.EQ.X2 RETURN.
      IERR = -2
      CALL XERROR ('DCHFEV -- INTERVAL ENDPOINTS EQUAL' &
                , 34, IERR, 1)
      RETURN
!------------- LAST LINE OF DCHFEV FOLLOWS -----------------------------
      END
      SUBROUTINE DINTRV(XT,LXT,X,ILO,ILEFT,MFLAG)
!***BEGIN PROLOGUE  DINTRV
!***DATE WRITTEN   800901   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  E3,K6
!***KEYWORDS  B-SPLINE,DATA FITTING,DOUBLE PRECISION,INTERPOLATION,
!             SPLINE
!***AUTHOR  AMOS, D. E., (SNLA)
!***PURPOSE  Computes the largest integer ILEFT in 1.LE.ILEFT.LE.LXT
!            such that XT(ILEFT).LE.X where XT(*) is a subdivision of
!            the X interval.
!***DESCRIPTION
!
!     Written by Carl de Boor and modified by D. E. Amos
!
!     Reference
!         SIAM J.  Numerical Analysis, 14, No. 3, June 1977, pp.441-472.
!
!     Abstract    **** a double precision routine ****
!         DINTRV is the INTERV routine of the reference.
!
!         DINTRV computes the largest integer ILEFT in 1 .LE. ILEFT .LE.
!         LXT such that XT(ILEFT) .LE. X where XT(*) is a subdivision of
!         the X interval.  Precisely,
!
!                      X .LT. XT(1)                1         -1
!         if  XT(I) .LE. X .LT. XT(I+1)  then  ILEFT=I  , MFLAG=0
!           XT(LXT) .LE. X                         LXT        1,
!
!         That is, when multiplicities are present in the break point
!         to the left of X, the largest index is taken for ILEFT.
!
!     Description of Arguments
!
!         Input      XT,X are double precision
!          XT      - XT is a knot or break point vector of length LXT
!          LXT     - length of the XT vector
!          X       - argument
!          ILO     - an initialization parameter which must be set
!                    to 1 the first time the spline array XT is
!                    processed by DINTRV.
!
!         Output
!          ILO     - ILO contains information for efficient process-
!                    ing after the initial call and ILO must not be
!                    changed by the user.  Distinct splines require
!                    distinct ILO parameters.
!          ILEFT   - largest integer satisfying XT(ILEFT) .LE. X
!          MFLAG   - signals when X lies out of bounds
!
!     Error Conditions
!         None
!***REFERENCES  C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,
!                 SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3,
!                 JUNE 1977, PP. 441-472.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  DINTRV
!
!
      INTEGER IHI, ILEFT, ILO, ISTEP, LXT, MFLAG, MIDDLE
      DOUBLE PRECISION X, XT
      DIMENSION XT(LXT)
!***FIRST EXECUTABLE STATEMENT  DINTRV
      IHI = ILO + 1
      IF (IHI.LT.LXT) GO TO 10
      IF (X.GE.XT(LXT)) GO TO 110
      IF (LXT.LE.1) GO TO 90
      ILO = LXT - 1
      IHI = LXT
!
   10 IF (X.GE.XT(IHI)) GO TO 40
      IF (X.GE.XT(ILO)) GO TO 100
!
! *** NOW X .LT. XT(IHI) . FIND LOWER BOUND
      ISTEP = 1
   20 IHI = ILO
      ILO = IHI - ISTEP
      IF (ILO.LE.1) GO TO 30
      IF (X.GE.XT(ILO)) GO TO 70
      ISTEP = ISTEP*2
      GO TO 20
   30 ILO = 1
      IF (X.LT.XT(1)) GO TO 90
      GO TO 70
! *** NOW X .GE. XT(ILO) . FIND UPPER BOUND
   40 ISTEP = 1
   50 ILO = IHI
      IHI = ILO + ISTEP
      IF (IHI.GE.LXT) GO TO 60
      IF (X.LT.XT(IHI)) GO TO 70
      ISTEP = ISTEP*2
      GO TO 50
   60 IF (X.GE.XT(LXT)) GO TO 110
      IHI = LXT
!
! *** NOW XT(ILO) .LE. X .LT. XT(IHI) . NARROW THE INTERVAL
   70 MIDDLE = (ILO+IHI)/2
      IF (MIDDLE.EQ.ILO) GO TO 100
!     NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1
      IF (X.LT.XT(MIDDLE)) GO TO 80
      ILO = MIDDLE
      GO TO 70
   80 IHI = MIDDLE
      GO TO 70
! *** SET OUTPUT AND RETURN
   90 MFLAG = -1
      ILEFT = 1
      RETURN
  100 MFLAG = 0
      ILEFT = ILO
      RETURN
  110 MFLAG = 1
      ILEFT = LXT
      RETURN
      END
      SUBROUTINE DPCHCE(IC,VC,N,X,H,SLOPE,D,INCFD,IERR)
!***BEGIN PROLOGUE  DPCHCE
!***REFER TO  DPCHIC
!***ROUTINES CALLED  DPCHDF,DPCHST,XERROR
!***REVISION DATE  870707   (YYMMDD)
!***DESCRIPTION
!
!         DPCHCE:  DPCHIC End Derivative Setter.
!
!    Called by DPCHIC to set end derivatives as requested by the user.
!    It must be called after interior derivative values have been set.
!                     -----
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the D-array.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!       PARAMETER  (INCFD = ...)
!       INTEGER  IC(2), N, IERR
!       DOUBLE PRECISION  VC(2), X(N), H(N), SLOPE(N), D(INCFD,N)
!
!       CALL  DPCHCE (IC, VC, N, X, H, SLOPE, D, INCFD, IERR)
!
!   Parameters:
!
!    IC -- (input) integer array of length 2 specifying desired
!          boundary conditions:
!          IC(1) = IBEG, desired condition at beginning of data.
!          IC(2) = IEND, desired condition at end of data.
!          ( see prologue to DPCHIC for details. )
!
!    VC -- (input) real*8 array of length 2 specifying desired boundary
!          values.  VC(1) need be set only if IC(1) = 2 or 3 .
!                   VC(2) need be set only if IC(2) = 2 or 3 .
!
!    N -- (input) number of data points.  (assumes N.GE.2)
!
!    X -- (input) real*8 array of independent variable values.  (the
!          elements of X are assumed to be strictly increasing.)
!
!    H -- (input) real*8 array of interval lengths.
!    SLOPE -- (input) real*8 array of data slopes.
!          If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!                 H(I) =  X(I+1)-X(I),
!             SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!
!    D -- (input) real*8 array of derivative values at the data points.
!          The value corresponding to X(I) must be stored in
!               D(1+(I-1)*INCFD),  I=1(1)N.
!         (output) the value of D at X(1) and/or X(N) is changed, if
!          necessary, to produce the requested boundary conditions.
!          no other entries in D are changed.
!
!    INCFD -- (input) increment between successive values in D.
!          This argument is provided primarily for 2-D applications.
!
!    IERR -- (output) error flag.
!          Normal return:
!             IERR = 0  (no errors).
!          Warning errors:
!             IERR = 1  if IBEG.LT.0 and D(1) had to be adjusted for
!                       monotonicity.
!             IERR = 2  if IEND.LT.0 and D(1+(N-1)*INCFD) had to be
!                       adjusted for monotonicity.
!             IERR = 3  if both of the above are true.
!
!    -------
!    WARNING:  This routine does no validity-checking of arguments.
!    -------
!
!  Fortran intrinsics used:  DABS, IABS.
!
!***END PROLOGUE  DPCHCE
!
! ----------------------------------------------------------------------
!
!  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
!                 Mathematics and Statistics Division,
!                 Lawrence Livermore National Laboratory.
!
!  Change record:
!    82-08-05   Converted to SLATEC library version.
!    87-07-07   Corrected XERROR calls for d.p. name(s).
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!    1. The function DPCHST(ARG1,ARG2)  is assumed to return zero if
!       either argument is zero, +1 if they are of the same sign, and
!       -1 if they are of opposite sign.
!    2. To produce a single precision version, simply:
!       a. Change DPCHCE to PCHCE wherever it occurs,
!       b. Change DPCHDF to PCHDF wherever it occurs,
!       c. Change DPCHST to PCHST wherever it occurs,
!       d. Change all references to the Fortran intrinsics to their
!          real equivalents,
!       e. Change the double precision declarations to real, and
!       f. Change the constants ZERO, HALF, ... to single precision.
!    3. One could reduce the number of arguments and amount of local
!       storage, at the expense of reduced code clarity, by passing in
!       the array WK (rather than splitting it into H and SLOPE) and
!       increasing its length enough to incorporate STEMP and XTEMP.
!    4. The two monotonicity checks only use the sufficient conditions.
!       Thus, it is possible (but unlikely) for a boundary condition to
!       be changed, even though the original interpolant was monotonic.
!       (At least the result is a continuous function of the data.)
!
!  DECLARE ARGUMENTS.
!
      INTEGER IC(2), N, INCFD, IERR
      DOUBLE PRECISION  VC(2), X(N), H(N), SLOPE(N), D(INCFD,N)
!
!  DELCARE LOCAL VARIABLES.
!
      INTEGER IBEG, IEND, IERF, INDEX, J, K
      DOUBLE PRECISION  HALF, STEMP(3), THREE, TWO, XTEMP(4), ZERO
      DOUBLE PRECISION  DPCHDF, DPCHST
!
!  INITIALIZE.
!
!      DATA  ZERO /0.D0/,  HALF/.5D0/,  TWO/2.D0/, THREE/3.D0/
      parameter (ZERO=0.D0, HALF=0.5D0, TWO=2.D0, THREE=3.D0)
!
!***FIRST EXECUTABLE STATEMENT  DPCHCE
      IBEG = IC(1)
      IEND = IC(2)
      IERR = 0
!
!  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
!
      IF ( IABS(IBEG).GT.N )  IBEG = 0
      IF ( IABS(IEND).GT.N )  IEND = 0
!
!  TREAT BEGINNING BOUNDARY CONDITION.
!
      IF (IBEG .EQ. 0)  GO TO 2000
      K = IABS(IBEG)
      IF (K .EQ. 1)  THEN
!       BOUNDARY VALUE PROVIDED.
         D(1,1) = VC(1)
      ELSE IF (K .EQ. 2)  THEN
!       BOUNDARY SECOND DERIVATIVE PROVIDED.
         D(1,1) = HALF*( (THREE*SLOPE(1) - D(1,2)) - HALF*VC(1)*H(1) )
      ELSE IF (K .LT. 5)  THEN
!       USE K-POINT DERIVATIVE FORMULA.
!       PICK UP FIRST K POINTS, IN REVERSE ORDER.
         DO 10  J = 1, K
            INDEX = K-J+1
!          INDEX RUNS FROM K DOWN TO 1.
            XTEMP(J) = X(INDEX)
            IF (J .LT. K)  STEMP(J) = SLOPE(INDEX-1)
   10    CONTINUE
!                -----------------------------
         D(1,1) = DPCHDF (K, XTEMP, STEMP, IERF)
!                -----------------------------
         IF (IERF .NE. 0)  GO TO 5001
      ELSE
!       USE 'NOT A KNOT' CONDITION.
         D(1,1) = ( THREE*(H(1)*SLOPE(2) + H(2)*SLOPE(1)) &
                  - TWO*(H(1)+H(2))*D(1,2) - H(1)*D(1,3) ) / H(2)
      ENDIF
!
      IF (IBEG .GT. 0)  GO TO 2000
!
!  CHECK D(1,1) FOR COMPATIBILITY WITH MONOTONICITY.
!
      IF (SLOPE(1) .EQ. ZERO)  THEN
         IF (D(1,1) .NE. ZERO)  THEN
            D(1,1) = ZERO
            IERR = IERR + 1
         ENDIF
      ELSE IF ( DPCHST(D(1,1),SLOPE(1)) .LT. ZERO)  THEN
         D(1,1) = ZERO
         IERR = IERR + 1
      ELSE IF ( DABS(D(1,1)) .GT. THREE*DABS(SLOPE(1)) )  THEN
         D(1,1) = THREE*SLOPE(1)
         IERR = IERR + 1
      ENDIF
!
!  TREAT END BOUNDARY CONDITION.
!
 2000 CONTINUE
      IF (IEND .EQ. 0)  GO TO 5000
      K = IABS(IEND)
      IF (K .EQ. 1)  THEN
!       BOUNDARY VALUE PROVIDED.
         D(1,N) = VC(2)
      ELSE IF (K .EQ. 2)  THEN
!       BOUNDARY SECOND DERIVATIVE PROVIDED.
         D(1,N) = HALF*( (THREE*SLOPE(N-1) - D(1,N-1)) +    &
                                                HALF*VC(2)*H(N-1) )
      ELSE IF (K .LT. 5)  THEN
!       USE K-POINT DERIVATIVE FORMULA.
!       PICK UP LAST K POINTS.
         DO 2010  J = 1, K
            INDEX = N-K+J
!          INDEX RUNS FROM N+1-K UP TO N.
            XTEMP(J) = X(INDEX)
            IF (J .LT. K)  STEMP(J) = SLOPE(INDEX)
 2010    CONTINUE
!                -----------------------------
         D(1,N) = DPCHDF (K, XTEMP, STEMP, IERF)
!                -----------------------------
         IF (IERF .NE. 0)  GO TO 5001
      ELSE
!       USE 'NOT A KNOT' CONDITION.
         D(1,N) = ( THREE*(H(N-1)*SLOPE(N-2) + H(N-2)*SLOPE(N-1))   &
                  - TWO*(H(N-1)+H(N-2))*D(1,N-1) - H(N-1)*D(1,N-2) )  &
                                                              / H(N-2)
      ENDIF
!
      IF (IEND .GT. 0)  GO TO 5000
!
!  CHECK D(1,N) FOR COMPATIBILITY WITH MONOTONICITY.
!
      IF (SLOPE(N-1) .EQ. ZERO)  THEN
         IF (D(1,N) .NE. ZERO)  THEN
            D(1,N) = ZERO
            IERR = IERR + 2
         ENDIF
      ELSE IF ( DPCHST(D(1,N),SLOPE(N-1)) .LT. ZERO)  THEN
         D(1,N) = ZERO
         IERR = IERR + 2
      ELSE IF ( DABS(D(1,N)) .GT. THREE*DABS(SLOPE(N-1)) )  THEN
         D(1,N) = THREE*SLOPE(N-1)
         IERR = IERR + 2
      ENDIF
!
!  NORMAL RETURN.
!
 5000 CONTINUE
      RETURN
!
!  ERROR RETURN.
!
 5001 CONTINUE
!    ERROR RETURN FROM DPCHDF.
!   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -1
      CALL XERROR ('DPCHCE -- ERROR RETURN FROM DPCHDF' &
                , 0, IERR, 1)
      RETURN
!------------- LAST LINE OF DPCHCE FOLLOWS -----------------------------
      END
      SUBROUTINE DPCHCI(N,H,SLOPE,D,INCFD)
!***BEGIN PROLOGUE  DPCHCI
!***REFER TO  DPCHIC
!***ROUTINES CALLED  DPCHST
!***REVISION DATE  870707   (YYMMDD)
!***DESCRIPTION
!
!         DPCHCI:  DPCHIC Initial Derivative Setter.
!
!    Called by DPCHIC to set derivatives needed to determine a monotone
!    piecewise cubic Hermite interpolant to the data.
!
!    Default boundary conditions are provided which are compatible
!    with monotonicity.  If the data are only piecewise monotonic, the
!    interpolant will have an extremum at each point where monotonicity
!    switches direction.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the D-array.
!
!    The resulting piecewise cubic Hermite function should be identical
!    (within roundoff error) to that produced by DPCHIM.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!       PARAMETER  (INCFD = ...)
!       INTEGER  N
!       DOUBLE PRECISION  H(N), SLOPE(N), D(INCFD,N)
!
!       CALL  DPCHCI (N, H, SLOPE, D, INCFD)
!
!   Parameters:
!
!    N -- (input) number of data points.
!          If N=2, simply does linear interpolation.
!
!    H -- (input) real*8 array of interval lengths.
!    SLOPE -- (input) real*8 array of data slopes.
!          If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!                 H(I) =  X(I+1)-X(I),
!             SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!
!    D -- (output) real*8 array of derivative values at data points.
!          If the data are monotonic, these values will determine a
!          a monotone cubic Hermite function.
!          The value corresponding to X(I) is stored in
!               D(1+(I-1)*INCFD),  I=1(1)N.
!          No other entries in D are changed.
!
!    INCFD -- (input) increment between successive values in D.
!          This argument is provided primarily for 2-D applications.
!
!    -------
!    WARNING:  This routine does no validity-checking of arguments.
!    -------
!
!  Fortran intrinsics used:  DABS, DMAX1, DMIN1.
!
!***END PROLOGUE  DPCHCI
!
! ----------------------------------------------------------------------
!
!  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
!                 Mathematics and Statistics Division,
!                 Lawrence Livermore National Laboratory.
!
!  Change record:
!    82-06-01   Modified end conditions to be continuous functions of
!               data when monotonicity switches in next interval.
!    82-06-02   1. Modified formulas so end conditions are less prone
!                  to over/underflow problems.
!               2. Minor modification to HSUM calculation.
!    82-08-05   Converted to SLATEC library version.
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!    1. The function  DPCHST(ARG1,ARG2)  is assumed to return zero if
!       either argument is zero, +1 if they are of the same sign, and
!       -1 if they are of opposite sign.
!    2. To produce a single precision version, simply:
!       a. Change DPCHCI to PCHCI wherever it occurs,
!       b. Change DPCHST to PCHST wherever it occurs,
!       c. Change all references to the Fortran intrinsics to their
!          single presision equivalents,
!       d. Change the double precision declarations to real, and
!       e. Change the constants ZERO and THREE to single precision.
!
!  DECLARE ARGUMENTS.
!
      INTEGER N, INCFD
      DOUBLE PRECISION  H(N), SLOPE(N), D(INCFD,N)
!
!  DECLARE LOCAL VARIBLES.
!
      INTEGER I, NLESS1
      DOUBLE PRECISION  DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, HSUM, &
           HSUMT3, THREE, W1, W2, ZERO
      DOUBLE PRECISION  DPCHST
!
!  INITIALIZE.
!
!      DATA  ZERO /0.D0/, THREE/3.D0/
       parameter (zero=0.D0, THREE=3.D0)
!***FIRST EXECUTABLE STATEMENT  DPCHCI
      NLESS1 = N - 1
      DEL1 = SLOPE(1)
!
!  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
!
      IF (NLESS1 .GT. 1)  GO TO 10
      D(1,1) = DEL1
      D(1,N) = DEL1
      GO TO 5000
!
!  NORMAL CASE  (N .GE. 3).
!
   10 CONTINUE
      DEL2 = SLOPE(2)
!
!  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!    SHAPE-PRESERVING.
!
      HSUM = H(1) + H(2)
      W1 = (H(1) + HSUM)/HSUM
      W2 = -H(1)/HSUM
      D(1,1) = W1*DEL1 + W2*DEL2
      IF ( DPCHST(D(1,1),DEL1) .LE. ZERO)  THEN
         D(1,1) = ZERO
      ELSE IF ( DPCHST(DEL1,DEL2) .LT. ZERO)  THEN
!       NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL1
         IF (DABS(D(1,1)) .GT. DABS(DMAX))  D(1,1) = DMAX
      ENDIF
!
!  LOOP THROUGH INTERIOR POINTS.
!
      DO 50  I = 2, NLESS1
         IF (I .EQ. 2)  GO TO 40
!
         HSUM = H(I-1) + H(I)
         DEL1 = DEL2
         DEL2 = SLOPE(I)
   40    CONTINUE
!
!       SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
!
         D(1,I) = ZERO
         IF ( DPCHST(DEL1,DEL2) .LE. ZERO)  GO TO 50
!
!       USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
!
         HSUMT3 = HSUM+HSUM+HSUM
         W1 = (HSUM + H(I-1))/HSUMT3
         W2 = (HSUM + H(I)  )/HSUMT3
         DMAX = DMAX1( DABS(DEL1), DABS(DEL2) )
         DMIN = DMIN1( DABS(DEL1), DABS(DEL2) )
         DRAT1 = DEL1/DMAX
         DRAT2 = DEL2/DMAX
         D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2)
!
   50 CONTINUE
!
!  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!    SHAPE-PRESERVING.
!
      W1 = -H(N-1)/HSUM
      W2 = (H(N-1) + HSUM)/HSUM
      D(1,N) = W1*DEL1 + W2*DEL2
      IF ( DPCHST(D(1,N),DEL2) .LE. ZERO)  THEN
         D(1,N) = ZERO
      ELSE IF ( DPCHST(DEL1,DEL2) .LT. ZERO)  THEN
!       NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL2
         IF (DABS(D(1,N)) .GT. DABS(DMAX))  D(1,N) = DMAX
      ENDIF
!
!  NORMAL RETURN.
!
 5000 CONTINUE
      RETURN
!------------- LAST LINE OF DPCHCI FOLLOWS -----------------------------
      END
      SUBROUTINE DPCHCS(SWITCH,N,H,SLOPE,D,INCFD,IERR)
!***BEGIN PROLOGUE  DPCHCS
!***REFER TO  DPCHIC
!***ROUTINES CALLED  DPCHST,PCHSW
!***REVISION DATE  870707   (YYMMDD)
!***DESCRIPTION
!
!        DPCHCS:  DPCHIC Monotonicity Switch Derivative Setter.
!
!    Called by  DPCHIC  to adjust the values of D in the vicinity of a
!    switch in direction of monotonicity, to produce a more "visually
!    pleasing" curve than that given by  DPCHIM .
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!       PARAMETER  (INCFD = ...)
!       INTEGER  N, IERR
!       DOUBLE PRECISION  SWITCH, H(N), SLOPE(N), D(INCFD,N)
!
!       CALL  DPCHCS (SWITCH, N, H, SLOPE, D, INCFD, IERR)
!
!   Parameters:
!
!    SWITCH -- (input) indicates the amount of control desired over
!          local excursions from data.
!
!    N -- (input) number of data points.  (assumes N.GT.2 .)
!
!    H -- (input) real*8 array of interval lengths.
!    SLOPE -- (input) real*8 array of data slopes.
!          If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!                 H(I) =  X(I+1)-X(I),
!             SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!
!    D -- (input) real*8 array of derivative values at the data points,
!          as determined by DPCHCI.
!         (output) derivatives in the vicinity of switches in direction
!          of monotonicity may be adjusted to produce a more "visually
!          pleasing" curve.
!          The value corresponding to X(I) is stored in
!               D(1+(I-1)*INCFD),  I=1(1)N.
!          No other entries in D are changed.
!
!    INCFD -- (input) increment between successive values in D.
!          This argument is provided primarily for 2-D applications.
!
!    IERR -- (output) error flag.  should be zero.
!          If negative, trouble in DPCHSW.  (should never happen.)
!
!    -------
!    WARNING:  This routine does no validity-checking of arguments.
!    -------
!
!  Fortran intrinsics used:  DABS, DMAX1, DMIN1.
!
!***END PROLOGUE  DPCHCS
!
! ----------------------------------------------------------------------
!
!  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
!                 Mathematics and Statistics Division,
!                 Lawrence Livermore National Laboratory.
!
!  Change record:
!    82-06-17   Redesigned to (1) fix  problem with lack of continuity
!               approaching a flat-topped peak (2) be cleaner and
!               easier to verify.
!               Eliminated subroutines PCHSA and PCHSX in the process.
!    82-06-22   1. Limited fact to not exceed one, so computed D is a
!                  convex combination of DPCHCI value and DPCHSD value.
!               2. Changed fudge from 1 to 4 (based on experiments).
!    82-06-23   Moved PCHSD to an inline function (eliminating MSWTYP).
!    82-08-05   Converted to SLATEC library version.
!    87-07-07   Corrected conversion to double precision.
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!    1. The function  DPCHST(ARG1,ARG2)  is assumed to return zero if
!       either argument is zero, +1 if they are of the same sign, and
!       -1 if they are of opposite sign.
!    2. To produce a single precision version, simply:
!       a. Change DPCHCS to PCHCS wherever it occurs,
!       b. Change DPCHSD to PCHSD wherever it occurs,
!       c. Change DPCHST to PCHST wherever it occurs,
!       d. Change DPCHSW to PCHSW wherever it occurs,
!       e. Change all references to the Fortran intrinsics to their
!          single precision equivalents,
!       f. Change the double precision declarations to real, and
!       g. Change the constants ZERO and ONE to single precision.
!
!  DECLARE ARGUMENTS.
!
      INTEGER N, INCFD, IERR
      DOUBLE PRECISION  SWITCH, H(N), SLOPE(N), D(INCFD,N)
!
!  DECLARE LOCAL VARIABLES.
!
      INTEGER I, INDX, K, NLESS1
      DOUBLE PRECISION  DEL(3), DEXT, DFLOC, DFMX, FACT, FUDGE, ONE,  &
           SLMAX, WTAVE(2), ZERO
      DOUBLE PRECISION  DPCHST
!
!  DEFINE INLINE FUNCTION FOR WEIGHTED AVERAGE OF SLOPES.
!
      DOUBLE PRECISION  DPCHSD, S1, S2, H1, H2
      DPCHSD(S1,S2,H1,H2) = (H2/(H1+H2))*S1 + (H1/(H1+H2))*S2
!
!  INITIALIZE.
!
!      DATA  ZERO /0.D0/,  ONE/1.D0/
!      DATA  FUDGE /4.D0/
      PARAMETER (ZERO=0.D0 , ONE=1.D0)
      FUDGE=4.D0
!***FIRST EXECUTABLE STATEMENT  DPCHCS
      IERR = 0
      NLESS1 = N - 1
!
!  LOOP OVER SEGMENTS.
!
      DO 900  I = 2, NLESS1
         IF ( DPCHST(SLOPE(I-1),SLOPE(I)) )  100, 300, 900
!            --------------------------
!
  100    CONTINUE
!
!....... SLOPE SWITCHES MONOTONICITY AT I-TH POINT .....................
!
!          DO NOT CHANGE D IF 'UP-DOWN-UP'.
            IF (I .GT. 2)  THEN
               IF ( DPCHST(SLOPE(I-2),SLOPE(I)) .GT. ZERO)  GO TO 900
!                  --------------------------
            ENDIF
            IF (I .LT. NLESS1)  THEN
               IF ( DPCHST(SLOPE(I+1),SLOPE(I-1)) .GT. ZERO)  GO TO 900
!                  ----------------------------
            ENDIF
!
!   ....... COMPUTE PROVISIONAL VALUE FOR D(1,I).
!
            DEXT = DPCHSD (SLOPE(I-1), SLOPE(I), H(I-1), H(I))
!
!   ....... DETERMINE WHICH INTERVAL CONTAINS THE EXTREMUM.
!
            IF ( DPCHST(DEXT, SLOPE(I-1)) )  200, 900, 250
!               -----------------------
!
  200       CONTINUE
!             DEXT AND SLOPE(I-1) HAVE OPPOSITE SIGNS --
!                       EXTREMUM IS IN (X(I-1),X(I)).
               K = I-1
!             SET UP TO COMPUTE NEW VALUES FOR D(1,I-1) AND D(1,I).
               WTAVE(2) = DEXT
               IF (K .GT. 1)  &
                 WTAVE(1) = DPCHSD (SLOPE(K-1), SLOPE(K), H(K-1), H(K))
               GO TO 400
!
  250       CONTINUE
!             DEXT AND SLOPE(I) HAVE OPPOSITE SIGNS --
!                       EXTREMUM IS IN (X(I),X(I+1)).
               K = I
!             SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
               WTAVE(1) = DEXT
               IF (K .LT. NLESS1)   &
                 WTAVE(2) = DPCHSD (SLOPE(K), SLOPE(K+1), H(K), H(K+1))
               GO TO 400
!
  300    CONTINUE
!
!....... AT LEAST ONE OF SLOPE(I-1) AND SLOPE(I) IS ZERO --
!                    CHECK FOR FLAT-TOPPED PEAK .......................
!
            IF (I .EQ. NLESS1)  GO TO 900
            IF ( DPCHST(SLOPE(I-1), SLOPE(I+1)) .GE. ZERO)  GO TO 900
!               -----------------------------
!
!          WE HAVE FLAT-TOPPED PEAK ON (X(I),X(I+1)).
            K = I
!          SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
            WTAVE(1) = DPCHSD (SLOPE(K-1), SLOPE(K), H(K-1), H(K))
            WTAVE(2) = DPCHSD (SLOPE(K), SLOPE(K+1), H(K), H(K+1))
!
  400    CONTINUE
!
!....... AT THIS POINT WE HAVE DETERMINED THAT THERE WILL BE AN EXTREMUM
!       ON (X(K),X(K+1)), WHERE K=I OR I-1, AND HAVE SET ARRAY WTAVE--
!          WTAVE(1) IS A WEIGHTED AVERAGE OF SLOPE(K-1) AND SLOPE(K),
!                   IF K.GT.1
!          WTAVE(2) IS A WEIGHTED AVERAGE OF SLOPE(K) AND SLOPE(K+1),
!                   IF K.LT.N-1
!
         SLMAX = DABS(SLOPE(K))
         IF (K .GT. 1)    SLMAX = DMAX1( SLMAX, DABS(SLOPE(K-1)) )
         IF (K.LT.NLESS1) SLMAX = DMAX1( SLMAX, DABS(SLOPE(K+1)) )
!
         IF (K .GT. 1)  DEL(1) = SLOPE(K-1) / SLMAX
         DEL(2) = SLOPE(K) / SLMAX
         IF (K.LT.NLESS1)  DEL(3) = SLOPE(K+1) / SLMAX
!
         IF ((K.GT.1) .AND. (K.LT.NLESS1))  THEN
!          NORMAL CASE -- EXTREMUM IS NOT IN A BOUNDARY INTERVAL.
            FACT = FUDGE* DABS(DEL(3)*(DEL(1)-DEL(2))*(WTAVE(2)/SLMAX))
            D(1,K) = D(1,K) + DMIN1(FACT,ONE)*(WTAVE(1) - D(1,K))
            FACT = FUDGE* DABS(DEL(1)*(DEL(3)-DEL(2))*(WTAVE(1)/SLMAX))
            D(1,K+1) = D(1,K+1) + DMIN1(FACT,ONE)*(WTAVE(2) - D(1,K+1))
         ELSE
!          SPECIAL CASE K=1 (WHICH CAN OCCUR ONLY IF I=2) OR
!                       K=NLESS1 (WHICH CAN OCCUR ONLY IF I=NLESS1).
            FACT = FUDGE* DABS(DEL(2))
            D(1,I) = DMIN1(FACT,ONE) * WTAVE(I-K+1)
!             NOTE THAT I-K+1 = 1 IF K=I  (=NLESS1),
!                       I-K+1 = 2 IF K=I-1(=1).
         ENDIF
!
!
!....... ADJUST IF NECESSARY TO LIMIT EXCURSIONS FROM DATA.
!
         IF (SWITCH .LE. ZERO)  GO TO 900
!
         DFLOC = H(K)*DABS(SLOPE(K))
         IF (K .GT. 1)    DFLOC = DMAX1(DFLOC, H(K-1)*DABS(SLOPE(K-1)))
         IF (K.LT.NLESS1) DFLOC = DMAX1(DFLOC, H(K+1)*DABS(SLOPE(K+1)))
         DFMX = SWITCH*DFLOC
         INDX = I-K+1
!       INDX = 1 IF K=I, 2 IF K=I-1.
!       ---------------------------------------------------------------
         CALL DPCHSW(DFMX, INDX, D(1,K), D(1,K+1), H(K), SLOPE(K), IERR)
!       ---------------------------------------------------------------
         IF (IERR .NE. 0)  RETURN
!
!....... END OF SEGMENT LOOP.
!
  900 CONTINUE
!
      RETURN
!------------- LAST LINE OF DPCHCS FOLLOWS -----------------------------
      END
      DOUBLE PRECISION FUNCTION DPCHDF(K,X,S,IERR)
!***BEGIN PROLOGUE  DPCHDF
!***REFER TO  DPCHCE,DPCHSP
!***ROUTINES CALLED  XERROR
!***REVISION DATE  870707   (YYMMDD)
!***DESCRIPTION
!
!         DPCHDF:   DPCHIP Finite Difference Formula
!
!    Uses a divided difference formulation to compute a K-point approx-
!    imation to the derivative at X(K) based on the data in X and S.
!
!    Called by  DPCHCE  and  DPCHSP  to compute 3- and 4-point boundary
!    derivative approximations.
!
! ----------------------------------------------------------------------
!
!    On input:
!       K      is the order of the desired derivative approximation.
!              K must be at least 3 (error return if not).
!       X      contains the K values of the independent variable.
!              X need not be ordered, but the values **MUST** be
!              distinct.  (Not checked here.)
!       S      contains the associated slope values:
!                 S(I) = (F(I+1)-F(I))/(X(I+1)-X(I)), I=1(1)K-1.
!              (Note that S need only be of length K-1.)
!
!    On return:
!       S      will be destroyed.
!       IERR   will be set to -1 if K.LT.2 .
!       DPCHDF  will be set to the desired derivative approximation if
!              IERR=0 or to zero if IERR=-1.
!
! ----------------------------------------------------------------------
!
!  Reference:  Carl de Boor, A Practical Guide to Splines, Springer-
!             Verlag (New York, 1978), pp. 10-16.
!
!***END PROLOGUE  DPCHDF
!
! ----------------------------------------------------------------------
!
!  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
!                 Mathematics and Statistics Division,
!                 Lawrence Livermore National Laboratory.
!
!  Change record:
!    82-08-05   Converted to SLATEC library version.
!    87-07-07   Corrected XERROR calls for d.p. name(s).
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!    To produce a single precision version, simply:
!       a. Change DPCHDF to PCHDF wherever it occurs,
!       b. Change the double precision declarations to real, and
!       c. Change the constant ZERO to single precision.
!
!  DECLARE ARGUMENTS.
!
      INTEGER K, IERR
      DOUBLE PRECISION  X(K), S(K)
!
!  DECLARE LOCAL VARIABLES.
!
      INTEGER I, J
      DOUBLE PRECISION  VALUE, ZERO
!      DATA  ZERO /0.D0/
      parameter (ZERO = 0.D0)
!
!  CHECK FOR LEGAL VALUE OF K.
!
!***FIRST EXECUTABLE STATEMENT  DPCHDF
      IF (K .LT. 3)  GO TO 5001
!
!  COMPUTE COEFFICIENTS OF INTERPOLATING POLYNOMIAL.
!
      DO 10  J = 2, K-1
         DO 9  I = 1, K-J
            S(I) = (S(I+1)-S(I))/(X(I+J)-X(I))
    9    CONTINUE
   10 CONTINUE
!
!  EVALUATE DERIVATIVE AT X(K).
!
      VALUE = S(1)
      DO 20  I = 2, K-1
         VALUE = S(I) + VALUE*(X(K)-X(I))
   20 CONTINUE
!
!  NORMAL RETURN.
!
      IERR = 0
      DPCHDF = VALUE
      RETURN
!
!  ERROR RETURN.
!
 5001 CONTINUE
!    K.LT.3 RETURN.
      IERR = -1
      CALL XERROR ('DPCHDF -- K LESS THAN THREE'  &
                , 0, IERR, 1)
      DPCHDF = ZERO
      RETURN
!------------- LAST LINE OF DPCHDF FOLLOWS -----------------------------
      END
    SUBROUTINE DPCHFE(N,X,F,D,INCFD,SKIP,NE,XE,FE,IERR)
!***BEGIN PROLOGUE  DPCHFE
!***DATE WRITTEN   811020   (YYMMDD)
!***REVISION DATE  870707   (YYMMDD)
!***CATEGORY NO.  E3
!***KEYWORDS  LIBRARY=SLATEC(PCHIP),
! TYPE=DOUBLE PRECISION(PCHFE-S DPCHFE-D),
! CUBIC HERMITE EVALUATION,HERMITE INTERPOLATION,
! PIECEWISE CUBIC EVALUATION
!***AUTHOR  FRITSCH, F. N., (LLNL)
! MATHEMATICS AND STATISTICS DIVISION
! LAWRENCE LIVERMORE NATIONAL LABORATORY
! P.O. BOX 808  (L-316)
! LIVERMORE, CA  94550
! FTS 532-4275, (415) 422-4275
!***PURPOSE  Evaluate a piecewise cubic Hermite function at an array of
! points.  May be used by itself for Hermite interpolation,
! or as an evaluator for DPCHIM or DPCHIC.
!***DESCRIPTION

! **** Double Precision version of PCHFE ****

! DPCHFE:  Piecewise Cubic Hermite Function Evaluator

! Evaluates the cubic Hermite function defined by  N, X, F, D  at
! the points  XE(J), J=1(1)NE.

! To provide compatibility with DPCHIM and DPCHIC, includes an
! increment between successive values of the F- and D-arrays.

! ----------------------------------------------------------------------

! Calling sequence:

! PARAMETER  (INCFD = ...)
! INTEGER  N, NE, IERR
! DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE)
! LOGICAL  SKIP

! CALL  DPCHFE (N, X, F, D, INCFD, SKIP, NE, XE, FE, IERR)

! Parameters:

! N -- (input) number of data points.  (Error return if N.LT.2 .)

! X -- (input) real*8 array of independent variable values.  The
! elements of X must be strictly increasing:
! X(I-1) .LT. X(I),  I = 2(1)N.
! (Error return if not.)

! F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
! the value corresponding to X(I).

! D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
! is the value corresponding to X(I).

! INCFD -- (input) increment between successive values in F and D.
! (Error return if  INCFD.LT.1 .)

! SKIP -- (input/output) logical variable which should be set to
! .TRUE. if the user wishes to skip checks for validity of
! preceding parameters, or to .FALSE. otherwise.
! This will save time in case these checks have already
! been performed (say, in DPCHIM or DPCHIC).
! SKIP will be set to .TRUE. on normal return.

! NE -- (input) number of evaluation points.  (Error return if
! NE.LT.1 .)

! XE -- (input) real*8 array of points at which the function is to
! be evaluated.

! NOTES:
! 1. The evaluation will be most efficient if the elements
! of XE are increasing relative to X;
! that is,   XE(J) .GE. X(I)
! implies    XE(K) .GE. X(I),  all K.GE.J .
! 2. If any of the XE are outside the interval [X(1),X(N)],
! values are extrapolated from the nearest extreme cubic,
! and a warning error is returned.

! FE -- (output) real*8 array of values of the cubic Hermite
! function defined by  N, X, F, D  at the points  XE.

! IERR -- (output) error flag.
! Normal return:
! IERR = 0  (no errors).
! Warning error:
! IERR.GT.0  means that extrapolation was performed at
! IERR points.
! "Recoverable" errors:
! IERR = -1  if N.LT.2 .
! IERR = -2  if INCFD.LT.1 .
! IERR = -3  if the X-array is not strictly increasing.
! IERR = -4  if NE.LT.1 .
! (The FE-array has not been changed in any of these cases.)
! NOTE:  The above errors are checked in the order listed,
! and following arguments have **NOT** been validated.

!***REFERENCES  (NONE)
!***ROUTINES CALLED  DCHFEV,XERROR
!***END PROLOGUE  DPCHFE

! ----------------------------------------------------------------------

! Change record:
! 82-08-03   Minor cosmetic changes for release 1.
! 87-07-07   Corrected XERROR calls for d.p. name(s).

! ----------------------------------------------------------------------

! Programming notes:

! 1. To produce a single precision version, simply:
! a. Change DPCHFE to PCHFE, and DCHFEV to CHFEV, wherever they
! occur,
! b. Change the double precision declaration to real,

! 2. Most of the coding between the call to DCHFEV and the end of
! the IR-loop could be eliminated if it were permissible to
! assume that XE is ordered relative to X.

! 3. DCHFEV does not assume that X1 is less than X2.  thus, it would
! be possible to write a version of DPCHFE that assumes a
! decreasing X-array by simply running the IR-loop backwards
! (and reversing the order of appropriate tests).

! 4. The present code has a minor bug, which I have decided is not
! worth the effort that would be required to fix it.
! If XE contains points in [X(N-1),X(N)], followed by points .LT.
! X(N-1), followed by points .GT.X(N), the extrapolation points
! will be counted (at least) twice in the total returned in IERR.

! DECLARE ARGUMENTS.

    INTEGER :: N, INCFD, NE, IERR
    real*8 ::  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE)
    LOGICAL ::  SKIP

! DECLARE LOCAL VARIABLES.

    INTEGER :: I, IERC, IR, J, JFIRST, NEXT(2), NJ

! VALIDITY-CHECK ARGUMENTS.

!***FIRST EXECUTABLE STATEMENT  DPCHFE
    IF (SKIP)  GO TO 5

    IF ( N < 2 )  GO TO 5001
    IF ( INCFD < 1 )  GO TO 5002
    DO 1  I = 2, N
        IF ( X(I) <= X(I-1) )  GO TO 5003
    1 ENDDO

! FUNCTION DEFINITION IS OK, GO ON.

    5 CONTINUE
    IF ( NE < 1 )  GO TO 5004
    IERR = 0
    SKIP = .TRUE.

! LOOP OVER INTERVALS.        (   INTERVAL INDEX IS  IL = IR-1  . )
! ( INTERVAL IS X(IL).LE.X.LT.X(IR) . )
    JFIRST = 1
    IR = 2
    10 CONTINUE

! SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.

    IF (JFIRST > NE)  GO TO 5000

! LOCATE ALL POINTS IN INTERVAL.

    DO 20  J = JFIRST, NE
        IF (XE(J) >= X(IR))  GO TO 30
    20 ENDDO
    J = NE + 1
    GO TO 40

! HAVE LOCATED FIRST POINT BEYOND INTERVAL.

    30 CONTINUE
    IF (IR == N)  J = NE + 1

    40 CONTINUE
    NJ = J - JFIRST

! SKIP EVALUATION IF NO POINTS IN INTERVAL.

    IF (NJ == 0)  GO TO 50

! EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .

! ----------------------------------------------------------------
    CALL DCHFEV (X(IR-1),X(IR), F(1,IR-1),F(1,IR), D(1,IR-1),D(1,IR) &
    ,NJ, XE(JFIRST), FE(JFIRST), NEXT, IERC)
! ----------------------------------------------------------------
    IF (IERC < 0)  GO TO 5005

    IF (NEXT(2) == 0)  GO TO 42
! IF (NEXT(2) .GT. 0)  THEN
! IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(2) TO THE
! RIGHT OF X(IR).

    IF (IR < N)  GO TO 41
! IF (IR .EQ. N)  THEN
! THESE ARE ACTUALLY EXTRAPOLATION POINTS.
    IERR = IERR + NEXT(2)
    GO TO 42
    41 CONTINUE
! ELSE
! WE SHOULD NEVER HAVE GOTTEN HERE.
    GO TO 5005
! ENDIF
! ENDIF
    42 CONTINUE

    IF (NEXT(1) == 0)  GO TO 49
! IF (NEXT(1) .GT. 0)  THEN
! IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE
! LEFT OF X(IR-1).

    IF (IR > 2)  GO TO 43
! IF (IR .EQ. 2)  THEN
! THESE ARE ACTUALLY EXTRAPOLATION POINTS.
    IERR = IERR + NEXT(1)
    GO TO 49
    43 CONTINUE
! ELSE
! XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
! EVALUATION INTERVAL.

! FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
    DO 44  I = JFIRST, J-1
        IF (XE(I) < X(IR-1))  GO TO 45
    44 ENDDO
! NOTE-- CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR
! IN DCHFEV.
    GO TO 5005

    45 CONTINUE
! RESET J.  (THIS WILL BE THE NEW JFIRST.)
    J = I

! NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
    DO 46  I = 1, IR-1
        IF (XE(J) < X(I)) GO TO 47
    46 ENDDO
! NB-- CAN NEVER DROP THROUGH HERE, SINCE XE(J).LT.X(IR-1).

    47 CONTINUE
! AT THIS POINT, EITHER  XE(J) .LT. X(1)
! OR      X(I-1) .LE. XE(J) .LT. X(I) .
! RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE
! CYCLING.
    IR = MAX0(1, I-1)
! ENDIF
! ENDIF
    49 CONTINUE

    JFIRST = J

! END OF IR-LOOP.

    50 CONTINUE
    IR = IR + 1
    IF (IR <= N)  GO TO 10

! NORMAL RETURN.

    5000 CONTINUE
    RETURN

! ERROR RETURNS.

    5001 CONTINUE
! N.LT.2 RETURN.
    IERR = -1
    CALL XERROR ('DPCHFE -- NUMBER OF DATA POINTS LESS THAN TWO' &
    , 0, IERR, 1)
    RETURN

    5002 CONTINUE
! INCFD.LT.1 RETURN.
    IERR = -2
    CALL XERROR ('DPCHFE -- INCREMENT LESS THAN ONE' &
    , 0, IERR, 1)
    RETURN

    5003 CONTINUE
! X-ARRAY NOT STRICTLY INCREASING.
    IERR = -3
    CALL XERROR ('DPCHFE -- X-ARRAY NOT STRICTLY INCREASING' &
    , 0, IERR, 1)
    RETURN

    5004 CONTINUE
! NE.LT.1 RETURN.
    IERR = -4
    CALL XERROR ('DPCHFE -- NUMBER OF EVALUATION POINTS LESS THAN ONE' &
    , 0, IERR, 1)
    RETURN

    5005 CONTINUE
! ERROR RETURN FROM DCHFEV.
! *** THIS CASE SHOULD NEVER OCCUR ***
    IERR = -5
    CALL XERROR ('DPCHFE -- ERROR RETURN FROM DCHFEV -- FATAL' &
    , 0, IERR, 2)
    RETURN
!------------- LAST LINE OF DPCHFE FOLLOWS -----------------------------
    end SUBROUTINE DPCHFE
    SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
!***BEGIN PROLOGUE  XERROR
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Processes an error (diagnostic) message.
!***DESCRIPTION
! Abstract
! XERROR processes a diagnostic message, in a manner
! determined by the value of LEVEL and the current value
! of the library error control flag, KONTRL.
! (See subroutine XSETF for details.)

! Description of Parameters
! --Input--
! MESSG - the Hollerith message to be processed, containing
! no more than 72 characters.
! NMESSG- the actual number of characters in MESSG.
! NERR  - the error number associated with this message.
! NERR must not be zero.
! LEVEL - error category.
! =2 means this is an unconditionally fatal error.
! =1 means this is a recoverable error.  (I.e., it is
! non-fatal if XSETF has been appropriately called.)
! =0 means this is a warning message only.
! =-1 means this is a warning message which is to be
! printed at most once, regardless of how many
! times this call is executed.

! Examples
! CALL XERROR('SMOOTH -- NUM WAS ZERO.',23,1,2)
! CALL XERROR('INTEG  -- LESS THAN FULL ACCURACY ACHIEVED.',
! 43,2,1)
! CALL XERROR('ROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL F
! 1ULLY COLLAPSED.',65,3,0)
! CALL XERROR('EXP    -- UNDERFLOWS BEING SET TO ZERO.',39,1,-1)

! Latest revision ---  19 MAR 1980
! Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
! HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
! 1982.
!***ROUTINES CALLED  XERRWV
!***END PROLOGUE  XERROR
    CHARACTER*(*) MESSG
!***FIRST EXECUTABLE STATEMENT  XERROR
    CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
    RETURN
    end SUBROUTINE XERROR
    SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
!***BEGIN PROLOGUE  XERRWV
!***DATE WRITTEN   800319   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Processes error message allowing 2 integer and two real
! values to be included in the message.
!***DESCRIPTION
! Abstract
! XERRWV processes a diagnostic message, in a manner
! determined by the value of LEVEL and the current value
! of the library error control flag, KONTRL.
! (See subroutine XSETF for details.)
! In addition, up to two integer values and two real
! values may be printed along with the message.

! Description of Parameters
! --Input--
! MESSG - the Hollerith message to be processed.
! NMESSG- the actual number of characters in MESSG.
! NERR  - the error number associated with this message.
! NERR must not be zero.
! LEVEL - error category.
! =2 means this is an unconditionally fatal error.
! =1 means this is a recoverable error.  (I.e., it is
! non-fatal if XSETF has been appropriately called.)
! =0 means this is a warning message only.
! =-1 means this is a warning message which is to be
! printed at most once, regardless of how many
! times this call is executed.
! NI    - number of integer values to be printed. (0 to 2)
! I1    - first integer value.
! I2    - second integer value.
! NR    - number of real values to be printed. (0 to 2)
! R1    - first real value.
! R2    - second real value.

! Examples
! CALL XERRWV('SMOOTH -- NUM (=I1) WAS ZERO.',29,1,2,
! 1   1,NUM,0,0,0.,0.)
! CALL XERRWV('QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM (
! 1R2).,54,77,1,0,0,0,2,ERRREQ,ERRMIN)

! Latest revision ---  19 MAR 1980
! Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
! HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
! 1982.
!***ROUTINES CALLED  FDUMP,I1MACH,J4SAVE,XERABT,XERCTL,XERPRT,XERSAV,
! XGETUA
!***END PROLOGUE  XERRWV
    CHARACTER*(*) MESSG
    CHARACTER(20) :: LFIRST
    CHARACTER(37) :: FORM
    DIMENSION LUN(5)
! GET FLAGS
!***FIRST EXECUTABLE STATEMENT  XERRWV
    LKNTRL = J4SAVE(2,0,.FALSE.)
    MAXMES = J4SAVE(4,0,.FALSE.)
! CHECK FOR VALID INPUT
    IF ((NMESSG > 0) .AND. (NERR /= 0) .AND. &
    (LEVEL >= (-1)) .AND. (LEVEL <= 2)) GO TO 10
    IF (LKNTRL > 0) CALL XERPRT('FATAL ERROR IN...',17)
    CALL XERPRT('XERROR -- INVALID INPUT',23)
    IF (LKNTRL > 0) CALL FDUMP
    IF (LKNTRL > 0) CALL XERPRT('JOB ABORT DUE TO FATAL ERROR.', &
    29)
    IF (LKNTRL > 0) CALL XERSAV(' ',0,0,0,KDUMMY)
    CALL XERABT('XERROR -- INVALID INPUT',23)
    RETURN
    10 CONTINUE
! RECORD MESSAGE
    JUNK = J4SAVE(1,NERR,.TRUE.)
    CALL XERSAV(MESSG,NMESSG,NERR,LEVEL,KOUNT)
! LET USER OVERRIDE
    LFIRST = MESSG
    LMESSG = NMESSG
    LERR = NERR
    LLEVEL = LEVEL
    CALL XERCTL(LFIRST,LMESSG,LERR,LLEVEL,LKNTRL)
! RESET TO ORIGINAL VALUES
    LMESSG = NMESSG
    LERR = NERR
    LLEVEL = LEVEL
    LKNTRL = MAX0(-2,MIN0(2,LKNTRL))
    MKNTRL = IABS(LKNTRL)
! DECIDE WHETHER TO PRINT MESSAGE
    IF ((LLEVEL < 2) .AND. (LKNTRL == 0)) GO TO 100
    IF (((LLEVEL == (-1)) .AND. (KOUNT > MIN0(1,MAXMES))) &
     .OR. ((LLEVEL == 0)   .AND. (KOUNT > MAXMES)) &
     .OR. ((LLEVEL == 1)   .AND. (KOUNT > MAXMES) .AND. (MKNTRL == 1)) &
     .OR. ((LLEVEL == 2)   .AND. (KOUNT > MAX0(1,MAXMES)))) GO TO 100
    IF (LKNTRL <= 0) GO TO 20
    CALL XERPRT(' ',1)
! INTRODUCTION
    IF (LLEVEL == (-1)) CALL XERPRT &
    ('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.',57)
    IF (LLEVEL == 0) CALL XERPRT('WARNING IN...',13)
    IF (LLEVEL == 1) CALL XERPRT &
    ('RECOVERABLE ERROR IN...',23)
    IF (LLEVEL == 2) CALL XERPRT('FATAL ERROR IN...',17)
    20 CONTINUE
! MESSAGE
    CALL XERPRT(MESSG,LMESSG)
    CALL XGETUA(LUN,NUNIT)
    ISIZEI = LOG10(FLOAT(I1MACH(9))) + 1.0
    ISIZEF = LOG10(FLOAT(I1MACH(10))**I1MACH(11)) + 1.0
    DO 50 KUNIT=1,NUNIT
        IUNIT = LUN(KUNIT)
        IF (IUNIT == 0) IUNIT = I1MACH(4)
        DO 22 I=1,MIN(NI,2)
            WRITE (FORM,21) I,ISIZEI
            21 FORMAT ('(11X,21HIN ABOVE MESSAGE, I',I1,'=,I',I2,')   ')
            IF (I == 1) WRITE (IUNIT,FORM) I1
            IF (I == 2) WRITE (IUNIT,FORM) I2
        22 ENDDO
        DO 24 I=1,MIN(NR,2)
            WRITE (FORM,23) I,ISIZEF+10,ISIZEF
            23 FORMAT ('(11X,21HIN ABOVE MESSAGE, R',I1,'=,E', &
            I2,'.',I2,')')
            IF (I == 1) WRITE (IUNIT,FORM) R1
            IF (I == 2) WRITE (IUNIT,FORM) R2
        24 ENDDO
        IF (LKNTRL <= 0) cycle
    ! ERROR NUMBER
        WRITE (IUNIT,30) LERR
        30 FORMAT (15H ERROR NUMBER =,I10)

    50 ENDDO
! TRACE-BACK
    IF (LKNTRL > 0) CALL FDUMP
    100 CONTINUE
    IFATAL = 0
    IF ((LLEVEL == 2) .OR. ((LLEVEL == 1) .AND. (MKNTRL == 2))) &
    IFATAL = 1
! QUIT HERE IF MESSAGE IS NOT FATAL
    IF (IFATAL <= 0) RETURN
    IF ((LKNTRL <= 0) .OR. (KOUNT > MAX0(1,MAXMES))) GO TO 120
! PRINT REASON FOR ABORT
    IF (LLEVEL == 1) CALL XERPRT &
    ('JOB ABORT DUE TO UNRECOVERED ERROR.',35)
    IF (LLEVEL == 2) CALL XERPRT &
    ('JOB ABORT DUE TO FATAL ERROR.',29)
! PRINT ERROR SUMMARY
    CALL XERSAV(' ',-1,0,0,KDUMMY)
    120 CONTINUE
! ABORT
    IF ((LLEVEL == 2) .AND. (KOUNT > MAX0(1,MAXMES))) LMESSG = 0
    CALL XERABT(MESSG,LMESSG)
    RETURN
    end SUBROUTINE XERRWV
    SUBROUTINE XERSAV(MESSG,NMESSG,NERR,LEVEL,ICOUNT)
!***BEGIN PROLOGUE  XERSAV
!***DATE WRITTEN   800319   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  Z
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Records that an error occurred.
!***DESCRIPTION
! Abstract
! Record that this error occurred.

! Description of Parameters
! --Input--
! MESSG, NMESSG, NERR, LEVEL are as in XERROR,
! except that when NMESSG=0 the tables will be
! dumped and cleared, and when NMESSG is less than zero the
! tables will be dumped and not cleared.
! --Output--
! ICOUNT will be the number of times this message has
! been seen, or zero if the table has overflowed and
! does not contain this message specifically.
! When NMESSG=0, ICOUNT will not be altered.

! Written by Ron Jones, with SLATEC Common Math Library Subcommittee
! Latest revision ---  19 Mar 1980
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
! HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
! 1982.
!***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
!***END PROLOGUE  XERSAV
    INTEGER :: LUN(5)
    CHARACTER*(*) MESSG
    CHARACTER(20) :: MESTAB(10),MES
    DIMENSION NERTAB(10),LEVTAB(10),KOUNT(10)
    SAVE MESTAB,NERTAB,LEVTAB,KOUNT,KOUNTX
! NEXT TWO DATA STATEMENTS ARE NECESSARY TO PROVIDE A BLANK
! ERROR TABLE INITIALLY
    DATA KOUNT(1),KOUNT(2),KOUNT(3),KOUNT(4),KOUNT(5), &
    KOUNT(6),KOUNT(7),KOUNT(8),KOUNT(9),KOUNT(10) &
    /0,0,0,0,0,0,0,0,0,0/
    DATA KOUNTX/0/
!***FIRST EXECUTABLE STATEMENT  XERSAV
    IF (NMESSG > 0) GO TO 80
! DUMP THE TABLE
    IF (KOUNT(1) == 0) RETURN
! PRINT TO EACH UNIT
    CALL XGETUA(LUN,NUNIT)
    DO 60 KUNIT=1,NUNIT
        IUNIT = LUN(KUNIT)
        IF (IUNIT == 0) IUNIT = I1MACH(4)
    ! PRINT TABLE HEADER
        WRITE (IUNIT,10)
        10 FORMAT (32H0          ERROR MESSAGE SUMMARY/ &
        51H MESSAGE START             NERR     LEVEL     COUNT)
    ! PRINT BODY OF TABLE
        DO 20 I=1,10
            IF (KOUNT(I) == 0) GO TO 30
            WRITE (IUNIT,15) MESTAB(I),NERTAB(I),LEVTAB(I),KOUNT(I)
            15 FORMAT (1X,A20,3I10)
        20 ENDDO
        30 CONTINUE
    ! PRINT NUMBER OF OTHER ERRORS
        IF (KOUNTX /= 0) WRITE (IUNIT,40) KOUNTX
        40 FORMAT (41H0OTHER ERRORS NOT INDIVIDUALLY TABULATED=,I10)
        WRITE (IUNIT,50)
        50 FORMAT (1X)
    60 ENDDO
    IF (NMESSG < 0) RETURN
! CLEAR THE ERROR TABLES
    DO 70 I=1,10
        KOUNT(I) = 0
    70 END DO
    KOUNTX = 0
    RETURN
    80 CONTINUE
! PROCESS A MESSAGE...
! SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
! OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
    MES = MESSG
    DO 90 I=1,10
        II = I
        IF (KOUNT(I) == 0) GO TO 110
        IF (MES /= MESTAB(I)) GO TO 90
        IF (NERR /= NERTAB(I)) GO TO 90
        IF (LEVEL /= LEVTAB(I)) GO TO 90
        GO TO 100
    90 ENDDO
! THREE POSSIBLE CASES...
! TABLE IS FULL
    KOUNTX = KOUNTX+1
    ICOUNT = 1
    RETURN
! MESSAGE FOUND IN TABLE
    100 KOUNT(II) = KOUNT(II) + 1
    ICOUNT = KOUNT(II)
    RETURN
! EMPTY SLOT FOUND FOR NEW MESSAGE
    110 MESTAB(II) = MES
    NERTAB(II) = NERR
    LEVTAB(II) = LEVEL
    KOUNT(II)  = 1
    ICOUNT = 1
    RETURN
    end SUBROUTINE XERSAV
    SUBROUTINE XGETUA(IUNITA,N)
!***BEGIN PROLOGUE  XGETUA
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Returns unit number(s) to which error messages are being
! sent.
!***DESCRIPTION
! Abstract
! XGETUA may be called to determine the unit number or numbers
! to which error messages are being sent.
! These unit numbers may have been set by a call to XSETUN,
! or a call to XSETUA, or may be a default value.

! Description of Parameters
! --Output--
! IUNIT - an array of one to five unit numbers, depending
! on the value of N.  A value of zero refers to the
! default unit, as defined by the I1MACH machine
! constant routine.  Only IUNIT(1),...,IUNIT(N) are
! defined by XGETUA.  The values of IUNIT(N+1),...,
! IUNIT(5) are not defined (for N .LT. 5) or altered
! in any way by XGETUA.
! N     - the number of units to which copies of the
! error messages are being sent.  N will be in the
! range from 1 to 5.

! Latest revision ---  19 MAR 1980
! Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
! HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
! 1982.
!***ROUTINES CALLED  J4SAVE
!***END PROLOGUE  XGETUA
    DIMENSION IUNITA(5)
!***FIRST EXECUTABLE STATEMENT  XGETUA
    N = J4SAVE(5,0,.FALSE.)
    DO 30 I=1,N
        INDEX = I+4
        IF (I == 1) INDEX = 3
        IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
    30 ENDDO
    RETURN
    end SUBROUTINE XGETUA


    SUBROUTINE FDUMP
!***BEGIN PROLOGUE  FDUMP
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  Z
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Symbolic dump (should be locally written).
!***DESCRIPTION
! ***Note*** Machine Dependent Routine
! FDUMP is intended to be replaced by a locally written
! version which produces a symbolic dump.  Failing this,
! it should be replaced by a version which prints the
! subprogram nesting list.  Note that this dump must be
! printed on each of up to five files, as indicated by the
! XGETUA routine.  See XSETUA and XGETUA for details.

! Written by Ron Jones, with SLATEC Common Math Library Subcommittee
! Latest revision ---  23 May 1979
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  FDUMP
!***FIRST EXECUTABLE STATEMENT  FDUMP
    RETURN
    end SUBROUTINE FDUMP
    FUNCTION J4SAVE(IWHICH,IVALUE,ISET)
!***BEGIN PROLOGUE  J4SAVE
!***REFER TO  XERROR
! Abstract
! J4SAVE saves and recalls several global variables needed
! by the library error handling routines.

! Description of Parameters
! --Input--
! IWHICH - Index of item desired.
! = 1 Refers to current error number.
! = 2 Refers to current error control flag.
! = 3 Refers to current unit number to which error
! messages are to be sent.  (0 means use standard.)
! = 4 Refers to the maximum number of times any
! message is to be printed (as set by XERMAX).
! = 5 Refers to the total number of units to which
! each error message is to be written.
! = 6 Refers to the 2nd unit for error messages
! = 7 Refers to the 3rd unit for error messages
! = 8 Refers to the 4th unit for error messages
! = 9 Refers to the 5th unit for error messages
! IVALUE - The value to be set for the IWHICH-th parameter,
! if ISET is .TRUE. .
! ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
! given the value, IVALUE.  If ISET=.FALSE., the
! IWHICH-th parameter will be unchanged, and IVALUE
! is a dummy parameter.
! --Output--
! The (old) value of the IWHICH-th parameter will be returned
! in the function value, J4SAVE.

! Written by Ron Jones, with SLATEC Common Math Library Subcommittee
! Adapted from Bell Laboratories PORT Library Error Handler
! Latest revision ---  23 MAY 1979
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
! HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
! 1982.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  J4SAVE
    LOGICAL :: ISET
    INTEGER :: IPARAM(9)
    SAVE IPARAM
    DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
    DATA IPARAM(5)/1/
    DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
!***FIRST EXECUTABLE STATEMENT  J4SAVE
    J4SAVE = IPARAM(IWHICH)
    IF (ISET) IPARAM(IWHICH) = IVALUE
    RETURN
    end FUNCTION J4SAVE
    SUBROUTINE XERABT(MESSG,NMESSG)
!***BEGIN PROLOGUE  XERABT
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Aborts program execution and prints error message.
!***DESCRIPTION
! Abstract
! ***Note*** machine dependent routine
! XERABT aborts the execution of the program.
! The error message causing the abort is given in the calling
! sequence, in case one needs it for printing on a dayfile,
! for example.

! Description of Parameters
! MESSG and NMESSG are as in XERROR, except that NMESSG may
! be zero, in which case no message is being supplied.

! Written by Ron Jones, with SLATEC Common Math Library Subcommittee
! Latest revision ---  19 MAR 1980
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
! HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
! 1982.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  XERABT
    CHARACTER*(*) MESSG
!***FIRST EXECUTABLE STATEMENT  XERABT
    STOP
    end SUBROUTINE XERABT
    SUBROUTINE XERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)
!***BEGIN PROLOGUE  XERCTL
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  R3C
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Allows user control over handling of individual errors.
!***DESCRIPTION
! Abstract
! Allows user control over handling of individual errors.
! Just after each message is recorded, but before it is
! processed any further (i.e., before it is printed or
! a decision to abort is made), a call is made to XERCTL.
! If the user has provided his own version of XERCTL, he
! can then override the value of KONTROL used in processing
! this message by redefining its value.
! KONTRL may be set to any value from -2 to 2.
! The meanings for KONTRL are the same as in XSETF, except
! that the value of KONTRL changes only for this message.
! If KONTRL is set to a value outside the range from -2 to 2,
! it will be moved back into that range.

! Description of Parameters

! --Input--
! MESSG1 - the first word (only) of the error message.
! NMESSG - same as in the call to XERROR or XERRWV.
! NERR   - same as in the call to XERROR or XERRWV.
! LEVEL  - same as in the call to XERROR or XERRWV.
! KONTRL - the current value of the control flag as set
! by a call to XSETF.

! --Output--
! KONTRL - the new value of KONTRL.  If KONTRL is not
! defined, it will remain at its original value.
! This changed value of control affects only
! the current occurrence of the current message.
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
! HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
! 1982.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  XERCTL
    CHARACTER(20) :: MESSG1
!***FIRST EXECUTABLE STATEMENT  XERCTL
    RETURN
    end SUBROUTINE XERCTL
    SUBROUTINE XERPRT(MESSG,NMESSG)
!***BEGIN PROLOGUE  XERPRT
!***DATE WRITTEN   790801   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  Z
!***KEYWORDS  ERROR,XERROR PACKAGE
!***AUTHOR  JONES, R. E., (SNLA)
!***PURPOSE  Prints error messages.
!***DESCRIPTION
! Abstract
! Print the Hollerith message in MESSG, of length NMESSG,
! on each file indicated by XGETUA.
! Latest revision ---  19 MAR 1980
!***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
! HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
! 1982.
!***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
!***END PROLOGUE  XERPRT
    INTEGER :: LUN(5)
    CHARACTER*(*) MESSG
! OBTAIN UNIT NUMBERS AND WRITE LINE TO EACH UNIT
!***FIRST EXECUTABLE STATEMENT  XERPRT
    CALL XGETUA(LUN,NUNIT)
    LENMES = LEN(MESSG)
    DO 20 KUNIT=1,NUNIT
        IUNIT = LUN(KUNIT)
        IF (IUNIT == 0) IUNIT = I1MACH(4)
        DO 10 ICHAR=1,LENMES,72
            LAST = MIN0(ICHAR+71 , LENMES)
            WRITE (IUNIT,'(1X,A)') MESSG(ICHAR:LAST)
        10 ENDDO
    20 ENDDO
    RETURN
    end SUBROUTINE XERPRT
      SUBROUTINE DPCHIC(IC,VC,SWITCH,N,X,F,D,INCFD,WK,NWK,IERR)
!***BEGIN PROLOGUE  DPCHIC
!***DATE WRITTEN   820218   (YYMMDD)
!***REVISION DATE  870707   (YYMMDD)
!***CATEGORY NO.  E1B
!***KEYWORDS  LIBRARY=SLATEC(PCHIP),
!            TYPE=DOUBLE PRECISION(PCHIC-S DPCHIC-D),
!            CUBIC HERMITE INTERPOLATION,MONOTONE INTERPOLATION,
!            PIECEWISE CUBIC INTERPOLATION,
!            SHAPE-PRESERVING INTERPOLATION
!***AUTHOR  FRITSCH, F. N., (LLNL)
!            MATHEMATICS AND STATISTICS DIVISION
!            LAWRENCE LIVERMORE NATIONAL LABORATORY
!            P.O. BOX 808  (L-316)
!            LIVERMORE, CA  94550
!            FTS 532-4275, (415) 422-4275
!***PURPOSE  Set derivatives needed to determine a piecewise monotone
!           piecewise cubic Hermite interpolant to given data.
!           User control is available over boundary conditions and/or
!           treatment of points where monotonicity switches direction.
!***DESCRIPTION
!
!      **** Double Precision version of PCHIC ****
!
!        DPCHIC:  Piecewise Cubic Hermite Interpolation Coefficients.
!
!    Sets derivatives needed to determine a piecewise monotone piece-
!    wise cubic interpolant to the data given in X and F satisfying the
!    boundary conditions specified by IC and VC.
!
!    The treatment of points where monotonicity switches direction is
!    controlled by argument SWITCH.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F- and D-arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by DPCHFE or DPCHFD.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!       PARAMETER  (INCFD = ...)
!       INTEGER  IC(2), N, NWK, IERR
!       DOUBLE PRECISION  VC(2), SWITCH, X(N), F(INCFD,N), D(INCFD,N),
!                         WK(NWK)
!
!       CALL DPCHIC (IC, VC, SWITCH, N, X, F, D, INCFD, WK, NWK, IERR)
!
!   Parameters:
!
!    IC -- (input) integer array of length 2 specifying desired
!          boundary conditions:
!          IC(1) = IBEG, desired condition at beginning of data.
!          IC(2) = IEND, desired condition at end of data.
!
!          IBEG = 0  for the default boundary condition (the same as
!                    used by DPCHIM).
!          If IBEG.NE.0, then its sign indicates whether the boundary
!                    derivative is to be adjusted, if necessary, to be
!                    compatible with monotonicity:
!             IBEG.GT.0  if no adjustment is to be performed.
!             IBEG.LT.0  if the derivative is to be adjusted for
!                    monotonicity.
!
!          Allowable values for the magnitude of IBEG are:
!          IBEG = 1  if first derivative at X(1) is given in VC(1).
!          IBEG = 2  if second derivative at X(1) is given in VC(1).
!          IBEG = 3  to use the 3-point difference formula for D(1).
!                    (Reverts to the default b.c. if N.LT.3 .)
!          IBEG = 4  to use the 4-point difference formula for D(1).
!                    (Reverts to the default b.c. if N.LT.4 .)
!          IBEG = 5  to set D(1) so that the second derivative is con-
!             tinuous at X(2). (Reverts to the default b.c. if N.LT.4.)
!             This option is somewhat analogous to the "not a knot"
!             boundary condition provided by DPCHSP.
!
!         NOTES (IBEG):
!          1. An error return is taken if ABS(IBEG).GT.5 .
!          2. Only in case  IBEG.LE.0  is it guaranteed that the
!             interpolant will be monotonic in the first interval.
!             If the returned value of D(1) lies between zero and
!             3*SLOPE(1), the interpolant will be monotonic.  This
!             is **NOT** checked if IBEG.GT.0 .
!          3. If IBEG.LT.0 and D(1) had to be changed to achieve mono-
!             tonicity, a warning error is returned.
!
!          IEND may take on the same values as IBEG, but applied to
!          derivative at X(N).  In case IEND = 1 or 2, the value is
!          given in VC(2).
!
!         NOTES (IEND):
!          1. An error return is taken if ABS(IEND).GT.5 .
!          2. Only in case  IEND.LE.0  is it guaranteed that the
!             interpolant will be monotonic in the last interval.
!             If the returned value of D(1+(N-1)*INCFD) lies between
!             zero and 3*SLOPE(N-1), the interpolant will be monotonic.
!             This is **NOT** checked if IEND.GT.0 .
!          3. If IEND.LT.0 and D(1+(N-1)*INCFD) had to be changed to
!             achieve monotonicity, a warning error is returned.
!
!    VC -- (input) real*8 array of length 2 specifying desired boundary
!          values, as indicated above.
!          VC(1) need be set only if IC(1) = 1 or 2 .
!          VC(2) need be set only if IC(2) = 1 or 2 .
!
!    SWITCH -- (input) indicates desired treatment of points where
!          direction of monotonicity switches:
!          Set SWITCH to zero if interpolant is required to be mono-
!          tonic in each interval, regardless of monotonicity of data.
!            NOTES:
!             1. This will cause D to be set to zero at all switch
!                points, thus forcing extrema there.
!             2. The result of using this option with the default boun-
!                dary conditions will be identical to using DPCHIM, but
!                will generally cost more compute time.
!                This option is provided only to facilitate comparison
!                of different switch and/or boundary conditions.
!          Set SWITCH nonzero to use a formula based on the 3-point
!             difference formula in the vicinity of switch points.
!          If SWITCH is positive, the interpolant on each interval
!             containing an extremum is controlled to not deviate from
!             the data by more than SWITCH*DFLOC, where DFLOC is the
!             maximum of the change of F on this interval and its two
!             immediate neighbors.
!          If SWITCH is negative, no such control is to be imposed.
!
!    N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!    X -- (input) real*8 array of independent variable values.  The
!          elements of X must be strictly increasing:
!               X(I-1) .LT. X(I),  I = 2(1)N.
!          (Error return if not.)
!
!    F -- (input) real*8 array of dependent variable values to be
!          interpolated.  F(1+(I-1)*INCFD) is value corresponding to
!          X(I).
!
!    D -- (output) real*8 array of derivative values at the data
!          points.  These values will determine a monotone cubic
!          Hermite function on each subinterval on which the data
!          are monotonic, except possibly adjacent to switches in
!          monotonicity. The value corresponding to X(I) is stored in
!               D(1+(I-1)*INCFD),  I=1(1)N.
!          No other entries in D are changed.
!
!    INCFD -- (input) increment between successive values in F and D.
!          This argument is provided primarily for 2-D applications.
!          (Error return if  INCFD.LT.1 .)
!
!    WK -- (scratch) real*8 array of working storage.  The user may
!          wish to know that the returned values are:
!             WK(I)     = H(I)     = X(I+1) - X(I) ;
!             WK(N-1+I) = SLOPE(I) = (F(1,I+1) - F(1,I)) / H(I)
!          for  I = 1(1)N-1.
!
!    NWK -- (input) length of work array.
!          (Error return if  NWK.LT.2*(N-1) .)
!
!    IERR -- (output) error flag.
!          Normal return:
!             IERR = 0  (no errors).
!          Warning errors:
!             IERR = 1  if IBEG.LT.0 and D(1) had to be adjusted for
!                       monotonicity.
!             IERR = 2  if IEND.LT.0 and D(1+(N-1)*INCFD) had to be
!                       adjusted for monotonicity.
!             IERR = 3  if both of the above are true.
!          "Recoverable" errors:
!             IERR = -1  if N.LT.2 .
!             IERR = -2  if INCFD.LT.1 .
!             IERR = -3  if the X-array is not strictly increasing.
!             IERR = -4  if ABS(IBEG).GT.5 .
!             IERR = -5  if ABS(IEND).GT.5 .
!             IERR = -6  if both of the above are true.
!             IERR = -7  if NWK.LT.2*(N-1) .
!            (The D-array has not been changed in any of these cases.)
!              NOTE:  The above errors are checked in the order listed,
!                  and following arguments have **NOT** been validated.
!
!***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE
!                CUBIC INTERPOLATION,' SIAM J.NUMER.ANAL. 17, 2 (APRIL
!                1980), 238-246.
!              2. F.N.FRITSCH AND J.BUTLAND, 'A METHOD FOR CONSTRUCTING
!                LOCAL MONOTONE PIECEWISE CUBIC INTERPOLANTS,' SIAM
!                J.SCI.STAT.COMPUT.5,2 (JUNE 1984), 300-304.
!              3. F.N.FRITSCH, 'PIECEWISE CUBIC INTERPOLATION PACKAGE,'
!                LLNL PREPRINT UCRL-87285 (JULY 1982).
!***ROUTINES CALLED  DPCHCE,DPCHCI,DPCHCS,XERROR
!***END PROLOGUE  DPCHIC
!
! ----------------------------------------------------------------------
!
!  Change record:
!    82-08-04   Converted to SLATEC library version.
!    87-07-07   Corrected XERROR calls for d.p. name(s).
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!    To produce a single precision version, simply:
!       a. Change DPCHIC to PCHIC wherever it occurs,
!       b. Change DPCHCE to PCHCE wherever it occurs,
!       c. Change DPCHCI to PCHCI wherever it occurs,
!       d. Change DPCHCS to PCHCS wherever it occurs,
!       e. Change the double precision declarations to real, and
!       f. Change the constant  ZERO  to single precision.
!
!  DECLARE ARGUMENTS.
!
      INTEGER IC(2), N, INCFD, NWK, IERR
      DOUBLE PRECISION  VC(2), SWITCH, X(N), F(INCFD,N), D(INCFD,N),  &
      WK(NWK)
!
!  DECLARE LOCAL VARIABLES.
!
      INTEGER I, IBEG, IEND, NLESS1
      DOUBLE PRECISION  ZERO
!      DATA  ZERO /0.D0/
      parameter (ZERO=0.D0)
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  DPCHIC
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
!
      IBEG = IC(1)
      IEND = IC(2)
      IERR = 0
      IF (IABS(IBEG) .GT. 5)  IERR = IERR - 1
      IF (IABS(IEND) .GT. 5)  IERR = IERR - 2
      IF (IERR .LT. 0)  GO TO 5004
!
!  FUNCTION DEFINITION IS OK -- GO ON.
!
      NLESS1 = N - 1
      IF ( NWK .LT. 2*NLESS1 )  GO TO 5007
!
!  SET UP H AND SLOPE ARRAYS.
!
      DO 20  I = 1, NLESS1
         WK(I) = X(I+1) - X(I)
         WK(NLESS1+I) = (F(1,I+1) - F(1,I)) / WK(I)
   20 CONTINUE
!
!  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
!
      IF (NLESS1 .GT. 1)  GO TO 1000
      D(1,1) = WK(2)
      D(1,N) = WK(2)
      GO TO 3000
!
!  NORMAL CASE  (N .GE. 3) .
!
 1000 CONTINUE
!
!  SET INTERIOR DERIVATIVES AND DEFAULT END CONDITIONS.
!
!    --------------------------------------
      CALL DPCHCI (N, WK(1), WK(N), D, INCFD)
!    --------------------------------------
!
!  SET DERIVATIVES AT POINTS WHERE MONOTONICITY SWITCHES DIRECTION.
!
      IF (SWITCH .EQ. ZERO)  GO TO 3000
!    ----------------------------------------------------
      CALL DPCHCS (SWITCH, N, WK(1), WK(N), D, INCFD, IERR)
!    ----------------------------------------------------
      IF (IERR .NE. 0)  GO TO 5008
!
!  SET END CONDITIONS.
!
 3000 CONTINUE
      IF ( (IBEG.EQ.0) .AND. (IEND.EQ.0) )  GO TO 5000
!    -------------------------------------------------------
      CALL DPCHCE (IC, VC, N, X, WK(1), WK(N), D, INCFD, IERR)
!    -------------------------------------------------------
      IF (IERR .LT. 0)  GO TO 5009
!
!  NORMAL RETURN.
!
 5000 CONTINUE
      RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!    N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('DPCHIC -- NUMBER OF DATA POINTS LESS THAN TWO' &
                , 0, IERR, 1)
      RETURN
!
 5002 CONTINUE
!    INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('DPCHIC -- INCREMENT LESS THAN ONE'  &
                , 0, IERR, 1)
      RETURN
!
 5003 CONTINUE
!    X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('DPCHIC -- X-ARRAY NOT STRICTLY INCREASING'  &
                , 0, IERR, 1)
      RETURN
!
 5004 CONTINUE
!    IC OUT OF RANGE RETURN.
      IERR = IERR - 3
      CALL XERROR ('DPCHIC -- IC OUT OF RANGE'   &
                , 0, IERR, 1)
      RETURN
!
 5007 CONTINUE
!    NWK .LT. 2*(N-1)  RETURN.
      IERR = -7
      CALL XERROR ('DPCHIC -- WORK ARRAY TOO SMALL'  &
                , 0, IERR, 1)
      RETURN
!
 5008 CONTINUE
!    ERROR RETURN FROM DPCHCS.
      IERR = -8
      CALL XERROR ('DPCHIC -- ERROR RETURN FROM DPCHCS' &
                , 0, IERR, 1)
      RETURN
!
 5009 CONTINUE
!    ERROR RETURN FROM DPCHCE.
!   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -9
      CALL XERROR ('DPCHIC -- ERROR RETURN FROM DPCHCE'  &
                , 0, IERR, 1)
      RETURN
!------------- LAST LINE OF DPCHIC FOLLOWS -----------------------------
      END
    SUBROUTINE DPCHSP(IC,VC,N,X,F,D,INCFD,WK,NWK,IERR)
!***BEGIN PROLOGUE  DPCHSP
!***DATE WRITTEN   820503   (YYMMDD)
!***REVISION DATE  870707   (YYMMDD)
!***CATEGORY NO.  E1B
!***KEYWORDS  LIBRARY=SLATEC(PCHIP),
! TYPE=DOUBLE PRECISION(PCHSP-S DPCHSP-D),
! CUBIC HERMITE INTERPOLATION,PIECEWISE CUBIC INTERPOLATION,
! SPLINE INTERPOLATION
!***AUTHOR  FRITSCH, F. N., (LLNL)
! MATHEMATICS AND STATISTICS DIVISION
! LAWRENCE LIVERMORE NATIONAL LABORATORY
! P.O. BOX 808  (L-316)
! LIVERMORE, CA  94550
! FTS 532-4275, (415) 422-4275
!***PURPOSE  Set derivatives needed to determine the Hermite represen-
! tation of the cubic spline interpolant to given data, with
! specified boundary conditions.
!***DESCRIPTION

! **** Double Precision version of PCHSP ****

! DPCHSP:   Piecewise Cubic Hermite Spline

! Computes the Hermite representation of the cubic spline inter-
! polant to the data given in X and F satisfying the boundary
! conditions specified by IC and VC.

! To facilitate two-dimensional applications, includes an increment
! between successive values of the F- and D-arrays.

! The resulting piecewise cubic Hermite function may be evaluated
! by DPCHFE or DPCHFD.

! NOTE:  This is a modified version of C. de Boor'S cubic spline
! routine CUBSPL.

! ----------------------------------------------------------------------

! Calling sequence:

! PARAMETER  (INCFD = ...)
! INTEGER  IC(2), N, NWK, IERR
! DOUBLE PRECISION  VC(2), X(N), F(INCFD,N), D(INCFD,N), WK(NWK)

! CALL  DPCHSP (IC, VC, N, X, F, D, INCFD, WK, NWK, IERR)

! Parameters:

! IC -- (input) integer array of length 2 specifying desired
! boundary conditions:
! IC(1) = IBEG, desired condition at beginning of data.
! IC(2) = IEND, desired condition at end of data.

! IBEG = 0  to set D(1) so that the third derivative is con-
! tinuous at X(2).  This is the "not a knot" condition
! provided by de Boor'S cubic spline routine CUBSPL.
! < This is the default boundary condition. >
! IBEG = 1  if first derivative at X(1) is given in VC(1).
! IBEG = 2  if second derivative at X(1) is given in VC(1).
! IBEG = 3  to use the 3-point difference formula for D(1).
! (Reverts to the default b.c. if N.LT.3 .)
! IBEG = 4  to use the 4-point difference formula for D(1).
! (Reverts to the default b.c. if N.LT.4 .)
! NOTES:
! 1. An error return is taken if IBEG is out of range.
! 2. For the "natural" boundary condition, use IBEG=2 and
! VC(1)=0.

! IEND may take on the same values as IBEG, but applied to
! derivative at X(N).  In case IEND = 1 or 2, the value is
! given in VC(2).

! NOTES:
! 1. An error return is taken if IEND is out of range.
! 2. For the "natural" boundary condition, use IEND=2 and
! VC(2)=0.

! VC -- (input) real*8 array of length 2 specifying desired boundary
! values, as indicated above.
! VC(1) need be set only if IC(1) = 1 or 2 .
! VC(2) need be set only if IC(2) = 1 or 2 .

! N -- (input) number of data points.  (Error return if N.LT.2 .)

! X -- (input) real*8 array of independent variable values.  The
! elements of X must be strictly increasing:
! X(I-1) .LT. X(I),  I = 2(1)N.
! (Error return if not.)

! F -- (input) real*8 array of dependent variable values to be
! interpolated.  F(1+(I-1)*INCFD) is value corresponding to
! X(I).

! D -- (output) real*8 array of derivative values at the data
! points.  These values will determine the cubic spline
! interpolant with the requested boundary conditions.
! The value corresponding to X(I) is stored in
! D(1+(I-1)*INCFD),  I=1(1)N.
! No other entries in D are changed.

! INCFD -- (input) increment between successive values in F and D.
! This argument is provided primarily for 2-D applications.
! (Error return if  INCFD.LT.1 .)

! WK -- (scratch) real*8 array of working storage.

! NWK -- (input) length of work array.
! (Error return if NWK.LT.2*N .)

! IERR -- (output) error flag.
! Normal return:
! IERR = 0  (no errors).
! "Recoverable" errors:
! IERR = -1  if N.LT.2 .
! IERR = -2  if INCFD.LT.1 .
! IERR = -3  if the X-array is not strictly increasing.
! IERR = -4  if IBEG.LT.0 or IBEG.GT.4 .
! IERR = -5  if IEND.LT.0 of IEND.GT.4 .
! IERR = -6  if both of the above are true.
! IERR = -7  if NWK is too small.
! NOTE:  The above errors are checked in the order listed,
! and following arguments have **NOT** been validated.
! (The D-array has not been changed in any of these cases.)
! IERR = -8  in case of trouble solving the linear system
! for the interior derivative values.
! (The D-array may have been changed in this case.)
! (             Do **NOT** use it!                )

!***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES, SPRINGER-
! VERLAG (NEW YORK, 1978), PP. 53-59.
!***ROUTINES CALLED  DPCHDF,XERROR
!***END PROLOGUE  DPCHSP

! ----------------------------------------------------------------------

! Change record:
! 82-08-04   Converted to SLATEC library version.
! 87-07-07   Corrected XERROR calls for d.p. name(s).

! ----------------------------------------------------------------------

! Programming notes:

! To produce a single precision version, simply:
! a. Change DPCHSP to PCHSP wherever it occurs,
! b. Change the double precision declarations to real, and
! c. Change the constants ZERO, HALF, ... to single precision.

! DECLARE ARGUMENTS.

    INTEGER :: IC(2), N, INCFD, NWK, IERR
    real*8 ::  VC(2), X(N), F(INCFD,N), D(INCFD,N), WK(2,N)

! DECLARE LOCAL VARIABLES.

    INTEGER :: IBEG, IEND, INDEX, J, NM1
    real*8 ::  G, HALF, ONE, STEMP(3), THREE, TWO, XTEMP(4), &
    ZERO
    real*8 ::  DPCHDF

!    DATA  ZERO /0.D0/, HALF/.5D0/, ONE/1.D0/, TWO/2.D0/, THREE/3.D0/
      parameter (ZERO=0.D0, HALF=0.5D0, ONE=1.D0, TWO=2.D0, THREE=3.D0)

! VALIDITY-CHECK ARGUMENTS.

!***FIRST EXECUTABLE STATEMENT  DPCHSP
    IF ( N < 2 )  GO TO 5001
    IF ( INCFD < 1 )  GO TO 5002
    DO 1  J = 2, N
        IF ( X(J) <= X(J-1) )  GO TO 5003
    1 ENDDO

    IBEG = IC(1)
    IEND = IC(2)
    IERR = 0
    IF ( (IBEG < 0) .OR. (IBEG > 4) )  IERR = IERR - 1
    IF ( (IEND < 0) .OR. (IEND > 4) )  IERR = IERR - 2
    IF ( IERR < 0 )  GO TO 5004

! FUNCTION DEFINITION IS OK -- GO ON.

    IF ( NWK < 2*N )  GO TO 5007

! COMPUTE FIRST DIFFERENCES OF X SEQUENCE AND STORE IN WK(1,.). ALSO,
! COMPUTE FIRST DIVIDED DIFFERENCE OF DATA AND STORE IN WK(2,.).
    DO 5  J=2,N
        WK(1,J) = X(J) - X(J-1)
        WK(2,J) = (F(1,J) - F(1,J-1))/WK(1,J)
    5 ENDDO

! SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.

    IF ( IBEG > N )  IBEG = 0
    IF ( IEND > N )  IEND = 0

! SET UP FOR BOUNDARY CONDITIONS.

    IF ( (IBEG == 1) .OR. (IBEG == 2) )  THEN
        D(1,1) = VC(1)
    ELSE IF (IBEG > 2)  THEN
    ! PICK UP FIRST IBEG POINTS, IN REVERSE ORDER.
        DO 10  J = 1, IBEG
            INDEX = IBEG-J+1
        ! INDEX RUNS FROM IBEG DOWN TO 1.
            XTEMP(J) = X(INDEX)
            IF (J < IBEG)  STEMP(J) = WK(2,INDEX)
        10 ENDDO
    ! --------------------------------
        D(1,1) = DPCHDF (IBEG, XTEMP, STEMP, IERR)
    ! --------------------------------
        IF (IERR /= 0)  GO TO 5009
        IBEG = 1
    ENDIF

    IF ( (IEND == 1) .OR. (IEND == 2) )  THEN
        D(1,N) = VC(2)
    ELSE IF (IEND > 2)  THEN
    ! PICK UP LAST IEND POINTS.
        DO 15  J = 1, IEND
            INDEX = N-IEND+J
        ! INDEX RUNS FROM N+1-IEND UP TO N.
            XTEMP(J) = X(INDEX)
            IF (J < IEND)  STEMP(J) = WK(2,INDEX+1)
        15 ENDDO
    ! --------------------------------
        D(1,N) = DPCHDF (IEND, XTEMP, STEMP, IERR)
    ! --------------------------------
        IF (IERR /= 0)  GO TO 5009
        IEND = 1
    ENDIF

! --------------------( BEGIN CODING FROM CUBSPL )--------------------

! **** A TRIDIAGONAL LINEAR SYSTEM FOR THE UNKNOWN SLOPES S(J) OF
! F  AT X(J), J=1,...,N, IS GENERATED AND THEN SOLVED BY GAUSS ELIM-
! INATION, WITH S(J) ENDING UP IN D(1,J), ALL J.
! WK(1,.) AND WK(2,.) ARE USED FOR TEMPORARY STORAGE.

! CONSTRUCT FIRST EQUATION FROM FIRST BOUNDARY CONDITION, OF THE FORM
! WK(2,1)*S(1) + WK(1,1)*S(2) = D(1,1)

    IF (IBEG == 0)  THEN
        IF (N == 2)  THEN
        ! NO CONDITION AT LEFT END AND N = 2.
            WK(2,1) = ONE
            WK(1,1) = ONE
            D(1,1) = TWO*WK(2,2)
        ELSE
        ! NOT-A-KNOT CONDITION AT LEFT END AND N .GT. 2.
            WK(2,1) = WK(1,3)
            WK(1,1) = WK(1,2) + WK(1,3)
            D(1,1) =((WK(1,2) + TWO*WK(1,1))*WK(2,2)*WK(1,3) &
            + WK(1,2)**2*WK(2,3)) / WK(1,1)
        ENDIF
    ELSE IF (IBEG == 1)  THEN
    ! SLOPE PRESCRIBED AT LEFT END.
        WK(2,1) = ONE
        WK(1,1) = ZERO
    ELSE
    ! SECOND DERIVATIVE PRESCRIBED AT LEFT END.
        WK(2,1) = TWO
        WK(1,1) = ONE
        D(1,1) = THREE*WK(2,2) - HALF*WK(1,2)*D(1,1)
    ENDIF

! IF THERE ARE INTERIOR KNOTS, GENERATE THE CORRESPONDING EQUATIONS AND
! CARRY OUT THE FORWARD PASS OF GAUSS ELIMINATION, AFTER WHICH THE J-TH
! EQUATION READS    WK(2,J)*S(J) + WK(1,J)*S(J+1) = D(1,J).

    NM1 = N-1
    IF (NM1 > 1)  THEN
        DO 20 J=2,NM1
            IF (WK(2,J-1) == ZERO)  GO TO 5008
            G = -WK(1,J+1)/WK(2,J-1)
            D(1,J) = G*D(1,J-1) &
            + THREE*(WK(1,J)*WK(2,J+1) + WK(1,J+1)*WK(2,J))
            WK(2,J) = G*WK(1,J-1) + TWO*(WK(1,J) + WK(1,J+1))
        20 ENDDO
    ENDIF

! CONSTRUCT LAST EQUATION FROM SECOND BOUNDARY CONDITION, OF THE FORM
! (-G*WK(2,N-1))*S(N-1) + WK(2,N)*S(N) = D(1,N)

! IF SLOPE IS PRESCRIBED AT RIGHT END, ONE CAN GO DIRECTLY TO BACK-
! SUBSTITUTION, SINCE ARRAYS HAPPEN TO BE SET UP JUST RIGHT FOR IT
! AT THIS POINT.
    IF (IEND == 1)  GO TO 30

    IF (IEND == 0)  THEN
        IF (N == 2 .AND. IBEG == 0)  THEN
        ! NOT-A-KNOT AT RIGHT ENDPOINT AND AT LEFT ENDPOINT AND N = 2.
            D(1,2) = WK(2,2)
            GO TO 30
        ELSE IF ((N == 2) .OR. (N == 3 .AND. IBEG == 0))  THEN
        ! EITHER (N=3 AND NOT-A-KNOT ALSO AT LEFT) OR (N=2 AND *NOT*
        ! NOT-A-KNOT AT LEFT END POINT).
            D(1,N) = TWO*WK(2,N)
            WK(2,N) = ONE
            IF (WK(2,N-1) == ZERO)  GO TO 5008
            G = -ONE/WK(2,N-1)
        ELSE
        ! NOT-A-KNOT AND N .GE. 3, AND EITHER N.GT.3 OR  ALSO NOT-A-
        ! KNOT AT LEFT END POINT.
            G = WK(1,N-1) + WK(1,N)
        ! DO NOT NEED TO CHECK FOLLOWING DENOMINATORS (X-DIFFERENCES).
            D(1,N) = ((WK(1,N)+TWO*G)*WK(2,N)*WK(1,N-1) &
            + WK(1,N)**2*(F(1,N-1)-F(1,N-2))/WK(1,N-1))/G
            IF (WK(2,N-1) == ZERO)  GO TO 5008
            G = -G/WK(2,N-1)
            WK(2,N) = WK(1,N-1)
        ENDIF
    ELSE
    ! SECOND DERIVATIVE PRESCRIBED AT RIGHT ENDPOINT.
        D(1,N) = THREE*WK(2,N) + HALF*WK(1,N)*D(1,N)
        WK(2,N) = TWO
        IF (WK(2,N-1) == ZERO)  GO TO 5008
        G = -ONE/WK(2,N-1)
    ENDIF

! COMPLETE FORWARD PASS OF GAUSS ELIMINATION.

    WK(2,N) = G*WK(1,N-1) + WK(2,N)
    IF (WK(2,N) == ZERO)   GO TO 5008
    D(1,N) = (G*D(1,N-1) + D(1,N))/WK(2,N)

! CARRY OUT BACK SUBSTITUTION

    30 CONTINUE
    DO 40 J=NM1,1,-1
        IF (WK(2,J) == ZERO)  GO TO 5008
        D(1,J) = (D(1,J) - WK(1,J)*D(1,J+1))/WK(2,J)
    40 ENDDO
! --------------------(  END  CODING FROM CUBSPL )--------------------

! NORMAL RETURN.

    RETURN

! ERROR RETURNS.

    5001 CONTINUE
! N.LT.2 RETURN.
    IERR = -1
    CALL XERROR ('DPCHSP -- NUMBER OF DATA POINTS LESS THAN TWO' &
    , 0, IERR, 1)
    RETURN

    5002 CONTINUE
! INCFD.LT.1 RETURN.
    IERR = -2
    CALL XERROR ('DPCHSP -- INCREMENT LESS THAN ONE' &
    , 0, IERR, 1)
    RETURN

    5003 CONTINUE
! X-ARRAY NOT STRICTLY INCREASING.
    IERR = -3
    CALL XERROR ('DPCHSP -- X-ARRAY NOT STRICTLY INCREASING' &
    , 0, IERR, 1)
    RETURN

    5004 CONTINUE
! IC OUT OF RANGE RETURN.
    IERR = IERR - 3
    CALL XERROR ('DPCHSP -- IC OUT OF RANGE' &
    , 0, IERR, 1)
    RETURN

    5007 CONTINUE
! NWK TOO SMALL RETURN.
    IERR = -7
    CALL XERROR ('DPCHSP -- WORK ARRAY TOO SMALL' &
    , 0, IERR, 1)
    RETURN

    5008 CONTINUE
! SINGULAR SYSTEM.
! *** THEORETICALLY, THIS CAN ONLY OCCUR IF SUCCESSIVE X-VALUES   ***
! *** ARE EQUAL, WHICH SHOULD ALREADY HAVE BEEN CAUGHT (IERR=-3). ***
    IERR = -8
    CALL XERROR ('DPCHSP -- SINGULAR LINEAR SYSTEM' &
    , 0, IERR, 1)
    RETURN

    5009 CONTINUE
! ERROR RETURN FROM DPCHDF.
! *** THIS CASE SHOULD NEVER OCCUR ***
    IERR = -9
    CALL XERROR ('DPCHSP -- ERROR RETURN FROM DPCHDF' &
    , 0, IERR, 1)
    RETURN
!------------- LAST LINE OF DPCHSP FOLLOWS -----------------------------
    end 
      DOUBLE PRECISION FUNCTION DPCHST(ARG1,ARG2)
!***BEGIN PROLOGUE  DPCHST
!***REFER TO  DPCHCE,DPCHCI,DPCHCS,DPCHIM
!***ROUTINES CALLED  (NONE)
!***REVISION DATE  870707   (YYMMDD)
!***DESCRIPTION
!
!        DPCHST:  DPCHIP Sign-Testing Routine.
!
!
!    Returns:
!       -1. if ARG1 and ARG2 are of opposite sign.
!        0. if either argument is zero.
!       +1. if ARG1 and ARG2 are of the same sign.
!
!    The object is to do this without multiplying ARG1*ARG2, to avoid
!    possible over/underflow problems.
!
!  Fortran intrinsics used:  SIGN.
!
!***END PROLOGUE  DPCHST
!
! ----------------------------------------------------------------------
!
!  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
!                 Mathematics and Statistics Division,
!                 Lawrence Livermore National Laboratory.
!
!  Change record:
!    82-08-05   Converted to SLATEC library version.
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!    To produce a single precision version, simply:
!       a. Change DPCHST to PCHST wherever it occurs,
!       b. Change all references to the Fortran intrinsics to their
!          single presision equivalents,
!       c. Change the double precision declarations to real, and
!       d. Change the constants  ZERO  and  ONE  to single precision.
!
!  DECLARE ARGUMENTS.
!
      DOUBLE PRECISION  ARG1, ARG2
!
!  DECLARE LOCAL VARIABLES.
!
      DOUBLE PRECISION  ONE, ZERO
!      DATA  ZERO /0.D0/,  ONE/1.D0/
      parameter (ZERO=0.D0, ONE=1.D0)
!
!  PERFORM THE TEST.
!
!***FIRST EXECUTABLE STATEMENT  DPCHST
      DPCHST = DSIGN(ONE,ARG1) * DSIGN(ONE,ARG2)
      IF ((ARG1.EQ.ZERO) .OR. (ARG2.EQ.ZERO))  DPCHST = ZERO
!
      RETURN
!------------- LAST LINE OF DPCHST FOLLOWS -----------------------------
      END
      SUBROUTINE DPCHSW(DFMAX,IEXTRM,D1,D2,H,SLOPE,IERR)
!***BEGIN PROLOGUE  DPCHSW
!***REFER TO  DPCHCS
!***ROUTINES CALLED  D1MACH,XERROR
!***REVISION DATE  870707   (YYMMDD)
!***DESCRIPTION
!
!        DPCHSW:  DPCHCS Switch Excursion Limiter.
!
!    Called by  DPCHCS  to adjust D1 and D2 if necessary to insure that
!    the extremum on this interval is not further than DFMAX from the
!    extreme data value.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!       INTEGER  IEXTRM, IERR
!       DOUBLE PRECISION  DFMAX, D1, D2, H, SLOPE
!
!       CALL  DPCHSW (DFMAX, IEXTRM, D1, D2, H, SLOPE, IERR)
!
!   Parameters:
!
!    DFMAX -- (input) maximum allowed difference between F(IEXTRM) and
!          the cubic determined by derivative values D1,D2.  (assumes
!          DFMAX.GT.0.)
!
!    IEXTRM -- (input) index of the extreme data value.  (assumes
!          IEXTRM = 1 or 2 .  Any value .NE.1 is treated as 2.)
!
!    D1,D2 -- (input) derivative values at the ends of the interval.
!          (Assumes D1*D2 .LE. 0.)
!         (output) may be modified if necessary to meet the restriction
!          imposed by DFMAX.
!
!    H -- (input) interval length.  (Assumes  H.GT.0.)
!
!    SLOPE -- (input) data SLOPE on the interval.
!
!    IERR -- (output) error flag.  should be zero.
!          If IERR=-1, assumption on D1 and D2 is not satisfied.
!          If IERR=-2, quadratic equation locating extremum has
!                      negative descriminant (should never occur).
!
!    -------
!    WARNING:  This routine does no validity-checking of arguments.
!    -------
!
!  Fortran intrinsics used:  DABS, DSIGN, DSQRT.
!
!***END PROLOGUE  DPCHSW
!
! ----------------------------------------------------------------------
!
!  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
!                 Mathematics and Statistics Division,
!                 Lawrence Livermore National Laboratory.
!
!  Change record:
!    82-08-05   Converted to SLATEC library version.
!    87-07-07   Replaced DATA statement for SMALL with a use of D1MACH.
!    87-07-07   Corrected XERROR calls for d.p. name(s).
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!    To produce a single precision version, simply:
!       a. Change DPCHSW to PCHSW wherever it occurs,
!       b. Change DPCHCS to PCHCS wherever it occurs,
!       c. Change D1MACH to R1MACH wherever it occurs,
!       d. Change all references to the Fortran intrinsics to their
!          single precision equivalents,
!       e. Change the double precision declarations to real, and
!       f. Change constants ZERO, ONE, TWO, ... to single precision.
!
!  DECLARE ARGUMENTS.
!
      INTEGER IEXTRM, IERR
      DOUBLE PRECISION  DFMAX, D1, D2, H, SLOPE

!  DECLARE LOCAL VARIABLES.
!
      DOUBLE PRECISION  CP, DMAX, FACT, LAMBDA, NU, ONE, PHI, RADCAL, &
                       RHO, SIGMA, SMALL, THAT, THIRD, THREE, TWO, ZERO
!      DATA  ZERO /0.D0/,  ONE /1.D0/,  TWO /2.D0/, THREE /3.D0/,   &
!           FACT /100.D0/
      parameter (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0)
!       THIRD SHOULD BE SLIGHTLY LESS THAN 1/3.
!      DATA  THIRD /0.33333D0/
      THIRD= 0.33333D0
      FACT=100.D0
!
!  NOTATION AND GENERAL REMARKS.
!
!    RHO IS THE RATIO OF THE DATA SLOPE TO THE DERIVATIVE BEING TESTED.
!    LAMBDA IS THE RATIO OF D2 TO D1.
!    THAT = T-HAT(RHO) IS THE NORMALIZED LOCATION OF THE EXTREMUM.
!    PHI IS THE NORMALIZED VALUE OF P(X)-F1 AT X = XHAT = X-HAT(RHO),
!          WHERE  THAT = (XHAT - X1)/H .
!       THAT IS, P(XHAT)-F1 = D*H*PHI,  WHERE D=D1 OR D2.
!    SIMILARLY,  P(XHAT)-F2 = D*H*(PHI-RHO) .
!
!    SMALL SHOULD BE A FEW ORDERS OF MAGNITUDE GREATER THAN MACHEPS.
!***FIRST EXECUTABLE STATEMENT  DPCHSW
      SMALL = FACT*D1MACH(4)
!
!  DO MAIN CALCULATION.
!
      IF (D1 .EQ. ZERO)  THEN
!
!       SPECIAL CASE -- D1.EQ.ZERO .
!
!         IF D2 IS ALSO ZERO, THIS ROUTINE SHOULD NOT HAVE BEEN CALLED.
         IF (D2 .EQ. ZERO)  GO TO 5001
!
         RHO = SLOPE/D2
!         EXTREMUM IS OUTSIDE INTERVAL WHEN RHO .GE. 1/3 .
         IF (RHO .GE. THIRD)  GO TO 5000
         THAT = (TWO*(THREE*RHO-ONE)) / (THREE*(TWO*RHO-ONE))
         PHI = THAT**2 * ((THREE*RHO-ONE)/THREE)
!
!         CONVERT TO DISTANCE FROM F2 IF IEXTRM.NE.1 .
         IF (IEXTRM .NE. 1)  PHI = PHI - RHO
!
!         TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
         DMAX = DFMAX / (H*DABS(PHI))
         IF (DABS(D2) .GT. DMAX)  D2 = DSIGN (DMAX, D2)
      ELSE
!
         RHO = SLOPE/D1
         LAMBDA = -D2/D1
         IF (D2 .EQ. ZERO)  THEN
!
!          SPECIAL CASE -- D2.EQ.ZERO .
!
!            EXTREMUM IS OUTSIDE INTERVAL WHEN RHO .GE. 1/3 .
            IF (RHO .GE. THIRD)  GO TO 5000
            CP = TWO - THREE*RHO
            NU = ONE - TWO*RHO
            THAT = ONE / (THREE*NU)
         ELSE
            IF (LAMBDA .LE. ZERO)  GO TO 5001
!
!          NORMAL CASE -- D1 AND D2 BOTH NONZERO, OPPOSITE SIGNS.
!
            NU = ONE - LAMBDA - TWO*RHO
            SIGMA = ONE - RHO
            CP = NU + SIGMA
            IF (DABS(NU) .GT. SMALL)  THEN
               RADCAL = (NU - (TWO*RHO+ONE))*NU + SIGMA**2
               IF (RADCAL .LT. ZERO)  GO TO 5002
               THAT = (CP - DSQRT(RADCAL)) / (THREE*NU)
            ELSE
               THAT = ONE/(TWO*SIGMA)
            ENDIF
         ENDIF
         PHI = THAT*((NU*THAT - CP)*THAT + ONE)
!
!         CONVERT TO DISTANCE FROM F2 IF IEXTRM.NE.1 .
         IF (IEXTRM .NE. 1)  PHI = PHI - RHO
!
!         TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
         DMAX = DFMAX / (H*DABS(PHI))
         IF (DABS(D1) .GT. DMAX)  THEN
            D1 = DSIGN (DMAX, D1)
            D2 = -LAMBDA*D1
         ENDIF
      ENDIF
!
!  NORMAL RETURN.
!
 5000 CONTINUE
      IERR = 0
      RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!    D1 AND D2 BOTH ZERO, OR BOTH NONZERO AND SAME SIGN.
      IERR = -1
      CALL XERROR ('DPCHSW -- D1 AND/OR D2 INVALID'  &
                , 0, IERR, 1)
      RETURN
!
 5002 CONTINUE
!    NEGATIVE VALUE OF RADICAL (SHOULD NEVER OCCUR).
      IERR = -2
      CALL XERROR ('DPCHSW -- NEGATIVE RADICAL'   &
                , 0, IERR, 1)
      RETURN
!------------- LAST LINE OF DPCHSW FOLLOWS -----------------------------
      END

