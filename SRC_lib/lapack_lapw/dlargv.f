      SUBROUTINE DLARGV( N, X, INCX, Y, INCY, C, INCC )
      INTEGER            INCC, INCX, INCY, N
      DOUBLE PRECISION   C( * ), X( * ), Y( * )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      INTEGER            I, IC, IX, IY
      DOUBLE PRECISION   TT, W, XI, YI
      INTRINSIC          ABS, MAX, SQRT
      IX = 1
      IY = 1
      IC = 1
      DO 10 I = 1, N
         XI = X( IX )
         YI = Y( IY )
         IF( XI.EQ.ZERO ) THEN
            C( IC ) = ZERO
            Y( IY ) = ONE
            X( IX ) = YI
         ELSE
            W = MAX( ABS( XI ), ABS( YI ) )
            XI = XI / W
            YI = YI / W
            TT = SQRT( XI*XI+YI*YI )
            C( IC ) = XI / TT
            Y( IY ) = YI / TT
            X( IX ) = W*TT
         END IF
         IX = IX + INCX
         IY = IY + INCY
         IC = IC + INCC
   10 CONTINUE
      RETURN
      END
