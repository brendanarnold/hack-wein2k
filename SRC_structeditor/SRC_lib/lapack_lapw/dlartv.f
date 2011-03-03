      SUBROUTINE DLARTV( N, X, INCX, Y, INCY, C, S, INCC )
      INTEGER            INCC, INCX, INCY, N
      DOUBLE PRECISION   C( * ), S( * ), X( * ), Y( * )
      INTEGER            I, IC, IX, IY
      DOUBLE PRECISION   XI, YI
      IX = 1
      IY = 1
      IC = 1
      DO 10 I = 1, N
         XI = X( IX )
         YI = Y( IY )
         X( IX ) = C( IC )*XI + S( IC )*YI
         Y( IY ) = C( IC )*YI - S( IC )*XI
         IX = IX + INCX
         IY = IY + INCY
         IC = IC + INCC
   10 CONTINUE
      RETURN
      END
