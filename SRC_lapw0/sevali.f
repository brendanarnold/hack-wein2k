!
       DOUBLE PRECISION FUNCTION SEVALI(I, N, X, Y, B, C, D)
       REAL*8 X(N), Y(N), B(N), C(N), D(N), THIRD, DX
       INTEGER I, N, J
       DATA THIRD/.333333333333D00/
!
       J=I
       DX=X(J+1)-X(J)
       SEVALI=DX*(Y(J)+DX*(.5D00*B(J)+DX*(THIRD*C(J)+DX*.25D00*D(J))))
!
       RETURN
       END
