!
       DOUBLE PRECISION FUNCTION SEVALIN(I, N, X, Y, B, C, D)
       REAL*8 X(N), Y(N), B(N), C(N), D(N), DX, DX3
       INTEGER I, N, J, J1
!
       J =I
       J1=I+1
       DX=X(J+1)-X(J)
       DX3=DX*DX*DX
       SEVALIN=.5D0*DX*(Y(J)+Y(J1))-DX3*(C(J)+C(J1))/12.D0
!
       RETURN
       END
