

      SUBROUTINE ORDER (X,M,N,INDX,AUX)
!
!  ORDERS ARRAY X(M,N) ACCORDING TO INDX(N).
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(M,N),AUX(N)
      INTEGER INDX(N)
      DO 3 J=1,M
         DO 1 I=1,N
            AUX(I)=X(J,INDX(I))
   1     CONTINUE
         DO 2 I=1,N
            X(J,I)=AUX(I)
   2     CONTINUE
   3  CONTINUE
      RETURN
      END
