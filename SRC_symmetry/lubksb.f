


      SUBROUTINE LUBKSB (A,N,NP,INDX,B)
!
! SOLVES A SYSTEM OF LINEAR EQUATIONS IN COMBINATION WITH 'LUDCMP'
! REF: W.PRESS ET AL. NUMERICAL RECIPIES. CAMBRIDGE U.P.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.D0,ONE=1.D0)
      DIMENSION A(NP,NP),INDX(N),B(N)
      II=0
      DO 12 I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)
         IF (II.NE.0) THEN
            DO 11 J=II,I-1
               SUM=SUM-A(I,J)*B(J)
  11        CONTINUE
         ELSE IF (SUM.NE.ZERO) THEN
            II=I
         ENDIF
         B(I)=SUM
  12  CONTINUE
      DO 14 I=N,1,-1
         SUM=B(I)
         IF (I.LT.N) THEN
            DO 13 J=I+1,N
               SUM=SUM-A(I,J)*B(J)
  13        CONTINUE
         ENDIF
         B(I)=SUM/A(I,I)
  14  CONTINUE
      RETURN
      END
