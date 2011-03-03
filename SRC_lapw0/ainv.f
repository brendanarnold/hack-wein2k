      SUBROUTINE AINV (A,N1,TOL,IERR,N2)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NMAX=200)
      CHARACTER*67 ERRMSG
      DIMENSION A(N2,N2),T(NMAX),TT(NMAX),NI(NMAX)
      DATA ZER/0.0D0/,ONE/1.0D0/,ZILCH/1.D-30/
      NC=N1
      IF(N2.GT.NMAX) GOTO 900
      IF(NC.LE.N2)GOTO 1
      GOTO 910
    1 ZERO=TOL
      IF (ZERO.LT.ZILCH) ZERO=ZILCH
      DO 20 I=1,NC
   20 NI(I)=0
      DO 14 K=1,NC
      TEMP=ZER
      DO 4 I=1,NC
      IF((NI(I).NE.0).OR.(ABS(A(I,I)).LE.TEMP)) GOTO 4
      TEMP=ABS(A(I,I))
      II=I
    4 CONTINUE
      IF (TEMP.LE.ZERO)THEN
      IERR=K
      RETURN
      ENDIF
      NI(II)=1
      T(II)=ONE
      TT(II)=ONE/A(II,II)
      A(II,II)=ZER
      DO 6 I=1,II-1
      T(I)=A(I,II)
      TT(I)=TT(II)*A(I,II)
      IF (NI(I).EQ.0) TT(I)=-TT(I)
      A(I,II)=ZER
    6 CONTINUE
      DO 8 I=II+1,NC
      IF (NI(I).EQ.0) THEN
      T(I)=A(II,I)
      ELSE
      T(I)=-A(II,I)
      ENDIF
      TT(I)=-A(II,I)*TT(II)
      A(II,I)=ZER
    8 CONTINUE
      DO 12 I=1,NC
      DO 12 J=I,NC
   12 A(I,J)=A(I,J)+T(I)*TT(J)
   14 CONTINUE
      DO 13 I=1,NC-1
      DO 13 J=I+1,NC
   13 A(J,I)=A(I,J)
      IERR=0
      RETURN
!
!        Error messages
!
  900 CALL OUTERR('AINV','dimension exceeded')
      STOP 'AINV - Error'
  910 WRITE (ERRMSG,9000) NC
      CALL OUTERR('AINV',ERRMSG)
      STOP 'AINV - Error'

!
 9000 FORMAT (I5,' IS TOO LARGE FOR AINV-(60)')
      END
