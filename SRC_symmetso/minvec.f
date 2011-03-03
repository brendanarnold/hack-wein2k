


      SUBROUTINE MINVEC(B0,BMIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (EPS=1.D-6)
      DOUBLE PRECISION &
         B0(3,3),BMIN(3,3),B(3,3),B2(3),BNEW(3),C(3,3),AUX(3,3)

!  FINDS THE LATTICE BASIS OF MINIMUN LENGTH, I.E. ANY OTHER BASIS
!  (NOT EQUIVALENT BY SYMMETRY) HAS ONE VECTOR LONGER.
!  WRITTEN BY J.MORENO AND J.SOLER (JSOLER AT EMDUAM11). AUGUST 1989.

      V0=ABS(VOLOFB(B0))
      IF (V0.LT.EPS) THEN
        WRITE(6,*) 'MINVEC: BASIS VECTORS ARE LINEARLY DEPENDENT'
        STOP
      ENDIF
      DO 20 I=1,3
        DO 10 J=1,3
          B(J,I)=B0(J,I)
  10    CONTINUE
        B2(I)=DOT(B(1,I),B(1,I),3)
  20  CONTINUE

  30  CONTINUE
      CALL ORDIX(B2,1,3,AUX)
      CALL ORDER(B2,1,3,AUX,AUX(1,2))
      CALL ORDER(B ,3,3,AUX,AUX(1,2))
      DO 50 I1=0,1
      DO 50 I2=-1,1
      DO 50 I3=-1,1
        IF (I1.EQ.0.AND.I2.NE.1) GO TO 50
        IF (I2.EQ.0.AND.I3.EQ.0) GO TO 50
        BNEW(1)=B(1,1)*I1+B(1,2)*I2+B(1,3)*I3
        BNEW(2)=B(2,1)*I1+B(2,2)*I2+B(2,3)*I3
        BNEW(3)=B(3,1)*I1+B(3,2)*I2+B(3,3)*I3
        BNEW2=DOT(BNEW,BNEW,3)
        DO 40 I=3,1,-1
          IF (BNEW2.GE.B2(I)) GO TO 50
            CALL VOLNEW(B,BNEW,I,VNEW)
            IF (ABS((VNEW-V0)/V0).LT.EPS) THEN
              B(1,I)=BNEW(1)
              B(2,I)=BNEW(2)
              B(3,I)=BNEW(3)
              B2(I)=BNEW2
              GO TO 30
            END IF
  40    CONTINUE
  50  CONTINUE

      IF (VOLOFB(B).LT.0.D0) THEN
        B(1,3)=-B(1,3)
        B(2,3)=-B(2,3)
        B(3,3)=-B(3,3)
      ENDIF
      CALL RECLAT(B0,AUX,0)
      DO 60 I=1,3
      DO 60 J=1,3
        C(J,I)=DNINT(DOT(AUX(1,J),B(1,I),3))
  60  CONTINUE
      DO 70 I=1,3
      DO 70 J=1,3
        B(J,I)=B0(J,1)*C(1,I)+B0(J,2)*C(2,I)+B0(J,3)*C(3,I)
  70  CONTINUE
      DO 80 I=1,3
      DO 80 J=1,3
        BMIN(J,I)=B(J,I)
  80  CONTINUE
      RETURN
      END
