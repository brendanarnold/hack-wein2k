      SUBROUTINE PGBSYM (A,TAU,IZ,NTAU,SYM,NS,alat)

!  FINDS THE SPACE GROUP OF A LATTICE WITH A BASE.
!  WRITTEN BY J.SOLER FROM AN ALGORITHM BY K.KUNC. JUN/91

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (TOL=1.D-4,ONE=1.D0,HALF=.5D0)
      DIMENSION A(3,3),SYM(4,3,48),TAU(3,NTAU),IZ(NTAU), &
         B(3,3),STAU(3),STAU1(3),X(3),Y(3),alat(3)
      logical tric,trici
      common /tricl/ tric,trici
      trici=.false.
      IF (NTAU.LE.0) RETURN
      CALL RECLAT (A,B,0)
      NS0=NS
      NS=0
      NSS=0
!     FOR EVERY SYMMETRY OPERATION ON INPUT
      DO 15 IS=1,NS0
        DO 3 I=1,3
          STAU1(I)=DOT(SYM(1,I,IS),TAU(1,1),3)
   3    CONTINUE
        DO 11 JTAU1=1,NTAU
!         ASSUME ATOM 1 TRANSFORMS INTO ATOM JTAU1 AND
!         FIND TRANSLATION VECTOR FOR THIS TO BE TRUE
          IF (IZ(JTAU1).NE.IZ(1)) GO TO 11
          DO 4 I=1,3
            X(I)=TAU(I,JTAU1)-STAU1(I)
   4      CONTINUE
!         MOVE TRANSLATION VECTOR TO FIRST UNIT CELL.
!         Y IS TRANSLATION VECTOR IN RECIPROCAL SPACE
          DO 5 I=1,3
            D=DOT(B(1,I),X,3)
            Y(I)=D-NINT(D)
   5      CONTINUE
!         FOR EACH ATOM ITAU OTHER THAN 1, CHECK IF ITAU TRANSFORMS
!         INTO SOME JTAU, USING ASSUMED TRANSLATION VECTOR
          DO 10 ITAU=2,NTAU
            DO 6 I=1,3
              STAU(I)=DOT(SYM(1,I,IS),TAU(1,ITAU),3)
   6        CONTINUE
            DO 9 JTAU=1,NTAU
              IF (IZ(JTAU).NE.IZ(ITAU)) GO TO 9
              DO 7 I=1,3
                X(I)=TAU(I,JTAU)-STAU(I)
   7          CONTINUE
              DO 8 I=1,3
                D=DOT(B(1,I),X,3)-Y(I)
                D=D-NINT(D)
                IF (ABS(D).GT.TOL) GO TO 9
   8          CONTINUE
!             ITAU DOES TRANSFORM INTO JTAU. CHECK ANOTHER ITAU
              GO TO 10
   9        CONTINUE
!           ITAU DOES NOT TRANSFORM INTO ANY JTAU.
!           TRY ANOTHER TRANSLATION VECTOR 
            GO TO 11
  10      CONTINUE
!         TRANSLATION VECTOR SUCCESFUL FOR EVERY ATOM.
!         ACCEPT SYMMETRY OPERATION
          GO TO 12
  11    CONTINUE
!         NO SUCCESFUL TRANSLATION VECTOR FOUND.
!         REJECT SYMMETRY OPERATION AND TRY NEXT ONE
          GO TO 15
  12    NS=NS+1
        DO 14 I=1,3
          DO 13 J=1,4
            SWAP=SYM(J,I,IS)
            SYM(J,I,IS)=SYM(J,I,NS)
            SYM(J,I,NS)=SYM(J,I,NSS+1)
            SYM(J,I,NSS+1)=SWAP
  13      CONTINUE
          SYM(4,I,NSS+1)=A(I,1)*Y(1)+A(I,2)*Y(2)+A(I,3)*Y(3)
  14    CONTINUE
!       IF TRANSLATION VECTOR IS ZERO, OPERATION IS SYMMORPHIC
        IF (ABS(Y(1))+ABS(Y(2))+ABS(Y(3)).LE.TOL) NSS=NSS+1
  15  CONTINUE

! PRINT SOME RESULTS
!*     WRITE(6,*)'PGBSYM: ORDER OF LATICE POINT GROUP (NO BASE)  =',NS0
      WRITE(6,*)'PGBSYM: ORDER OF LATTICE SPACE GROUP (WITH BASE) =',NS
!*     WRITE(6,*)'PGBSYM: ORDER OF ORIGIN POINT GROUP (WITH BASE) = ',NSS
      IF (NSS.EQ.NS) WRITE(6,*) 'PGBSYM: SPACE GROUP IS SYMMORPHIC'
      IF (NSS.NE.NS) WRITE(6,*) 'PGBSYM: NON-SYMMORPHIC SPACE GROUP OR', &
                                  ' NON-STANDARD ORIGIN OF COORDINATES'
      DO 16 IS=1,NS
        D=ABS(SYM(1,1,IS)+ONE)+ABS(SYM(2,1,IS))+ABS(SYM(3,1,IS)) &
         +ABS(SYM(1,2,IS))+ABS(SYM(2,2,IS)+ONE)+ABS(SYM(3,2,IS)) &
         +ABS(SYM(1,3,IS))+ABS(SYM(2,3,IS))+ABS(SYM(3,3,IS)+ONE)
        IF (D.LT.TOL) THEN
          WRITE(6,*) 'PGBSYM: SPACE GROUP CONTAINS INVERSION'
          trici=.true.
          D=ABS(SYM(4,1,IS))+ABS(SYM(4,2,IS))+ABS(SYM(4,3,IS))
        if((ns0.eq.48).or.(ns0.eq.16).or.(ns0.eq.8)) then
          IF (D.GT.TOL) WRITE(6,'(8X,A,3F13.8)') &
      'BUT ATOMS SHOULD BE SHIFTED BY',(-HALF*SYM(4,I,IS)/alat(i),I=1,3)
          GO TO 99
        else 
          IF (D.GT.TOL) WRITE(6,'(8X,A,3F13.8,/,8x,a)') &
      'BUT ATOMS SHOULD BE SHIFTED BY',(-HALF*SYM(4,I,IS),I=1,3) &
       ,'(NOTE: You must convert carthesian to internal coordinates)'
          GO TO 99
        endif
        ENDIF
  16  CONTINUE
      WRITE(6,*) 'PGBSYM: SPACE GROUP DOES NOT CONTAIN INVERSION'
  99  CONTINUE
      END
