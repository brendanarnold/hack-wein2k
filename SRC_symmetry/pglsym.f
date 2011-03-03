      SUBROUTINE PGLSYM (A,SYM,NSYM)

!*  FINDS THE CRYSTAL SYSTEM AND THE POINT GROUP OF THE LATTICE.
!*  WRITTEN BY JOSE M. SOLER (JSOLER AT EMDUAM11) 22/4/90

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NRD=28,TOL=1.D-4,ZERO=0.D0,ONE=1.D0)
      DIMENSION A(3,3),SYM(4,3,48),R(3,NRD),R2(NRD), &
         B(3,3),B2(3),BNEW(3,3),BINV(3,3),U(3,3),AUX(3,3)
      logical tric,trici
      common /tricl/ tric,trici
      INTEGER NS(7)
      CHARACTER*12 NAME(7)
      DATA (NS(I),NAME(I),I=1,7)              /  2,'TRICLINIC   ', &
         4,'MONOCLINIC  ',   8,'ORTHORHOMBIC',  16,'TETRAGONAL  ', &
        48,'CUBIC       ',  12,'TRIGONAL    ',  24,'HEXAGONAL   '/

!* FIND MINIMUN LATTICE VECTORS B AND THEIR STARS R
! this has been disabled to treat supercells
      do i=1,3
      do j=1,3
      b(i,j)=a(i,j)
      enddo
      enddo
!      CALL MINVEC(A,B)
      B2(1)=B(1,1)**2+B(2,1)**2+B(3,1)**2
      B2(2)=B(1,2)**2+B(2,2)**2+B(3,2)**2
      B2(3)=B(1,3)**2+B(2,3)**2+B(3,3)**2
      NR=0
      DO 10 I1=-1,1
      DO 10 I2=-1,1
      DO 10 I3=-1,1
        R(1,NR+1)=B(1,1)*I1+B(1,2)*I2+B(1,3)*I3
        R(2,NR+1)=B(2,1)*I1+B(2,2)*I2+B(2,3)*I3
        R(3,NR+1)=B(3,1)*I1+B(3,2)*I2+B(3,3)*I3
        R2(NR+1)=R(1,NR+1)**2+R(2,NR+1)**2+R(3,NR+1)**2
        IF (ABS(R2(NR+1)-B2(1)).LT.TOL .OR. &
            ABS(R2(NR+1)-B2(2)).LT.TOL .OR. &
            ABS(R2(NR+1)-B2(3)).LT.TOL) NR=NR+1
   10 CONTINUE

!* FIND THE INVERSE OF MATRIX B
      CALL MATINV (B,3,3,BINV,AUX)

!* FIND SYMMETRY OPERATIONS OF THE LATTICE POINT GROUP
      CALL PUT (ZERO,SYM,4*3*48)
      NSYM=0
!*     ITERATE OVER POSSIBLE SETS OF TRANSFORMED BASIS VECTORS
      DO 70 I1=1,NR
      DO 70 I2=1,NR
      DO 70 I3=1,NR
        IF (ABS(R2(I1)-B2(1)).GT.TOL) GOTO 70
        IF (ABS(R2(I2)-B2(2)).GT.TOL) GOTO 70
        IF (ABS(R2(I3)-B2(3)).GT.TOL) GOTO 70
          DO 20 J=1,3
            BNEW(J,1)=R(J,I1)
            BNEW(J,2)=R(J,I2)
            BNEW(J,3)=R(J,I3)
   20     CONTINUE

!*         FIND THE TRANSFORMATION MATRIX U FOR U*B=BNEW
          DO 30 I=1,3
          DO 30 J=1,3
            U(J,I)=BNEW(J,1)*BINV(1,I)+ &
                   BNEW(J,2)*BINV(2,I)+ &
                   BNEW(J,3)*BINV(3,I)
   30     CONTINUE

!*         CHECK IF THE MATRIX IS UNITARY AND STORE IT
          DO 50 I=1,3
            DO 40 J=1,I
              D=U(1,J)*U(1,I)+U(2,J)*U(2,I)+U(3,J)*U(3,I)
              IF (J.EQ.I) D=D-ONE
              IF (ABS(D).GT.TOL) GOTO 70
   40       CONTINUE
   50     CONTINUE
          NSYM=NSYM+1
          DO 60 I=1,3
          DO 60 J=1,3
!           SYM(J,I,NSYM)=U(J,I)
            SYM(J,I,NSYM)=U(I,J)
   60     CONTINUE
   70 CONTINUE

!* PRINT THE NAME OF THE CRYSTAL SYSTEM
            tric=.false.
      DO 80 IHG=1,7
         IF (NSYM.EQ.NS(IHG)) THEN
            WRITE (6,*) 'PGLSYM: THE CRYSTAL SYSTEM IS ',NAME(IHG)
            WRITE (6,*) 'PGLSYM: ORDER OF LATTICE POINT GROUP', &
                        ' (NO BASE) =  ',NSYM
            if(ihg.eq.1) tric=.true.
 !           if(ihg.eq.2) tric=.true.
            RETURN
         ENDIF
   80 CONTINUE
      WRITE (6,*) 'PLGSYM: CRYSTAL SYSTEM NOT FOUND. NR,NSYM=',NR,NSYM
      END
