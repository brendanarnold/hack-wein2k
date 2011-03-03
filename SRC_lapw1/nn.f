      SUBROUTINE NN(NAT)
!
!
!
      use char, only  : LATTIC
      use orth, only  : ORTHO
      use struk, only : POS, ALAT, ALPHA, RMT, MULT
      IMPLICIT NONE
      INTEGER  NAT
      INCLUDE 'param.inc'
!
!
      INTEGER            MM, NNN, N1, N2, N3, NC
      INTEGER            K, I3, JAT, L, JATOM, INDEX, M, I2, I1, J
      DOUBLE PRECISION   DIST
      DOUBLE PRECISION   PI, SINGAM, COSGAM, SUMRAD
      CHARACTER*67       ERRMSG
      PARAMETER          (NNN=10000)
      INTEGER            NR(NNN), NNAT(NNN)
      DOUBLE PRECISION   DIF(3),XX(3),PP(3),P(3),DISTS(NNN),PNN(3,NNN), HELP(3), BRNN(3,3)
!
!-----------------------------------------------------------------------
!
!
!
!     INDEX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE
!     NONEQUIVALENT ATOMS
!
!
!.....SET UP LATTICE, AND IF REQUIRED ROTATION MATRICES
      CALL DIRLAT (NAT,alpha,brnn)
      pi=4.d0*atan(1.d0)
      cosgam=cos(alpha(3))
      singam=sin(alpha(3))
      INDEX=0
!
      INDEX=0
      DO 200 JATOM=1,NAT
      DO 190 M=1,MULT(JATOM)
      INDEX=INDEX+1
      DO 150 J=1,3
  150 XX(J)=POS(J,INDEX)
      NC=0
          DO 180 I1=-2,2
          DO 180 I2=-2,2
          DO 180 I3=-2,2
          IF(ortho) THEN
            P(1)=I1*BRnn(1,1)+I2*BRnn(2,1)+I3*BRnn(3,1)
            P(2)=I1*BRnn(1,2)+I2*BRnn(2,2)+I3*BRnn(3,2)
            P(3)=I1*BRnn(1,3)+I2*BRnn(2,3)+I3*BRnn(3,3)
            ELSE
            P(1)=I1
            P(2)=I2
            P(3)=I3
            IF(LATTIC(1:3).eq.'CXZ') THEN
              P(1)=I1*0.5d0+i3*0.5d0
              P(2)=I2
              P(3)=-I1*0.5d0+i3*0.5d0
            END IF
           ENDIF
          K=0
      DO 120 JAT=1,NAT
      DO 110 MM=1,MULT(JAT)
      K=K+1
      DIST=0.d0
      DO 100 L=1,3
      PP(L)=POS(L,K)+P(L)
  100 DIF(L)=XX(L)-PP(L)
      IF (.not.ortho) THEN
        help(1)=dif(1)
        help(2)=dif(2)
        help(3)=dif(3)
      if(lattic(1:1).eq.'R') then
        dif(1)=help(1)*BRnn(1,1)+help(2)*BRnn(2,1)+help(3)*BRnn(3,1)
        dif(2)=help(1)*BRnn(1,2)+help(2)*BRnn(2,2)+help(3)*BRnn(3,2)
        dif(3)=help(1)*BRnn(1,3)+help(2)*BRnn(2,3)+help(3)*BRnn(3,3)
      elseif(lattic(1:3).eq.'CXZ') then
        dif(1)=help(1)*singam
        dif(2)=(help(1)*cosgam*alat(1)+help(2)*alat(2))/alat(2)
        dif(3)=help(3)
      else
        dif(1)=(help(1)*BRnn(1,1)*ALAT(1)+help(2)*BRnn(2,1)*ALAT(2) &
               +help(3)*BRnn(3,1)*ALAT(3))/ALAT(1)
        dif(2)=(help(1)*BRnn(1,2)*ALAT(1)+help(2)*BRnn(2,2)*ALAT(2) &
               +help(3)*BRnn(3,2)*ALAT(3))/ALAT(2)
        dif(3)=(help(1)*BRnn(1,3)*ALAT(1)+help(2)*BRnn(2,3)*ALAT(2) &
               +help(3)*BRnn(3,3)*ALAT(3))/ALAT(3)
      endif
      ENDIF
      DO 103 L=1,3
  103 DIST=DIST+DIF(L)*DIF(L)*ALAT(L)*ALAT(L)
      DIST=SQRT(DIST)
      IF(DIST.GT.40.01d0) GO TO 110
      IF(DIST.LT..001) GO TO 110
      NC=NC+1
      if(nc.gt.nnn) goto 900
      DISTS(NC)=DIST
      NNAT(NC)=JAT
      DO 105 L=1,3
  105 PNN(L,NC)=PP(L)
!     WRITE(6,2)JAT,NAME(JAT),PP(1),PP(2),PP(3),DIST
!   2 FORMAT(' TO ATOM:',I2,2X,A10,' AT',3F8.4,
!    * ' IS ',F10.5,' A.U.')
  110 CONTINUE
  120 CONTINUE
  180 CONTINUE
      CALL ORD2(DISTS,NR,NC)
      N1=1
      N2=NR(N1)
      N3=NNAT(N2)
      SUMRAD=RMT(JATOM)+RMT(N3)
      IF(M.EQ.1) THEN
      IF(SUMRAD.GT.DISTS(N1)) GOTO 910
      ENDIF
!
      DO 185 N1=1,NC
      N2=NR(N1)
      N3=NNAT(N2)
      SUMRAD=RMT(JATOM)+RMT(N3)
      IF(SUMRAD.GT.DISTS(N1)) GOTO 910
!      if(dists(n1).lt.8.d0*dists(1)) then
!      write(66,3) N3,NAME(N3),(PNN(L,N2),L=1,3),DISTS(N1)
!      end if
!    3 FORMAT(' ATOM:',I2,2X,A10,' AT',3F8.4,
!     * ' IS ',F10.5,' A.U.')
  185 CONTINUE
  190 CONTINUE
  200 CONTINUE
!

      RETURN

!        Error messages
!
  900 CALL OUTERR('NN','nnn too small')
      STOP 'NN - Error'
  910 CALL OUTERR('NN','overlapping spheres')
      WRITE (ERRMSG,9000) JATOM,RMT(JATOM),N3,RMT(N3)
      CALL OUTERR('NN',ERRMSG)
      WRITE (ERRMSG,9010) SUMRAD, DISTS(N1)
      CALL OUTERR('NN',ERRMSG)
      STOP 'NN - Error'
 9000 FORMAT('RMT(',I2,')=',F7.5,' AND RMT(',I2,')=',F7.5)
 9010 FORMAT('SUMS TO',F8.5,' GT NNN-DIST=',F8.5)
!
!        End of 'NN'
!
      END
!
      SUBROUTINE DIRLAT (NAT,alpha,brnn)
!
!     LATGEN GENERATES THE BRAVAIS MATRIX BRnn(3,3), WHICH TRANSFORMS
!     A RECIPROCAL LATTICE VECTOR OF A SPECIAL COORDINATE SYSTEM (IN
!     UNITS OF 2 PI/ALAT(I)) INTO CARTESIAN SYSTEM
!
      use char, only : LATTIC
      IMPLICIT NONE

      DOUBLE PRECISION  ALPHA(3)
      DOUBLE PRECISION  BRnn(3,3)
      INTEGER           NAT

      DOUBLE PRECISION PI, GAMMA, BETA, ALPH, COSG1, GAMMA1, SINAB, COSAB
!-----------------------------------------------------------------------
!
      pi=4.d0*atan(1.d0)
      gamma=alpha(3)
      beta=alpha(2)
      alph=alpha(1)
      cosg1=(cos(gamma)-cos(alph)*cos(beta))/sin(alph)/sin(beta)
      gamma1=acos(cosg1)
      IF(LATTIC(1:1).EQ.'H') GOTO 10
      IF(LATTIC(1:1).EQ.'S') GOTO 20
      IF(LATTIC(1:1).EQ.'P') GOTO 20
      IF(LATTIC(1:1).EQ.'B') GOTO 30
      IF(LATTIC(1:1).EQ.'F') GOTO 40
      IF(LATTIC(1:1).EQ.'C') GOTO 50
      IF(LATTIC(1:1).EQ.'R') GOTO 80
      GOTO 900
!
!.....HEXAGONAL CASE
 10   BRnn(1,1)=SQRT(3.)/2.
      BRnn(1,2)=-.5d0
      BRnn(1,3)=0.0d0
      BRnn(2,1)=0.0d0
      BRnn(2,2)=1.0d0
      BRnn(2,3)=0.0d0
      BRnn(3,1)=0.0d0
      BRnn(3,2)=0.0d0
      BRnn(3,3)=1.d0
      GOTO 100
!
!.....PRIMITIVE LATTICE CASE
 20   continue

      BRnn(1,1)=1.0d0*sin(gamma1)*sin(beta)
      BRnn(1,2)=1.0d0*cos(gamma1)*sin(beta)
      BRnn(1,3)=1.0d0*cos(beta)
      BRnn(2,1)=0.0d0
      BRnn(2,2)=1.0d0*sin(alph)
      BRnn(2,3)=1.0d0*cos(alph)
      BRnn(3,1)=0.0d0
      BRnn(3,2)=0.0d0
      BRnn(3,3)=1.0d0
      GOTO 100
!
!.....BC CASE (DIRECT LATTICE)
 30   CONTINUE
      BRnn(1,1)=-0.5d0
      BRnn(1,2)=0.5d0
      BRnn(1,3)=0.5d0
      BRnn(2,1)=0.5d0
      BRnn(2,2)=-0.5d0
      BRnn(2,3)=0.5d0
      BRnn(3,1)=0.5d0
      BRnn(3,2)=0.5d0
      BRnn(3,3)=-0.5d0
      GOTO 100
!
!.....FC CASE (DIRECT LATTICE)
 40   CONTINUE
      BRnn(1,1)=0.0d0
      BRnn(1,2)=0.5d0
      BRnn(1,3)=0.5d0
      BRnn(2,1)=0.5d0
      BRnn(2,2)=0.0d0
      BRnn(2,3)=0.5d0
      BRnn(3,1)=0.5d0
      BRnn(3,2)=0.5d0
      BRnn(3,3)=0.0d0
      GOTO 100
!
!.....CXY  CASE (DIRECT LATTICE)
 50   CONTINUE
      IF(LATTIC(2:3).EQ.'XZ') GOTO 60
      IF(LATTIC(2:3).EQ.'YZ') GOTO 70
      BRnn(1,1)=0.5d0
      BRnn(1,2)=-0.5d0
      BRnn(1,3)=0.0d0
      BRnn(2,1)=0.5d0
      BRnn(2,2)=0.5d0
      BRnn(2,3)=0.0d0
      BRnn(3,1)=0.0d0
      BRnn(3,2)=0.0d0
      BRnn(3,3)=1.0d0
      GOTO 100
!
!.....CXZ  CASE (DIRECT LATTICE)
 60   CONTINUE
!.....CXZ ORTHOROMBIC CASE
      IF(ABS(ALPHA(3)-PI/2.0D0).LT.0.0001) then
         BRnn(1,1)=0.5d0
         BRnn(1,2)=0.0d0
         BRnn(1,3)=-0.5d0
         BRnn(2,1)=0.0d0
         BRnn(2,2)=1.0d0
         BRnn(2,3)=0.0d0
         BRnn(3,1)=0.5d0
         BRnn(3,2)=0.0d0
         BRnn(3,3)=0.5d0
         GOTO 100
      ELSE
!.....CXZ MONOCLINIC CASE
         write(6,*) 'gamma not equal 90'
         SINAB=SIN(ALPHA(3))
         COSAB=COS(ALPHA(3))
!
         BRNN(1,1)=0.5d0*sinab
         BRNN(1,2)=0.5d0*cosab
         BRNN(1,3)=-0.5d0
         BRNN(2,1)=0.0d0
         BRNN(2,2)=1.0d0
         BRNN(2,3)=0.0d0
         BRNN(3,1)=0.5d0*sinab
         BRNN(3,2)=0.5d0*cosab
         BRNN(3,3)=0.5d0
         GOTO 100
      ENDIF
!
!.....CYZ  CASE (DIRECT LATTICE)
 70   CONTINUE
      BRnn(1,1)=1.0d0
      BRnn(1,2)=0.0d0
      BRnn(1,3)=0.0d0
      BRnn(2,1)=0.0d0
      BRnn(2,2)=0.5d0
      BRnn(2,3)=0.5d0
      BRnn(3,1)=0.0d0
      BRnn(3,2)=-0.5d0
      BRnn(3,3)=0.5d0
      GOTO 100
!.....RHOMBOHEDRAL CASE
 80   BRnn(1,1)=1/2.d0/sqrt(3.d0)
      BRnn(1,2)=-1/2.d0
      BRnn(1,3)=1/3.d0
      BRnn(2,1)=1/2.d0/SQRT(3.d0)
      BRnn(2,2)=1/2.D0
      BRnn(2,3)=1/3.d0
      BRnn(3,1)=-1/SQRT(3.d0)
      BRnn(3,2)=0.d0
      BRnn(3,3)=1/3.d0
      GOTO 100
!
 100  CONTINUE
!      write(6,*) 'Bravais Matrix:'
!      write(6,999) brnn
 999  format(3f15.5)
!
      RETURN
!
!        Error messages
!
  900 CALL OUTERR('NN','wrong lattice.')
      STOP 'NN - Error'
!
!        End of 'LATDIR'
!
      END
      SUBROUTINE ORD2(A,NR,IMAX)
!     ORDERS ELEMENTS IN ARRAY A INCREASING IN SIZE
!       REORDERS CORRESPONDING INDICES (NR)
      IMPLICIT NONE
      DOUBLE PRECISION A(*)
      INTEGER NR(*)
      INTEGER IMAX

      LOGICAL CONT
      INTEGER I, NHLP
      DOUBLE PRECISION HLP

      DO 50 I=1,IMAX
   50 NR(I)=I
  100 I=1
      CONT=.FALSE.
  110 I=I+1
      IF(A(I).LT.A(I-1))THEN
!       INTERCHANGE I AND (I-1) ELEMENT IN ARRAYS A AND NR
        HLP=A(I)
        A(I)=A(I-1)
        A(I-1)=HLP
        NHLP=NR(I)
        NR(I)=NR(I-1)
        NR(I-1)=NHLP
        CONT=.TRUE.
      ENDIF
      IF(I.LT.IMAX)GO TO 110
      IF(CONT) GO TO 100
      RETURN
      END
