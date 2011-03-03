      SUBROUTINE radint(JATOM,is)
      use potnlc
      use intu
!ad
!ad  calculation of radial integrals:
!ad  U   = u
!ad  DU  = du/dr
!ad  UP  = du/dE
!ad  DUP = dUP/dr
!ad 
!ad  for LO's: 
!ad  UL  = u_lo
!ad  DUL = du_lo/dr
!ad
!ad  the radial functions are calculated in ATPAR and LOMAIN
!ad
      INCLUDE 'param.inc'
!
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL  lapw(0:lmax2),loor(0:lomax),lloor(0:lmax2) 
!ad
      COMMON /XA/     R(NRAD),BK(3)
!      COMMON /POTNLC/ VR(NRAD),RZERO(NATO),DX(NATO),JRI(NATO)
      COMMON /RADFU/  U(NRAD,LMAX1),UP(NRAD,LMAX1), &
                      DU(NRAD,LMAX1),DUP(NRAD,LMAX1)    
      COMMON /UHELP/  UA(NRAD), UB(NRAD), UAB(NRAD), UBA(NRAD), &
                      UAA(NRAD),UBB(NRAD)
      common /loabc/   alo(0:lomax,nloat),blo(0:lomax,nloat), &
                       clo(0:lomax,nloat), &
                       elo(0:lomax,nloat),plo(0:lomax), &
                       dplo(0:lomax),pelo(0:lomax), &
                       dpelo(0:lomax),peilo(0:lomax),      &
                       pi12lo(0:lomax),pe12lo(0:lomax), &
                       UL(nrad,lomax+1),DUL(nrad,lomax+1)
      common /lolog/  nlo,nlov,nlon,lapw,ilo(0:lomax),loor,lloor
!      COMMON /INTU/   Duu1(LMAX1,NATO,2),Duu2(LMAX1,NATO,2), &
!                      Duup1(LMAX1,NATO,2),Duup2(LMAX1,NATO,2), &
!                      Dupu1(LMAX1,NATO,2),Dupu2(LMAX1,NATO,2), &
!                      Dupup1(LMAX1,NATO,2),Dupup2(LMAX1,NATO,2), &
!                      Ruu(LMAX1,NATO,2),Ruup(LMAX1,NATO,2), &
!                      Rupu(LMAX1,NATO,2),Rupup(LMAX1,NATO,2)
!      COMMON /INTUL/  Duul1(LOMAX1,NATO,2), Duul2(LOMAX1,NATO,2), &
!                      Dulup1(LOMAX1,NATO,2),Dulup2(LOMAX1,NATO,2), &
!                      Dupul1(LOMAX1,NATO,2),Dupul2(LOMAX1,NATO,2), &
!                      Dulu1(LOMAX1,NATO,2), Dulu2(LOMAX1,NATO,2), &
!                      Dulul1(LOMAX1,NATO,2),Dulul2(LOMAX1,NATO,2), &
!                      Ruul(LOMAX1,NATO,2),  Rulu(LOMAX1,NATO,2), &
!                      Rupul(LOMAX1,NATO,2), Rulup(LOMAX1,NATO,2), &
!                      Rulul(LOMAX1,NATO,2)
  
!ad
!ad..................calculation of radial integrals....................
!ad
      TWO=2.d0
      IMAX=JRI(JATOM)
      DO 2 I=1,IMAX
      R(I)=Rnot(JATOM)*EXP((I-1)*DX(JATOM))
    2 continue
!ad
!ad..................RUU(L)..RUUP(L)..RUPU(L)..RUPUP(L).................
!ad
      DO 103 L=1,LMAX
      DO 105 I=1,IMAX
        SR=SQRT(R(I))
        UA(I) =U(I,L)/SR
        UB(I) =U(I,L+1)/SR
        UAB(I)=UP(I,L)/SR
        UBA(I)=UP(I,L+1)/SR
  105 CONTINUE
      CALL RINT13(.FALSE.,UA,UA,UB,UB,RUU(L,JATOM,is),JATOM)
      CALL RINT13(.FALSE.,UB,UB,UAB,UAB,RUUP(L,JATOM,is),JATOM)
      CALL RINT13(.FALSE.,UBA,UBA,UA,UA,RUPU(L,JATOM,is),JATOM)
      CALL RINT13(.FALSE.,UBA,UBA,UAB,UAB,RUPUP(L,JATOM,is),JATOM)
  103 CONTINUE
!
!ad
!ad...............DUU1(L)..DUUP1(L)..DUPU1(L)..DUPUP1(L)................
!ad
      DO 203 L=1,LMAX
      DO 205 I=1,IMAX
        UA(I) =U(I,L+1)
        UB(I) =DU(I,L)*TWO
        UAB(I)=DUP(I,L)*TWO
        UBA(I)=UP(I,L+1)
  205 CONTINUE
      CALL RINT13(.FALSE.,UA,UA,UB,UB,DUU1(L,JATOM,is),JATOM)
      CALL RINT13(.FALSE.,UA,UA,UAB,UAB,DUUP1(L,JATOM,is),JATOM)
      CALL RINT13(.FALSE.,UBA,UBA,UB,UB,DUPU1(L,JATOM,is),JATOM)
      CALL RINT13(.FALSE.,UBA,UBA,UAB,UAB,DUPUP1(L,JATOM,is),JATOM)
  203 CONTINUE

!ad
!ad...............DUU2(L)..DUUP2(L)..DUPU2(L)..DUPUP2(L)................
!ad
      DO 303 L=1,LMAX
      DO 305 I=1,IMAX
        UA(I) =U(I,L)
        UB(I) =DU(I,L+1)*TWO
        UAB(I)=DUP(I,L+1)*TWO
        UBA(I)=UP(I,L)
  305 CONTINUE
      CALL RINT13(.FALSE.,UA,UA,UB,UB,DUU2(L,JATOM,is),JATOM)
      CALL RINT13(.FALSE.,UA,UA,UAB,UAB,DUUP2(L,JATOM,is),JATOM)
      CALL RINT13(.FALSE.,UBA,UBA,UB,UB,DUPU2(L,JATOM,is),JATOM)
      CALL RINT13(.FALSE.,UBA,UBA,UAB,UAB,DUPUP2(L,JATOM,is),JATOM)
  303 CONTINUE

!ad
!ad..............................for LO's...............................
!ad
      DO 503 L=0,LOMAX
      L1=L+1
      L11=L1+1
      IF (.not.loor(L)) THEN
        DUUL1(L1,JATOM,is) =  0
        DUPUL1(L1,JATOM,is) = 0
        DULU2(L1,JATOM,is)  = 0
        DULUP2(L1,JATOM,is) = 0
        RUUL(L1,JATOM,is)   = 0
        RUPUL(L1,JATOM,is)  = 0
      ELSE 
        DO 505 I=1,IMAX
        SR=SQRT(R(I))
        UA(I) =U (I,L11)/SR
        UB(I) =UP(I,L11)/SR
        UAB(I)=UL(I,L1) /SR
  505   CONTINUE
!ad
!ad............DUUL1..DUPUL1..DULU2..DULUP2..RUUL..RUPUL................
!ad
        CALL RINT13(.FALSE.,U(1,L11),U(1,L11),DUL(1,L1),DUL(1,L1) &
                  ,DUUL1(L1,JATOM,is),JATOM)
        DUUL1(L1,JATOM,is)=TWO*DUUL1(L1,JATOM,is)
        CALL RINT13(.FALSE.,UP(1,L11),UP(1,L11),DUL(1,L1),DUL(1,L1) &
                  ,DUPUL1(L1,JATOM,is),JATOM)
        DUPUL1(L1,JATOM,is)=TWO*DUPUL1(L1,JATOM,is)
        CALL RINT13(.FALSE.,UL(1,L1),UL(1,L1),DU(1,L11),DU(1,L11) &
                  ,DULU2(L1,JATOM,is),JATOM)
        DULU2(L1,JATOM,is)=TWO*DULU2(L1,JATOM,is)
        CALL RINT13(.FALSE.,UL(1,L1),UL(1,L1),DUP(1,L11),DUP(1,L11) &
                  ,DULUP2(L1,JATOM,is),JATOM)
        DULUP2(L1,JATOM,is)=TWO*DULUP2(L1,JATOM,is)
        CALL RINT13(.FALSE.,UA,UA,UAB,UAB &
                  ,RUUL(L1,JATOM,is),JATOM)
        CALL RINT13(.FALSE.,UB,UB,UAB,UAB &
                  ,RUPUL(L1,JATOM,is),JATOM)
      END IF
      IF (L+1.le.LOMAX) THEN
      IF (.not.loor(L+1)) THEN
        DULU1(L1,JATOM,is)  = 0
        DULUP1(L1,JATOM,is) = 0
        DUUL2(L1,JATOM,is)  = 0
        DUPUL2(L1,JATOM,is) = 0
        DULUL1(L1,JATOM,is) = 0
        DULUL2(L1,JATOM,is) = 0
        RULU(L1,JATOM,is)   = 0
        RULUP(L1,JATOM,is)  = 0
        RULUL(L1,JATOM,is)  = 0
      ELSE
        DO 506 I=1,IMAX
              SR=SQRT(R(I))
           UA(I)=U(I,L1)/SR
           UB(I)=UP(I,L1)/SR
          UAB(I)=UL(I,L11)/SR
          UBA(I)=UL(I,L1)/SR
  506   CONTINUE
!ad
!ad..........DULU1..DULUP1..DUUL2..DUPUL2..DULUP2..DULUL1..DULUL2.......
!ad
        CALL RINT13(.FALSE.,UL(1,L11),UL(1,L11),DU(1,L1),DU(1,L1) &
                  ,DULU1(L1,JATOM,is),JATOM)
        DULU1(L1,JATOM,is)=TWO*DULU1(L1,JATOM,is)
        CALL RINT13(.FALSE.,UL(1,L11),UL(1,L11),DUP(1,L1),DUP(1,L1) &
                  ,DULUP1(L1,JATOM,is),JATOM)
        DULUP1(L1,JATOM,is)=TWO*DULUP1(L1,JATOM,is)
        CALL RINT13(.FALSE.,U(1,L1),U(1,L1),DUL(1,L11),DUL(1,L11) &
                  ,DUUL2(L1,JATOM,is),JATOM)
        DUUL2(L1,JATOM,is)=TWO*DUUL2(L1,JATOM,is)
        CALL RINT13(.FALSE.,UP(1,L1),UP(1,L1),DUL(1,L11),DUL(1,L11) &
                  ,DUPUL2(L1,JATOM,is),JATOM)
        DUPUL2(L1,JATOM,is)=TWO*DUPUL2(L1,JATOM,is)
        CALL RINT13(.FALSE.,UL(1,L11),UL(1,L11),DUL(1,L1),DUL(1,L1) &
                  ,DULUL1(L1,JATOM,is),JATOM)
        DULUL1(L1,JATOM,is)=TWO*DULUL1(L1,JATOM,is)
        CALL RINT13(.FALSE.,UL(1,L1),UL(1,L1),DUL(1,L11),DUL(1,L11) &
                  ,DULUL2(L1,JATOM,is),JATOM)
        DULUL2(L1,JATOM,is)=TWO*DULUL2(L1,JATOM,is)
!ad
!ad...........................RULU..RULUP..RULUL........................
!ad
        CALL RINT13(.FALSE.,UAB,UAB,UA,UA &
                  ,RULU(L1,JATOM,is),JATOM)
        CALL RINT13(.FALSE.,UAB,UAB,UB,UB &
                  ,RULUP(L1,JATOM,is),JATOM)
        CALL RINT13(.FALSE.,UAB,UAB,UBA,UBA &
                  ,RULUL(L1,JATOM,is),JATOM)
      END IF
      END IF

  503 CONTINUE

  
!
!....... END CALCULATION OF RADIAL INTEGRALS ...............
      RETURN 
      END
