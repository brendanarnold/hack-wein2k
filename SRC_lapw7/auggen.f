      SUBROUTINE AUGGEN(REL,NAT,WHPSI)
      use struct
      use radgrd
      use lolog
      use loabc
      use atspdt
      use radfu
!     last changes: 19.07.00 ba (radfu2 stuff)
!                   15.08.00 ub (merging and pollishing)
!                   29.08.00 ub (pollishing)
!                   01.11.00 ub (rescaling of small rel. comp.)
!
! augmentation functions psi_a,lm(r,E) := u_a,l(|r|,E) Y_lm(r/|r|)
!
! This program requires LOMAX <= LMAX7 !
! ----------------------------------------------------------------
! Input:
! REL    --  .TRUE. for relativistic augmentation
! NAT    --  number of groups of symmetry equivalent atoms
!
! COMMON /STRUCT/
! ZZ     -- the nuclear charge of the atoms
! RMT    -- the muffin tin radius of the atoms
!
! Output:
! COMMON /ATSPDT/ 
! P  (l,a) :           u_a,l(Rmt_a,E_a,l)
! PE (l,a) :      d/dE u_a,l(Rmt_a,E_a,l)
! DP (l,a) :      d/dr u_a,l(Rmt_a,E_a,l)
! DPE(l,a) : d/dr d/dE u_a,l(Rmt_a,E_a,l)
!
! COMMON /RADFU/ (large or small relativistic component)
! RRAD(r,l,a) :      u(s)_a,l (r,E _a,l)
! RADE(r,l,a) : d/dE u(s)_a,l (r,E _a,l)
! RADL(r,l,a) :      u(s)_a,l (r,E'_a,l)
!
! COMMON /LOABC/ 
! ALO(l,a) : A_a,l
! BLO(l,a) : B_a,l
! CLO(l,a) : C_a,l
!
! COMMON /LOLOG/
! NLO       : total number of local orbitals of atom type JATOM
! LOOR(l,a) : .TRUE. if there is a local orbital for atom type a
!                    and angular momentum l
! LORB(a)   : .TRUE. if there are local orbitals for atom type a
!
! Data Passing: (from OUTWIN to AUGGEN)
! COMMON /WORK/
! A(r) : r * u _l(r,E)     for the current potential and energy
! B(r) : r * us_l(r,E) * c the corresponding 2nd rel. component
!
! c = 274.074 au (Rydberg) which is set to 10^+11 au (Rydberg)
! for non-relativistic calculations (in accordance to RINT13)
!
! ----------------------------------------------------------------
! SPIN-ORBITLESS RELATIVISTIC EQUATIONS USED [ L IS 7 (LMAX) ]
! REQUIRED SUBROUTINES ARE OUTWIN AND RINT13
! ----------------------------------------------------------------
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
!
      CHARACTER WHPSI*5
      LOGICAL   REL, LARGE
!
      DIMENSION AE(NRAD),BE(NRAD),VR(NRAD), &
                E(0:LMAX7),ELO(0:LOMAX),PEI(0:LMAX7)
      DIMENSION RRAD1(NRAD,0:LMAX7),RRAD2(NRAD,0:LMAX7), &
                RADE1(NRAD,0:LMAX7),RADE2(NRAD,0:LMAX7), &
                RADL1(NRAD,0:LOMAX),RADL2(NRAD,0:LOMAX)
!
!      COMMON /STRUCT/ POS(3,NDIF),ZZ(NATO),RMT(NATO),MULT(NATO),
!     *                IATNR(NDIF),ANAME(NATO)
!      COMMON /RADGRD/ RM(NRAD,NATO),RNOT(NATO),DX(NATO),JRI(NATO)
!
!      COMMON /LOABC/  ALO(0:LOMAX,NATO),BLO(0:LOMAX,NATO),
!     *                CLO(0:LOMAX,NATO) 
!      COMMON /ATSPDT/ P (0:LMAX7,NATO),DP (0:LMAX7,NATO),
!     *                PE(0:LMAX7,NATO),DPE(0:LMAX7,NATO)
!      COMMON /RADFU/  RRAD(NRAD,0:LMAX7,NATO),RADE(NRAD,0:LMAX7,NATO), &
!                      RADL(NRAD,0:LOMAX,NATO)
      COMMON /WORK1/   A(NRAD),B(NRAD)
!      COMMON /LOLOG/  NLO,LOOR(0:LOMAX,NATO),LORB(NATO)
!---------------------------------------------------------------------
      CFAC = 1.0D0 / 274.074D0
      IF(.NOT.REL) CFAC = 1.0D-11
      LARGE = WHPSI .EQ. 'LARGE'
      IF(LARGE) CFAC = 1.0D0
!
!     << skip header of *.vsp file >>
      READ(18,1000)

      WRITE(6,2000)
      NLO = 0
      DO 10 JATOM=1,NAT
        IMAX = JRI(JATOM)
        ZNUC1 = Znuc (JATOM)
        RMT2 = RMT(JATOM)*RMT(JATOM)
        WRITE(6,2010)JATOM,RMT(JATOM)
!
!       << read atomic potential r * V(r) and convert into Rydberg >>
        READ(18,1010)
        READ(18,1020)(VR(J),J=1,IMAX)
        READ(18,1030)
        DO 20 J=1,IMAX
          VR(J) = 0.5D0 * VR(J)
   20   CONTINUE
!
!       << read augmentation energies and check for local orbitals >>
        READ(10) E
        READ(10) ELO
        LORB(JATOM) = .FALSE.
        DO 30 L=0,LOMAX
          IF(ELO(L).LT.995D0)THEN
            LOOR(L,JATOM) = .TRUE.
            LORB(JATOM) = .TRUE.
            NLO = NLO + (2*L+1)*MULT(JATOM)
          ELSE
            LOOR(L,JATOM) = .FALSE.
          ENDIF
   30   CONTINUE
!
!       << regular augmentation functions >>
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        WRITE(6,2020)
        DELE=2.0D-3
        DELEI=0.25D0/DELE
        DO 40 l=0,LMAX7
            if(e(l).gt.150.d0) then
              write(6,*) ' APW+lo basis not supported in LAPW7'
              write(6,*) ' Rerun lapw1 with LAPW-basis only'
              stop ' Rerun LAPW1 with LAPW-basis only '
            endif
          FL=L
          EI=E(L)/2.0d0
!
!         << compute d/dE u_l(r,E) by finite differences >>
          E1=EI-DELE
          CALL OUTWIN(REL,VR,RM(1,JATOM),DX(JATOM),IMAX,E1, &
                      FL,UVB,DUVB,NODES,ZNUC1)
          CALL RINT13(REL,A,B,A,B,RNORM,JATOM)
          RNORM = 1.0D0/SQRT(RNORM)
          DO 50 M=1,IMAX
            AE(M) = RNORM * A(M)
            BE(M) = RNORM * B(M)
   50     CONTINUE
          UVB  = RNORM * UVB
          DUVB = RNORM * DUVB

          E1=EI+DELE
          CALL OUTWIN(REL,VR,RM(1,JATOM),DX(JATOM),IMAX,E1, &
                      FL,UVE,DUVE,NODES,ZNUC1)
          CALL RINT13(REL,A,B,A,B,RNORM,JATOM)
          RNORM = 1.0D0/SQRT(RNORM)
          DO 60 M=1,IMAX
            AE(M) = DELEI*(RNORM*A(M)-AE(M))
            BE(M) = DELEI*(RNORM*B(M)-BE(M))
   60     CONTINUE
          UVE  = DELEI*(RNORM*UVE -UVB )
          DUVE = DELEI*(RNORM*DUVE-DUVB)
!
!         << now compute u_l(r,E) >>
          CALL OUTWIN(REL,VR,RM(1,JATOM),DX(JATOM),IMAX,EI, &
                      FL,UV,DUV,NODES,ZNUC1)
          CALL RINT13(REL,A,B,A,B,RNORM,JATOM)
          RNORM = 1.0D0/SQRT(RNORM)
          DO 70 M=1,IMAX
            A(M) = RNORM*A(M)
            B(M) = RNORM*B(M)
   70     CONTINUE
          P (L,JATOM) = RNORM*UV
          DP(L,JATOM) = RNORM*DUV
!
!         << insure orthogonality of d/dE u_l(r,E) on u_l(r,E) >>
          CALL RINT13(REL,A,B,AE,BE,CROSS,JATOM)
          DO 80 M=1,IMAX
            AE(M) = (AE(M)-CROSS*A(M))
            BE(M) = (BE(M)-CROSS*B(M))
   80     CONTINUE
          PE (L,JATOM) = UVE -CROSS*P (L,JATOM)
          DPE(L,JATOM) = DUVE-CROSS*DP(L,JATOM)
!
!         << store augmentation functions >>
          DO 90 I=1,IMAX
            RRAD1(I,L) = A(I)
            RADE1(I,L) = AE(I)
            RRAD2(I,L) = B(I)
            RADE2(I,L) = BE(I)
   90     CONTINUE
! 
!         << evaluate non-trivial radial integrals >>
          CALL RINT13(REL,AE,BE,AE,BE,PEI(L),JATOM)

          WRITE(6,2030) L,E(L),P(L,JATOM),DP(L,JATOM), &
                        PE(L,JATOM),DPE(L,JATOM),CROSS,PEI(L)
   40   CONTINUE
!
!       << local orbitals >>
!       ~~~~~~~~~~~~~~~~~~~~
        IF(.NOT.LORB(JATOM))GOTO 15

        WRITE(6,2040)
        DELE=2.0D-3
        DELEI=0.25D0/DELE
        DO 100 L=0,LOMAX
        IF(.NOT.LOOR(L,JATOM)) GOTO 100
          FL=L
          EI=ELO(L)/2.d0
!
!         << compute u_l(r,E') >>
          CALL OUTWIN(REL,VR,RM(1,JATOM),DX(JATOM),IMAX,EI, &
                      FL,UV,DUV,NODES,ZNUC1)
          CALL RINT13(REL,A,B,A,B,RNORM,JATOM)
          RNORM = 1.0d0/SQRT(RNORM)
          DO 110 M=1,IMAX
            A(M) = RNORM*A(M)
            B(M) = RNORM*B(M)
  110     CONTINUE
          PLO  = RNORM*UV
          DPLO = RNORM*DUV
!
!         << store augmentation functions >>
          DO 120 I=1,IMAX
            RADL1(I,L) = A(I)
            RADL2(I,L) = B(I)
  120     CONTINUE
! 
!         << evaluate non-trivial radial integrals >>
          CALL RINT13(REL,RRAD1(1,L),RRAD2(1,L),A,B,PI12LO,JATOM)
          CALL RINT13(REL,RADE1(1,L),RADE2(1,L),A,B,PE12LO,JATOM)

          WRITE(6,2050)L,ELO(L),PLO,DPLO,PI12LO,PE12LO
! 
!         << the augmentation coefficients A_a,l, B_a,l and C_a,l >>
! ------------------------------------------------------------------
! A_a,l = C_a,l * A'  and  B_a,l = C_a,l * B'  with
! 
! A' := [ d/dr d/dE u_a,l(Rmt,E_a,l)      u_a,l(Rmt,E'_a,l) -
!              d/dE u_a,l(Rmt,E_a,l) d/dr u_a,l(Rmt,E'_a,l) ] * Rmt^2
!
! B' := [           u_a,l(Rmt,E_a,l) d/dr u_a,l(Rmt,E'_a,l) -
!              d/dr u_a,l(Rmt,E_a,l)      u_a,l(Rmt,E'_a,l) ] * Rmt^2
! and
!
! C_a,l = [ 1 + A'^2
!             + B'^2 < d/dE u_l(r,E_l) | d/dE u_l(r,E _l) >
!             + 2 A' <      u_l(r,E_l) |      u_l(r,E'_l) >
!             + 2 B' < d/dE u_l(r,E_l) |      u_l(r,E'_l) > ]^(-1/2)
! ------------------------------------------------------------------
          XAC = ( PLO *DPE(L,JATOM) - DPLO*PE (L,JATOM) ) * RMT2
          XBC = ( DPLO*P  (L,JATOM) - PLO *DP (L,JATOM) ) * RMT2
          XCC = XAC*( XAC        + 2.0D0*PI12LO ) &
              + XBC*( XBC*PEI(L) + 2.0D0*PE12LO ) + 1.0D0
          CLO(L,JATOM) = 1.0D0/MAX( SQRT(XCC) , 0.005D0 )
          ALO(L,JATOM) = XAC * CLO(L,JATOM)
          BLO(L,JATOM) = XBC * CLO(L,JATOM)
  100   CONTINUE
!
        WRITE(6,2060)
        DO 130 L=0,LOMAX
        IF(.NOT.LOOR(L,JATOM)) GOTO 130
          WRITE(6,2070)L,ALO(L,JATOM),BLO(L,JATOM),CLO(L,JATOM)
  130   CONTINUE

   15   CONTINUE
!   
!       << re-scale radial augmentation functions  >>
!       << with VR just being a working array here >>
        DO 135 I=1,IMAX
          VR(I) = CFAC / RM(I,JATOM)
  135   CONTINUE
        DO 140 L=0,LMAX7
           IF(LARGE)THEN
              DO 145 I=1,IMAX
                 RRAD(I,L,JATOM) = RRAD1(I,L) * VR(I)
                 RADE(I,L,JATOM) = RADE1(I,L) * VR(I)
  145         CONTINUE
           ELSE
              DO 146 I=1,IMAX
                 RRAD(I,L,JATOM) = RRAD2(I,L) * VR(I)
                 RADE(I,L,JATOM) = RADE2(I,L) * VR(I)
  146         CONTINUE
           ENDIF
  140   CONTINUE
        IF(LORB(JATOM))THEN
          DO 150 L=0,LOMAX
            IF(LOOR(L,JATOM))THEN
               IF(LARGE)THEN
                  DO 155 I=1,IMAX
                     RADL(I,L,JATOM) = RADL1(I,L) * VR(I)
  155             CONTINUE
               ELSE
                  DO 156 I=1,IMAX
                     RADL(I,L,JATOM) = RADL2(I,L) * VR(I)
  156             CONTINUE
               ENDIF
            ENDIF
  150     CONTINUE
        ENDIF

   10 CONTINUE
      RETURN
!
 1000 FORMAT(//)
 1010 FORMAT(/////)
 1020 FORMAT(3X,4E19.12)
 1030 FORMAT(/////)
!
 2000 FORMAT(/' AUGMENTATION FUNCTIONS' &
             /' ----------------------' &
             /' psi_lm(r) = u_l(|r|) * Y_lm(r/|r|)')
 2010 FORMAT(/' augmentation functions for atom',i4,' (Rmt=',F6.3,')')
 2020 FORMAT(/' regular augmentation functions at r = Rmt' &
             /'  L    E(L)      u(r,E)     u''(r,E)   dE u(r,E)' &
             ,'  dE u''(r,E)  <u|dE u>  |dE u|^2')
 2030 FORMAT(I3,F8.3,1P,4E12.4,2E10.2)
 2040 FORMAT(/' local orbitals at r = Rmt' &
             /'  L   E''(L)     u(r,E'')    u''(r,E'')' &
             ,'  <u(E'')|u(E)>  <u(E'')|dE u(E)>')
 2050 FORMAT(I3,F8.3,1P,2E12.4,E14.2,E17.2)
 2060 FORMAT(/' local orbital coefficients:' &
             ,' A_l u_l(r,E) + B_l d/dE u_l(r,E) + C_l u_l(r,E'')' &
             /'  L   A_L           B_L           C_L')
 2070 FORMAT(I3,1P,3E14.6)
!
      END
