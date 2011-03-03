      SUBROUTINE ATPAR(NT,REL,NAT,LNSMAX,ZZ,force)
!
      use loabc, only : ALO, BLO, CLO, ELO, PLO, DPLO, PELO, DPELO, PEILO, PI12LO, PE12LO, init_loabc
      use radial, only: A,B,AP,BP,AE,BE,VLM,A1,B1,AE1,BE1,A1LO,B1LO, AE1LO, BE1LO, init_radial, end_radial
      use lolog, only : nlo, loor, lapw, ilo, init_lolog
      use atspdt, only: INS, NL, LMMAX, LQIND, LM, LQNS, DP, DPE, E, P, PE, PEI, VNS1, VNS2, VNS3, GFAC, init_atspdt
      use char, only  : LATTIC, NAME
      use loint, only : VNS1LO, VNS2LO, VNS3LO, init_loint
      use potnlc, only: JRJ, DH, RO, VR
      use rotmat, only: ROTIJ, ROTLOC
      use struk, only : POS, ALAT, ALPHA, RMT, V, PIA, VI, IATNR, MULT
      use parallel, only: MYID
      IMPLICIT NONE
      INCLUDE 'param.inc'
!
!        Arguments
!
      INTEGER            LNSMAX, NAT, NT
      DOUBLE PRECISION   ZZ(NAT), TRAPLEVEL
      parameter (traplevel=1D-3)
      LOGICAL            REL, DOTRAP
!
!     ..................................................................
!
!        ATPAR calculates the solutions u(l) of the radial Schroedinger
!        Equation, the energy derivatives ue(l), the radial non muffin-
!        tin integrals uvu, uevu, and the Gaunt-coefficients.
!        Spin-orbitless relativistic equations are used.
!
!     ..................................................................
!
!        Local Scalars
!
      INTEGER            I, ICOUNT, IDUMMY, INDEX, J, JATOM, JC, JCOL
      INTEGER            JRI, JROW, L, L0, L0BEFO, LALT, LATOM, LL
      INTEGER            LLBEFO, LLL, LLMM, LM1, LM2, LMBEFO, LMX, LP
      INTEGER            LPBEFO, LQX, M, MINU, MM, MP, NLR, NODEL
      INTEGER            NODES, NODEU, IAPW, JLO
      DOUBLE PRECISION   AINT, CFEIN1, CFEIN2, CROSS, DE, DELE, DELEI
      DOUBLE PRECISION   DUV, DUVB, DUVE, DX, E1, EI, FL, OVLP, RNOT
      DOUBLE PRECISION   TRX, TRY, UV, UVB, UVE
      COMPLEX*16         IMAG, IMAG1
      CHARACTER*4        EMAIN
      CHARACTER*67       ERRMSG
!fb
!     variables for [yu91]-force formalism
      logical     force
!fe
!
!        External Functions
!
      DOUBLE PRECISION   GAUNT1
      EXTERNAL           GAUNT1
!
!        External Subroutines
!
      EXTERNAL           ABC,CBCOMB,FORFHS,OUTERR,OUTWIN,RINT13,SELECT
!
!        Intrinsic Functions
!
      INTRINSIC          ABS, DBLE, MOD, SQRT
!
!        Data statements
!
      DATA DELE /2.0D-3/
      DATA IMAG /(0.0D+0,1.0D+0)/
!
!fb
!      force=.false.
!fe
      call init_loabc(LOMAX, NLOAT, NAT)
      call init_radial(NRAD, LOMAX, NSLMAX, NLOAT, LMMX)
      call init_lolog(LOMAX, NAT, LMAX)
      call init_atspdt(NAT, LMMX, NGAU, LMAX, NSLMAX)
      call init_loint(LOMAX+1, LMMX, NSLMAX, NAT)
!        initialize constants
!        CFEIN1 and CFEIN2 are used in RINT13 to calculate the
!        integrals of radialfunctions and are derivates of the
!        Feinstruktur-constant (in Hartrees)
!
      IF (REL) THEN
         CFEIN1 = 1.0D+0
         CFEIN2 = 1.d0/137.0359895d0**2
      ELSE
         CFEIN1 = 1.0D+0
         CFEIN2 = 1.0D-22
      ENDIF
!
!        read total spherical potential V(0,0) of TAPE18=VSP
!        norm of V(0,0)=V00(R)*R/SQRT(4.0D+0*PI)
!
      READ (18,5070)
      DO 20 JATOM = 1, NAT
         READ (18,5010)
         READ (18,5020) IDUMMY
         READ (18,5060)
         READ (18,5040) (VR(J,JATOM),J=1,JRJ(JATOM))
         READ (18,5060)
         READ (18,5050)
         DO J = 1, JRJ(JATOM)
            VR(J,JATOM) = VR(J,JATOM)/2.0D+0
         ENDDO
   20 CONTINUE
!
      READ (19,5050,END=30,IOSTAT=INS)
      REWIND 19
      READ (19,5070)
!fb
!     force-output only 
!     if non-spherical potential exists
      IF(myid.eq.0) then
      if(force)  write(71,'( &
          &" non-spherical Matrixelements x Gaunt-Coefficient (d=du/dE):",/, &
          &" l0 ll lp m0 mm mp     uVu dVd",/, &
          &"                       uVd dVu")')
      endif
!fe
   30 CONTINUE
!
!        start loop over atoms in unitcell
!
      LATOM = 0
      INDEX = 0
      NLO = 0
      DO 340 JATOM = 1, NAT
         do l=0,lmax-1
            lapw(l,jatom)=.true.
         enddo
         RNOT = RO(JATOM)
         LATOM = LATOM + 1
         JRI = JRJ(JATOM)
         DX = DH(JATOM)
         iapw=0
         READ (5,*,err=41) EI, NLR,iapw
         if(iapw.eq.1) then
            do l=0,lmax-1
               lapw(l,jatom)=.false.
            enddo
         endif
 41      continue
!         WRITE (6,6040) IATNR(JATOM), (POS(I,LATOM),I=1,3),
!     ,        MULT(JATOM)
         LATOM = LATOM + MULT(JATOM) - 1
!         WRITE (6,6080) NAME(JATOM)
         DO 40 JCOL = 1, 3
!            WRITE (6,6090) (ROTLOC(JROW,JCOL,JATOM),JROW=1,3)
 40      CONTINUE
         DO 60 M = 1, MULT(JATOM)
            INDEX = INDEX + 1
!            WRITE (6,6100) M, (POS(JC,INDEX),JC=1,3)
            DO 50 JCOL = 1, 3
!               WRITE (6,6110) (ROTIJ(JROW,JCOL,INDEX),JROW=1,3)
 50         CONTINUE
 60      CONTINUE
         if(iapw.eq.1) then
!            WRITE (6,6000) NAME(JATOM), EI,'  APW'
            IF(myid.eq.0) WRITE (21,6000) NAME(JATOM), EI,'  APW'
         else
!            WRITE (6,6000) NAME(JATOM), EI,' LAPW'
            IF(myid.eq.0) WRITE (21,6000) NAME(JATOM), EI,' LAPW'
         endif
         DO 70 L = 1, LMAX
            IF ((L-1) .LE. (LOMAX)) THEN
               LOOR(L-1,JATOM) = .FALSE.
               ilo(l-1,jatom)=0
               do i=1,nloat
                  ELO(L-1,i,JATOM) = 997.0D+0
                  if(lapw(l-1,jatom)) ELO(L-1,i,JATOM) = 999.0D+0
               enddo
            ENDIF
            E(L,JATOM) = EI
 70      CONTINUE
         L = -1
         DO 80 J = 1, NLR
            LALT = L
            iapw=0
            READ (5,5000,err=42) L, EI, DE, EMAIN, iapw
 42         continue
            IF (DE .GT. 1.0D-5) THEN
               CALL SELECT(L,EI,DE,EMAIN,REL, &
                    VR(1,JATOM),RNOT,DX,JRI,ZZ(JATOM))
            ELSE
               IF(myid.eq.0) WRITE (6,6010) L, EI
               IF(myid.eq.0) WRITE (21,6010) L, EI
            ENDIF
            IF (L .EQ. LALT) THEN
               ilo(l,jatom)=ilo(l,jatom)+1
               ELO(L,ilo(l,jatom),JATOM) = EI
               NLO = NLO + ((2*L+1))*MULT(JATOM)
               IF(myid.eq.0) write(6,*)  '            LOCAL ORBITAL'
               IF(myid.eq.0) write(21,*) '            LOCAL ORBITAL'
               loor(l,jatom)=.true.
            ELSE
               if(iapw.eq.1) then
                  lapw(l,jatom)=.false.
                  NLO = NLO + ((2*L+1))*MULT(JATOM)
                  IF(myid.eq.0) write(6,*) '            APW+lo'
                  IF(myid.eq.0) write(21,*) '            APW+lo'
                  ilo(l,jatom)=1
                  elo(l,1,jatom)=EI
                  elo(l,2,jatom)=997.0D+0
               else
                  IF(myid.eq.0) write(6,*) '            LAPW'
                  IF(myid.eq.0) write(21,*) '            LAPW'
               endif
               E(L+1,JATOM) = EI
            ENDIF
 80      CONTINUE
         
         IF(myid.eq.0) WRITE (6,6020)
         DELEI = .250D+0/DELE
         DO 150 J = 1, NT
!  Experimental Trap location
            DOTRAP=.true.
1500        CONTINUE
            L = J - 1
            FL = L
            EI = E(J,JATOM)/2.D0
!
!        calculate energy-derivative by finite difference
!        DELE is the up and downward energy-shift in Hartrees
!
            E1 = EI - DELE
            CALL OUTWIN(REL,VR(1,JATOM),RNOT,DX,JRI,E1,FL,UVB,DUVB, &
                 NODEL,ZZ(JATOM))
            CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM)
            TRX = 1.0D+0/SQRT(OVLP)
            DO 100 M = 1, JRI
               AE(M) = TRX*A(M)
               BE(M) = TRX*B(M)
 100        CONTINUE
            UVB = TRX*UVB
            DUVB = TRX*DUVB
            E1 = EI + DELE
            CALL OUTWIN(REL,VR(1,JATOM),RNOT,DX,JRI,E1,FL,UVE,DUVE, &
                 NODEU,ZZ(JATOM))
            CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM)
            TRX = 1.0D+0/SQRT(OVLP)
            UVE = DELEI*(TRX*UVE - UVB)
            DUVE = DELEI*(TRX*DUVE - DUVB)
            DO 110 M = 1, JRI
               AE(M) = DELEI*(TRX*A(M) - AE(M))
               BE(M) = DELEI*(TRX*B(M) - BE(M))
 110        CONTINUE
!     
!     calculate function at EI
!
            CALL OUTWIN(REL,VR(1,JATOM),RNOT,DX,JRI,EI,FL,UV,DUV,NODES &
                 ,ZZ(JATOM))
            CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM)
            TRX = 1.0D+0/SQRT(OVLP)
            P(J,JATOM) = TRX*UV
            DP(J,JATOM) = TRX*DUV
            DO 120 M = 1, JRI
               A(M) = TRX*A(M)
               B(M) = TRX*B(M)
 120        CONTINUE
!
!     insure orthogonalization
!
            CALL RINT13(CFEIN1,CFEIN2,A,B,AE,BE,CROSS,JATOM)
            TRY = -CROSS
            IF (TRY .LT. -0.05D+0) WRITE (6,6070) L, TRY, OVLP
            DO 130 M = 1, JRI
               AE(M) = AE(M) + TRY*A(M)
               BE(M) = BE(M) + TRY*B(M)
 130        CONTINUE
            PE(J,JATOM) = UVE + TRY*P(J,JATOM)
            DPE(J,JATOM) = DUVE + TRY*DP(J,JATOM)
            CALL RINT13(CFEIN1,CFEIN2,AE,BE,AE,BE,PEI(J,JATOM),JATOM)
            IF (J .LE. NSLMAX) THEN
! gm remove again
!               write(6,*) 'atom,j',jatom,j
               DO 140 M = 1, JRI
                  A1(M,J) = A(M)
                  B1(M,J) = B(M)
                  AE1(M,J) = AE(M)
                  BE1(M,J) = BE(M)
!                  rpip=rnot*(exp(dx*(m-1.d0)))
!                  write(6,*) rpip,a1(m,j)/rpip,ae1(m,j)/rpip
 140           CONTINUE
            ENDIF

! Trap P(J,JATOM) too small
        IF((Abs(P(J,JATOM)) .lt. traplevel).and.DOTRAP)then
                DOTRAP=.false.
! An alternative might be to set the level to LAPW....?
                E(J,JATOM)=E(J,JATOM)-0.05
                IF(myid.eq.0) write(21,1501)E(J,JATOM)
                IF(myid.eq.0) write(6,1501)E(J,JATOM)
1501            format(':WARN  : P(J,JATOM) almost zero. Shifted Energy by -0.05 down to ',F10.4)
                goto 1500
        endif
            IF(myid.eq.0) WRITE (6,6030) L, P(J,JATOM), DP(J,JATOM), PE(J,JATOM), &
                 DPE(J,JATOM), PEI(J,JATOM), NODEL, NODES, &
                 NODEU
 150     CONTINUE
!     
!     and now for local orbitals
!     
         IF(myid.eq.0) write(6,6021)
         DO 220 L = 0, LOMAX
            IF (LOOR(L,JATOM)) THEN
               J = L + 1
               FL = L
               EI = ELO(L,ilo(l,jatom),JATOM)/2.D0
!     
!     calculate function at EI
!     
               CALL OUTWIN(REL,VR(1,JATOM),RNOT,DX,JRI,EI,FL,UV,DUV, &
                    NODES,ZZ(JATOM))
               CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM)
               TRX = 1.0D+0/SQRT(OVLP)
               PLO(L,JATOM) = TRX*UV
               DPLO(L,JATOM) = TRX*DUV
               DO 180 M = 1, JRI
                  A(M) = TRX*A(M)
                  B(M) = TRX*B(M)
 180           CONTINUE

               CALL RINT13(CFEIN1,CFEIN2,A1(1,J),B1(1,J),A,B, &
                    PI12LO(L,JATOM),JATOM)
               CALL RINT13(CFEIN1,CFEIN2,AE1(1,J),BE1(1,J),A,B, &
                    PE12LO(L,JATOM),JATOM)
               IF (J .LE. NSLMAX) THEN
                  DO 200 M = 1, JRI
                     A1LO(M,L) = A(M)
                     B1LO(M,L) = B(M)
 200              CONTINUE
               ENDIF
               IF(myid.eq.0) WRITE (6,6031) L, PLO(L,JATOM), DPLO(L,JATOM), &
                   PI12LO(L,JATOM),PE12LO(L,JATOM), &
                   NODEL, NODES, NODEU
            endif
 220     CONTINUE
         do l=0,lomax
            do jlo=1,ilo(l,jatom)
               CALL ABC(JATOM,L,jlo,lapw(l,jatom))
            enddo
         enddo

         IF (INS .EQ. 0) THEN
! fb
            IF(myid.eq.0) then
            if (force) write(71,'(i3," .Atom ")') jatom
            endif
!fe
!
!        read total nonspherical potential of TAPE19=VNS
!        norm of VLM=VLM*1
!
            READ (19,5010)
            READ (19,5020) LMMAX(JATOM)
            LM2 = 0
            DO 230 LM1 = 1, LMMAX(JATOM)
               READ (19,5030) LLL, MM
               LL = ABS(LLL)
               LM2 = LM2 + 1
               IF (LM2 .GT. LMMX) THEN
!
!        Error: LMMX too small - stop execution
!
                  GOTO 900
               ENDIF
               READ (19,5040) (VLM(J,LM2),J=1,JRJ(JATOM))
!     
               IF (LL .GT. LNSMAX*2) THEN
                  LM2 = LM2 - 1
               ELSE
                  LM(1,LM2,JATOM) = LLL
                  LM(2,LM2,JATOM) = MM
               ENDIF
               READ (19,5060)
 230        CONTINUE
            READ (19,5050)
!
!        calculate possible nonspherical contributions to H
!
            IF(myid.eq.0) WRITE (6,*) ' POSSIBLE NONSPHERICAL CONTRIBUTIONS TO H: '
            IF(myid.eq.0) WRITE (6,*) '(L0,LL,LP,M, MM,MP, GNT,  U V U,    UE V UE,', &
                 '   U V UE     )'
            LQIND(JATOM) = 0
            LLMM = LM2
!     
            DO 300 L0 = 0, LNSMAX
               DO 290 LP = 0, LNSMAX
                  DO 280 LMX = 1, LLMM
                     LL = ABS(LM(1,LMX,JATOM))
                     MM = LM(2,LMX,JATOM)
!     
!     which gaunts are not zero ?
!
                     IF (MOD((L0+LP+LL),2) .EQ. 1) GOTO 280
                     IF ((L0+LP-LL) .LT. 0) GOTO 280
                     IF ((L0-LP+LL) .LT. 0) GOTO 280
                     IF ((LP-L0+LL) .LT. 0) GOTO 280
!     
 240                 CONTINUE
                     DO 270 M = -L0, L0
                        DO 260 MP = -LP, LP
                           IF (-M+MM+MP) 260, 250, 260
 250                       CONTINUE
                           LQIND(JATOM) = LQIND(JATOM) + 1
                           IF (LQIND(JATOM) .GT. NGAU) THEN
!
!     Error: more then NGAU gaunt-factors - stop execution
!
                              GOTO 910
                           ENDIF
                           GFAC(LQIND(JATOM),JATOM) = &
                                GAUNT1(L0,LL,LP,M,MM,MP)
!                           GFAC(LQIND(JATOM),JATOM) =
!     &                        GAUNT1(L0,LL,LP,M,MM,MP)*IMAG**(LP-L0)
                           IF (IATNR(JATOM) .GT. 0) THEN
			      IF (LL.EQ.3.OR.LL.EQ.7.OR.LL.EQ.9) THEN
				 IF(MM.LT.0) THEN
                                    GFAC(LQIND(JATOM),JATOM) = &
                                      -GFAC(LQIND(JATOM),JATOM)*IMAG
				 ELSE
                                    GFAC(LQIND(JATOM),JATOM) = &
                                      GFAC(LQIND(JATOM),JATOM)*IMAG
				 ENDIF
                              ENDIF
                           ELSE
                              IF (MM .NE. 0) THEN
                                 MINU = 1
                                 IMAG1 = (1.0D+0,0.0D+0)
                                 IF (LM(1,LMX,JATOM) .LT. 0) THEN
                                    IMAG1 = -IMAG
                                    MINU = -1
                                 ENDIF
                                 IF (MOD(MM,2) .EQ. 1) THEN
                                    IMAG1 = -IMAG1
                                    MINU = -MINU
                                 ENDIF
                                 IF (MM .GT. 0) MINU = 1
                                 GFAC(LQIND(JATOM),JATOM) = &
                                      GFAC(LQIND(JATOM),JATOM)*IMAG1* &
                                      DBLE(MINU)/SQRT(2.0D+0)
                              ENDIF
                           ENDIF
                           LQNS(1,LQIND(JATOM),JATOM) = L0 + 1
                           LQNS(2,LQIND(JATOM),JATOM) = LL
                           LQNS(3,LQIND(JATOM),JATOM) = LP + 1
                           LQNS(4,LQIND(JATOM),JATOM) = LMX
                           LQNS(5,LQIND(JATOM),JATOM) = M
                           LQNS(6,LQIND(JATOM),JATOM) = MP
!                           IF(myid.eq.0) WRITE(6,6050) L0, LL, LP, M, MM, MP,
!     ,                                   GFAC(LQIND(JATOM),JATOM)
 260                    CONTINUE
 270                 CONTINUE
!
                     IF (MM.GT.0) THEN
                        MM = -MM
                        GOTO 240
                     ENDIF
 280              CONTINUE
 290           CONTINUE
 300        CONTINUE
!
            IF ((IATNR(JATOM) .GT. 0) .AND. (LMMAX(JATOM) .GT. 1)) THEN
!
!        combine the radial total potential coefficients according to
!        cubic harmonic
!
               CALL CBCOMB(JRI,LLMM,VLM,LM,JATOM)
            ENDIF
!
!        integrate product of radialfunctions (and derivatives) and
!        total potential coefficients (combinations) from 0 to R-MT
!
            ICOUNT = 0
            L0BEFO = 999
            LLBEFO = 999
            LPBEFO = 999
            LMBEFO = 999
            DO 330 LQX = 1, LQIND(JATOM)
               L0 = LQNS(1,LQX,JATOM) - 1
               LL = LQNS(2,LQX,JATOM)
               LP = LQNS(3,LQX,JATOM) - 1
               LM1 = LQNS(4,LQX,JATOM)
               IF ((L0 .EQ. L0BEFO) .AND. (LL .EQ. LLBEFO) .AND. &
                   (LP .EQ. LPBEFO) .AND. (LM1 .EQ. LMBEFO)) THEN
                  M = LQNS(5,LQX,JATOM)
                  MP = LQNS(6,LQX,JATOM)
!                  WRITE (6,6050) L0, LL, LP, M, LM(2,LM1,JATOM), MP,
!     ,                 GFAC(LQX,JATOM)
               ELSE
!     
                  DO 310 M = 1, JRI
!                     RI = RNOT*EXP(DX*(M-1))
                     A(M) = A1(M,L0+1)*VLM(M,LM1)
                     B(M) = B1(M,L0+1)*VLM(M,LM1)
                     AE(M) = AE1(M,L0+1)*VLM(M,LM1)
                     BE(M) = BE1(M,L0+1)*VLM(M,LM1)
 310              CONTINUE
                  ICOUNT = ICOUNT + 1
!
                  CALL RINT13(CFEIN1,CFEIN2,A,B,A1(1,LP+1),B1(1,LP+1), &
                       AINT,JATOM)
                  VNS1(L0+1,LM1,LP+1,JATOM) = AINT
                  CALL RINT13(CFEIN1,CFEIN2,AE,BE,AE1(1,LP+1), &
                       BE1(1,LP+1),AINT,JATOM)
                  VNS2(L0+1,LM1,LP+1,JATOM) = AINT
                  CALL RINT13(CFEIN1,CFEIN2,A,B,AE1(1,LP+1),BE1(1,LP+1), &
                       AINT,JATOM)
                  VNS3(L0+1,LM1,LP+1,JATOM) = AINT
!
!     integrals needed for local orbitals
!
                  IF (L0 .LE. LOMAX) THEN
                     VNS1LO(L0+1,LM1,LP+1,JATOM) = 0.0D0
                     VNS2LO(L0+1,LM1,LP+1,JATOM) = 0.0D0
                     if(loor(l0,jatom)) THEN
!                     if (lapw(jatom).and.ilo(l0,jatom).eq.1) THEN
                        DO M = 1, JRI
                           A(M) = A1LO(M,L0)*VLM(M,LM1)
                           B(M) = B1LO(M,L0)*VLM(M,LM1)
                        ENDDO
                     else
                        goto 320
                     endif
                     CALL RINT13(CFEIN1,CFEIN2,A,B,A1(1,LP+1), &
                          B1(1,LP+1),AINT,JATOM)
                     VNS1LO(L0+1,LM1,LP+1,JATOM) = AINT
                     CALL RINT13(CFEIN1,CFEIN2,A,B,AE1(1,LP+1), &
                          BE1(1,LP+1),AINT,JATOM)
                     VNS2LO(L0+1,LM1,LP+1,JATOM) = AINT
                     IF (LP .LE. LOMAX) THEN
!     IF (LOOR(LP,JATOM)) THEN
                        AINT = 0.0d0
                        if(loor(lp,jatom)) &
                             CALL RINT13(CFEIN1,CFEIN2,A,B,A1LO(1,LP), &
                             B1LO(1,LP),AINT,JATOM)
                        VNS3LO(L0+1,LM1,LP+1,JATOM) = AINT
                     ENDIF
                     M = LQNS(5,LQX,JATOM)
                     MP = LQNS(6,LQX,JATOM)
!                     WRITE (6,6060) L0, LL, LP, M, LM(2,LM1,JATOM), MP,
!     ,                    GFAC(LQX,JATOM),
!     ,                    VNS1LO(L0+1,LM1,LP+1,JATOM),
!     ,                    VNS2LO(L0+1,LM1,LP+1,JATOM),
!     ,                    VNS3LO(L0+1,LM1,LP+1,JATOM), LM1
                  ENDIF
 320              continue
!
!        end of lo-integrals
!
                  L0BEFO = L0
                  LLBEFO = LL
                  LPBEFO = LP
                  LMBEFO = LM1
                  M = LQNS(5,LQX,JATOM)
                  MP = LQNS(6,LQX,JATOM)
!                  WRITE (6,6050) L0, LL, LP, M, LM(2,LM1,JATOM), MP, &
!                                GFAC(LQX,JATOM),&
!                                VNS1(L0+1,LM1,LP+1,JATOM),&
!                                VNS2(L0+1,LM1,LP+1,JATOM),&
!                                VNS3(L0+1,LM1,LP+1,JATOM), LM1
               ENDIF
  330       CONTINUE
!fb
            if (force) then
                call forfhs(jatom) 
!
!               final line of [yu91]-output
!               final line for each atom = 6 x 0 and 8 x 0.0d0
                IF(myid.eq.0) write(71,'(6i3,4e19.12,/,6(3x),4e19.12)') &
                      (0, i=1,6),(0.d0,i=1,8)
            endif
!fe
            IF(myid.eq.0) WRITE (6,*) '   NUMBER OF RADIAL INTEGRALS FOR ATOM', &
                        JATOM,' = ',ICOUNT
         ENDIF
!
!        end of loop over atoms in unitcell
!
  340 CONTINUE
      
      CALL END_RADIAL
!
      RETURN
!
  900 CONTINUE
      WRITE (ERRMSG,9000) 'LMMX too small', LMMX, LM2
      CALL OUTERR('ATPAR',ERRMSG)
      STOP 'ATPAR - Error'
!
  901 CONTINUE
      WRITE (ERRMSG,9001) 'LOMAX too small', LOMAX, L
      CALL OUTERR('ATPAR',ERRMSG)
      STOP 'ATPAR - Error'
!
  902 CONTINUE
      WRITE (ERRMSG,9002) 'NLOAT too small', NLOAT,l,jatom
      CALL OUTERR('ATPAR',ERRMSG)
      STOP 'ATPAR - Error'
!
  910 CONTINUE
      CALL OUTERR('ATPAR','more than NGAU gaunts')
      WRITE (ERRMSG,9010) '  NGAU,   L0,   LP,   LL,    M,   MP,   MM'
      CALL OUTERR('ATPAR',ERRMSG)
      WRITE (ERRMSG,9020) NGAU, L0, LP, LL, M, MP, MM
      CALL OUTERR('ATPAR',ERRMSG)
      STOP 'ATPAR - Error'
!
 5000 FORMAT (1X,I1,2F10.5,A4,i2)
 5010 FORMAT (3X)
 5020 FORMAT (15X,I3,/,/)
 5030 FORMAT (15X,I3,5X,I2,/)
 5040 FORMAT (3X,4e19.12)
 5050 FORMAT (/,/,/)
 5060 FORMAT (/)
 5070 FORMAT (/,/)
 6000 FORMAT (/,10X,'ATOMIC SPHERE DEPENDENT PARAMETERS FOR ATOM  ',A10, &
              /,10X,'OVERALL ENERGY PARAMETER IS',F10.4 &
              /,10X,'OVERALL BASIS SET ON ATOM IS',A5)
 6010 FORMAT (10X,'E(',I2,2H)=,F10.4)
 6020 FORMAT (/,10X,'POTENTIAL PARAMETERS ',/,11X,'L',5X,'U(R)',10X, &
              'U''(R)',9X,'DU/DE',8X,'DU''/DE',6X,'NORM-U''')
 6021 FORMAT (/,10X,'LOCAL ORBITAL POTENTIAL PARAMETERS',/, &
           11X,'L',5X,'U(R)',10X, &
           'U''(R)',5X,'NORM U1U2',8X,'NORM UE1U2')
 6030 FORMAT (10X,I2,5E14.6,5X,3I2)
 6031 FORMAT (10X,I2,4E14.6,5X,3I2)
 6040 FORMAT (/,/,/,' AT.NR.:',I2,5X,'POSITION: ',3F8.3,5X, &
              'MULTIPLICITY:',I5)
 6050 FORMAT (3X,6I3,5F15.8,I3)
 6060 FORMAT ('lo ',6I3,5F15.8,I3)
 6070 FORMAT (10X,'FOR L=',I2,' CORRECTION=',E14.5,' OVERLAP=',E14.5)
 6080 FORMAT (/,/,3X,'NOT EQUIV ATOM ',A10,'  LOCAL ROTATION MATRIX')
 6090 FORMAT (30X,3F10.5)
 6100 FORMAT (13X,'EQUIV ATOM ',I3,3X,'POSITION: ',3F8.3)
 6110 FORMAT (30X,3F10.5)
 9000 FORMAT (A,' (LMMX =',I5,', LM2 =',I5,')')
 9001 FORMAT (A,' (LOMAX =',I5,', L  =',I5,')')
 9002 FORMAT (A,' (NLOAT =',I3,', L =',I3,' JATOM =',i3,')')
 9010 FORMAT (A)
 9020 FORMAT (I6,',',6(I5,','))
!
!        End of 'ATPAR'
!
      END






