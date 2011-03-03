      SUBROUTINE HAMILT(LTOP,NAT,NV,U,UP,DU,DUP,MULT,R, &
                        POS,V,VOLUME,CHN,E,INS)
!
#ifdef Parallel
      use matrices, only : H, S, HSROWS, &
                           HPANEL, SPANEL, KZZ, XK, YK, ZK, RK
#else
      use matrices, only : HS, HSDIAG, HSROWS, &
                           HPANEL, SPANEL, KZZ, XK, YK, ZK, RK
#endif
      use albl, only : AL, BL, init_albl
      use loabc, only : ALO, BLO, CLO, ELO, PLO, DPLO, PELO, &
                        DPELO, PEILO, PI12LO, PE12LO
      use parallel, only : ICTXTALL, DESCHS, DESCZ, NPE, MYID,  &
                           NPROW, NPCOL, MYROWHS, MYCOLHS, &
                           BLOCKSIZE, LDHS, LCOLHS, LCOLZ, &
                           INIT_PARALLEL, DESCPANEL, BARRIER
      use lapw_timer, only : time_albl, time_loop260, time_phase,   &
                        time_legendre, time_step_function, time_h,   &
        	       time_us, time_loop210, time_loop230,   &
      		       time_loop240, time_overlap, time_lo,   &
      		       START_TIMER, STOP_TIMER, READ_CPU_TIME, &
                       init_all_timer
      use lolog, only : nlo, loor, lapw, ilo
      use rotmat, only: ROTIJ, ROTLOC
      use struk, only : NDF, multmax
      use out, only   : WARP
!      use out, only   : WARP1, WARP2, WARP3
      IMPLICIT NONE
      INCLUDE 'param.inc'
!
!        Parameters
!
      DOUBLE PRECISION   ONE, TOL, UTOP, ZERO, HALF
      DOUBLE PRECISION   TWO, THREE, FOUR
      PARAMETER          (ONE=1.0D+0, TOL=1.0D-8)
      PARAMETER          (TWO=2.0D+0, THREE=3.0D0, FOUR=4.0D0)
      PARAMETER          (HALF=0.5D0)
      PARAMETER          (UTOP=0.9999999D+0, ZERO=0.0D+0)
      DOUBLE COMPLEX         CIMAG
      PARAMETER          (CIMAG = (0.0D0,1.0D0))
      DOUBLE COMPLEX         CZERO
      PARAMETER          (CZERO=(0.0D0,0.0D0))
!
!        Arguments
!
      INTEGER            INS, LTOP, NAT, NV
      INTEGER            MULT(NAT)
      DOUBLE PRECISION   VOLUME
      DOUBLE PRECISION   R(NAT)
!      DOUBLE PRECISION   XK(HSROWS+1), YK(HSROWS+1), ZK(HSROWS+1)
      DOUBLE PRECISION   CHN(0:LMAX-1,NAT), DU(0:LMAX-1,NAT)
      DOUBLE PRECISION   DUP(0:LMAX-1,NAT), E(0:LMAX-1,NAT)
      DOUBLE PRECISION   POS(3,NDF)
      DOUBLE PRECISION   U(0:LMAX-1,NAT), UP(0:LMAX-1,NAT), V(NAT)

!        YL(lm,:,i) - temporary storage of spherical harmonics for
!                     local orbitals
!
      DOUBLE COMPLEX, allocatable :: YL(:,:,:)
!
!        Locals
!
      DOUBLE PRECISION   V3
      INTEGER            jlo,jlop
      DOUBLE PRECISION   atmp1, tmp3, stmp, htmp, btmp1
      DOUBLE PRECISION   akinlor,akinlorup,akinloiup,akinbr,akinbi
      INTEGER            I, INDATM, J, JEQ, JEQO, JNEQ, K, L
      INTEGER            L0M0, MO1, MO2, N, NNLO, IATM, IHELP, IOUTER
      INTEGER            IMIN, IMAX, II
      DOUBLE PRECISION   ARGX, ARGY, ARGZ, ATMP, BESR, BTMP
      DOUBLE PRECISION   ARGXX, ARGYY, ARGZZ
      DOUBLE PRECISION   C1R, C1I, C2R, C2I, C3R, C3I
      DOUBLE PRECISION   C11R, C11I, C12R, C12I, C13R, C13I
      DOUBLE PRECISION   C1RUP, C1IUP, C2RUP, C2IUP, C3RUP, C3IUP
      DOUBLE PRECISION   C11RUP,C11IUP,C12RUP,C12IUP,C13RUP,C13IUP
      DOUBLE PRECISION   C6R, C6I, C7R, C7I, C8R, C8I
      DOUBLE PRECISION   CABR, CABI, CABCR, CABCI, CLR, CLI
      DOUBLE PRECISION   CLYLR, CLYLI, CLYLNR, CLYLNI
      DOUBLE PRECISION   DLL1, PI4R2V
      DOUBLE PRECISION   DTIME1, DTIME2, DTIMH, DTIMLG, DTIMMA, DTIMPH
      DOUBLE PRECISION   DTIME3, dtime240, dtime241, dtime242
      DOUBLE PRECISION   dtime260, dtime211, dtime231, dtime210,dtime230
      DOUBLE PRECISION   DTIMS, DTIMUS, PI, Q, RKN, TMP, X2

      DOUBLE PRECISION, allocatable ::   DEN(:), DFJ(:), FJ(:), FK(:)
      DOUBLE PRECISION, allocatable ::   PPLX(:), PPLY(:), ROTV1(:)
      DOUBLE PRECISION, allocatable ::   ROTV2(:), SLEN(:), TMP1(:)
      DOUBLE PRECISION, allocatable ::   TMP2(:), VEC(:), VIFPR4(:)
      DOUBLE PRECISION, allocatable ::   X(:), XL(:)
      DOUBLE PRECISION, allocatable ::   P(:,:)
      DOUBLE PRECISION, allocatable ::   help1(:), help_cos(:)
!_COMPLEX      DOUBLE COMPLEX, allocatable ::     help_exp(:)
      DOUBLE PRECISION, allocatable ::   help_sin(:), tmp_y(:)
#if defined (INTEL_VML)
!_COMPLEX      DOUBLE PRECISION, allocatable ::   help_tmpsin(:), help_tmpcos(:)
      DOUBLE PRECISION, allocatable ::   tmp_ytmp(:)
#endif
!_REAL      DOUBLE PRECISION ::   PHS
!_REAL      DOUBLE PRECISION, allocatable ::   PHASE(:,:)
!_COMPLEX      DOUBLE COMPLEX ::     PHS
!_COMPLEX      DOUBLE COMPLEX, allocatable ::     PHASE(:,:)
      DOUBLE COMPLEX, allocatable ::         PHSC(:)
!
!        External Functions
!
!_REAL      DOUBLE PRECISION   USTPHX
!_COMPLEX      DOUBLE COMPLEX         USTPHX
      DOUBLE COMPLEX     WARPIN
      EXTERNAL           USTPHX, WARPIN
!
!        External Subroutines
!
      EXTERNAL           CPUTIM, DVBES1, ROTATE, SPHBES, YLM
!
!        Intrinsic Functions
!
      INTRINSIC          ABS, ATAN, DCONJG, COS, DBLE, DCMPLX, EXP, SIN
      INTRINSIC          SQRT

      CALL BARRIER

      allocate( DEN(HSROWS), DFJ(0:LMAX-1), FJ(0:LMAX-1), FK(3) )
      allocate( PPLX(HSROWS), PPLY(HSROWS), ROTV1(3) )
      allocate( ROTV2(3), SLEN(HSROWS), TMP1(HSROWS) )
      allocate( TMP2(HSROWS), VEC(3), VIFPR4(NAT) )
      allocate( X(HSROWS), XL(LMAX-2) )
      allocate( P(HSROWS,0:LMAX-1) )
      allocate( help1(HSROWS), help_cos(HSROWS), help_sin(HSROWS) )
#if defined (INTEL_VML)
!_COMPLEX      allocate( help_tmpcos(HSROWS), help_tmpsin(HSROWS) )
      allocate( tmp_ytmp(HSROWS) )
#endif
      allocate( tmp_y(HSROWS) )
!_COMPLEX      allocate( help_exp(HSROWS) )
      allocate( PHASE(HSROWS,NAT) )
      allocate( PHSC(HSROWS) )

!      allocate( YL(0:(LOMAX+1)*(LOMAX+1)-1,HSROWS,NDF) )
      allocate( YL(0:(LOMAX+1)*(LOMAX+1)-1,HSROWS,multmax) )
!
      PI = FOUR*ATAN(ONE)
      DO 10 I = 1, LMAX - 2
         XL(I) = DBLE(I)/DBLE(I+1)
   10 CONTINUE
      DTIMPH = ZERO
      DTIMUS = ZERO
      DTIMLG = ZERO
      DTIMMA = ZERO
      DTIMS = ZERO
      DTIMH = ZERO
!
!        precompute al(kn) and bl(kn)
!
!      allocate( AL(HSROWS, 0:LMAX-1, NDF) )
!      allocate( BL(HSROWS, 0:LMAX-1,NDF) )
!      call init_albl(hsrows,lmax,ndf)
      call init_albl(hsrows,ltop,nat)
      CALL START_TIMER(time_albl)
      DO 40 N = 1, NV
         RKN = RK(N)
         DO 30 JNEQ = 1, NAT
            CALL SPHBES(LTOP-1,R(JNEQ)*RKN,FJ)
            CALL DVBES1(FJ,DFJ,RKN,R(JNEQ),LTOP)
            DO 20 L = 0, LTOP - 1
               if(lapw(l,jneq)) then
               AL(N,L,JNEQ) = RKN*DFJ(L)*UP(L,JNEQ) - FJ(L)*DUP(L,JNEQ)
               BL(N,L,JNEQ) = FJ(L)*DU(L,JNEQ) - RKN*DFJ(L)*U(L,JNEQ)
               else
               AL(N,L,JNEQ) = FJ(L)/U(L,JNEQ)/R(JNEQ)**2
               BL(N,L,JNEQ) = ZERO
               endif
   20       CONTINUE
   30    CONTINUE
   40 CONTINUE
      CALL STOP_TIMER(time_albl)
      WRITE (6,*) 'Time for al,bl: ', READ_CPU_TIME(time_albl)
!
!        compute Overlap-matrix
!
      DO 50 JNEQ = 1, NAT
         VIFPR4(JNEQ) = FOUR*PI*VOLUME*R(JNEQ)**4
   50 CONTINUE


      CALL START_TIMER(time_loop260)

      DO IOUTER = 0, NV/(NPE*BLOCKSIZE)
	IMIN = (IOUTER * NPE + MYID) * BLOCKSIZE + 1
	IMAX = MIN(NV, (IOUTER * NPE + MYID + 1) * BLOCKSIZE)
	IHELP = 0
	imin1: DO I = IMIN, IMAX
	  IHELP = IHELP + 1
!      DO 260 I = 1, NV
!
!        calculate phase-factors
!
         CALL START_TIMER(time_phase)
         INDATM = 0
         DO 90 JNEQ = 1, NAT
            DO 60 J = 1, I
!_REAL               PHASE(j,JNEQ) = ZERO
!_COMPLEX               PHASE(j,JNEQ) = CZERO
   60       CONTINUE
            DO 80 JEQ = 1, MULT(JNEQ)
               INDATM = INDATM + 1

                  ARGXX = POS(1,INDATM)*KZZ(1,I)
                  ARGXX = ARGXX + POS(2,INDATM)*KZZ(2,I)
                  ARGXX = ARGXX + POS(3,INDATM)*KZZ(3,I)

               DO 70 J = 1, I
                  ARGX = POS(1,INDATM)*KZZ(1,J) - ARGXX
                  ARGX = POS(2,INDATM)*KZZ(2,J) + ARGX
                  ARGX = POS(3,INDATM)*KZZ(3,J) + ARGX
              help1(j)=TWO*PI*(ARGX)


   70          CONTINUE

!MASS-Lib http://www.rs6000.ibm.com/resource/technology/MASS
!MASS-Lib help_cos(j)-cos(help1(j)), j=1,i
! or vmllib from Intel
!
#if defined (INTEL_VML)
!_REAL                call vdcos(i,help1,help_cos)
!_COMPLEX             call vdcos(i,help1,help_tmpcos)
!_COMPLEX             call vdsin(i,help1,help_tmpsin)
!_COMPLEX             do j = 1, i
!_COMPLEX                help_exp(j) = dcmplx(help_tmpcos(j),help_tmpsin(j))
!_COMPLEX             end do
#else
!_REAL                call vcos(help_cos,help1,i)
!_COMPLEX             call vcosisin(help_exp,help1,i)
#endif

              do 71 j = 1,i
!_REAL                  PHASE(j,JNEQ) = PHASE(j,JNEQ) + help_cos(j)
!_COMPLEX                   PHASE(j,JNEQ) = PHASE(j,JNEQ) + help_exp(j)
   71           continue

   80       CONTINUE
   90    CONTINUE
         CALL STOP_TIMER(time_phase)
         DTIMPH = DTIMPH + READ_CPU_TIME(time_phase)
	 CALL START_TIMER(time_legendre)
!
!        precompute Legendre-Polynomials
!
         DO 100 J = 1, I
            X(J) = XK(I)*XK(J) + YK(I)*YK(J) + ZK(I)*ZK(J)
            DEN(J) = RK(I)*RK(J)
  100    CONTINUE
         DO 110 J = 1, I
            IF (DEN(J) .GT. TOL) X(J) = X(J)/DEN(J)
  110    CONTINUE
         DO 120 J = 1, I
            PPLX(J) = X(J)*X(J)
            PPLY(J) = PPLX(J) - ONE
            P(J,0) = ONE
            P(J,1) = X(J)
  120    CONTINUE
         DO 140 K = 2, LTOP - 2
            DO 130 J = 1, I
               IF (ABS(X(J)) .LT. UTOP) THEN
                  P(J,K) = PPLX(J) + XL(K-1)*PPLY(J)
                  PPLX(J) = X(J)*P(J,K)
                  PPLY(J) = PPLX(J) - P(J,K-1)
                  P(J,LTOP-1) = PPLX(J) + XL(LTOP-2)*PPLY(J)
               ELSE
                  P(J,K) = X(J)*P(J,K-1)
               ENDIF
  130       CONTINUE
  140    CONTINUE
         DO 150 J = 1, I
            IF (ABS(X(J)) .LT. UTOP) THEN
               P(J,LTOP-1) = PPLX(J) + XL(LTOP-2)*PPLY(J)
            ELSE
               P(J,LTOP-1) = X(J)*P(J,LTOP-2)
            ENDIF
  150    CONTINUE
         CALL STOP_TIMER(time_legendre)
         DTIMLG = DTIMLG + READ_CPU_TIME(time_legendre)
	 CALL START_TIMER(time_us)
	 CALL START_TIMER(time_step_function)
!
!        initialize Overlap-matrixelements with step-function values
!
         IF (INS .EQ. 0) THEN
            DO 160 J = 1, I
!_REAL               SPANEL(J,IHELP) = ZERO
!_COMPLEX               SPANEL(J,IHELP) = CZERO
               SLEN(J) = SQRT((XK(J)-XK(I))*(XK(J)-XK(I))+ &
                              (YK(J)-YK(I))*(YK(J)-YK(I))+ &
                              (ZK(J)-ZK(I))*(ZK(J)-ZK(I)))
  160       CONTINUE
            DO 180 JNEQ = 1, NAT
              j1: do j=1,i-1
                help1(j) = slen(j)*r(jneq)
                tmp_y(j)   = help1(j)*help1(j)*help1(j)
              end do j1

!MASS-Lib http://www.rs6000.ibm.com/resource/technology/MASS
!MASS-Lib help_sin(j)=sin(help1(j)), help_cos(j)-cos(help1(j))

#if defined (INTEL_VML)
              call vdsincos(i-1,help1,help_sin,help_cos)
#else
              call vsincos(help_sin,help_cos,help1,i-1)
#endif

!MASS-Lib http://www.rs6000.ibm.com/resource/technology/MASS
!MASS-Lib: tmp_y(j)=1/tmp_y(j), j=1,i-1

#if defined (INTEL_VML)
              call vdinv(i-1, tmp_y, tmp_ytmp)
              call dcopy(i-1, tmp_ytmp, 1, tmp_y, 1)
#else
              call vrec(tmp_y,tmp_y,i-1)
#endif

              v3=THREE*V(JNEQ)
              DO 170 J = 1, I-1
                BESR = (help1(j)*help_COS(j)-help_SIN(j))*tmp_y(j)
                SPANEL(J,IHELP) = SPANEL(J,IHELP)+v3*PHASE(j,JNEQ)*BESR
  170         CONTINUE
!
              SPANEL(I,IHELP) = SPANEL(I,IHELP) + ONE/DBLE(NAT) - &
                           V(JNEQ)*DBLE(MULT(JNEQ))

  180       CONTINUE
            CALL STOP_TIMER(time_step_function)
            DTIMS = DTIMS + READ_CPU_TIME(time_step_function)
	    CALL START_TIMER(time_h)
            DO 190 J = 1, I
               X2 = XK(J)*XK(i) + YK(J)*YK(i) + ZK(J)*ZK(i)
!_REAL               HPANEL(J,IHELP) = X2* SPANEL(J,IHELP) + &
!_REAL                            DBLE(WARPIN(KZZ(1,J)-KZZ(1,I), &
!_REAL                                        KZZ(2,J)-KZZ(2,I), &
!_REAL                                        KZZ(3,J)-KZZ(3,I)))
!_COMPLEX               HPANEL(J,IHELP) = X2* SPANEL(J,IHELP) + &
!_COMPLEX                                 WARPIN(KZZ(1,J)-KZZ(1,I), &
!_COMPLEX                                        KZZ(2,J)-KZZ(2,I), &
!_COMPLEX                                        KZZ(3,J)-KZZ(3,I))
  190       CONTINUE
	    CALL STOP_TIMER(time_h)
            DTIMH = DTIMH + READ_CPU_TIME(time_h)
         ELSE
            DO 200 J = 1, I
               FK(1) = XK(J) - XK(I)
               FK(2) = YK(J) - YK(I)
               FK(3) = ZK(J) - ZK(I)
               X2 = XK(J)*XK(i) + YK(J)*YK(i) + ZK(J)*ZK(i)
               SPANEL(J,IHELP) = USTPHX(FK,PHASE,j,NAT,V,R,MULT)
               HPANEL(J,IHELP) = X2* SPANEL(J,IHELP)
!		       'ustphx: ', USTPHX(FK,PHASE,j,NAT,V,R,MULT), &
!		       'x2: ', X2, &
!		       'FK: ', FK, &
!		       'PHASE: ', PHASE, &
!		       'J: ', j, &
!		       'NAT: ', NAT, &
!		       'V: ', V, &
!		       'R: ', R, &
!		       'MULT: ', MULT
  200       CONTINUE
         ENDIF
!
         CALL STOP_TIMER(time_us)
         DTIMUS = DTIMUS + READ_CPU_TIME(time_us)
	 CALL START_TIMER(time_overlap)
         DO 250 JNEQ = 1, NAT
!
!        calculate Overlap-matrixelement
!

      CALL START_TIMER(time_loop210)

            DO 210 J = 1, I
               TMP1(J) = ZERO
               TMP2(J) = ZERO
  210       CONTINUE
      CALL STOP_TIMER(time_loop210)
      dtime210=dtime210+READ_CPU_TIME(time_loop210)
      CALL START_TIMER(time_loop230)

            DO 230 L = 0, LTOP - 1
           if(.not.lapw(l,jneq)) then
              ATMP = AL(I,L,JNEQ)
              atmp1 = ATMP*U(L,JNEQ)*DU(L,JNEQ)*R(JNEQ)**2
              DLL1 = DBLE(L+L+1)
              j2: DO J = 1, I
                 TMP3 = DLL1*P(J,L)
                 STMP = ATMP*AL(J,L,JNEQ)
                 STMP = STMP*TMP3
                 TMP2(J) = TMP2(J) + STMP
                 HTMP = STMP*E(L,JNEQ)
                 TMP1(J) = TMP1(J) + HTMP
!  surface part of kinetic energy 
                 tmp1(j) = tmp1(j)+TMP3*ATMP1*AL(J,L,JNEQ)
              end do j2
        else
!....for lapw all BL lines must be kept
           ATMP = AL(I,L,JNEQ)
           BTMP = BL(I,L,JNEQ)*CHN(L,JNEQ)
           atmp1 = R(JNEQ)**2* (ATMP*U(L,JNEQ)*DU(L,JNEQ) + &
                    BL(I,L,JNEQ)*UP(L,JNEQ)*DU(L,JNEQ))
           btmp1 = R(JNEQ)**2*(BL(I,L,JNEQ)*UP(L,JNEQ)*DUP(L,JNEQ) + &
                    AL(I,L,JNEQ)*U(L,JNEQ)*DUp(L,JNEQ))
           DLL1 = DBLE(L+L+1)
           j3: DO J = 1, I
              TMP3 = DLL1*P(J,L)
              STMP = ATMP*AL(J,L,JNEQ)
              STMP = STMP + BTMP*BL(J,L,JNEQ)
              STMP = STMP*TMP3
              TMP2(J) = TMP2(J) + STMP
              HTMP = ATMP*BL(J,L,JNEQ)
              HTMP = HTMP*TMP3
              HTMP = HTMP + STMP*E(L,JNEQ)
              TMP1(J) = TMP1(J) + HTMP
!  surface part of kinetic energy 
              tmp1(j) = tmp1(j)+TMP3*( &
                   ATMP1*AL(J,L,JNEQ)+ &
                   Btmp1*BL(J,L,JNEQ))
           end do j3
        endif
  230       CONTINUE

      CALL STOP_TIMER(time_loop230)
      dtime230=dtime230+READ_CPU_TIME(time_loop230)
      CALL START_TIMER(time_loop240)
            DO 240 J = 1, I
               TMP1(J) = TMP1(J)*VIFPR4(JNEQ)
               TMP2(J) = TMP2(J)*VIFPR4(JNEQ)
               SPANEL(J,IHELP) = SPANEL(J,IHELP) + PHASE(J,JNEQ)*TMP2(J)
               HPANEL(J,IHELP) = HPANEL(J,IHELP) + PHASE(J,JNEQ)*TMP1(J)
  240       CONTINUE

      CALL STOP_TIMER(time_loop240)
      dtime240=dtime240+READ_CPU_TIME(time_loop240)
  250    CONTINUE
         CALL STOP_TIMER(time_overlap)
         DTIMMA = DTIMMA + READ_CPU_TIME(time_overlap)
	END DO imin1

	IHELP = 0
	IMIN = IOUTER * NPE * BLOCKSIZE + 1
	IMAX = MIN(NV, (IOUTER + 1) * NPE * BLOCKSIZE)
#ifdef Parallel
!_REAL          CALL PDGEMR2D(IMAX,IMAX-IMIN+1, HPANEL,1,1, DESCPANEL, &
!_REAL                             H, 1,IMIN, DESCHS, ICTXTALL)
!_REAL          CALL PDGEMR2D(IMAX,IMAX-IMIN+1, SPANEL,1,1, DESCPANEL, &
!_REAL                             S, 1,IMIN, DESCHS, ICTXTALL)
!_COMPLEX          CALL PZGEMR2D(IMAX,IMAX-IMIN+1, HPANEL,1,1, DESCPANEL, &
!_COMPLEX                             H, 1,IMIN, DESCHS, ICTXTALL)
!_COMPLEX          CALL PZGEMR2D(IMAX,IMAX-IMIN+1, SPANEL,1,1, DESCPANEL, &
!_COMPLEX                             S, 1,IMIN, DESCHS, ICTXTALL)
        IHELP=IMAX-IMIN+1
#else
	imin2: DO I = IMIN, IMAX
	  IHELP = IHELP + 1
!_REAL        CALL DCOPY(I-1, HPANEL(1,IHELP), 1, HS(I, 1), HSROWS )
!_REAL        CALL DCOPY(I, SPANEL(1,IHELP), 1, HS(1,I), 1 )
!_COMPLEX        CALL ZCOPY(I-1, HPANEL(1,IHELP), 1, HS(I, 1), HSROWS )
!_COMPLEX        CALL ZCOPY(I, SPANEL(1,IHELP), 1, HS(1,I), 1 )
          HSDIAG(I)=HPANEL(I,IHELP)
	END DO imin2
#endif
      END DO
!  260 CONTINUE

      CALL STOP_TIMER(time_loop260)

      CALL BARRIER
      WRITE (6,*) 'nlo', NLO
      IF (NLO .NE. 0) THEN
!
!------ calculate LO-Overlap-matrixelements ---------------------------
!
         CALL START_TIMER(time_lo)
         INDATM = 0
         IATM = 0
         NNLO = 0
!
!        initialize matrixelements
!
!_REAL       HPANEL(1:HSROWS, 1:BLOCKSIZE)=ZERO
!_REAL       SPANEL(1:HSROWS, 1:BLOCKSIZE)=ZERO
!_COMPLEX       HPANEL(1:HSROWS, 1:BLOCKSIZE)=CZERO
!_COMPLEX       SPANEL(1:HSROWS, 1:BLOCKSIZE)=CZERO

	 I = NV + 1
	 IMIN = I
	 IMAX = MIN(NV + NLO,  &
                    ((IMIN-1)/(BLOCKSIZE*NPE)+1)*BLOCKSIZE*NPE)
         DO 340 JNEQ = 1, NAT
            DO 272 JEQ = 1, MULT(JNEQ)
               IATM = IATM + 1
!
!        precalculate spherical-harmonics for Alm, Blm
!
               DO 271 J = 1, NV
                  VEC(1) = XK(J)
                  VEC(2) = YK(J)
                  VEC(3) = ZK(J)
                  CALL ROTATE(VEC,ROTIJ(1,1,IATM),ROTV1)
                  CALL ROTATE(ROTV1,ROTLOC(1,1,JNEQ),ROTV2)
                  CALL YLM(ROTV2,LOMAX,YL(0,J,JEQ))
  271          CONTINUE
  272       CONTINUE
            PI4R2V = FOUR*PI*SQRT(VOLUME)*R(JNEQ)*R(JNEQ)
            CLR = ONE
            CLI = ZERO
            DO 330 L = 0, LOMAX
               do 331 jlo=1,ilo(l,jneq)
!               IF (LOOR(L,JNEQ)) THEN
!        precalculate factors of a,b,c for LO's
!        C1,C2,C3 are needed for the Overlap-matrix update
!        C11,C12,C13 are needed for the Hamilton-matrix update
!
          C1R  = ALO(L,jlo,JNEQ) + PI12LO(L,JNEQ)*CLO(L,jlo,JNEQ)
                  C1I  = ZERO
                  C2R  = BLO(L,jlo,JNEQ)*CHN(L,JNEQ) + &
                         PE12LO(L,JNEQ)*CLO(L,jlo,JNEQ)
                  C2I  = ZERO
          C11R = ALO(L,jlo,JNEQ)*E(L,JNEQ) + HALF*(BLO(L,jlo,JNEQ) + &
                 CLO(L,jlo,JNEQ)*PI12LO(L,JNEQ)*(ELO(L,jlo,JNEQ) + &
                         E(L,JNEQ)))
                  C11I = ZERO
                  C12R = BLO(L,jlo,JNEQ)*CHN(L,JNEQ)*E(L,JNEQ) + &
                       HALF*(ALO(L,jlo,JNEQ) + CLO(L,jlo,JNEQ)* &
                       PI12LO(L,JNEQ) + CLO(L,jlo,JNEQ)*PE12LO(L,JNEQ)* &
                       (ELO(L,jlo,JNEQ) + E(L,JNEQ)))
                  C12I = ZERO
          C3R  = CLO(L,jlo,JNEQ) +  &
                       PE12LO(L,JNEQ)*BLO(L,jlo,JNEQ) + &
                       ALO(L,jlo,JNEQ)*PI12LO(L,JNEQ)
          C3I  = ZERO
          C13R = CLO(L,jlo,JNEQ)*ELO(L,jlo,JNEQ)+ &
                       HALF*(BLO(L,jlo,JNEQ)*PI12LO(L,JNEQ) + &
                       (ELO(L,jlo,JNEQ) + E(L,JNEQ))*(ALO(L,jlo,JNEQ)* &
              PI12LO(L,JNEQ)+PE12LO(L,JNEQ)*BLO(L,jlo,JNEQ)))
                  C13I = ZERO
!
! Precalculate surface kinetic energy factor
!
          akinloR=HALF*r(jneq)**2*(alo(l,jlo,jneq)*DU(L,JNEQ)+ &
                       blo(l,jlo,jneq)*DUP(L,JNEQ)+ &
                       CLO(l,jlo,jneq)*DPLO(l,jneq))

                  DO 320 JEQO = 1, MULT(JNEQ)
                     DO 310 MO1 = -L, L
                        NNLO = NNLO + 1
!
!    compute only local columns
!
		    IHELP=MOD(I,BLOCKSIZE)
		    IF(IHELP.EQ.0) IHELP=BLOCKSIZE
                      DO 290 JEQ = 1, MULT(JNEQ)
                        INDATM = INDATM + 1
!                        IF(MOD((I+BLOCKSIZE*NPE-1)/BLOCKSIZE, 
!     +                         NPE).EQ.MYID) THEN
                           DO 275 J = 1, NV + NNLO
!
!        determine phase factor
!
                              ARGX = POS(1,INDATM)* &
                                     (KZZ(1,NV+NNLO)-KZZ(1,J))
                              ARGY = POS(2,INDATM)* &
                                     (KZZ(2,NV+NNLO)-KZZ(2,J))
                              ARGZ = POS(3,INDATM)* &
                                     (KZZ(3,NV+NNLO)-KZZ(3,J))
                              PHSC(J) = EXP(CMPLX(ZERO,TWO)*PI &
                                            *(ARGX+ARGY+ARGZ))
                              PHSC(J) = DCONJG(PHSC(J))
  275                      CONTINUE
!
!        calculate spherical-harmonics for Alm, Blm
!
                           VEC(1) = XK(NV+NNLO)
                           VEC(2) = YK(NV+NNLO)
                           VEC(3) = ZK(NV+NNLO)
                           CALL ROTATE(VEC,ROTIJ(1,1,INDATM),ROTV1)
                           CALL ROTATE(ROTV1,ROTLOC(1,1,JNEQ),ROTV2)
                           CALL YLM(ROTV2,L,YL(0,NV+NNLO,JEQ))
                           DO 280 MO2 = -L, L
                              L0M0 = L*(L+1) + MO2
!
!        CLYLN = (CIMAG**L)*DCONJG(YL(L0M0,NV+NNLO,JEQ))
!
                              CLYLNR = CLR*DBLE(YL(L0M0,NV+NNLO,JEQ)) + &
                                       CLI*DIMAG(YL(L0M0,NV+NNLO,JEQ))
                              CLYLNI = CLI*DBLE(YL(L0M0,NV+NNLO,JEQ)) - &
                                       CLR*DIMAG(YL(L0M0,NV+NNLO,JEQ))
!
!        update C1,C2,C3 for the Overlap-matrix 
!        C11,C12,C13 for the Hamilton-matrix 
!
                              C1RUP = C1R*CLYLNR*PI4R2V
                              C1IUP = C1R*CLYLNI*PI4R2V
                              C2RUP = C2R*CLYLNR*PI4R2V
                              C2IUP = C2R*CLYLNI*PI4R2V
                              C3RUP = C3R*CLYLNR*PI4R2V
                              C3IUP = C3R*CLYLNI*PI4R2V
                              C11RUP = C11R*CLYLNR*PI4R2V
                              C11IUP = C11R*CLYLNI*PI4R2V
                              C12RUP = C12R*CLYLNR*PI4R2V
                              C12IUP = C12R*CLYLNI*PI4R2V
                              C13RUP = C13R*CLYLNR*PI4R2V
                              C13IUP = C13R*CLYLNI*PI4R2V
                              akinlorup = akinloR*CLYLNR*PI4R2V
                              akinloiup = akinloR*CLYLNI*PI4R2V

                              DO 277 J = 1, NV
!
!        CLYL = (CIMAG**L)*DCONJG(YL(L0M0,J,JEQ))
!
                                 CLYLR = CLR*DBLE(YL(L0M0,J,JEQ)) + &
                                         CLI*DIMAG(YL(L0M0,J,JEQ))
                                 CLYLI = CLI*DBLE(YL(L0M0,J,JEQ)) - &
                                         CLR*DIMAG(YL(L0M0,J,JEQ))
!
!        C6 ... Alm, C7 ... Blm
!
                                 C6R = CLYLR*AL(J,L,JNEQ)*PI4R2V
                                 C6I = CLYLI*AL(J,L,JNEQ)*PI4R2V
                                 C7R = CLYLR*BL(J,L,JNEQ)*PI4R2V
                                 C7I = CLYLI*BL(J,L,JNEQ)*PI4R2V
!
!        CAB = DCONJG(C6)*C1 + DCONJG(C7)*C2
!
                                 CABR =  C6R*C1RUP + C7R*C2RUP + &
                                         C6I*C1IUP + C7I*C2IUP
                                 CABI = -C6I*C1RUP - C7I*C2RUP + &
                                         C6R*C1IUP + C7R*C2IUP
!
!        update Overlap-matrix
!
                        IF(MOD((I+BLOCKSIZE*NPE-1)/BLOCKSIZE,  &
                               NPE).EQ.MYID) THEN
                                SPANEL(J,IHELP) = SPANEL(J,IHELP) + &
                                             (CABR*DBLE(PHSC(J)) - &
                                             CABI*DIMAG(PHSC(J)))
!_COMPLEX                                SPANEL(J,IHELP) = SPANEL(J,IHELP) +  &
!_COMPLEX                                    CIMAG *  (CABR*DIMAG(PHSC(J)) + &
!_COMPLEX                                             CABI*DBLE(PHSC(J)))
                        END IF
!
!        CAB = DCONJG(C6)*C11 + DCONJG(C7)*C12
!
                                 CABR =  C6R*C11RUP + C7R*C12RUP + &
                                         C6I*C11IUP + C7I*C12IUP
                                 CABI = -C6I*C11RUP - C7I*C12RUP + &
                                         C6R*C11IUP + C7R*C12IUP
                     akinbr = (c6r*u(l,jneq)+c7r*up(l,jneq)) * akinlorup &
                             +(c6i*u(l,jneq)+c7i*up(l,jneq)) * akinloiup
                     akinbi = (c6r*u(l,jneq)+c7r*up(l,jneq)) * akinloiup &
                             -(c6i*u(l,jneq)+c7i*up(l,jneq)) * akinlorup
!
!        update Hamilton-matrix
!
                        IF(MOD((I+BLOCKSIZE*NPE-1)/BLOCKSIZE,  &
                               NPE).EQ.MYID) THEN
                                 HPANEL(J,IHELP) = HPANEL(J,IHELP) + &
                                              (CABR*DBLE(PHSC(J)) - &
                                               CABI*DIMAG(PHSC(J)))
!_COMPLEX                                 HPANEL(J,IHELP) = HPANEL(J,IHELP) + &
!_COMPLEX                                       CIMAG *(CABR*DIMAG(PHSC(J)) + &
!_COMPLEX                                               CABI*DBLE(PHSC(J)))
                                 HPANEL(J,IHELP) = HPANEL(J,IHELP) + &
                                              (akinbr*DBLE(PHSC(J)) - &
                                               akinbi*DIMAG(PHSC(J)))
!_COMPLEX                                 HPANEL(J,IHELP) = HPANEL(J,IHELP) +  &
!_COMPLEX                                     CIMAG *  (akinbr*DIMAG(PHSC(J)) + &
!_COMPLEX                                               akinbi*DBLE(PHSC(J)))
                        END IF
  277                         CONTINUE
                        do  279 jlop=1,jlo
                              DO 278 J = NV + NNLO - (MO1+L) - &
             (2*L+1)*(JEQO-1) - (jlo-jlop)*(2*l+1)*mult(jneq), NV + NNLO &
             -(jlo-jlop)*(MO1+L+1)-(jlo-jlop)*(JEQO-1)*(2*l+1)
!
!        CLYL = (CIMAG**L)*DCONJG(YL(L0M0,J,JEQ))
!
                                 CLYLR = CLR*DBLE(YL(L0M0,J,JEQ)) + &
                                         CLI*DIMAG(YL(L0M0,J,JEQ))
                                 CLYLI = CLI*DBLE(YL(L0M0,J,JEQ)) - &
                                         CLR*DIMAG(YL(L0M0,J,JEQ))
!
!        C6 ... CONJG(ALOlm), C7 ... CONJG(BLOlm), C8 ... CONJG(CLOlm)
!
                                 C6R = CLYLR*ALO(L,jlop,JNEQ)*PI4R2V
                                 C6I = CLYLI*ALO(L,jlop,JNEQ)*PI4R2V
                                 C7R = CLYLR*BLO(L,jlop,JNEQ)*PI4R2V
                                 C7I = CLYLI*BLO(L,jlop,JNEQ)*PI4R2V
                                 C8R = CLYLR*CLO(L,jlop,JNEQ)*PI4R2V
                                 C8I = CLYLI*CLO(L,jlop,JNEQ)*PI4R2V
!                                                         
!        CABC = DCONJG(C6)*C1 + DCONJG(C7)*C2 + DCONJG(C8)*C3
!
                                 CABCR =  C6R*C1RUP + C6I*C1IUP + &
                                          C7R*C2RUP + C7I*C2IUP +    &
                                          C8R*C3RUP + C8I*C3IUP    
                                 CABCI = -C6I*C1RUP + C6R*C1IUP - &
                                          C7I*C2RUP + C7R*C2IUP - &
                                          C8I*C3RUP + C8R*C3IUP 
!
!        update Overlap-matrix
!
                        IF(MOD((I+BLOCKSIZE*NPE-1)/BLOCKSIZE,  &
                               NPE).EQ.MYID) THEN
                                 SPANEL(J,IHELP) = SPANEL(J,IHELP) + &
                                              (CABCR*DBLE(PHSC(J)) - &
                                               CABCI*DIMAG(PHSC(J)))
!_COMPLEX                                 SPANEL(J,IHELP) = SPANEL(J,IHELP) +  &
!_COMPLEX                                     CIMAG *  (CABCR*DIMAG(PHSC(J)) + &
!_COMPLEX                                               CABCI*DBLE(PHSC(J)))
                        END IF
!
!        CABC = DCONJG(C6)*C11 + DCONJG(C7)*C12 + DCONJG(C8)*C13
!
                                 CABCR =  C6R*C11RUP + C6I*C11IUP + &
                                          C7R*C12RUP + C7I*C12IUP +    &
                                          C8R*C13RUP + C8I*C13IUP    
                                 CABCI = -C6I*C11RUP + C6R*C11IUP - &
                                          C7I*C12RUP + C7R*C12IUP - &
                                          C8I*C13RUP + C8R*C13IUP 
!
!        update Hamilton-matrix
!
                        IF(MOD((I+BLOCKSIZE*NPE-1)/BLOCKSIZE,  &
                               NPE).EQ.MYID) THEN
                                    HPANEL(J,IHELP) = HPANEL(J,IHELP) + &
                                              (CABCR*DBLE(PHSC(J)) - &
                                              CABCI*DIMAG(PHSC(J)))
!_COMPLEX                                 HPANEL(J,IHELP) = HPANEL(J,IHELP) +  &
!_COMPLEX                                      CIMAG * (CABCR*DIMAG(PHSC(J)) + &
!_COMPLEX                                               CABCI*DBLE(PHSC(J)))
                        END IF
  278                         CONTINUE
  279                       CONTINUE
  280                      CONTINUE
!                          END IF
  290                   CONTINUE
! Achtung: PGI Compiler: Write-Statements um BARRIER sind notwendig!
!                        write(6,*) 'Vor Barrier (lo)'
!                        CALL BARRIER
!                        write(6,*) 'Nach Barrier (lo)'
			IF(I.EQ.IMAX) THEN
#ifdef Parallel
			    IHELP = MOD(IMIN,BLOCKSIZE*NPE)
			    IF(IHELP.EQ.0) IHELP=BLOCKSIZE*NPE
!_REAL          CALL PDGEMR2D(IMAX,IMAX-IMIN+1, HPANEL,1,IHELP, DESCPANEL, &
!_REAL                             H, 1,IMIN, DESCHS, ICTXTALL)
!_REAL          CALL PDGEMR2D(IMAX,IMAX-IMIN+1, SPANEL,1,IHELP, DESCPANEL, &
!_REAL                             S, 1,IMIN, DESCHS, ICTXTALL)
!_COMPLEX          CALL PZGEMR2D(IMAX,IMAX-IMIN+1, HPANEL,1,IHELP, DESCPANEL, &
!_COMPLEX                             H, 1,IMIN, DESCHS, ICTXTALL)
!_COMPLEX          CALL PZGEMR2D(IMAX,IMAX-IMIN+1, SPANEL,1,IHELP, DESCPANEL, &
!_COMPLEX                             S, 1,IMIN, DESCHS, ICTXTALL)
			    IHELP = MOD(IMAX,BLOCKSIZE*NPE)
			    IF(IHELP.EQ.0) IHELP=BLOCKSIZE*NPE
#else
			  DO II = IMIN, IMAX
			    IHELP = MOD(II,BLOCKSIZE*NPE)
			    IF(IHELP.EQ.0) IHELP=BLOCKSIZE*NPE
!_REAL        CALL DCOPY(II-1, HPANEL(1,IHELP), 1, HS(II, 1), HSROWS )
!_REAL        CALL DCOPY(II, SPANEL(1,IHELP), 1, HS(1,II), 1 )
!_COMPLEX        CALL ZCOPY(II-1, HPANEL(1,IHELP), 1, HS(II, 1), HSROWS )
!_COMPLEX        CALL ZCOPY(II, SPANEL(1,IHELP), 1, HS(1,II), 1 )
                            HSDIAG(II)=HPANEL(II,IHELP)
			  END DO
#endif
!_REAL                      HPANEL(1:HSROWS,1:BLOCKSIZE)=ZERO
!_REAL                      SPANEL(1:HSROWS,1:BLOCKSIZE)=ZERO
!_COMPLEX                      HPANEL(1:HSROWS,1:BLOCKSIZE)=CZERO
!_COMPLEX                      SPANEL(1:HSROWS,1:BLOCKSIZE)=CZERO
	                  IMIN = IMAX+1
	                  IMAX = MIN(NV + NLO,  &
                             ((IMIN-1)/(BLOCKSIZE*NPE)+1)*BLOCKSIZE*NPE)
			END IF
			I = I + 1
                        INDATM = INDATM - MULT(JNEQ)
  310                CONTINUE
  320             CONTINUE
  331           CONTINUE
  330       CONTINUE
            TMP = CLR
            CLR = -CLI
            CLI = TMP
            INDATM = INDATM + MULT(JNEQ)
  340    CONTINUE
         nevertrue: if(imax.gt.imin) THEN
#ifdef Parallel
	   IHELP = MOD(IMIN,BLOCKSIZE*NPE)
	   IF(IHELP.EQ.0) IHELP=BLOCKSIZE*NPE
!_REAL           CALL PDGEMR2D(IMAX,IMAX-IMIN+1, HPANEL,1,IHELP, DESCPANEL, &
!_REAL                             H, 1,IMIN, DESCHS, ICTXTALL)
!_REAL           CALL PDGEMR2D(IMAX,IMAX-IMIN+1, SPANEL,1,IHELP, DESCPANEL, &
!_REAL                             S, 1,IMIN, DESCHS, ICTXTALL)
!_COMPLEX           CALL PZGEMR2D(IMAX,IMAX-IMIN+1, HPANEL,1,IHELP, DESCPANEL, &
!_COMPLEX                             H, 1,IMIN, DESCHS, ICTXTALL)
!_COMPLEX           CALL PZGEMR2D(IMAX,IMAX-IMIN+1, SPANEL,1,IHELP, DESCPANEL, &
!_COMPLEX                             S, 1,IMIN, DESCHS, ICTXTALL)
#else
	   DO II = IMIN, IMAX
	     IHELP = MOD(II,BLOCKSIZE*NPE)
	     IF(IHELP.EQ.0) IHELP=BLOCKSIZE*NPE
!_REAL        CALL DCOPY(II-1, HPANEL(1,IHELP), 1, HS(II, 1), HSROWS )
!_REAL        CALL DCOPY(II, SPANEL(1,IHELP), 1, HS(1,II), 1 )
!_COMPLEX        CALL ZCOPY(II-1, HPANEL(1,IHELP), 1, HS(II, 1), HSROWS )
!_COMPLEX        CALL ZCOPY(II, SPANEL(1,IHELP), 1, HS(1,II), 1 )
             HSDIAG(II)=HPANEL(II,IHELP)
	   END DO
#endif
         END IF nevertrue
	 CALL STOP_TIMER(time_lo)
         WRITE (6,*) 'Time for los: ',READ_CPU_TIME(time_lo)
      ENDIF

      deallocate( DEN, DFJ, FJ, FK )
      deallocate( PPLX, PPLY, ROTV1 )
      deallocate( ROTV2, SLEN, TMP1 )
      deallocate( TMP2, VEC, VIFPR4 )
      deallocate( X, XL )
      deallocate( P )
      deallocate( help1, help_cos, help_sin )
#if defined (INTEL_VML)
!_COMPLEX      deallocate( help_tmpcos, help_tmpsin )
      deallocate( tmp_ytmp )
#endif
      deallocate( tmp_y )
!_COMPLEX      deallocate( help_exp )
      deallocate( PHASE )
      deallocate( PHSC )
      deallocate( YL )

      CALL BARRIER

      WRITE (6,*) 'Time for LOOP 260:      ',READ_CPU_TIME(time_loop260)
!
      RETURN
!
!         End of 'HAMILT'
!
      END
