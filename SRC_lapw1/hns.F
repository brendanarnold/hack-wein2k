      SUBROUTINE  HNS(NV,NAT,NNSL,POS,VI,R,MULT, &
                      LQIND,LQNS,GFAC,ROTLOC,ROTIJ,VNS1,VNS2, &
                      VNS3)
!
#ifdef Parallel
      use matrices, only : H, HSROWS, HPANEL, SPANEL, PAR_SYR2K, KZZ, XK, YK, ZK
#else 
      use matrices, only : HS, HSDIAG, HSROWS, HPANEL, SPANEL, KZZ, XK, YK, ZK
#endif
      use parallel, only : BARRIER, myid, ictxtall, npe, blocksize, ictxtallr, &
			   descr, lld, ldhs, lcolhs
      use albl, only : AL, BL
      use loabc, only : ALO, BLO, CLO, ELO, PLO, DPLO, PELO, &
                        DPELO, PEILO, PI12LO, PE12LO
      use lolog, only : nlo, loor, lapw, ilo
      use loint, only : VNS1LO, VNS2LO, VNS3LO
      use struk, only : NDF
      use lapw_timer, only : READ_CPU_TIME, READ_WALL_TIME, time_pwpw, time_pwlo,  &
                             time_lolo,time_hns30,  START_TIMER, STOP_TIMER
      IMPLICIT NONE
      INCLUDE 'param.inc'
!
!        Arguments
!
      INTEGER            NNSL, NV, NAT
      INTEGER            MULT(NAT), LQIND(NAT)
      INTEGER            LQNS(6,NGAU,NAT)
      INTEGER            jlo
      DOUBLE PRECISION   VI
      DOUBLE PRECISION   R(NAT)
      DOUBLE PRECISION   ROTLOC(3,3,NAT), ROTIJ(3,3,NDF)
      DOUBLE PRECISION   POS(3,NDF)
!                , XK(HSROWS+1), YK(HSROWS+1), ZK(HSROWS+1)
      DOUBLE PRECISION   VNS1(0:NSLMAX-1,LMMX,0:NSLMAX-1,NAT)
      DOUBLE PRECISION   VNS2(0:NSLMAX-1,LMMX,0:NSLMAX-1,NAT)
      DOUBLE PRECISION   VNS3(0:NSLMAX-1,LMMX,0:NSLMAX-1,NAT)
      DOUBLE COMPLEX         GFAC(NGAU,NAT)
!
!
!        YL(lm,:,i) - temporary storage of spherical harmonics for
!                     local orbitals
!
      DOUBLE COMPLEX, allocatable :: YL(:,:)
!
      DOUBLE PRECISION, ALLOCATABLE, TARGET :: A1R(:,:), A1I(:,:)
      DOUBLE PRECISION, ALLOCATABLE, TARGET :: A2R(:,:), A2I(:,:)
      DOUBLE PRECISION, POINTER :: B1R(:,:), B1I(:,:)
      DOUBLE PRECISION, POINTER :: B2R(:,:), B2I(:,:)
      DOUBLE PRECISION, POINTER :: C1R(:,:), C1I(:,:)
      DOUBLE PRECISION, POINTER :: C2R(:,:), C2I(:,:)
      DOUBLE PRECISION   TMP(HSROWS)
!_COMPLEX      DOUBLE COMPLEX, ALLOCATABLE :: A1(:,:), A2(:,:)
!_COMPLEX      DOUBLE COMPLEX A1x,A2x,B1x,B2x,C1x,C2x
!_COMPLEX      DOUBLE PRECISION   HTMPR(HSROWS), HTMPI(HSROWS)
      DOUBLE PRECISION, ALLOCATABLE :: ALMR(:,:), ALMI(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: BLMR(:,:), BLMI(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: CLMR(:,:), CLMI(:,:)
!
!        Local Parameters
!
      DOUBLE COMPLEX         CIMAG
      PARAMETER          (CIMAG  = (0.0D0,1.0D0))
      DOUBLE PRECISION   ONE, ZERO, HALF, TWO, FOUR
      PARAMETER          (ONE=1.0D0, ZERO=0.0D0, HALF=0.5D0, TWO=2.0D0, FOUR=4.0D0)
      DOUBLE COMPLEX     CONE, CZERO, CHALF
      PARAMETER          (CONE=(1.0D0,0.0D0), CZERO=(0.0D0,0.0D0))
      PARAMETER          (CHALF=(0.5D0,0.0D0))
!
!        Locals
!
      INTEGER            I, I1, I2, J, JNEQ, JEQ, INDATM
      INTEGER            L, LINDEX
      INTEGER            L0M0, LPMP, L0, LP, M0, MP, LMX
      INTEGER            L0IDX, LMXIDX, LPIDX, M0IDX, MPIDX
      INTEGER            LMXPTR(NGAU+2), LPPTR(NGAU+2)
      INTEGER            L0PTR(NGAU+2), M0PTR(NGAU+2), MPPTR(NGAU+2)
      INTEGER            L0FLD(NGAU), LMXFLD(NGAU), LPFLD(NGAU)
      INTEGER            M0FLD(NGAU), MPFLD(NGAU)
      INTEGER            IOUTER,IMIN,IMAX,IHELP,II
      INTEGER            JOUTER,JMIN,JMAX,JHELP,JJ
      INTEGER            NSL2
!      INTEGER, EXTERNAL :: NUMROC
      PARAMETER          (NSL2=NSLMAX*NSLMAX)
      DOUBLE PRECISION   ARGX, ARGY, ARGZ, RHOATM, PI
      DOUBLE PRECISION   VEC(3), ROTV1(3), ROTV2(3)
      DOUBLE PRECISION   TUUR, TUUI, TUDR, TUDI
      DOUBLE PRECISION   TDUR, TDUI, TDDR, TDDI
!
!        integrals for local orbitals
!
      DOUBLE PRECISION   TUUR21, TUUI21, TUDR21, TUDI21, TUUR12, TUUI12
      DOUBLE PRECISION   TDUR12, TDUI12, TUUR22, TUUI22
!
      DOUBLE COMPLEX         CIL
      DOUBLE COMPLEX         PHS(HSROWS), GNTFLD(NGAU)
      LOGICAL            LOORL0
!
!        External Subroutines
!
      EXTERNAL           LMSORT, LOHNS, ROTATE, YLM
!
!        Intrinsic Functions
!
      INTRINSIC          DCONJG, EXP, SQRT, ATAN, DBLE, DIMAG
!
      CALL BARRIER
#ifdef Parallel
      allocate( A1R(lld, NSL2), A1I(lld, NSL2) )
      allocate( B1R(lld, NSL2), B1I(lld, NSL2) )
      allocate( C1R(lld, NSL2), C1I(lld, NSL2) )
      allocate( A2R(lld, NSL2), A2I(lld, NSL2) )
      allocate( B2R(lld, NSL2), B2I(lld, NSL2) )
      allocate( C2R(lld, NSL2), C2I(lld, NSL2) )

      allocate( ALMR(lld, 0:NSL2-1) )
      allocate( ALMI(lld, 0:NSL2-1) )
      allocate( BLMR(lld, 0:NSL2-1) )
      allocate( BLMI(lld, 0:NSL2-1) )
      allocate( CLMR(lld, 0:NSL2-1) )
      allocate( CLMI(lld, 0:NSL2-1) )

      allocate( YL(0:NSL2-1,lld) )
#else
      allocate( A1R(HSROWS, NSL2*3), A1I(HSROWS, NSL2*3) )
      allocate( A2R(HSROWS, NSL2*3), A2I(HSROWS, NSL2*3) )
      B1R=>A1R(1:HSROWS,  NSL2+1:2*NSL2)
      B2R=>A2R(1:HSROWS,  NSL2+1:2*NSL2)
      B1I=>A1I(1:HSROWS,  NSL2+1:2*NSL2)
      B2I=>A2I(1:HSROWS,  NSL2+1:2*NSL2)
      C1R=>A1R(1:HSROWS,2*NSL2+1:3*NSL2)
      C2R=>A2R(1:HSROWS,2*NSL2+1:3*NSL2)
      C1I=>A1I(1:HSROWS,2*NSL2+1:3*NSL2)
      C2I=>A2I(1:HSROWS,2*NSL2+1:3*NSL2)

      allocate( ALMR(HSROWS, 0:NSL2-1) )
      allocate( ALMI(HSROWS, 0:NSL2-1) )
      allocate( BLMR(HSROWS, 0:NSL2-1) )
      allocate( BLMI(HSROWS, 0:NSL2-1) )
      allocate( CLMR(HSROWS, 0:NSL2-1) )
      allocate( CLMI(HSROWS, 0:NSL2-1) )

      allocate( YL(0:NSL2-1,HSROWS) )
#endif
#ifndef Parallel
!_COMPLEX      allocate( A1(HSROWS, NSL2) )
!_COMPLEX      allocate( A2(HSROWS, NSL2) )
#endif


      PI = FOUR*ATAN(ONE)
      INDATM = 0
      DO 110 JNEQ = 1, NAT
         CALL LMSORT(NAT,NGAU,JNEQ,LQIND,LQNS,GFAC,L0PTR, &
                     L0FLD,M0PTR,M0FLD,LPPTR,LPFLD,MPPTR,MPFLD, &
                     LMXPTR,LMXFLD,GNTFLD)
         RHOATM = FOUR*PI*SQRT(VI)*R(JNEQ)*R(JNEQ)
         DO 120 JEQ = 1, MULT(JNEQ)
            call start_timer(time_pwpw)
            INDATM = INDATM + 1
!RED        DO 170 I = 1, NV+NLO
            IHELP = 0
            DO IOUTER = 0, (NV+NLO)/(NPE*BLOCKSIZE)
             IMIN = (IOUTER * NPE + MYID) * BLOCKSIZE + 1
             IMAX = MIN(NV+NLO, (IOUTER * NPE + MYID + 1) * BLOCKSIZE)
             DO I = IMIN, IMAX
               IHELP = IHELP + 1

!
!        calculate phase-factor for Alm, Blm
!
               ARGX = POS(1,INDATM)*KZZ(1,I)
               ARGY = POS(2,INDATM)*KZZ(2,I)
               ARGZ = POS(3,INDATM)*KZZ(3,I)
               PHS(I) = EXP(CIMAG*TWO*PI*(ARGX+ARGY+ARGZ))
!
!        calculate spherical-harmonics for Alm, Blm
!
               VEC(1) = XK(I)
               VEC(2) = YK(I)
               VEC(3) = ZK(I)
               CALL ROTATE (VEC,ROTIJ(1,1,INDATM),ROTV1)
               CALL ROTATE (ROTV1,ROTLOC(1,1,JNEQ),ROTV2)
               CALL YLM (ROTV2,NNSL-1,YL(0,IHELP))
!RED  170       CONTINUE
             END DO          ! DO I =
           END DO            ! DO IOUTER =
!
!        precompute CONJG(Alm), CONJG(Blm)
!
            CIL = CONE
            L0M0 = 0
            DO 140 L0 = 0, NNSL - 1
               DO 150 M0 = -L0, L0
!RED              DO 160 I = 1, NV
                  IHELP = 0
                  DO IOUTER = 0, (NV+NLO)/(NPE*BLOCKSIZE)
                   IMIN = (IOUTER * NPE + MYID) * BLOCKSIZE + 1
                   IMAX = MIN(NV+NLO, (IOUTER * NPE + MYID + 1) * BLOCKSIZE)
                   DO I = IMIN, IMAX
                     IHELP = IHELP + 1
                     if (I <= NV) then
                        ALMR(IHELP,L0M0) = RHOATM*DBLE(CIL*PHS(I)* &
                                       DCONJG(YL(L0M0,IHELP)))
                        ALMI(IHELP,L0M0) = RHOATM*DIMAG(CIL*PHS(I)* &
                                       DCONJG(YL(L0M0,IHELP)))
                        BLMR(IHELP,L0M0) = BL(I,L0,JNEQ)*ALMR(IHELP,L0M0)
                        BLMI(IHELP,L0M0) = BL(I,L0,JNEQ)*ALMI(IHELP,L0M0)
                        ALMR(IHELP,L0M0) = AL(I,L0,JNEQ)*ALMR(IHELP,L0M0)
                        ALMI(IHELP,L0M0) = AL(I,L0,JNEQ)*ALMI(IHELP,L0M0)
                        CLMR(IHELP,L0M0) = ZERO
                        CLMI(IHELP,L0M0) = ZERO
                     else
                        ALMR(IHELP,L0M0) = ZERO
                        ALMI(IHELP,L0M0) = ZERO
                        BLMR(IHELP,L0M0) = ZERO
                        BLMI(IHELP,L0M0) = ZERO
                        CLMR(IHELP,L0M0) = ZERO
                        CLMI(IHELP,L0M0) = ZERO
                     endif
                    END DO          ! DO I =
                  END DO            ! DO IOUTER =

                  IF (L0 .LE. LOMAX) THEN
!                     IF (LOOR(L0,JNEQ)) THEN
                      do jlo=1,ilo(l0,jneq)
                       CALL LOHNS (JNEQ,MULT,I1,I2,L0,jlo)   
!RED                    DO 166 I = NV+I1,NV+I2
                        IHELP = 0
                        DO IOUTER = 0, (NV+NLO)/(NPE*BLOCKSIZE)
                         IMIN = (IOUTER * NPE + MYID) * BLOCKSIZE + 1
                         IMAX = MIN(NV+NLO, (IOUTER * NPE + MYID + 1) * BLOCKSIZE)
                         DO I = IMIN, IMAX
                           IHELP = IHELP + 1
                           if ((i>=NV+I1).and.(i<=NV+I2)) then
                              ALMR(IHELP,L0M0) = RHOATM*DBLE(CIL*PHS(I)* &
                                             DCONJG(YL(L0M0,IHELP)))
                              ALMI(IHELP,L0M0) = RHOATM*DIMAG(CIL*PHS(I)* &
                                             DCONJG(YL(L0M0,IHELP)))
                              BLMR(IHELP,L0M0) = BLO(L0,jlo,JNEQ)* &
                                             ALMR(IHELP,L0M0)
                              BLMI(IHELP,L0M0) = BLO(L0,jlo,JNEQ)* &
                                             ALMI(IHELP,L0M0)
                              CLMR(IHELP,L0M0) = CLO(L0,jlo,JNEQ)* &
                                             ALMR(IHELP,L0M0)
                              CLMI(IHELP,L0M0) = CLO(L0,jlo,JNEQ)* &
                                             ALMI(IHELP,L0M0)
                              ALMR(IHELP,L0M0) = ALO(L0,jlo,JNEQ)* &
                                             ALMR(IHELP,L0M0)
                              ALMI(IHELP,L0M0) = ALO(L0,jlo,JNEQ)* &
                                             ALMI(IHELP,L0M0)
                           end if
                         END DO          ! DO I =
                        END DO            ! DO IOUTER =
                     ENDdo
                  ENDIF
                  L0M0 = L0M0 + 1
  150          CONTINUE
               CIL = CIL*CIMAG
  140       CONTINUE
           call stop_timer(time_pwpw)
           call start_timer(time_pwlo)
!
!        determine nonspherical contribution to Hamilton-matrix elements
!
!	************* initialize Loop-Index
            LINDEX=0
            DO 90 L0IDX = L0PTR(1), L0PTR(2) - 1
               L0 = L0FLD(L0IDX)    
               DO 100 M0IDX = M0PTR(L0IDX), M0PTR(L0IDX+1) - 1
!       ************ increment Loop-Index 
	          LINDEX=LINDEX+1
                  M0 = M0FLD(M0IDX)
                  L0M0 = L0*(L0 + 1) + M0
#ifdef Parallel
                  DO 80 J = 1, lld
#else
                  DO 80 J = 1, NV+NLO
#endif
                     A1R(J,LINDEX) = ALMR(J,L0M0)
                     A1I(J,LINDEX) = ALMI(J,L0M0)
                     B1R(J,LINDEX) = BLMR(J,L0M0)
                     B1I(J,LINDEX) = BLMI(J,L0M0)
                     C1R(J,LINDEX) = CLMR(J,L0M0) 
                     C1I(J,LINDEX) = CLMI(J,L0M0) 
                     A2R(J,LINDEX) = ZERO
                     A2I(J,LINDEX) = ZERO
                     B2R(J,LINDEX) = ZERO
                     B2I(J,LINDEX) = ZERO
                     C2R(J,LINDEX) = ZERO
                     C2I(J,LINDEX) = ZERO
   80             CONTINUE
                  DO 40 LPIDX = LPPTR(M0IDX), LPPTR(M0IDX+1) - 1
                     LP = LPFLD(LPIDX)
                     DO 50 MPIDX = MPPTR(LPIDX), MPPTR(LPIDX+1) - 1
                        MP = MPFLD(MPIDX)
                        LPMP = LP*(LP + 1) + MP
                        TUUR = ZERO
                        TUUI = ZERO
                        TDUR = ZERO
                        TDUI = ZERO
                        TUDR = ZERO
                        TUDI = ZERO
                        TDDR = ZERO
                        TDDI = ZERO
!
                        TUUR21 = ZERO
                        TUUI21 = ZERO
                        TUDR21 = ZERO
                        TUDI21 = ZERO
                        TUUR12 = ZERO
                        TUUI12 = ZERO
                        TDUR12 = ZERO
                        TDUI12 = ZERO
                        TUUR22 = ZERO
                        TUUI22 = ZERO
!
                        DO 10 LMXIDX = LMXPTR(MPIDX), LMXPTR(MPIDX+1)-1
                           LMX = LMXFLD(LMXIDX)
                           TUUR = TUUR + DBLE(GNTFLD(LMXIDX))* &
                                  VNS1(LP,LMX,L0,JNEQ)
                           TUUI = TUUI + DIMAG(GNTFLD(LMXIDX))* &
                                  VNS1(LP,LMX,L0,JNEQ)
!
                           TDUR = TDUR + DBLE(GNTFLD(LMXIDX))* &
                                  VNS3(L0,LMX,LP,JNEQ)
                           TDUI = TDUI + DIMAG(GNTFLD(LMXIDX))* &
                                  VNS3(L0,LMX,LP,JNEQ)
!
                           TUDR = TUDR + DBLE(GNTFLD(LMXIDX))* &
                                  VNS3(LP,LMX,L0,JNEQ)
                           TUDI = TUDI + DIMAG(GNTFLD(LMXIDX))* &
                                  VNS3(LP,LMX,L0,JNEQ)
!
                           TDDR = TDDR + DBLE(GNTFLD(LMXIDX))* &
                                  VNS2(LP,LMX,L0,JNEQ)
                           TDDI = TDDI + DIMAG(GNTFLD(LMXIDX))* &
                                  VNS2(LP,LMX,L0,JNEQ)
! 
!        part for local orbitals
!
                           LOORL0 = .FALSE.
                           IF (L0 .LE. LOMAX) THEN
                              LOORL0 = LOOR(L0,JNEQ)
                              IF (LOORL0) THEN
                                 TUUR12 = TUUR12+DBLE(GNTFLD(LMXIDX))* &
                                          VNS1LO(L0+1,LMX,LP+1,JNEQ)
                                 TUUI12 = TUUI12+DIMAG(GNTFLD(LMXIDX))* &
                                          VNS1LO(L0+1,LMX,LP+1,JNEQ)
!
                                 TDUR12 = TDUR12+DBLE(GNTFLD(LMXIDX))* &
                                          VNS2LO(L0+1,LMX,LP+1,JNEQ)
                                 TDUI12 = TDUI12+DIMAG(GNTFLD(LMXIDX))* &
                                          VNS2LO(L0+1,LMX,LP+1,JNEQ)
                              ENDIF
                           ENDIF
                           IF (LP .LE. LOMAX) THEN
                              IF (LOOR(LP,JNEQ)) THEN
                                 IF (LOORL0) THEN
                                    TUUR22 = TUUR22 + &
                                             DBLE(GNTFLD(LMXIDX))* &
                                             VNS3LO(LP+1,LMX,L0+1,JNEQ)
                                    TUUI22 = TUUI22 + &
                                             DIMAG(GNTFLD(LMXIDX))* &
                                             VNS3LO(LP+1,LMX,L0+1,JNEQ)
                                 ENDIF
                                 TUUR21 = TUUR21 + DBLE(GNTFLD(LMXIDX))* &
                                          VNS1LO(LP+1,LMX,L0+1,JNEQ)
                                 TUUI21 = TUUI21+DIMAG(GNTFLD(LMXIDX))* &
                                          VNS1LO(LP+1,LMX,L0+1,JNEQ)
!
                                 TUDR21 = TUDR21 + DBLE(GNTFLD(LMXIDX))* &
                                          VNS2LO(LP+1,LMX,L0+1,JNEQ)
                                 TUDI21 = TUDI21+DIMAG(GNTFLD(LMXIDX))* &
                                          VNS2LO(LP+1,LMX,L0+1,JNEQ)
                              ENDIF
                           ENDIF
!
!        end of local orbital part
!
   10                   CONTINUE
           call start_timer(time_hns30)
#ifdef Parallel
                        DO 30 I = 1, lld
#else
                        DO 30 I = 1, NV+NLO
#endif
                           TMP(I) = TUUR*ALMR(I,LPMP)
                           TMP(I) = TMP(I) + TUUI*ALMI(I,LPMP)
                           TMP(I) = TMP(I) + TDUR*BLMR(I,LPMP)
                           TMP(I) = TMP(I) + TDUI*BLMI(I,LPMP) 
                           TMP(I) = TMP(I) + TUUR21*CLMR(I,LPMP)
                           TMP(I) = TMP(I) + TUUI21*CLMI(I,LPMP) 
                           A2R(I,LINDEX)  = A2R(I,LINDEX) + TMP(I)
!
                           TMP(I) = TUUI*ALMR(I,LPMP)
                           TMP(I) = TMP(I) - TUUR*ALMI(I,LPMP)
                           TMP(I) = TMP(I) - TDUR*BLMI(I,LPMP)
                           TMP(I) = TMP(I) + TDUI*BLMR(I,LPMP)
                           TMP(I) = TMP(I) - TUUR21*CLMI(I,LPMP)
                           TMP(I) = TMP(I) + TUUI21*CLMR(I,LPMP)
                           A2I(I,LINDEX)  = A2I(I,LINDEX) + TMP(I)
!
                           TMP(I) = TUDR*ALMR(I,LPMP)
                           TMP(I) = TMP(I) + TUDI*ALMI(I,LPMP)
                           TMP(I) = TMP(I) + TDDR*BLMR(I,LPMP)
                           TMP(I) = TMP(I) + TDDI*BLMI(I,LPMP)
                           TMP(I) = TMP(I) + TUDR21*CLMR(I,LPMP)
                           TMP(I) = TMP(I) + TUDI21*CLMI(I,LPMP) 
                           B2R(I,LINDEX) = B2R(I,LINDEX) + TMP(I)
!
                           TMP(I) = TUDI*ALMR(I,LPMP)
                           TMP(I) = TMP(I) - TUDR*ALMI(I,LPMP)
                           TMP(I) = TMP(I) - TDDR*BLMI(I,LPMP)
                           TMP(I) = TMP(I) + TDDI*BLMR(I,LPMP)
                           TMP(I) = TMP(I) - TUDR21*CLMI(I,LPMP)
                           TMP(I) = TMP(I) + TUDI21*CLMR(I,LPMP) 
                           B2I(I,LINDEX) = B2I(I,LINDEX) + TMP(I)
!
                           TMP(I) = TUUR12*ALMR(I,LPMP)
                           TMP(I) = TMP(I) + TUUI12*ALMI(I,LPMP)
                           TMP(I) = TMP(I) + TDUR12*BLMR(I,LPMP)
                           TMP(I) = TMP(I) + TDUI12*BLMI(I,LPMP)
                           TMP(I) = TMP(I) + TUUR22*CLMR(I,LPMP)
                           TMP(I) = TMP(I) + TUUI22*CLMI(I,LPMP) 
                           C2R(I,LINDEX) = C2R(I,LINDEX) + TMP(I)
!
                           TMP(I) = TUUI12*ALMR(I,LPMP)
                           TMP(I) = TMP(I) - TUUR12*ALMI(I,LPMP)
                           TMP(I) = TMP(I) - TDUR12*BLMI(I,LPMP)
                           TMP(I) = TMP(I) + TDUI12*BLMR(I,LPMP)
                           TMP(I) = TMP(I) - TUUR22*CLMI(I,LPMP)
                           TMP(I) = TMP(I) + TUUI22*CLMR(I,LPMP) 
                           C2I(I,LINDEX) = C2I(I,LINDEX) + TMP(I)
   30                   CONTINUE
           call stop_timer(time_hns30)
   50                CONTINUE
   40             CONTINUE
!            *** speed improvement	
  100          CONTINUE
   90       CONTINUE
           call stop_timer(time_pwlo)
           call start_timer(time_lolo)
#ifdef Parallel
           CALL PAR_SYR2K(1,NV,NLO,A1R,A1I,A2R,A2I,LINDEX,H)
           CALL PAR_SYR2K(1,NV,NLO,B1R,B1I,B2R,B2I,LINDEX,H)
           CALL PAR_SYR2K(1,NV,NLO,C1R,C1I,C2R,C2I,LINDEX,H)
!          CALL PAR_SYR2K(2,NV,NLO,C1R,C1I,C2R,C2I,LINDEX,H)
#else
!_COMPLEX             A1(1:NV+NLO,1:LINDEX)=A1R(1:NV+NLO,1:LINDEX)+CIMAG*A1I(1:NV+NLO,1:LINDEX)
!_COMPLEX             A2(1:NV+NLO,1:LINDEX)=A2R(1:NV+NLO,1:LINDEX)+CIMAG*A2I(1:NV+NLO,1:LINDEX)
!               Note: not 'C', but 'T' for second matrix
!_COMPLEX             CALL ZHER2M('L','N','T',NV+NLO,LINDEX, &
!_COMPLEX                         CONE,A2,HSROWS,A1,HSROWS, &
!_COMPLEX                         CONE, HS,HSROWS)
!                  Reuse A1 and A2 for B1 and B2
!_COMPLEX             A1(1:NV+NLO,1:LINDEX)=B1R(1:NV+NLO,1:LINDEX)+CIMAG*B1I(1:NV+NLO,1:LINDEX)
!_COMPLEX             A2(1:NV+NLO,1:LINDEX)=B2R(1:NV+NLO,1:LINDEX)+CIMAG*B2I(1:NV+NLO,1:LINDEX)
!_COMPLEX             CALL ZHER2M('L','N','T',NV+NLO,LINDEX, &
!_COMPLEX                         CONE,A2,HSROWS,A1,HSROWS, &
!_COMPLEX                         CONE, HS,HSROWS)
!                  Reuse A1 and A2 for C1 and C2
!_COMPLEX             A1(NV+1:NV+NLO,1:LINDEX)=C1R(NV+1:NV+NLO,1:LINDEX)+CIMAG*C1I(NV+1:NV+NLO,1:LINDEX)
!_COMPLEX             A2(NV+1:NV+NLO,1:LINDEX)=C2R(NV+1:NV+NLO,1:LINDEX)+CIMAG*C2I(NV+1:NV+NLO,1:LINDEX)
!_COMPLEX             IF(NLO.GT.0) &
!_COMPLEX	       CALL ZHER2M('L','N','T',NLO,LINDEX, &
!_COMPLEX                         CONE,A2(NV+1,1),HSROWS,A1(NV+1,1),HSROWS, &
!_COMPLEX                         CONE, HS(NV+1,NV+1),HSROWS)
!_REAL             CALL DSYR2M('L','N','T',NV+NLO,LINDEX, &
!_REAL                         ONE,A2R,HSROWS,A1R,HSROWS, &
!_REAL                         ONE, HS,HSROWS)
!_REAL             CALL DSYR2M('L','N','T',NV+NLO,LINDEX, &
!_REAL                         -ONE,A2I,HSROWS,A1I,HSROWS, &
!_REAL                         ONE, HS,HSROWS)
!_REAL             CALL DSYR2M('L','N','T',NV+NLO,LINDEX, &
!_REAL                         ONE,B2R,HSROWS,B1R,HSROWS, &
!_REAL                         ONE, HS,HSROWS)
!_REAL             CALL DSYR2M('L','N','T',NV+NLO,LINDEX, &
!_REAL                         -ONE,B2I,HSROWS,B1I,HSROWS, &
!_REAL                         ONE, HS,HSROWS)
!_REAL             IF(NLO.GT.0) THEN
!_REAL	            CALL DSYR2M('L','N','T',NLO,LINDEX, &
!_REAL                         ONE,C2R(NV+1,1),HSROWS,C1R(NV+1,1),HSROWS, &
!_REAL                         ONE, HS(NV+1,NV+1),HSROWS)
!_REAL	            CALL DSYR2M('L','N','T',NLO,LINDEX, &
!_REAL                         -ONE,C2I(NV+1,1),HSROWS,C1I(NV+1,1),HSROWS, &
!_REAL                         ONE, HS(NV+1,NV+1),HSROWS)
!_REAL             END IF
#endif
           call stop_timer(time_lolo)
!      write(6,6000) READ_CPU_TIME(time_pwpw), READ_CPU_TIME(time_pwlo), READ_CPU_TIME(time_lolo), READ_CPU_TIME(time_hns30)
!      write(6,6001) READ_WALL_TIME(time_pwpw), READ_WALL_TIME(time_pwlo), READ_WALL_TIME(time_lolo), READ_WALL_TIME(time_hns30)
 6000 format(5X,"HNS: Loop 140: ",F8.1," Loop 90: ",F8.1," matrix update: ",F8.1, "loop 30:",f8.1,"(cpu time)")
 6001 format(5X,"HNS: Loop 140: ",F8.1," Loop 90: ",F8.1," matrix update: ",F8.1, "loop 30:",f8.1,"(wall time)")
  120    CONTINUE
  110 CONTINUE
      write(6,6001) READ_WALL_TIME(time_pwpw), READ_WALL_TIME(time_pwlo), READ_WALL_TIME(time_lolo), READ_WALL_TIME(time_hns30)

      deallocate( A1R, A1I )
      deallocate( A2R, A2I )

#ifdef Parallel
      deallocate( B1R, B1I )
      deallocate( B2R, B2I )
      deallocate( C1R, C1I )
      deallocate( C2R, C2I )
#endif

#ifndef Parallel
!_COMPLEX      deallocate( A1 )
!_COMPLEX      deallocate( A2 )
#endif

      deallocate( ALMR )
      deallocate( ALMI )
      deallocate( BLMR )
      deallocate( BLMI )
      deallocate( CLMR )
      deallocate( CLMI )
      deallocate( YL )
 999  CONTINUE
      CALL BARRIER
      RETURN
!
!        End of 'HNS'
!
      END
