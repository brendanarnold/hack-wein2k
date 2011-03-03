      SUBROUTINE DRMNG(D, FX, G, IV, LIV, LV, N, V, X,HESS)
      implicit real*8 (a-h,o-z)
!
!  ***  CARRY OUT  DMNG (UNCONSTRAINED MINIMIZATION) ITERATIONS, USING
!  ***  DOUBLE-DOGLEG/BFGS STEPS.
!
!  ***  PARAMETER DECLARATIONS  ***
!
      INTEGER LIV, LV, N
      INTEGER IV(LIV)
      DOUBLE PRECISION D(N), FX, G(N), V(LV), X(N),HESS(*)
!
!--------------------------  PARAMETER USAGE  --------------------------
!
! D.... SCALE VECTOR.
! FX... FUNCTION VALUE.
! G.... GRADIENT VECTOR.
! IV... INTEGER VALUE ARRAY.
! LIV.. LENGTH OF IV (AT LEAST 60).
! LV... LENGTH OF V (AT LEAST 71 + N*(N+13)/2).
! N.... NUMBER OF VARIABLES (COMPONENTS IN X AND G).
! V.... FLOATING-POINT VALUE ARRAY.
! X.... VECTOR OF PARAMETERS TO BE OPTIMIZED.
!
!  ***  DISCUSSION  ***
!
!        PARAMETERS IV, N, V, AND X ARE THE SAME AS THE CORRESPONDING
!     ONES TO  DMNG (WHICH SEE), EXCEPT THAT V CAN BE SHORTER (SINCE
!     THE PART OF V THAT  DMNG USES FOR STORING G IS NOT NEEDED).
!     MOREOVER, COMPARED WITH  DMNG, IV(1) MAY HAVE THE TWO ADDITIONAL
!     OUTPUT VALUES 1 AND 2, WHICH ARE EXPLAINED BELOW, AS IS THE USE
!     OF IV(TOOBIG) AND IV(NFGCAL).  THE VALUE IV(G), WHICH IS AN
!     OUTPUT VALUE FROM  DMNG (AND  DMNF), IS NOT REFERENCED BY
!     DRMNG OR THE SUBROUTINES IT CALLS.
!        FX AND G NEED NOT HAVE BEEN INITIALIZED WHEN DRMNG IS CALLED
!     WITH IV(1) = 12, 13, OR 14.
!
! IV(1) = 1 MEANS THE CALLER SHOULD SET FX TO F(X), THE FUNCTION VALUE
!             AT X, AND CALL DRMNG AGAIN, HAVING CHANGED NONE OF THE
!             OTHER PARAMETERS.  AN EXCEPTION OCCURS IF F(X) CANNOT BE
!             (E.G. IF OVERFLOW WOULD OCCUR), WHICH MAY HAPPEN BECAUSE
!             OF AN OVERSIZED STEP.  IN THIS CASE THE CALLER SHOULD SET
!             IV(TOOBIG) = IV(2) TO 1, WHICH WILL CAUSE DRMNG TO IG-
!             NORE FX AND TRY A SMALLER STEP.  THE PARAMETER NF THAT
!              DMNG PASSES TO CALCF (FOR POSSIBLE USE BY CALCG) IS A
!             COPY OF IV(NFCALL) = IV(6).
! IV(1) = 2 MEANS THE CALLER SHOULD SET G TO G(X), THE GRADIENT VECTOR
!             OF F AT X, AND CALL DRMNG AGAIN, HAVING CHANGED NONE OF
!             THE OTHER PARAMETERS EXCEPT POSSIBLY THE SCALE VECTOR D
!             WHEN IV(DTYPE) = 0.  THE PARAMETER NF THAT  DMNG PASSES
!             TO CALCG IS IV(NFGCAL) = IV(7).  IF G(X) CANNOT BE
!             EVALUATED, THEN THE CALLER MAY SET IV(TOOBIG) TO 0, IN
!             WHICH CASE DRMNG WILL RETURN WITH IV(1) = 65.
!.
!  ***  GENERAL  ***
!
!     CODED BY DAVID M. GAY (DECEMBER 1979).  REVISED SEPT. 1982.
!     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH SUPPORTED
!     IN PART BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS
!     MCS-7600324 AND MCS-7906671.
!
!        (SEE  DMNG FOR REFERENCES.)
!
!+++++++++++++++++++++++++++  DECLARATIONS  ++++++++++++++++++++++++++++
!
!  ***  LOCAL VARIABLES  ***
!
      INTEGER DG1, DUMMY, G01, I, K, L, LSTGST, NWTST1, RSTRST, STEP1, &
              TEMP1, W, X01, Z
      DOUBLE PRECISION T
!
!     ***  CONSTANTS  ***
!
      DOUBLE PRECISION HALF, NEGONE, ONE, ONEP2, ZERO
!
!  ***  NO INTRINSIC FUNCTIONS  ***
!
!  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***
!
      LOGICAL STOPX
      DOUBLE PRECISION DD7TPR, DRLDST, DV2NRM
      EXTERNAL DA7SST,DD7DOG,DIVSET, DD7TPR,DITSUM, DL7ITV, DL7IVM, &
               DL7TVM, DL7UPD,DL7VML,DPARCK, DRLDST, STOPX,DV2AXY, &
              DV7CPY, DV7SCP, DV7VMP, DV2NRM, DW7ZBF
!
! DA7SST.... ASSESSES CANDIDATE STEP.
! DD7DOG.... COMPUTES DOUBLE-DOGLEG (CANDIDATE) STEP.
! DIVSET.... SUPPLIES DEFAULT IV AND V INPUT COMPONENTS.
! DD7TPR... RETURNS INNER PRODUCT OF TWO VECTORS.
! DITSUM.... PRINTS ITERATION SUMMARY AND INFO ON INITIAL AND FINAL X.
! DL7ITV... MULTIPLIES INVERSE TRANSPOSE OF LOWER TRIANGLE TIMES VECTOR.
! DL7IVM... MULTIPLIES INVERSE OF LOWER TRIANGLE TIMES VECTOR.
! DL7TVM... MULTIPLIES TRANSPOSE OF LOWER TRIANGLE TIMES VECTOR.
! LUPDT.... UPDATES CHOLESKY FACTOR OF HESSIAN APPROXIMATION.
! DL7VML.... MULTIPLIES LOWER TRIANGLE TIMES VECTOR.
! DPARCK.... CHECKS VALIDITY OF INPUT IV AND V VALUES.
! DRLDST... COMPUTES V(RELDX) = RELATIVE STEP SIZE.
! STOPX.... RETURNS .TRUE. IF THE BREAK KEY HAS BEEN PRESSED.
! DV2AXY.... COMPUTES SCALAR TIMES ONE VECTOR PLUS ANOTHER.
! DV7CPY.... COPIES ONE VECTOR TO ANOTHER.
! DV7SCP... SETS ALL ELEMENTS OF A VECTOR TO A SCALAR.
! DV7VMP... MULTIPLIES VECTOR BY VECTOR RAISED TO POWER (COMPONENTWISE).
! DV2NRM... RETURNS THE 2-NORM OF A VECTOR.
! DW7ZBF... COMPUTES W AND Z FOR DL7UPD CORRESPONDING TO BFGS UPDATE.
!
!  ***  SUBSCRIPTS FOR IV AND V  ***
!
      INTEGER CNVCOD, DG, DGNORM, DINIT, DSTNRM, DST0, F, F0, FDIF, &
              GTHG, GTSTEP, G0, INCFAC, INITH, IRC, KAGQT, LMAT, LMAX0, &
              LMAXS, MODE, MODEL, MXFCAL, MXITER, NEXTV, NFCALL, NFGCAL, &
              NGCALL, NITER, NREDUC, NWTSTP, PREDUC, RADFAC, RADINC, &
              RADIUS, RAD0, RELDX, RESTOR, STEP, STGLIM, STLSTG, TOOBIG, &
              TUNER4, TUNER5, VNEED, XIRC, X0
!
!  ***  IV SUBSCRIPT VALUES  ***
!
!/6
!     DATA CNVCOD/55/, DG/37/, G0/48/, INITH/25/, IRC/29/, KAGQT/33/,
!    1     MODE/35/, MODEL/5/, MXFCAL/17/, MXITER/18/, NFCALL/6/,
!    2     NFGCAL/7/, NGCALL/30/, NITER/31/, NWTSTP/34/, RADINC/8/,
!    3     RESTOR/9/, STEP/40/, STGLIM/11/, STLSTG/41/, TOOBIG/2/,
!    4     VNEED/4/, XIRC/13/, X0/43/
!/7
      PARAMETER (CNVCOD=55, DG=37, G0=48, INITH=25, IRC=29, KAGQT=33, &
                 MODE=35, MODEL=5, MXFCAL=17, MXITER=18, NFCALL=6, &
                 NFGCAL=7, NGCALL=30, NITER=31, NWTSTP=34, RADINC=8, &
                 RESTOR=9, STEP=40, STGLIM=11, STLSTG=41, TOOBIG=2, &
                 VNEED=4, XIRC=13, X0=43)
!/
!
!  ***  V SUBSCRIPT VALUES  ***
!
!/6
!     DATA DGNORM/1/, DINIT/38/, DSTNRM/2/, DST0/3/, F/10/, F0/13/,
!    1     FDIF/11/, GTHG/44/, GTSTEP/4/, INCFAC/23/, LMAT/42/,
!    2     LMAX0/35/, LMAXS/36/, NEXTV/47/, NREDUC/6/, PREDUC/7/,
!    3     RADFAC/16/, RADIUS/8/, RAD0/9/, RELDX/17/, TUNER4/29/,
!    4     TUNER5/30/
!/7
      PARAMETER (DGNORM=1, DINIT=38, DSTNRM=2, DST0=3, F=10, F0=13, &
                 FDIF=11, GTHG=44, GTSTEP=4, INCFAC=23, LMAT=42, &
                 LMAX0=35, LMAXS=36, NEXTV=47, NREDUC=6, PREDUC=7, &
                 RADFAC=16, RADIUS=8, RAD0=9, RELDX=17, TUNER4=29, &
                 TUNER5=30)
!/
!
!/6
!     DATA HALF/0.5D+0/, NEGONE/-1.D+0/, ONE/1.D+0/, ONEP2/1.2D+0/,
!    1     ZERO/0.D+0/
!/7
      PARAMETER (HALF=0.5D+0, NEGONE=-1.D+0, ONE=1.D+0, ONEP2=1.2D+0, &
                 ZERO=0.D+0)
!/
!
!+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
!
!        write(6,*)'Deb line 1',V(iv(LMAT)),iv(LMAT)
      I = IV(1)
!        write(6,*)'Deb jump ',i
      IF (I .EQ. 1) GO TO 50
      IF (I .EQ. 2) GO TO 60
!
!  ***  CHECK VALIDITY OF IV AND V INPUT VALUES  ***
!
      IF (IV(1) .EQ. 0) CALL DIVSET(2, IV, LIV, LV, V)
      IF (IV(1) .EQ. 12 .OR. IV(1) .EQ. 13) &
           IV(VNEED) = IV(VNEED) + N*(N+13)/2
!       These two lines seem to be needed -- why?
        tmp=V(iv(LMAT))
      CALL DPARCK(2, D, IV, LIV, LV, N, V)
        V(iv(LMAT))=tmp
!        write(6,*)'Debug allocated ',IV(LMAT),V(IV(LMAT))
      I = IV(1) - 2
!        write(6,*)'Debug jump 2',i
      IF (I .GT. 12) GO TO 999
      GO TO (190, 190, 190, 190, 190, 190, 120, 90, 120, 10, 10, 20), I
!
!  ***  STORAGE ALLOCATION  ***
!
 10   L = IV(LMAT)
      IV(X0) = L + N*(N+1)/2
      IV(STEP) = IV(X0) + N
      IV(STLSTG) = IV(STEP) + N
      IV(G0) = IV(STLSTG) + N
      IV(NWTSTP) = IV(G0) + N
      IV(DG) = IV(NWTSTP) + N
      IV(NEXTV) = IV(DG) + N
      IF (IV(1) .NE. 13) GO TO 20
         IV(1) = 14
         GO TO 999
!
!  ***  INITIALIZATION  ***
!
 20   IV(NITER) = 0
      IV(NFCALL) = 1
      IV(NGCALL) = 1
      IV(NFGCAL) = 1
      IV(MODE) = -1
      IV(MODEL) = 1
      IV(STGLIM) = 1
      IV(TOOBIG) = 0
      IV(CNVCOD) = 0
      IV(RADINC) = 0
      V(RAD0) = ZERO
!        write(6,*)'Debug odd 1',V(72)
      IF (V(DINIT) .GE. ZERO) CALL DV7SCP(N, D, V(DINIT))
!        write(6,*)'Debug odd 1',V(72)
      IF (IV(INITH) .NE. 1) GO TO 40
!        write(6,*)'Debug X setting hessian to D'
!
!     ***  SET THE INITIAL HESSIAN APPROXIMATION TO DIAG(D)**-2  ***
!
         L = IV(LMAT)
         CALL DV7SCP(N*(N+1)/2, V(L), ZERO)
         K = L - 1
         KK=0
         DO 30 I = 1, N
              K = K + I
              T = D(I)
              IF (T .LE. ZERO) T = ONE
              KK = KK + I
              HESS(KK) = T
              V(K) = T
 30           CONTINUE
!        do I = 1, (N*(N+1))/2
!                HESS(I)=V(L+I-1)
!        enddo
!
!  ***  COMPUTE INITIAL FUNCTION VALUE  ***
!
 40   IV(1) = 1
        L=iv(LMAT)
!        write(6,*)'Deb 1 ',v(L)
      GO TO 999
!
 50   V(F) = FX
      IF (IV(MODE) .GE. 0) GO TO 190
      V(F0) = FX
      IV(1) = 2
      IF (IV(TOOBIG) .EQ. 0) GO TO 999
         IV(1) = 63
         GO TO 350
!
!  ***  MAKE SURE GRADIENT COULD BE COMPUTED  ***
!
 60   IF (IV(TOOBIG) .EQ. 0) GO TO 70
         IV(1) = 65
         GO TO 350
!
 70   DG1 = IV(DG)
      CALL DV7VMP(N, V(DG1), G, D, -1)
      V(DGNORM) = DV2NRM(N, V(DG1))
!
      IF (IV(CNVCOD) .NE. 0) GO TO 340
      IF (IV(MODE) .EQ. 0) GO TO 300
!
!  ***  ALLOW FIRST STEP TO HAVE SCALED 2-NORM AT MOST V(LMAX0)  ***
!
      V(RADIUS) = V(LMAX0)
!
      IV(MODE) = 0
!
!
!-----------------------------  MAIN LOOP  -----------------------------
!
!
!  ***  PRINT ITERATION SUMMARY, CHECK ITERATION LIMIT  ***
!
 80   CALL DITSUM(D, G, IV, LIV, LV, N, V, X)
 90   K = IV(NITER)
!        write(6,*)'Deb 2',V(IV(LMAT))
      IF (K .LT. IV(MXITER)) GO TO 100
         IV(1) = 10
         GO TO 350
!
!  ***  UPDATE RADIUS  ***
!
 100  IV(NITER) = K + 1
      IF (K .GT. 0) V(RADIUS) = V(RADFAC) * V(DSTNRM)
!
!  ***  INITIALIZE FOR START OF NEXT ITERATION  ***
!
      G01 = IV(G0)
      X01 = IV(X0)
      V(F0) = V(F)
      IV(IRC) = 4
      IV(KAGQT) = -1
!
!     ***  COPY X TO X0, G TO G0  ***
!
      CALL DV7CPY(N, V(X01), X)
      CALL DV7CPY(N, V(G01), G)
!
!  ***  CHECK STOPX AND FUNCTION EVALUATION LIMIT  ***
!
 110  IF (.NOT. STOPX(DUMMY)) GO TO 130
         IV(1) = 11
         GO TO 140
!
!     ***  COME HERE WHEN RESTARTING AFTER FUNC. EVAL. LIMIT OR STOPX.
!
 120  IF (V(F) .GE. V(F0)) GO TO 130
         V(RADFAC) = ONE
         K = IV(NITER)
         GO TO 100
!
 130  IF (IV(NFCALL) .LT. IV(MXFCAL)) GO TO 150
         IV(1) = 9
 140     IF (V(F) .GE. V(F0)) GO TO 350
!
!        ***  IN CASE OF STOPX OR FUNCTION EVALUATION LIMIT WITH
!        ***  IMPROVED V(F), EVALUATE THE GRADIENT AT X.
!
              IV(CNVCOD) = IV(1)
              GO TO 290
!
!. . . . . . . . . . . . .  COMPUTE CANDIDATE STEP  . . . . . . . . . .
!
 150  STEP1 = IV(STEP)
      DG1 = IV(DG)
      NWTST1 = IV(NWTSTP)
      IF (IV(KAGQT) .GE. 0) GO TO 160
         L = IV(LMAT)
!        write(6,*)'Deb 3 ',V(L),L
         CALL DL7IVM(N, V(NWTST1), V(L), G)
!        write(6,*)'deb 4 ',V(L),L
         V(NREDUC) = HALF * DD7TPR(N, V(NWTST1), V(NWTST1))
!        write(6,*)'Deb 5 ',V(L),L
         CALL DL7ITV(N, V(NWTST1), V(L), V(NWTST1))
!        write(6,*)'Deb 6 ',V(L),V(NWTST1),D(1),L
         CALL DV7VMP(N, V(STEP1), V(NWTST1), D, 1)
         V(DST0) = DV2NRM(N, V(STEP1))
!        write(6,*)'Deb 7 ',V(DSTO),N
         CALL DV7VMP(N, V(DG1), V(DG1), D, -1)
         CALL DL7TVM(N, V(STEP1), V(L), V(DG1))
         V(GTHG) = DV2NRM(N, V(STEP1))
         IV(KAGQT) = 0
 160  continue
!      print*,'v(NWTST1)',IV(KAGQT),n,v(NWTST1:NWTST1+n-1)
      CALL DD7DOG(V(DG1), LV, N, V(NWTST1), V(STEP1), V(1))
      IF (IV(IRC) .NE. 6) GO TO 170
         IF (IV(RESTOR) .NE. 2) GO TO 190
         RSTRST = 2
         GO TO 200
!
!  ***  CHECK WHETHER EVALUATING F(X0 + STEP) LOOKS WORTHWHILE  ***
!
 170  IV(TOOBIG) = 0
      IF (V(DSTNRM) .LE. ZERO) GO TO 190
      IF (IV(IRC) .NE. 5) GO TO 180
      IF (V(RADFAC) .LE. ONE) GO TO 180
      IF (V(PREDUC) .GT. ONEP2 * V(FDIF)) GO TO 180
         IF (IV(RESTOR) .NE. 2) GO TO 190
         RSTRST = 0
         GO TO 200
!
!  ***  COMPUTE F(X0 + STEP)  ***
!
 180  X01 = IV(X0)
      STEP1 = IV(STEP)
      CALL DV2AXY(N, X, ONE, V(STEP1), V(X01))
      IV(NFCALL) = IV(NFCALL) + 1
      IV(1) = 1
      GO TO 999
!
!. . . . . . . . . . . . .  ASSESS CANDIDATE STEP  . . . . . . . . . . .
!
 190  RSTRST = 3
 200  X01 = IV(X0)
      V(RELDX) = DRLDST(N, D, X, V(X01))
      CALL DA7SST(IV, LIV, LV, V)
      STEP1 = IV(STEP)
      LSTGST = IV(STLSTG)
      I = IV(RESTOR) + 1
      GO TO (240, 210, 220, 230), I
 210  CALL DV7CPY(N, X, V(X01))
      GO TO 240
 220   CALL DV7CPY(N, V(LSTGST), V(STEP1))
       GO TO 240
 230     CALL DV7CPY(N, V(STEP1), V(LSTGST))
         CALL DV2AXY(N, X, ONE, V(STEP1), V(X01))
         V(RELDX) = DRLDST(N, D, X, V(X01))
         IV(RESTOR) = RSTRST
!
 240  K = IV(IRC)
      GO TO (250,280,280,280,250,260,270,270,270,270,270,270,330,300), K
!
!     ***  RECOMPUTE STEP WITH CHANGED RADIUS  ***
!
 250     V(RADIUS) = V(RADFAC) * V(DSTNRM)
         GO TO 110
!
!  ***  COMPUTE STEP OF LENGTH V(LMAXS) FOR SINGULAR CONVERGENCE TEST.
!
 260  V(RADIUS) = V(LMAXS)
      GO TO 150
!
!  ***  CONVERGENCE OR FALSE CONVERGENCE  ***
!
 270  IV(CNVCOD) = K - 4
      IF (V(F) .GE. V(F0)) GO TO 340
         IF (IV(XIRC) .EQ. 14) GO TO 340
              IV(XIRC) = 14
!
!. . . . . . . . . . . .  PROCESS ACCEPTABLE STEP  . . . . . . . . . . .
!
 280  IF (IV(IRC) .NE. 3) GO TO 290
         STEP1 = IV(STEP)
         TEMP1 = IV(STLSTG)
!
!     ***  SET  TEMP1 = HESSIAN * STEP  FOR USE IN GRADIENT TESTS  ***
!
         L = IV(LMAT)
         CALL DL7TVM(N, V(TEMP1), V(L), V(STEP1))
         CALL DL7VML(N, V(TEMP1), V(L), V(TEMP1))
!
!  ***  COMPUTE GRADIENT  ***
!
 290  IV(NGCALL) = IV(NGCALL) + 1
      IV(1) = 2
      GO TO 999
!
!  ***  INITIALIZATIONS -- G0 = G - G0, ETC.  ***
!
 300  G01 = IV(G0)
      CALL DV2AXY(N, V(G01), NEGONE, V(G01), G)
      STEP1 = IV(STEP)
      TEMP1 = IV(STLSTG)
      IF (IV(IRC) .NE. 3) GO TO 320
!
!  ***  SET V(RADFAC) BY GRADIENT TESTS  ***
!
!     ***  SET  TEMP1 = DIAG(D)**-1 * (HESSIAN*STEP + (G(X0)-G(X)))  ***
!
         CALL DV2AXY(N, V(TEMP1), NEGONE, V(G01), V(TEMP1))
         CALL DV7VMP(N, V(TEMP1), V(TEMP1), D, -1)
!
!        ***  DO GRADIENT TESTS  ***
!
         IF (DV2NRM(N, V(TEMP1)) .LE. V(DGNORM) * V(TUNER4)) &
                        GO TO 310
              IF (DD7TPR(N, G, V(STEP1)) &
                        .GE. V(GTSTEP) * V(TUNER5))  GO TO 320
 310               V(RADFAC) = V(INCFAC)
!
!  ***  UPDATE H, LOOP  ***
!
 320  W = IV(NWTSTP)
      Z = IV(X0)
      L = IV(LMAT)
      CALL DW7ZBF(V(L), N, V(STEP1), V(W), V(G01), V(Z))
!
!     ** USE THE N-VECTORS STARTING AT V(STEP1) AND V(G01) FOR SCRATCH..
      CALL DL7UPD(V(TEMP1), V(STEP1), V(L), V(G01), V(L), N, V(W), V(Z))
      DO J=1,(N*(N+1))/2
        HESS(J)=V(J+L-1)
      ENDDO
      HESS(1)=V(L)
!      open(unit=81,file='.min.tmp')
!      write(81,*)(HESS(J),J=1,(N*(N+1))/2)
!      close(unit=81)
      IV(1) = 2
      GO TO 80
!
!. . . . . . . . . . . . . .  MISC. DETAILS  . . . . . . . . . . . . . .
!
!  ***  BAD PARAMETERS TO ASSESS  ***
!
 330  IV(1) = 64
      GO TO 350
!
!  ***  PRINT SUMMARY OF FINAL ITERATION AND OTHER REQUESTED ITEMS  ***
!
 340  IV(1) = IV(CNVCOD)
      IV(CNVCOD) = 0
 350  CALL DITSUM(D, G, IV, LIV, LV, N, V, X)
!
 999  RETURN
!
!  ***  LAST LINE OF DRMNG FOLLOWS  ***
      END
