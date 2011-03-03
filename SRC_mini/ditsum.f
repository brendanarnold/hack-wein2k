      SUBROUTINE DITSUM(D, G, IV, LIV, LV, P, V, X)
      implicit real*8 (a-h,o-z)
!
!  ***  PRINT ITERATION SUMMARY FOR ***SOL (VERSION 2.3)  ***
!
!  ***  PARAMETER DECLARATIONS  ***
!
      INTEGER LIV, LV, P
      INTEGER IV(LIV)
      DOUBLE PRECISION D(P), G(P), V(LV), X(P)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  ***  LOCAL VARIABLES  ***
!
      INTEGER ALG, I, IV1, M, NF, NG, OL, PU
!/6S
!     REAL MODEL1(6), MODEL2(6)
!/7S
      CHARACTER*4 MODEL1(6), MODEL2(6)
!/
      DOUBLE PRECISION NRELDF, OLDF, PRELDF, RELDF, ZERO
!
!  ***  NO EXTERNAL FUNCTIONS OR SUBROUTINES  ***
!
!  ***  SUBSCRIPTS FOR IV AND V  ***
!
      INTEGER ALGSAV, DSTNRM, F, FDIF, F0, NEEDHD, NFCALL, NFCOV, NGCOV, &
              NGCALL, NITER, NREDUC, OUTLEV, PREDUC, PRNTIT, PRUNIT, &
              RELDX, SOLPRT, STATPR, STPPAR, SUSED, X0PRT
!
!  ***  IV SUBSCRIPT VALUES  ***
!
!/6
!     DATA ALGSAV/51/, NEEDHD/36/, NFCALL/6/, NFCOV/52/, NGCALL/30/,
!    1     NGCOV/53/, NITER/31/, OUTLEV/19/, PRNTIT/39/, PRUNIT/21/,
!    2     SOLPRT/22/, STATPR/23/, SUSED/64/, X0PRT/24/
!/7
      PARAMETER (ALGSAV=51, NEEDHD=36, NFCALL=6, NFCOV=52, NGCALL=30, &
                 NGCOV=53, NITER=31, OUTLEV=19, PRNTIT=39, PRUNIT=21, &
                 SOLPRT=22, STATPR=23, SUSED=64, X0PRT=24)
!/
!
!  ***  V SUBSCRIPT VALUES  ***
!
!/6
!     DATA DSTNRM/2/, F/10/, F0/13/, FDIF/11/, NREDUC/6/, PREDUC/7/,
!    1     RELDX/17/, STPPAR/5/
!/7
      PARAMETER (DSTNRM=2, F=10, F0=13, FDIF=11, NREDUC=6, PREDUC=7, &
                 RELDX=17, STPPAR=5)
!/
!
!/6
!     DATA ZERO/0.D+0/
!/7
      PARAMETER (ZERO=0.D+0)
!/
!/6S
!     DATA MODEL1(1)/4H    /, MODEL1(2)/4H    /, MODEL1(3)/4H    /,
!    1     MODEL1(4)/4H    /, MODEL1(5)/4H  G /, MODEL1(6)/4H  S /,
!    2     MODEL2(1)/4H G  /, MODEL2(2)/4H S  /, MODEL2(3)/4HG-S /,
!    3     MODEL2(4)/4HS-G /, MODEL2(5)/4H-S-G/, MODEL2(6)/4H-G-S/
!/7S
      DATA MODEL1/'    ','    ','    ','    ','  G ','  S '/, &
           MODEL2/' G  ',' S  ','G-S ','S-G ','-S-G','-G-S'/
	SAVE MODEL1, MODEL2
!/
!
!-------------------------------  BODY  --------------------------------
!
      PU = IV(PRUNIT)
      IF (PU .EQ. 0) GO TO 999
      IV1 = IV(1)
      IF (IV1 .GT. 62) IV1 = IV1 - 51
      OL = IV(OUTLEV)
      ALG = MOD(IV(ALGSAV)-1,2) + 1
      IF (IV1 .LT. 2 .OR. IV1 .GT. 15) GO TO 370
      IF (IV1 .GE. 12) GO TO 120
      IF (IV1 .EQ. 2 .AND. IV(NITER) .EQ. 0) GO TO 390
      IF (OL .EQ. 0) GO TO 120
      IF (IV1 .GE. 10 .AND. IV(PRNTIT) .EQ. 0) GO TO 120
      IF (IV1 .GT. 2) GO TO 10
         IV(PRNTIT) = IV(PRNTIT) + 1
         IF (IV(PRNTIT) .LT. IABS(OL)) GO TO 999
 10   NF = IV(NFCALL) - IABS(IV(NFCOV))
      IV(PRNTIT) = 0
      RELDF = ZERO
      PRELDF = ZERO
      OLDF = DMAX1(DABS(V(F0)), DABS(V(F)))
      IF (OLDF .LE. ZERO) GO TO 20
         RELDF = V(FDIF) / OLDF
         PRELDF = V(PREDUC) / OLDF
 20   IF (OL .GT. 0) GO TO 60
!
!        ***  PRINT SHORT SUMMARY LINE  ***
!
         IF (IV(NEEDHD) .EQ. 1 .AND. ALG .EQ. 1) WRITE(PU,30)
 30   FORMAT(/':MIN ',10H   IT   NF,6X,1HE,7X,5HRELDF,3X,6HPRELDF,3X, &
             5HRELDX,2X,13HMODEL  STPPAR)
         IF (IV(NEEDHD) .EQ. 1 .AND. ALG .EQ. 2) WRITE(PU,40)
 40   FORMAT(/,':MIN ',10H    IT  NF,6X,1HE,7X,5HRELDF,3X,6HPRELDF,3X, &
             5HRELDX,2X,6HSTPPAR)
         IV(NEEDHD) = 0
         IF (ALG .EQ. 2) GO TO 50
         M = IV(SUSED)
         WRITE(PU,100) IV(NITER), NF, V(F), RELDF, PRELDF, V(RELDX), &
                       MODEL1(M), MODEL2(M), V(STPPAR)
         GO TO 120
!
 50      WRITE(PU,110) IV(NITER), NF, V(F), RELDF, PRELDF, V(RELDX), &
                       V(STPPAR)
         GO TO 120
!
!     ***  PRINT LONG SUMMARY LINE  ***
!
 60   IF (IV(NEEDHD) .EQ. 1 .AND. ALG .EQ. 1) WRITE(PU,70)
 70   FORMAT(/':MIN ',6HIT  NF,6X,1HE,7X,5HRELDF,2X,6HPRELDF,2X, &
             5HRELDX,2X,13HMODEL  STPPAR,2X,6HD*STEP,2X,7HNPRELDF)
      IF (IV(NEEDHD) .EQ. 1 .AND. ALG .EQ. 2) WRITE(PU,80)
 80   FORMAT(/,':MIN ',6HIT  NF,5X,1HE,7X,5HRELDF,3X,6HPRELDF,4X, &
             5HRELDX,3X,6HSTPPAR,3X,6HD*STEP,3X,7HNPRELDF)
!        write(6,*)'Alg is ',alg
      IV(NEEDHD) = 0
      NRELDF = ZERO
      IF (OLDF .GT. ZERO) NRELDF = V(NREDUC) / OLDF
      IF (ALG .EQ. 2) GO TO 90
      M = IV(SUSED)
      WRITE(PU,100) IV(NITER), NF, V(F), RELDF, PRELDF, V(RELDX), &
                   MODEL1(M), MODEL2(M), V(STPPAR), V(DSTNRM), NRELDF, &
                   V(16)
      GO TO 120
!
 90   WRITE(PU,110) IV(NITER), NF, V(F), RELDF, PRELDF, &
                   V(RELDX), V(STPPAR), V(DSTNRM), NRELDF
 100  FORMAT(':MIN ',I2,I3,D10.4,2D8.1,D8.1,A3,A4,2D8.1,2D9.2)
 110  FORMAT(':MIN ',I2,I3,D11.4,2D9.1,3D9.1,2D9.1)
!
 120  IF (IV1 .LE. 2) GO TO 999
      I = IV(STATPR)
      IF (I .EQ. (-1)) GO TO 460
      IF (I + IV1 .LT. 0) GO TO 460
      GO TO (999, 999, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, &
             330, 350, 500),  IV1
!
 130  WRITE(PU,140)
 140  FORMAT(/26H ***** X-CONVERGENCE *****)
      GO TO 430
!
 150  WRITE(PU,160)
 160  FORMAT(/42H ***** RELATIVE FUNCTION CONVERGENCE *****)
      GO TO 430
!
 170  WRITE(PU,180)
 180  FORMAT(/49H ***** X- AND RELATIVE FUNCTION CONVERGENCE *****)
      GO TO 430
!
 190  WRITE(PU,200)
 200  FORMAT(/42H ***** ABSOLUTE FUNCTION CONVERGENCE *****)
      GO TO 430
!
 210  WRITE(PU,220)
 220  FORMAT(/33H ***** SINGULAR CONVERGENCE *****)
      GO TO 430
!
 230  WRITE(PU,240)
 240  FORMAT(/30H ***** FALSE CONVERGENCE *****)
      GO TO 430
!
 250  WRITE(PU,260)
 260  FORMAT(/38H ***** FUNCTION EVALUATION LIMIT *****)
      GO TO 430
!
 270  WRITE(PU,280)
 280  FORMAT(/28H ***** ITERATION LIMIT *****)
      GO TO 430
!
 290  WRITE(PU,300)
 300  FORMAT(/18H ***** STOPX *****)
      GO TO 430
!
 310  WRITE(PU,320)
 320  FORMAT(/44H ***** INITIAL F(X) CANNOT BE COMPUTED *****)
!
      GO TO 390
!
 330  WRITE(PU,340)
 340  FORMAT(/37H ***** BAD PARAMETERS TO ASSESS *****)
      GO TO 999
!
 350  WRITE(PU,360)
 360  FORMAT(/43H ***** GRADIENT COULD NOT BE COMPUTED *****)
      IF (IV(NITER) .GT. 0) GO TO 460
      GO TO 390
!
 370  WRITE(PU,380) IV(1)
 380  FORMAT(/14H ***** IV(1) =,I5,6H *****)
      GO TO 999
!
!  ***  INITIAL CALL ON DITSUM  ***
!
 390  IF (IV(X0PRT) .NE. 0) WRITE(PU,400) (I, X(I), D(I), I = 1, P)
 400  FORMAT(/23H     I     INITIAL X(I),8X,4HD(I)//(1X,I5,D17.6,D14.3))
!     *** THE FOLLOWING ARE TO AVOID UNDEFINED VARIABLES WHEN THE
!     *** FUNCTION EVALUATION LIMIT IS 1...
      V(DSTNRM) = ZERO
      V(FDIF) = ZERO
      V(NREDUC) = ZERO
      V(PREDUC) = ZERO
      V(RELDX) = ZERO
      IF (IV1 .GE. 12) GO TO 999
      IV(NEEDHD) = 0
      IV(PRNTIT) = 0
      IF (OL .EQ. 0) GO TO 999
      IF (OL .LT. 0 .AND. ALG .EQ. 1) WRITE(PU,30)
      IF (OL .LT. 0 .AND. ALG .EQ. 2) WRITE(PU,40)
      IF (OL .GT. 0 .AND. ALG .EQ. 1) WRITE(PU,70)
      IF (OL .GT. 0 .AND. ALG .EQ. 2) WRITE(PU,80)
      IF (ALG .EQ. 1) WRITE(PU,410) IV(NFCALL), V(F)
      IF (ALG .EQ. 2) WRITE(PU,420) IV(NFCALL), V(F)
 410  FORMAT(/4H   0,I5,D10.4)
 420  FORMAT(/4H   0,I5,D11.4)
      GO TO 999
!
!  ***  PRINT VARIOUS INFORMATION REQUESTED ON SOLUTION  ***
!
 430  IV(NEEDHD) = 1
      IF (IV(STATPR) .LE. 0) GO TO 460
         OLDF = DMAX1(DABS(V(F0)), DABS(V(F)))
         PRELDF = ZERO
         NRELDF = ZERO
         IF (OLDF .LE. ZERO) GO TO 440
              PRELDF = V(PREDUC) / OLDF
              NRELDF = V(NREDUC) / OLDF
 440     NF = IV(NFCALL) - IV(NFCOV)
         NG = IV(NGCALL) - IV(NGCOV)
         WRITE(PU,450) V(F), V(RELDX), NF, NG, PRELDF, NRELDF
 450  FORMAT(/9H FUNCTION,D17.6,8H   RELDX,D17.3/12H FUNC. EVALS, &
         I8,9X,11HGRAD. EVALS,I8/7H PRELDF,D16.3,6X,7HNPRELDF,D15.3)
!
 460  IF (IV(SOLPRT) .EQ. 0) GO TO 999
         IV(NEEDHD) = 1
         IF (IV(ALGSAV) .GT. 2) GO TO 999
         WRITE(PU,470)
 470  FORMAT(/':MIN ',20H   I      FINAL X(I),8X,4HD(I),10X,4HG(I)/)
         DO 480 I = 1, P
 480          WRITE(PU,490) I, X(I), D(I), G(I)
 490     FORMAT(':Min ',1X,I5,D16.6,2D14.3)
      GO TO 999
!
 500  WRITE(PU,510)
 510  FORMAT(/24H INCONSISTENT DIMENSIONS)
 999  RETURN
!  ***  LAST CARD OF DITSUM FOLLOWS  ***
      END
