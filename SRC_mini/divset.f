      SUBROUTINE DIVSET(ALG, IV, LIV, LV, V)
      implicit real*8 (a-h,o-z)
!
!  ***  SUPPLY ***SOL (VERSION 2.3) DEFAULT VALUES TO IV AND V  ***
!
!  ***  ALG = 1 MEANS REGRESSION CONSTANTS.
!  ***  ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS.
!
      INTEGER LIV, LV
      INTEGER ALG, IV(LIV)
      DOUBLE PRECISION V(LV)
!
      EXTERNAL I7MDCN,DV7DFL
      INTEGER I7MDCN
! I7MDCN... RETURNS MACHINE-DEPENDENT INTEGER CONSTANTS.
! DV7DFL.... PROVIDES DEFAULT VALUES TO V.
!
      INTEGER ALG1, MIV, MV
      INTEGER MINIV(4), MINV(4)
!
!  ***  SUBSCRIPTS FOR IV  ***
!
      INTEGER ALGSAV, COVPRT, COVREQ, DRADPR, DTYPE, HC, IERR, INITH, &
              INITS, IPIVOT, IVNEED, LASTIV, LASTV, LMAT, MXFCAL, &
              MXITER, NFCOV, NGCOV, NVDFLT, NVSAVE, OUTLEV, PARPRT, &
              PARSAV, PERM, PRUNIT, QRTYP, RDREQ, RMAT, SOLPRT, STATPR, &
              VNEED, VSAVE, X0PRT
!
!  ***  IV SUBSCRIPT VALUES  ***
!
!/6
!     DATA ALGSAV/51/, COVPRT/14/, COVREQ/15/, DRADPR/101/, DTYPE/16/,
!    1     HC/71/, IERR/75/, INITH/25/, INITS/25/, IPIVOT/76/,
!    2     IVNEED/3/, LASTIV/44/, LASTV/45/, LMAT/42/, MXFCAL/17/,
!    3     MXITER/18/, NFCOV/52/, NGCOV/53/, NVDFLT/50/, NVSAVE/9/,
!    4     OUTLEV/19/, PARPRT/20/, PARSAV/49/, PERM/58/, PRUNIT/21/,
!    5     QRTYP/80/, RDREQ/57/, RMAT/78/, SOLPRT/22/, STATPR/23/,
!    6     VNEED/4/, VSAVE/60/, X0PRT/24/
!/7
      PARAMETER (ALGSAV=51, COVPRT=14, COVREQ=15, DRADPR=101, DTYPE=16, &
                 HC=71, IERR=75, INITH=25, INITS=25, IPIVOT=76, &
                 IVNEED=3, LASTIV=44, LASTV=45, LMAT=42, MXFCAL=17, &
                 MXITER=18, NFCOV=52, NGCOV=53, NVDFLT=50, NVSAVE=9, &
                 OUTLEV=19, PARPRT=20, PARSAV=49, PERM=58, PRUNIT=21, &
                 QRTYP=80, RDREQ=57, RMAT=78, SOLPRT=22, STATPR=23, &
                 VNEED=4, VSAVE=60, X0PRT=24)
!/
      DATA MINIV(1)/82/, MINIV(2)/59/, MINIV(3)/103/, MINIV(4)/103/, &
           MINV(1)/98/, MINV(2)/71/, MINV(3)/101/, MINV(4)/85/
      save MINV
!
!-------------------------------  BODY  --------------------------------
!
      IF (PRUNIT .LE. LIV) IV(PRUNIT) = I7MDCN(1)
      IF (ALGSAV .LE. LIV) IV(ALGSAV) = ALG
      IF (ALG .LT. 1 .OR. ALG .GT. 4) GO TO 40
      MIV = MINIV(ALG)
      IF (LIV .LT. MIV) GO TO 20
      MV = MINV(ALG)
      IF (LV .LT. MV) GO TO 30
      ALG1 = MOD(ALG-1,2) + 1
      CALL DV7DFL(ALG1, LV, V)
      IV(1) = 12
      IF (ALG .GT. 2) IV(DRADPR) = 1
      IV(IVNEED) = 0
      IV(LASTIV) = MIV
      IV(LASTV) = MV
      IV(LMAT) = MV + 1
      IV(MXFCAL) = 200
      IV(MXITER) = 150
      IV(OUTLEV) = 1
      IV(PARPRT) = 1
      IV(PERM) = MIV + 1
      IV(SOLPRT) = 1
      IV(STATPR) = 1
      IV(VNEED) = 0
      IV(X0PRT) = 1
!
      IF (ALG1 .GE. 2) GO TO 10
!
!  ***  REGRESSION  VALUES
!
      IV(COVPRT) = 3
      IV(COVREQ) = 1
      IV(DTYPE) = 1
      IV(HC) = 0
      IV(IERR) = 0
      IV(INITS) = 0
      IV(IPIVOT) = 0
      IV(NVDFLT) = 32
      IV(VSAVE) = 58
      IF (ALG .GT. 2) IV(VSAVE) = IV(VSAVE) + 3
      IV(PARSAV) = IV(VSAVE) + NVSAVE
      IV(QRTYP) = 1
      IV(RDREQ) = 3
      IV(RMAT) = 0
      GO TO 999
!
!  ***  GENERAL OPTIMIZATION VALUES
!
 10   IV(DTYPE) = 0
      IV(INITH) = 1
      IV(NFCOV) = 0
      IV(NGCOV) = 0
      IV(NVDFLT) = 25
      IV(PARSAV) = 47
      IF (ALG .GT. 2) IV(PARSAV) = 61
      GO TO 999
!
 20   IV(1) = 15
      GO TO 999
!
 30   IV(1) = 16
      GO TO 999
!
 40   IV(1) = 67
!
 999  RETURN
!  ***  LAST CARD OF DIVSET FOLLOWS  ***
      END