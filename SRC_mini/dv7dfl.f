      SUBROUTINE DV7DFL(ALG, LV, V)
      implicit real*8 (a-h,o-z)
!
!  ***  SUPPLY ***SOL (VERSION 2.3) DEFAULT VALUES TO V  ***
!
!  ***  ALG = 1 MEANS REGRESSION CONSTANTS.
!  ***  ALG = 2 MEANS GENERAL UNCONSTRAINED OPTIMIZATION CONSTANTS.
!
      INTEGER ALG, LV
      DOUBLE PRECISION V(LV)
!
      EXTERNAL DR7MDC
      DOUBLE PRECISION DR7MDC
! DR7MDC... RETURNS MACHINE-DEPENDENT CONSTANTS
!
      DOUBLE PRECISION MACHEP, MEPCRT, ONE, SQTEPS, THREE
!
!  ***  SUBSCRIPTS FOR V  ***
!
      INTEGER AFCTOL, BIAS, COSMIN, DECFAC, DELTA0, DFAC, DINIT, DLTFDC, &
              DLTFDJ, DTINIT, D0INIT, EPSLON, ETA0, FUZZ, HUBERC, &
              INCFAC, LMAX0, LMAXS, PHMNFC, PHMXFC, RDFCMN, RDFCMX, &
              RFCTOL, RLIMIT, RSPTOL, SCTOL, SIGMIN, TUNER1, TUNER2, &
              TUNER3, TUNER4, TUNER5, XCTOL, XFTOL
!
!/6
!     DATA ONE/1.D+0/, THREE/3.D+0/
!/7
      PARAMETER (ONE=1.D+0, THREE=3.D+0)
!/
!
!  ***  V SUBSCRIPT VALUES  ***
!
!/6
!     DATA AFCTOL/31/, BIAS/43/, COSMIN/47/, DECFAC/22/, DELTA0/44/,
!    1     DFAC/41/, DINIT/38/, DLTFDC/42/, DLTFDJ/43/, DTINIT/39/,
!    2     D0INIT/40/, EPSLON/19/, ETA0/42/, FUZZ/45/, HUBERC/48/,
!    3     INCFAC/23/, LMAX0/35/, LMAXS/36/, PHMNFC/20/, PHMXFC/21/,
!    4     RDFCMN/24/, RDFCMX/25/, RFCTOL/32/, RLIMIT/46/, RSPTOL/49/,
!    5     SCTOL/37/, SIGMIN/50/, TUNER1/26/, TUNER2/27/, TUNER3/28/,
!    6     TUNER4/29/, TUNER5/30/, XCTOL/33/, XFTOL/34/
!/7
      PARAMETER (AFCTOL=31, BIAS=43, COSMIN=47, DECFAC=22, DELTA0=44, &
                 DFAC=41, DINIT=38, DLTFDC=42, DLTFDJ=43, DTINIT=39, &
                 D0INIT=40, EPSLON=19, ETA0=42, FUZZ=45, HUBERC=48, &
                 INCFAC=23, LMAX0=35, LMAXS=36, PHMNFC=20, PHMXFC=21, &
                 RDFCMN=24, RDFCMX=25, RFCTOL=32, RLIMIT=46, RSPTOL=49, &
                 SCTOL=37, SIGMIN=50, TUNER1=26, TUNER2=27, TUNER3=28, &
                 TUNER4=29, TUNER5=30, XCTOL=33, XFTOL=34)
!/
!
!-------------------------------  BODY  --------------------------------
!
      MACHEP = DR7MDC(3)
      V(AFCTOL) = 1.D-20
      IF (MACHEP .GT. 1.D-10) V(AFCTOL) = MACHEP**2
      V(DECFAC) = 0.5D+0
      SQTEPS = DR7MDC(4)
!     Tweak here ?
      V(DFAC) = 0.6D+0
!     This acts as a lower-bound for scaling
!     Making this relatively small seems to be better
!      V(DFAC) = 0.25d+0
!     But... this also effects the updating!
      V(DELTA0) = SQTEPS
      V(DTINIT) = 1.D-6
      MEPCRT = MACHEP ** (ONE/THREE)
      V(D0INIT) = 1.D+0
      V(EPSLON) = 0.1D+0
      V(INCFAC) = 2.D+0
      V(LMAX0) = 1.D+0
      V(LMAXS) = 1.D+0
      V(PHMNFC) = -0.1D+0
      V(PHMXFC) = 0.1D+0
      V(RDFCMN) = 0.1D+0
      V(RDFCMX) = 4.D+0
      V(RFCTOL) = DMAX1(1.D-10, MEPCRT**2)
      V(SCTOL) = V(RFCTOL)
      V(TUNER1) = 0.1D+0
      V(TUNER2) = 1.D-4
      V(TUNER3) = 0.75D+0
      V(TUNER4) = 0.5D+0
      V(TUNER5) = 0.75D+0
      V(XCTOL) = SQTEPS
      V(XFTOL) = 1.D+2 * MACHEP
!
      IF (ALG .GE. 2) GO TO 10
!
!  ***  REGRESSION  VALUES
!
      V(COSMIN) = DMAX1(1.D-6, 1.D+2 * MACHEP)
      V(DINIT) = 0.D+0
      V(DLTFDC) = MEPCRT
      V(DLTFDJ) = SQTEPS
      V(FUZZ) = 1.5D+0
      V(HUBERC) = 0.7D+0
      V(RLIMIT) = DR7MDC(5)
      V(RSPTOL) = 1.D-3
      V(SIGMIN) = 1.D-4
      GO TO 999
!
!  ***  GENERAL OPTIMIZATION VALUES
!
 10   V(BIAS) = 0.8D+0
      V(DINIT) = -1.0D+0
      V(ETA0) = 1.0D+3 * MACHEP
!
 999  RETURN
!  ***  LAST CARD OF DV7DFL FOLLOWS  ***
      END
