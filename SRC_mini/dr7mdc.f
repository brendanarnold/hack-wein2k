      DOUBLE PRECISION FUNCTION DR7MDC(K)
      implicit real*8 (a-h,o-z)
!
!  ***  RETURN MACHINE DEPENDENT CONSTANTS USED BY NL2SOL  ***
!
      INTEGER K
!
!  ***  THE CONSTANT RETURNED DEPENDS ON K...
!
!  ***        K = 1... SMALLEST POS. ETA SUCH THAT -ETA EXISTS.
!  ***        K = 2... SQUARE ROOT OF ETA.
!  ***        K = 3... UNIT ROUNDOFF = SMALLEST POS. NO. MACHEP SUCH
!  ***                 THAT 1 + MACHEP .GT. 1 .AND. 1 - MACHEP .LT. 1.
!  ***        K = 4... SQUARE ROOT OF MACHEP.
!  ***        K = 5... SQUARE ROOT OF BIG (SEE K = 6).
!  ***        K = 6... LARGEST MACHINE NO. BIG SUCH THAT -BIG EXISTS.
!
      DOUBLE PRECISION BIG, ETA, MACHEP
!/+
      DOUBLE PRECISION DSQRT
!/
!
      DOUBLE PRECISION D1MACH, ZERO
      EXTERNAL D1MACH
      DATA BIG/0.D+0/, ETA/0.D+0/, MACHEP/0.D+0/, ZERO/0.D+0/
      save big,eta,machep,zero
      IF (BIG .GT. ZERO) GO TO 1
         BIG = D1MACH(2)
         ETA = D1MACH(1)
         MACHEP = D1MACH(4)
 1    CONTINUE
!
!-------------------------------  BODY  --------------------------------
!
      GO TO (10, 20, 30, 40, 50, 60), K
!
 10   DR7MDC = ETA
      GO TO 999
!
 20   DR7MDC = DSQRT(256.D+0*ETA)/16.D+0
      GO TO 999
!
 30   DR7MDC = MACHEP
      GO TO 999
!
 40   DR7MDC = DSQRT(MACHEP)
      GO TO 999
!
 50   DR7MDC = DSQRT(BIG/256.D+0)*16.D+0
      GO TO 999
!
 60   DR7MDC = BIG
!
 999  RETURN
!  ***  LAST CARD OF DR7MDC FOLLOWS  ***
      END
