      SUBROUTINE DD7DOG(DIG, LV, N, NWTSTP, STEP, V)
      implicit real*8 (a-h,o-z)
!
!  ***  COMPUTE DOUBLE DOGLEG STEP  ***
!
!  ***  PARAMETER DECLARATIONS  ***
!
      INTEGER LV, N
      DOUBLE PRECISION DIG(N), NWTSTP(N), STEP(N), V(LV)
!
!  ***  PURPOSE  ***
!
!        THIS SUBROUTINE COMPUTES A CANDIDATE STEP (FOR USE IN AN UNCON-
!     STRAINED MINIMIZATION CODE) BY THE DOUBLE DOGLEG ALGORITHM OF
!     DENNIS AND MEI (REF. 1), WHICH IS A VARIATION ON POWELL*S DOGLEG
!     SCHEME (REF. 2, P. 95).
!
!--------------------------  PARAMETER USAGE  --------------------------
!
!    DIG (INPUT) DIAG(D)**-2 * G -- SEE ALGORITHM NOTES.
!      G (INPUT) THE CURRENT GRADIENT VECTOR.
!     LV (INPUT) LENGTH OF V.
!      N (INPUT) NUMBER OF COMPONENTS IN  DIG, G, NWTSTP,  AND  STEP.
! NWTSTP (INPUT) NEGATIVE NEWTON STEP -- SEE ALGORITHM NOTES.
!   STEP (OUTPUT) THE COMPUTED STEP.
!      V (I/O) VALUES ARRAY, THE FOLLOWING COMPONENTS OF WHICH ARE
!             USED HERE...
! V(BIAS)   (INPUT) BIAS FOR RELAXED NEWTON STEP, WHICH IS V(BIAS) OF
!             THE WAY FROM THE FULL NEWTON TO THE FULLY RELAXED NEWTON
!             STEP.  RECOMMENDED VALUE = 0.8 .
! V(DGNORM) (INPUT) 2-NORM OF DIAG(D)**-1 * G -- SEE ALGORITHM NOTES.
! V(DSTNRM) (OUTPUT) 2-NORM OF DIAG(D) * STEP, WHICH IS V(RADIUS)
!             UNLESS V(STPPAR) = 0 -- SEE ALGORITHM NOTES.
! V(DST0) (INPUT) 2-NORM OF DIAG(D) * NWTSTP -- SEE ALGORITHM NOTES.
! V(GRDFAC) (OUTPUT) THE COEFFICIENT OF  DIG  IN THE STEP RETURNED --
!             STEP(I) = V(GRDFAC)*DIG(I) + V(NWTFAC)*NWTSTP(I).
! V(GTHG)   (INPUT) SQUARE-ROOT OF (DIG**T) * (HESSIAN) * DIG -- SEE
!             ALGORITHM NOTES.
! V(GTSTEP) (OUTPUT) INNER PRODUCT BETWEEN G AND STEP.
! V(NREDUC) (OUTPUT) FUNCTION REDUCTION PREDICTED FOR THE FULL NEWTON
!             STEP.
! V(NWTFAC) (OUTPUT) THE COEFFICIENT OF  NWTSTP  IN THE STEP RETURNED --
!             SEE V(GRDFAC) ABOVE.
! V(PREDUC) (OUTPUT) FUNCTION REDUCTION PREDICTED FOR THE STEP RETURNED.
! V(RADIUS) (INPUT) THE TRUST REGION RADIUS.  D TIMES THE STEP RETURNED
!             HAS 2-NORM V(RADIUS) UNLESS V(STPPAR) = 0.
! V(STPPAR) (OUTPUT) CODE TELLING HOW STEP WAS COMPUTED... 0 MEANS A
!             FULL NEWTON STEP.  BETWEEN 0 AND 1 MEANS V(STPPAR) OF THE
!             WAY FROM THE NEWTON TO THE RELAXED NEWTON STEP.  BETWEEN
!             1 AND 2 MEANS A TRUE DOUBLE DOGLEG STEP, V(STPPAR) - 1 OF
!             THE WAY FROM THE RELAXED NEWTON TO THE CAUCHY STEP.
!             GREATER THAN 2 MEANS 1 / (V(STPPAR) - 1) TIMES THE CAUCHY
!             STEP.
!
!-------------------------------  NOTES  -------------------------------
!
!  ***  ALGORITHM NOTES  ***
!
!        LET  G  AND  H  BE THE CURRENT GRADIENT AND HESSIAN APPROXIMA-
!     TION RESPECTIVELY AND LET D BE THE CURRENT SCALE VECTOR.  THIS
!     ROUTINE ASSUMES DIG = DIAG(D)**-2 * G  AND  NWTSTP = H**-1 * G.
!     THE STEP COMPUTED IS THE SAME ONE WOULD GET BY REPLACING G AND H
!     BY  DIAG(D)**-1 * G  AND  DIAG(D)**-1 * H * DIAG(D)**-1,
!     COMPUTING STEP, AND TRANSLATING STEP BACK TO THE ORIGINAL
!     VARIABLES, I.E., PREMULTIPLYING IT BY DIAG(D)**-1.
!
!  ***  REFERENCES  ***
!
! 1.  DENNIS, J.E., AND MEI, H.H.W. (1979), TWO NEW UNCONSTRAINED OPTI-
!             MIZATION ALGORITHMS WHICH USE FUNCTION AND GRADIENT
!             VALUES, J. OPTIM. THEORY APPLIC. 28, PP. 453-482.
! 2. POWELL, M.J.D. (1970), A HYBRID METHOD FOR NON-LINEAR EQUATIONS,
!             IN NUMERICAL METHODS FOR NON-LINEAR EQUATIONS, EDITED BY
!             P. RABINOWITZ, GORDON AND BREACH, LONDON.
!
!  ***  GENERAL  ***
!
!     CODED BY DAVID M. GAY.
!     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH SUPPORTED
!     BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS MCS-7600324 AND
!     MCS-7906671.
!
!------------------------  EXTERNAL QUANTITIES  ------------------------
!
!  ***  INTRINSIC FUNCTIONS  ***
!/+
      DOUBLE PRECISION DSQRT
!/
!--------------------------  LOCAL VARIABLES  --------------------------
!
      INTEGER I
      DOUBLE PRECISION CFACT, CNORM, CTRNWT, GHINVG, FEMNSQ, GNORM, &
                       NWTNRM, RELAX, RLAMBD, T, T1, T2
      DOUBLE PRECISION HALF, ONE, TWO, ZERO
!
!  ***  V SUBSCRIPTS  ***
!
      INTEGER BIAS, DGNORM, DSTNRM, DST0, GRDFAC, GTHG, GTSTEP, &
              NREDUC, NWTFAC, PREDUC, RADIUS, STPPAR
!
!  ***  DATA INITIALIZATIONS  ***
!
!/6
!     DATA HALF/0.5D+0/, ONE/1.D+0/, TWO/2.D+0/, ZERO/0.D+0/
!/7
      PARAMETER (HALF=0.5D+0, ONE=1.D+0, TWO=2.D+0, ZERO=0.D+0)
!/
!
!/6
!     DATA BIAS/43/, DGNORM/1/, DSTNRM/2/, DST0/3/, GRDFAC/45/,
!    1     GTHG/44/, GTSTEP/4/, NREDUC/6/, NWTFAC/46/, PREDUC/7/,
!    2     RADIUS/8/, STPPAR/5/
!/7
      PARAMETER (BIAS=43, DGNORM=1, DSTNRM=2, DST0=3, GRDFAC=45, &
                 GTHG=44, GTSTEP=4, NREDUC=6, NWTFAC=46, PREDUC=7, &
                 RADIUS=8, STPPAR=5)
!/
!
!+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
!
!      print*, 'nwtstp=',n,NWTSTP
!      print*, 'v=', V
      NWTNRM = V(DST0)
      RLAMBD = ONE
      IF (NWTNRM .GT. ZERO) RLAMBD = V(RADIUS) / NWTNRM
      GNORM = V(DGNORM)
      GHINVG = TWO * V(NREDUC)
      V(GRDFAC) = ZERO
      V(NWTFAC) = ZERO
      IF (RLAMBD .LT. ONE) GO TO 30
!
!        ***  THE NEWTON STEP IS INSIDE THE TRUST REGION  ***
!
         V(STPPAR) = ZERO
         V(DSTNRM) = NWTNRM
         V(GTSTEP) = -GHINVG
         V(PREDUC) = V(NREDUC)
         V(NWTFAC) = -ONE
         DO 20 I = 1, N
 20           STEP(I) = -NWTSTP(I)
         GO TO 999
!
 30   V(DSTNRM) = V(RADIUS)
      CFACT = (GNORM / V(GTHG))**2
!     ***  CAUCHY STEP = -CFACT * G.
      CNORM = GNORM * CFACT
      RELAX = ONE - V(BIAS) * (ONE - GNORM*CNORM/GHINVG)
      IF (RLAMBD .LT. RELAX) GO TO 50
!
!        ***  STEP IS BETWEEN RELAXED NEWTON AND FULL NEWTON STEPS  ***
!
         V(STPPAR)  =  ONE  -  (RLAMBD - RELAX) / (ONE - RELAX)
         T = -RLAMBD
         V(GTSTEP) = T * GHINVG
         V(PREDUC) = RLAMBD * (ONE - HALF*RLAMBD) * GHINVG
         V(NWTFAC) = T
         DO 40 I = 1, N
 40           STEP(I) = T * NWTSTP(I)
         GO TO 999
!
 50   IF (CNORM .LT. V(RADIUS)) GO TO 70
!
!        ***  THE CAUCHY STEP LIES OUTSIDE THE TRUST REGION --
!        ***  STEP = SCALED CAUCHY STEP  ***
!
         T = -V(RADIUS) / GNORM
         V(GRDFAC) = T
         V(STPPAR) = ONE  +  CNORM / V(RADIUS)
         V(GTSTEP) = -V(RADIUS) * GNORM
      V(PREDUC) = V(RADIUS)*(GNORM - HALF*V(RADIUS)*(V(GTHG)/GNORM)**2)
         DO 60 I = 1, N
 60           STEP(I) = T * DIG(I)
         GO TO 999
!
!     ***  COMPUTE DOGLEG STEP BETWEEN CAUCHY AND RELAXED NEWTON  ***
!     ***  FEMUR = RELAXED NEWTON STEP MINUS CAUCHY STEP  ***
!
 70   CTRNWT = CFACT * RELAX * GHINVG / GNORM
!     *** CTRNWT = INNER PROD. OF CAUCHY AND RELAXED NEWTON STEPS,
!     *** SCALED BY GNORM**-1.
      T1 = CTRNWT - GNORM*CFACT**2
!     ***  T1 = INNER PROD. OF FEMUR AND CAUCHY STEP, SCALED BY
!     ***  GNORM**-1.
      T2 = V(RADIUS)*(V(RADIUS)/GNORM) - GNORM*CFACT**2
      T = RELAX * NWTNRM
      FEMNSQ = (T/GNORM)*T - CTRNWT - T1
!     ***  FEMNSQ = SQUARE OF 2-NORM OF FEMUR, SCALED BY GNORM**-1.
      T = T2 / (T1 + DSQRT(T1**2 + FEMNSQ*T2))
!     ***  DOGLEG STEP  =  CAUCHY STEP  +  T * FEMUR.
      T1 = (T - ONE) * CFACT
      V(GRDFAC) = T1
      T2 = -T * RELAX
      V(NWTFAC) = T2
      V(STPPAR) = TWO - T
      V(GTSTEP) = T1*GNORM**2 + T2*GHINVG
      V(PREDUC) = -T1*GNORM * ((T2 + ONE)*GNORM) &
                       - T2 * (ONE + HALF*T2)*GHINVG &
                        - HALF * (V(GTHG)*T1)**2
      DO 80 I = 1, N
 80      STEP(I) = T1*DIG(I) + T2*NWTSTP(I)
!
 999  if(v(stppar).eq.0)then
          write(6,661)':DD7DOG Newton Step, radius',V(DSTNRM)
661     format(a,e10.3)
662     format(a,f10.3,' Radius ',e10.3)
      else if(V(stppar).le.1)then
          write(6,662)':DD7DOG Relaxed Newton ',V(STPPAR),v(DSTNRM)
      else if(V(STPPAR).le.2.)then
          write(6,662)':DD7DOG Dogleg Step',V(STPPAR)-1,v(DSTNRM)
      else if(V(STPPAR).gt.2.)then
          write(6,662)':DD7DOG Cauchy Step',1./(v(stppar)-1),v(dstnrm)
      endif
      RETURN
!  ***  LAST LINE OF DD7DOG FOLLOWS  ***
      END
