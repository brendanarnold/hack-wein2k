      SUBROUTINE EXCH17(D,S,U,V,X,EX,VX)
!  GGA91 EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
!  INPUT D : DENSITY
!  INPUT S:  ABS(GRAD D)/(2*KF*D)
!  INPUT U:  (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KF)**3)
!  INPUT V: (LAPLACIAN D)/(D*(2*KF)**2)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA A1,A2,A3,A4/0.19645D0,0.27430D0,0.15084D0,100.D0/
      DATA AA,AA1,AA2/0.8074d0,0.545663D0,0.56270D0/
      DATA AX,A,B1/-0.7385588D0,7.7956D0,0.004D0/
      DATA THRD,THRD4/0.333333333333D0,1.33333333333D0/
      FAC = AX*D**THRD
      S2 = S*S
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
      P0 = 1.D0/DSQRT(1.D0+A*A*S2)
      P1 = DLOG(A*S+1.D0/P0)
      P2 = DEXP(-A4*S2)
      P3 = 1.D0/(1.D0+A1*S*P1+B1*S4)
      P4 = 1.D0+A1*S*P1+(A2-A3*P2)*S2
      F = P3*P4
      EX = FAC*F
!  LOCAL EXCHANGE OPTION
!     EX = FAC
!  ENERGY DONE. NOW THE POTENTIAL:
!eq8      S2 = S*S
!eq8      P0 = 1.D0/DSQRT(1.D0+AA*AA*S2)
!eq8      P1 = DLOG(AA*S+1.D0/P0)
!eq8      P3 = 1.D0/(1.D0+AA1*S*P1)
!eq8
!eq8      P2 = DEXP(-A4*S2)
!eq8      P4 = (AA2-A3*P2)*S2
!eq8
!eq8      F = P3*P4
!eq8      F = AX*D**THRD*2.0d0**THRD4*F
!eq7
!eq7      X2 = X*X
!eq7      P0 = 1.D0/DSQRT(1.D0+X2)
!eq7      P1 = DLOG(X+1.D0/P0)
!eq7      P3 = 1.D0/(1.D0+6.0d0*0.323d0*X*P1)
!eq7
!eq7      P4 = X2
!eq7
!eq7      F = P3*P4
!eq7      F = -2.0d0*0.323d0*D**THRD*F
!eq7
!eq7_8      VX = FAC*THRD4+F
!
      X2 = X*X
      P0 = 1.D0/DSQRT(1.D0+X2)
      P1 = DLOG(X+1.D0/P0)
      P3 = 1.D0/(1.D0+3.0d0*0.05d0*X*P1)
!
      P4 = X2
!
      F = P3*P4
      F = -0.05*D**THRD*F
      VX = FAC*THRD4+F
!  LOCAL EXCHANGE OPTION:
!     VX = FAC*THRD4
      RETURN
      END
