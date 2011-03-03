      SUBROUTINE FXHPBE(RHO,S,TAU,EX)
      IMPLICIT REAL*8(A-H,O-Z)      
      PARAMETER(THRD=1.D0/3.D0)
      PARAMETER(PI=3.14159265358979323846264338327950D0)
      PARAMETER(AX=-0.738558766382022405884230032680836D0)
!....old      PARAMETER(DD=0.089651D0,XK=0.804D0)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! CONSTRUCT LDA EXCHANGE ENERGY DENSITY
!      EXUNIF = AX*RHO**(THRD)
!      S2 = S*S
!      QT=3.D0/(2.D0*(3.D0*PI**2*RHO)**(2.D0/3.D0)*RHO)*TAU-
!     *9.D0/20.D0-S2/12.D0
!      XX=10.D0/81.D0*S2+DD*S2*S2-73.D0/405.D0*S2*QT+
!     *146.D0/2025.D0*QT**2+3.D0/(2.*XK)*(10.D0/81.D0*S2)**2
!      FX = 1.D0+XK-XK/SQRT(1.D0+2.D0*XX/XK)
!      EX = EXUNIF*FX
!      RETURN
!      END
!      SUBROUTINE FXHPBE(RHO,S,TAU,EX)
!      IMPLICIT REAL*8(A-H,O-Z)      
!      PARAMETER(THRD=1.D0/3.D0)
!      PARAMETER(PI=3.14159265358979323846264338327950D0)
!      PARAMETER(AX=-0.738558766382022405884230032680836D0)
      PARAMETER(DD=0.113D0,XK=0.804D0)
!----------------------------------------------------------------------
      EXUNIF = AX*RHO**(THRD)
      S2 = S*S
      QT=3.D0/(2.D0*(3.D0*PI**2*RHO)**(2.D0/3.D0)*RHO)*TAU- &
      9.D0/20.D0-S2/12.D0
      XX=10.D0/81.D0*S2+DD*S2*S2-73.D0/405.D0*S2*QT+ &
      146.D0/2025.D0*QT**2+1.D0/(XK)*(10.D0/81.D0*S2)**2
      FX = 1.D0+XK-XK/(1.D0+XX/XK)
      EX = EXUNIF*FX
      RETURN
      END
