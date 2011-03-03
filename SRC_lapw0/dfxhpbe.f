      SUBROUTINE DFXHPBE(fu,gmagu,tauup,fd,gmagd,taudn,exu,exd)
      IMPLICIT REAL*8(A-H,O-Z)      
      PARAMETER(CON13=1.D0/3.D0)
      PARAMETER(PI=3.14159265358979323846264338327950D0)
!      PARAMETER(CON=0.6203504908994000166680068120477D0)
!      PARAMETER(ALPHA=1.919158292677513006624820326246D0)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! RHO=DENSITY
! ZETA=RELATIVE SPIN POLARIZATION
!      PI=4.*ATAN(1.)
      PI32=3.D0*PI**2
      CON=(3.D0/(4.D0*PI))**CON13
      ALPHA=(9.D0*PI/4.D0)**CON13
! UP DENSITY
      RHO2=2.d0*fu
      IF(RHO2.GT.1.D-18)THEN
        FK=(PI32*RHO2)**CON13
        S=2.d0*Gmagu/(2.D0*FK*RHO2)
        CALL FXHPBE(RHO2,S,2.d0*TAUup,EXU)
      ELSE
       EXU=0.D0
      ENDIF
! REPEAT FOR DOWN
      RHO2=2.d0*fd
      IF(RHO2.GT.1D-18)THEN
        FK=(PI32*RHO2)**CON13
        S=2.d0*Gmagd/(2.D0*FK*RHO2)
        CALL FXHPBE(RHO2,S,2.d0*TAUdn,EXD)
      ELSE
       EXD=0.D0
      ENDIF
      RETURN
      END
