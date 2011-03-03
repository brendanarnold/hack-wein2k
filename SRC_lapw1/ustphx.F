!_REAL      DOUBLE PRECISION FUNCTION USTPHX (SK,PHSAVE,j,NAT,V,R,MULT)
!_COMPLEX      COMPLEX*16 FUNCTION USTPHX(SK,PHSAVE,j,NAT,V,R,MULT)
      USE matrices, only: HSROWS
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INTEGER            NAT,  MULT(NAT)
      DOUBLE PRECISION   R(NAT),  V(NAT),  SK(3)
!_REAL      DOUBLE PRECISION   PHSAVE(HSROWS,NAT)
!_COMPLEX      COMPLEX*16         PHSAVE(HSROWS,NAT)
!
!     For rectangular ('cubic case') and hexagonal lattices
!
      DOUBLE PRECISION   TINY
!_REAL      DOUBLE PRECISION   DONE,  DZERO,  DTHREE
!_COMPLEX      COMPLEX*16         DONE,  DZERO,  DTHREE
!_REAL      PARAMETER          (DZERO  = 0.0D+0)
!_REAL      PARAMETER          (DONE   = 1.0D+0)
!_REAL      PARAMETER          (DTHREE = 3.0D+0)
!_COMPLEX      PARAMETER          (DZERO  = (0.0D+0,0.0D+0))
!_COMPLEX      PARAMETER          (DONE   = (1.0D+0,0.0D+0))
!_COMPLEX      PARAMETER          (DTHREE = (3.0D+0,0.0D+0))
      PARAMETER          (TINY   = 1.0D-10)
!
      INTEGER            JATOM, J
      DOUBLE PRECISION   Q,  BESR,  S
!_REAL      DOUBLE PRECISION   U
!_COMPLEX      COMPLEX*16         U
!
      INTRINSIC          DBLE,  SIN,  COS
!
      S = SQRT(SK(1)*SK(1)+SK(2)*SK(2)+SK(3)*SK(3))
      IF (S .GT. TINY) THEN
!
!        Define stepfunction U(K) for K .NE. 0
!
         U = DZERO
         DO 10 JATOM = 1, NAT
            Q = S*R(JATOM)
            BESR = (Q*COS(Q) - SIN(Q))/(Q*Q*Q)
            U = U + V(JATOM)*PHSAVE(j,JATOM)*BESR
 10      CONTINUE
         USTPHX = DTHREE*U
      ELSE
!
!        Define stepfunction U(K) for K = 0
!
         U = DZERO
         DO 20 JATOM = 1, NAT
            U = U + V(JATOM)*DBLE(MULT(JATOM))
 20      CONTINUE
         USTPHX = DONE - U
      ENDIF
!
      RETURN
!
!        End of 'USTPHX'
!
      END
