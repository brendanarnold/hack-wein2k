      SUBROUTINE LOHNS(JNEQ,MULT,I1,I2,L0,jlo)
!
      use lolog, only : NLO, ILO
      use comi, only  : NAT
      IMPLICIT NONE
      INCLUDE 'param.inc'
!
!        Arguments
!
      INTEGER            I1, I2, JNEQ, L0
      INTEGER            MULT(NAT)
!
!     ..................................................................
!
!        LOHNS calculates indices of loop in hns
!
!     ..................................................................
!
!        Local Scalars
!
      INTEGER            I, L
      INTEGER            JLO, JLO1
!
      I2 = 1
      DO 20 I = 1, JNEQ - 1
         DO 10 L = 0, LOMAX
!            I1 = I2
!            IF (LOOR(L,I)) I2 = I2 + (2*L+1)*MULT(I)
             do jlo1=1,ilo(l,i)
              I2 = I2 + (2*L+1)*MULT(I)
             enddo
   10    CONTINUE
   20 CONTINUE
      DO 30 L = 0, L0-1
!         I1 = I2
!         IF (LOOR(L,JNEQ)) I2 = I2 + (2*L+1)*MULT(JNEQ)
             do jlo1=1,ilo(l,jneq)
              I2 = I2 + (2*L+1)*MULT(jneq)
             enddo
   30 CONTINUE

      do jlo1=1,jlo
         I1 = I2
         I2 = I2 + (2*L+1)*MULT(jneq)
      enddo
!
      I2 = I2 - 1
!
      RETURN
!
!        End of 'LOHNS'
!
      END
