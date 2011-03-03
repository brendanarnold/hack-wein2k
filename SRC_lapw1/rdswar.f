      SUBROUTINE RDSWAR
!
      use out, only      : NKK
      use totpot, only   : JA, JB, JC, POTK, init_totpot
      IMPLICIT NONE
      INCLUDE 'param.inc'
  CHARACTER*19     nkktext
!
!     ..................................................................
!
!        read interstitial total potential (for warpin) from tape
! 
!     ..................................................................
!
!        Local Scalars
!
      INTEGER            JK
!
      READ(19,1000)
!      READ(19,1010) NKK
  read(19,'(/,a19)') nkktext
  read(nkktext,'(9x,i10)',err=6767) nkk
  goto 6768
 6767 read(nkktext,'(13x,i6)') nkk
 6768 continue
      call init_totpot(NKK)
      DO 10 JK = 1, NKK
         READ(19,1020) JA(JK), JB(JK), JC(JK), POTK(JK)
   10 CONTINUE
      REWIND 19
!
      RETURN
!
 1000 FORMAT(3X)
 1010 FORMAT(/,13X,I6)
 1020 FORMAT(3X,3I5,2E19.12)
!
!        End of 'RDSWAR'
!
      END
