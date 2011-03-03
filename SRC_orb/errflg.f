      SUBROUTINE ERRFLG(FNAME,MSG)
      CHARACTER*(*)      FNAME, MSG
!
      OPEN (99,FILE=FNAME,ERR=900)
      WRITE (99,9000) MSG
      CLOSE (99)
      OPEN (99,FILE=FNAME,ERR=900)
!
      RETURN
!
  900 STOP 'ERRFLG - couldn''t open errorflag-file.'
!
 9000 FORMAT (A)
!
      END
