      program mini
!
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*11       STATUS,FORM
      CHARACTER*67       ERRMSG
      CHARACTER*80       DEFFN, ERRFN, FNAME
!
!
      CALL GTFNAM(DEFFN,ERRFN)
      CALL ERRFLG(ERRFN,'Error in MINI')
      OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
   10 CONTINUE
         READ (1,*,END=20,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL
         OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
      GOTO 10
   20 CONTINUE
      CLOSE (1)
      call haupt(errfn)
      CALL ERRCLR(ERRFN)
      stop '>>  (mini) arrived at end -> exit'
!
!        error handling
!
  910 INFO = 1
!
!        'lapw2.def' couldn't be opened
!
      WRITE (ERRMSG,9000) FNAME
      CALL OUTERR('MINI',ERRMSG)
      GOTO 999
  920 INFO = 2
!
!        file FNAME couldn't be opened
!
      WRITE (ERRMSG,9010) IUNIT
      CALL OUTERR('MINI',ERRMSG)
      WRITE (ERRMSG,9020) FNAME
      CALL OUTERR('MINI',ERRMSG)
      WRITE (ERRMSG,9030) STATUS, FORM
      CALL OUTERR('MINI',ERRMSG)
      GOTO 999
 960  INFO = 7
!
!        Error reading file 'lapw2.def'
!
      WRITE (ERRMSG,9040) FNAME
      CALL OUTERR('MINI',ERRMSG)
      GOTO 999
  999 STOP 'MINI - Error'
 9000 FORMAT('can''t open definition file ',A40)
 9010 FORMAT('can''t open unit: ',I2)
 9020 FORMAT('       filename: ',A50)
 9030 FORMAT('         status: ',A,'  form: ',A)
 9040 FORMAT('Error reading file: ',A47)
      end
