!$hp9000_800 intrinsics on
      PROGRAM LAPW7
!     last changes: 29.08.00 ub (removing CALL to PRTRH)
!
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*80      FNAME
      CHARACTER*11      STATUS,FORM
!
!     << establish file conections >>
      iarg=iargc()
      CALL GETARG(iarg,FNAME)
!      CALL GETARG(2,FNAME)
!      IF(FNAME.EQ.'      ') CALL GETARG(1,FNAME)
      OPEN(1,FILE=FNAME,STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=8002)
      GOTO 8003
!
 8000 WRITE(*,*) ' ERROR IN OPENING LAPW7.DEF !!!!'
      STOP 'LAPW7.DEF'
!
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, &
                 '  FORM:',FORM
      STOP 'OPEN FAILED'
!
 8001 CONTINUE
!
!     << start wave function evaluation >>
      CALL MAIN2
      STOP
      END
