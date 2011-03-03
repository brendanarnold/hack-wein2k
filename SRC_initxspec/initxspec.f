      PROGRAM INIXSPEC
!
!     Prepare input files for the calculation
!     of X-ray spectra
!
!     (C)1996 by Joachim Luitz
!
      IMPLICIT REAL*8 (A-H,O-Z)

      EXTERNAL OUTERR

      CHARACTER *80 DEFFN,ERRFN
      CHARACTER *67 ERRMSG,DESCR
      CHARACTER *80 QTLL
      CHARACTER *2 DUMMY,TESTS
      CHARACTER *80 FNAME
      CHARACTER *11 STATUS,FORM

      INTEGER IUNIT,IRECL
      INTEGER i,j,k
      INTEGER pos,natom,lc

      LOGICAL FOUND
      
      call GTFNAM(DEFFN,ERRFN)
      CALL ERRFLG(ERRFN,'Error in INIT_XSPEC')
      OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
 8003 CONTINUE
         READ (1,*,END=8001,ERR=911) IUNIT,FNAME,STATUS,FORM,IRECL
!        IF(IUNIT.NE.6)then
! IF(IUNIT.NE.32)then
         OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM, &
               ERR=912)
! ENDIF
!        ENDIF
      GOTO 8003
 8001 CONTINUE

!     READ(20,3220)NATO

!     get input data from *.inxs

      READ(5,'(a)')DESCR
      READ(5,*)NATOM
      READ(5,*)NC
      READ(5,*)LC
!      READ(5,*)KAPPA
      READ(5,*)SPLIT,XINT1,XINT2
      READ(5,*)EMIN,DE,EMAX

      READ(30,3269)EF
 2222 CONTINUE
      READ(30,3270)IATOM,QTLL
      IF(NATOM.NE.IATOM) GOTO 2222

      WRITE(8,'(a)')'Autocreate '//DESCR
      WRITE(8,'(4f11.7)')(EMIN/13.6058d0+EF),(DE/13.6058d0),(EMAX/13.6058d0+EF),0.0d0

!     How many partial DOS do we need?
      IF (LC.EQ.0)then
         write(8,'(i5)') 2
      else
         write(8,'(i5)') 3
      endif

!     which part.DOS do we need? begin with lc+1

      ll=lc+1
 
!     check for position of ll in QTL-File
      DUMMY=','//char(48+ll)
!      call STRLN(QTLL,k)
      k=len_trim(qtll)
      do i=1,k-1
         TESTS=QTLL(i:i+1)
         IF(TESTS.EQ.DUMMY)then
            FOUND=.TRUE.
            pos=1
            DO j=1,i
               if(QTLL(j:j).EQ.',')then
                  pos=pos+1
               endif
            enddo
          endif
      enddo
      IF(FOUND)then
         write(8,3272)NATOM,pos,'l+1'
      else
         goto 913
      endif

      if(lc.gt.0)then
         FOUND=.FALSE.
         ll=lc-1
!        check for position of ll in QTL-File
         DUMMY=','//char(48+ll)
!         call STRLN(QTLL,k)
         k=len_trim(qtll)
         do i=1,k-1
            TESTS=QTLL(i:i+1)
            IF(TESTS.EQ.DUMMY)then
               FOUND=.TRUE.
               pos=1
               DO j=1,i
                  if(QTLL(j:j).EQ.',')then
                     pos=pos+1
                  endif
               enddo
            endif
         enddo
         IF(FOUND)then
            write(8,3272)NATOM,pos,'l-1'
         else
            goto 913
         endif
      endif
      write(8,3272)0,1,'tot'
      CALL ERRCLR(ERRFN)
      STOP 'INIT_XSPEC - done'

! errors

!     xspec.def couldn't be opened
 910  WRITE(ERRMSG,9000) FNAME
      CALL OUTERR('INIT_XSPEC',ERRMSG)
      GOTO 999

!     xspec.def is corrupt
 911  WRITE(ERRMSG,9001)FNAME
      CALL OUTERR('INIT_XSPEC',ERRMSG)
      GOTO 999

!     Error opening a unit
 912  WRITE (ERRMSG,9010) IUNIT
      CALL OUTERR('INIT_XSPEC',ERRMSG)
      WRITE (ERRMSG,9020) FNAME
      CALL OUTERR('INIT_XSPEC',ERRMSG)
      WRITE (ERRMSG,9030) STATUS, FORM
      CALL OUTERR('INIT_XSPEC',ERRMSG)
      GOTO 999

!     L not found in QTL file
 913  WRITE(ERRMSG,9090)LL
      CALL OUTERR('INIT_XSPEC',ERRMSG)
      GOTO 999

 999  STOP 'INIT_XSPEC - Error'

! formats
 3220 FORMAT(/,28x,i2)
 3269 FORMAT(//,56x,f10.5/)
 3270 FORMAT (7x,i2,22x,60a)
 3271 FORMAT(2a)
 3272 FORMAT(3x,i2,3x,i2,3x,a6)
 9000 FORMAT('can''t open definition file ',A40)
 9001 FORMAT('can''t read definition file ',A40)
 9010 FORMAT('can''t open unit: ',I2)
 9020 FORMAT('       filename: ',A50)
 9030 FORMAT('         status: ',A,'  form: ',A)
 9090 FORMAT('l=',I5,' not found in QTL-File!')
      END



