      PROGRAM INITELNES
!
!     Prepare input files for the calculation
!     of X-ray spectra
!
!     (C)1996 by Joachim Luitz
!     new version for ELNES by Pierre-Henri Louf
!
      IMPLICIT REAL*8 (A-H,O-Z)

      EXTERNAL OUTERR

      CHARACTER *80 DEFFN,ERRFN
      CHARACTER *67 ERRMSG
      CHARACTER *80 QTLL,DESCR
      CHARACTER *2 DUMMY,TESTS
      CHARACTER *80 FNAME
      CHARACTER *11 STATUS,FORM

      INTEGER IUNIT,IRECL
      INTEGER i,j,k
      INTEGER pos,natom,lc

      DOUBLE PRECISION Q, QPr, CosAlpha, DeltaE, SPLIT, PreV
      INTEGER LambdaMax
      DOUBLE PRECISION Energy
      DOUBLE PRECISION ThetaX, ThetaY
      DOUBLE PRECISION TwoThetaB
      DOUBLE PRECISION DeltaThX, DeltaThY
      INTEGER NStep
      CHARACTER Choice, HyperFine
      
      LOGICAL FOUND
      INTEGER Em
      INTEGER LMax

      INTEGER OurISPLIT

!     Special 99 case:
      CHARACTER*3 CharDOS(16)


!     for the moment I set LMax to 3 (f DOS)
      LMax = 3
!
      call GTFNAM(DEFFN,ERRFN)
      CALL ERRFLG(ERRFN,'Error in INIT_ELNES')
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
           

!     get input data from *.innes (number 5)
      READ(5,'(a)')DESCR
      READ(5,*)NATOM
      READ(5,*)NC
      READ(5,*)LC
      READ(5,*)EMIN,DE,EMAX
      PRINT *, EMIN,DE,EMAX
!      READ(5,*) DeltaE, SPLIT, PreV
!      READ(5,*) Energy
!      READ(5,*) ThetaX, ThetaY
!      READ(5,*) TwoThetaB
!      READ(5,*) DeltaThX, DeltaThY, NStep
!      READ(5,*) LambdaMax
!      READ(5,'(A)') Choice
!      READ(5,'(A)') HyperFine
!     we don't need more in case.innes

      CALL GetSPLIT(NATOM, OurISPLIT)

      READ(30,3269)EF
 2222 CONTINUE
      READ(30,3270)IATOM,QTLL
      IF(NATOM.NE.IATOM) GOTO 2222
      
      WRITE(8,*)'Autocreate '//DESCR
      WRITE(8,*)(EMIN/13.6058+EF),(DE/13.6058),(EMAX/13.6058+EF),0.0
      
!     
      IF (OurISPLIT.EQ.99) THEN
         WRITE(CharDOS(1),'(3A)') '0 0'
         WRITE(CharDOS(2),'(3A)') '1-1'
         WRITE(CharDOS(3),'(3A)') '1 0'
         WRITE(CharDOS(4),'(3A)') '1 1'
         WRITE(CharDOS(5),'(3A)') '2-2'
         WRITE(CharDOS(6),'(3A)') '2-1'
         WRITE(CharDOS(7),'(3A)') '2 0'
         WRITE(CharDOS(8),'(3A)') '2 1'
         WRITE(CharDOS(9),'(3A)') '2 2'
         WRITE(CharDOS(10),'(3A)') '3-3'
         WRITE(CharDOS(11),'(3A)') '3-2'
         WRITE(CharDOS(12),'(3A)') '3-1'
         WRITE(CharDOS(13),'(3A)') '3 0'
         WRITE(CharDOS(14),'(3A)') '3 1'
         WRITE(CharDOS(15),'(3A)') '3 2'
         WRITE(CharDOS(16),'(3A)') '3 3'

         WRITE(8,*) 273
         J = 0
         DO LL=1, 16
            DO LLP=1, LL
!     real part
               J = J+1
               WRITE (8,3272) NATOM,(J+100),CharDOS(LL)//CharDOS(LLP)
!     imaginary part
               J= J+1
               WRITE (8,3272) NATOM,(J+100),CharDOS(LL)//CharDOS(LLP)
            ENDDO
         ENDDO
!     number 273:
         WRITE(8,3272)0,1,'tot'
      ELSE
!     the atom is "normal" in the QTL file: 13 columns

                  
!      IF (HyperFine.EQ.'N') THEN
!     How many partial DOS do we need?
!     We need all the DOS from l=0 to LMax (LMax+1) AND total-DOS
!         write(8,*) LMax+2
!      ELSE IF (HyperFine.EQ.'H') THEN
!         write(8,*) LMax*LMax+2
!      ENDIF

!     we need: 0 1(-1 0 1) 2(-2 -1 0 1 2) 3 and total : 13 DOS cases.
!     They are given on 2 files case.dos1ev and case.dos2ev.
         WRITE(8,*) 13

!ph  **************************************** PH Louf 17 11 97
!ph for ll=lc +- 1 => ll=0 -> LMax
         do ll=0, LMax
            DUMMY=','//char(48+ll)
!     48 is the ascii position of '0'
            call STRLN(QTLL,k)
!     we get the length of the line in k (cf. strlen.f)
!     the do ... is the same as in initxspec.f
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
               WRITE(8,3272) NATOM, pos, 'l='//char(48+ll)
               IF ((ll.NE.0).AND.(ll.NE.LMAX)) THEN
                  DO Em = ll, 1, -1
                     write(8,3272)NATOM,pos+ll-Em+1,  &
                          char(48+ll)//' -'//char(48+Em)
                  ENDDO
                  WRITE (8,3272) NATOM, pos+ll+1, char(48+ll)//'  0'
                  DO Em = 1, ll
                     write(8,3272)NATOM,pos+ll+Em+1, &
                          char(48+ll)//' +'//char(48+Em)
                  ENDDO
               ENDIF
            ELSE
!     IF (HyperFine.EQ.'N') THEN
!     Normal case: write l-total DOS.
!               write(8,3272)NATOM,pos,'l='//char(48+ll)
!            ELSE IF (HyperFine.EQ.'H') THEN
!     Case HyperFine: lm DOS to write.
!               IF ((ll.EQ.0).OR.(ll.EQ.LMAX)) THEN
!                  write(8,3272)NATOM,pos,'l='//char(48+ll)
!               ELSE
!                  DO Em = 1, 2*ll+1
!                     write(8,3272)NATOM,pos+Em,'l='//char(48+ll)
!                  ENDDO
!               ENDIF
!            ENDIF
!         else
               goto 913
            endif
            FOUND=.FALSE.
         enddo
         
         write(8,3272)0,1,'tot'
      ENDIF
      CALL ERRCLR(ERRFN)
      STOP 'INIT_ELNES - done'

! errors

!     elnes.def couldn't be opened
 910  WRITE(ERRMSG,9000) FNAME
      CALL OUTERR('INIT_ELNES',ERRMSG)
      GOTO 999

!     elnes.def is corrupt
 911  WRITE(ERRMSG,9001)FNAME
      CALL OUTERR('INIT_ELNES',ERRMSG)
      GOTO 999

!     Error opening a unit
 912  WRITE (ERRMSG,9010) IUNIT
      CALL OUTERR('INIT_ELNES',ERRMSG)
      WRITE (ERRMSG,9020) FNAME
      CALL OUTERR('INIT_ELNES',ERRMSG)
      WRITE (ERRMSG,9030) STATUS, FORM
      CALL OUTERR('INIT_ELNES',ERRMSG)
      GOTO 999

!     L not found in QTL file
 913  WRITE(ERRMSG,9090)LL
      CALL OUTERR('INIT_ELNES',ERRMSG)
      GOTO 999

 999  STOP 'INIT_ELNES - Error'

! formats
 3220 FORMAT(/,28x,i2)
 3269 FORMAT(//,56x,f10.5/)
 3270 FORMAT (7x,i2,22x,60a)
 3271 FORMAT(2a)
 3272 FORMAT(3x,i2,2x,i3,3x,a6)
 24   FORMAT(3X,I2,2X,I3,3X,A6) 

 9000 FORMAT('can''t open definition file ',A40)
 9001 FORMAT('can''t read definition file ',A40)
 9010 FORMAT('can''t open unit: ',I2)
 9020 FORMAT('       filename: ',A50)
 9030 FORMAT('         status: ',A,'  form: ',A)
 9090 FORMAT('l=',I5,' not found in QTL-File!')
      END



