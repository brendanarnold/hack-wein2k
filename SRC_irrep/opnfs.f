      SUBROUTINE OPNFS(DEFFN,ERRFN,FL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'param.inc'
!
      CHARACTER*80    DEFFN,ERRFN,FNAME,TNAME
      CHARACTER*11    STATUS,FORM
      CHARACTER*1     TJUNK
      LOGICAL         FL(FLMAX), FLG1
!**********************************************************************
!
      iin=0
      DO 5 II=1,FLMAX
 5       FL(II)=.FALSE.
!
      CALL GTFNAM(DEFFN,ERRFN)
      OPEN(1,FILE=DEFFN,STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL


!        check for in1c file for complex case
      
!
      DO 10 I = LEN(fname), 1, -1
        IF (fname(I:I) .EQ. '.') THEN
          IF (fname(i+1:i+8 ) .eq. 'vectorso')    fl(2)=.true.
          IF (fname(i+1:i+10) .eq. 'vectorsodn')  fl(2)=.true.
          IF (fname(i+1:i+8 ) .eq. 'vectorup')    fl(4)=.true.
          IF (fname(i+1:i+10) .eq. 'vectorsoup')  fl(4)=.true.
!
          IF (fname(i+1:i+6 ) .eq. 'struct') GOTO 20
          GOTO 11
        ENDIF
   10 CONTINUE
      GOTO 11
!
!.....test if case.in1 or case.in1c is found
 20   flg1=.false.

      tname = fname(1:i)//'in1'
      OPEN(IUNIT,FILE=TNAME,STATUS='OLD',FORM='FORMATTED',ERR=21)
      READ(IUNIT,FMT='(A1)',ERR=21,END=21) tjunk
      flg1=.true.
 21   CLOSE (IUNIT) 
      tname = fname(1:i)//'in1s'
      OPEN(IUNIT,FILE=TNAME,STATUS='OLD',FORM='FORMATTED',ERR=22)
      READ(IUNIT,FMT='(A1)',ERR=22,END=22) tjunk
      flg1=.true.
 22   CLOSE (IUNIT) 
      tname = fname(1:i)//'in1c'
      OPEN(IUNIT,FILE=TNAME,STATUS='OLD',FORM='FORMATTED',ERR=23)
      READ(IUNIT,FMT='(A1)',ERR=23,END=23) tjunk
      fl(1)=.true.
 23   CLOSE (IUNIT)
      tname = fname(1:i)//'in1cs'
      OPEN(IUNIT,FILE=TNAME,STATUS='OLD',FORM='FORMATTED',ERR=24)
      READ(IUNIT,FMT='(A1)',ERR=24,END=24) tjunk
      fl(1)=.true.
 24   CLOSE (IUNIT)
      
!
!
 11   OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=8002)
      GOTO 8003



 8000 WRITE(*,*) ' ERROR IN OPENING KSTCLASS.DEF !!!!'
      STOP 'KSTCLASS.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, &
      '  FORM:',FORM
      STOP 'OPEN FAILED'
 8001 CONTINUE
      IF( (.not.flg1).and.(.not.fl(1)) ) THEN
          WRITE(*,*) ' ERROR IN OPENING: NEEDS case.in1 or case.in1c !!!!'
          STOP 'case.in1(c)'
      ENDIF
      IF( flg1.and.fl(1) ) THEN
         WRITE(*,*)
         WRITE(*,*) 'WARNING !!!   Both case.in1[s] and case.in1c[s] exist.'
         WRITE(*,*) &
           'Make sure that the program handles the eigenfunctions correctly.'
         WRITE(*,*) 'Here, IRREP assumes complex valued eigenfunctions.'
         WRITE(*,*) 'If the eigenfunctions are real valued, please remove case.in1c[s].'
         WRITE(*,*)
     
      ENDIF     
      CLOSE (1)
!
      END  
