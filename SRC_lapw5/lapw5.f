!$hp9000_800 intrinsics on
      PROGRAM LAPW5                                                     
!     PROGRAM LAPW5(INPUT,OUTPUT,POTE,SIGMA,CLM,PLOT,TAPE9=CLM,         
!    *TAPE12=SIGMA,TAPE8=POTE,TAPE5=INPUT,TAPE6=OUTPUT,                 
!    *TAPE20=PLOT,CLM1,TAPE11=CLM1,RHO,TAPE10=RHO,tape8=struct)                      
!
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*80      FNAME
      CHARACTER*11      STATUS,FORM 
      iarg=iargc()
      if(iarg.ne.1) goto 900
      CALL GETARG(1,fname)
!      call getarg(2,fname)
!      if(fname.eq.'      ') call getarg(1,fname)
      OPEN(1,FILE=fname,STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM, &
      ERR=8002)
      GOTO 8003
 8000 WRITE(*,*) ' ERROR IN OPENING LAPW5.DEF !!!!'
      STOP 'LAPW5.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, &
      '  FORM:',FORM
      STOP 'OPEN FAILED'
 8001 CONTINUE
      CALL MAIN1                                                         
      CALL PRTRH                                                        
      STOP                                                              
 900 STOP 'GTFNAM - Exactly one commandline argument has to be given.'
      END                                                               
