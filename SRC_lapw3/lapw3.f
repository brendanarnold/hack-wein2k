!$hp9000_800 intrinsics on
      PROGRAM LAPW3                                                     
      use atomgrid
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
!                                                                       
!                                                                       
      CHARACTER*80     FNAME                                      
      CHARACTER*11     STATUS,FORM                                      
      CHARACTER*4 LATTIC,ISYM,IREL,cform                                
      CHARACTER*10 TITLE(8)                                             
      CHARACTER*10     ANAME                                     
      LOGICAL REL                                                       
      COMMON /CHAR/   TITLE,LATTIC                   
      COMMON /COM/ XWT,EMIN,EF,ELECN,REL,MODUS,NBAND,NAT                
      COMMON /STRUK/ AA,BB,CC,alpha(3),PIA(3),VOL,ndif 
      COMMON /GENER/ GX1,GY1,GZ1,GX2,GY2,GZ2,GX3,GY3,GZ3,BR2(3,3)       
!      COMMON /CUBWK1/ DUM(999,ncom)                                     
      common /CUBWK1/ DUM(NRADINT,NCOM)
      COMMON /POTNLC/ DUMMY2(nrad)
!      DATA PI/3.1415926535898d0/                                        
      PI=4.D0*DATAN(1.D0)
!
      iarg=iargc()
      if(iarg.ne.1) STOP 'Exactly one commandline argument must be given'
      call getarg(1,fname)
      OPEN(1,FILE=fname,STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM, &
      ERR=8002)
      GOTO 8003
 8000 WRITE(*,*) ' ERROR IN OPENING LAPW3.DEF !!!!'
      STOP 'LAPW3.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, &
      '  FORM:',FORM
      STOP 'OPEN FAILED'
 8001 CONTINUE
      CLOSE (1)
!
      SQRT3=SQRT(3.d0)                                                    
      READ(20,1000) TITLE                                               
      READ(20,1001) LATTIC,NAT,cform,IREL  
      nato=nat                             
      allocate ( POS(nato*48,3))
      allocate ( RMT(nato),V(nato),Rnot(nato),DX(nato))
      allocate ( IATNR(nato),MULT(nato),ISPLIT(nato),jri(nato))   
      allocate ( ROTLOC(3,3,NATO),ROTIJ(3,3,Nato*48),TAUIJ(3,Nato*48))
      DO 9575 IND=1,nato                                                
      RMT(IND)=0.0d0                                                      
      V(IND)=0.0d0                                                        
 9575 MULT(IND)=0                                                       
      WRITE(66,1004) TITLE,LATTIC,IREL                                  
      WRITE(66,1006) NAT                                                
        assign 2021 to iform1
        assign 2071 to iform2
 2021 FORMAT(3X,4E19.12)                                                
 2071 FORMAT(3X,3I5,2E19.12)                                            
!                                                                       
      REL=.FALSE.                                                       
      IF(IREL.EQ.'RELA') REL=.TRUE.                                     
      READ(20,1007) AA,BB,CC,alpha                                            
      if(alpha(1).eq.0.d0) alpha(1)=90.d0
      if(alpha(2).eq.0.d0) alpha(2)=90.d0
      if(alpha(3).eq.0.d0) alpha(3)=90.d0
      WRITE(66,1008) AA,BB,CC,alpha                                           
!
      call latgen(nat)
! 
      CALL ATPAR(NAT)                                                   
!
        CALL FOURIR                                              
      STOP                                                              
 1000 FORMAT(8A10)                                                      
 1001 FORMAT(A4,23X,I3,1x,a4,/,13X,A4,18X,A4)                           
 1004 FORMAT(8A10,//,10X,A4,' LATTIC ASSUMED',//,10X,2A4,            &
      '-CALCULATION',//)                                                
 1006 FORMAT(10X,'NUMBER OF ATOMS:',I3)                                 
 1007 FORMAT(6F10.7)                                                    
 1008 FORMAT('  LATTICE CONSTANTS  ',6F10.5)                        
      END                                                               
