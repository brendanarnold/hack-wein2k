!$hp9000_800 intrinsics on
      PROGRAM HFSD                                                      
!                                                                       
!               H A R T R E E   F O C K   DIRAC SLATER                  
!                                                                       
!     RELHF1 CALCULATES THE ATOMIC EIGENVALUES OF A SPHERICAL TOTAL     
!     POTENTIAL                                                         
!     J P DESCLAUX   CEA PARIS 1969                                     
!     LAST UPDATE : 26. FEB 87    FOR NAS, IN INSLD TO READ FILE POTE   
!                                                                       
!         5            ...   INPUT                                      
!         6            ...   OUTPUT                                     
!         7            ...   CLMCORE                                    
!         8            ...   VSP                                        
!         9 .. 19                                                       
!         20           ...   POTE                                       
!         21           ...   SCFNEW                                     
!         22           ...   SCFOLD                                     
!                                                                       
!     PROGRAM HFSD (INPUT,OUTPUT,POTE,VSP,CLMCORE,SCFNEW,TAPE10,        
!    *    SCFOLD,   TAPE5=INPUT,TAPE6=OUTPUT,TAPE7=CLMCORE,TAPE8=VSP,   
!    *              TAPE20=POTE,TAPE21=SCFNEW,TAPE22=SCFOLD)            
!                                                                       
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
!fb
!     parameter for calculation of [yu91]-force
      logical force
      parameter (force=.true.)
!fe
      CHARACTER *10 BAR                                                 
      CHARACTER *11 STATUS,FORM                                         
      CHARACTER *77 RMARGN                                              
      CHARACTER*67       ERRMSG
      CHARACTER*80       DEFFN, ERRFN
      CHARACTER *80 FNAME                                               
      CHARACTER*4 IREL,cform                                       
      LOGICAL OLD,RELA                                                  
      COMMON        DEN(30),DQ1(30),DFL(30),NQN(30),NQL(30),NK(30),      &
                    NMAX(30),NEL(30),NORB                               
      COMMON /DIRA/ DV(NRAD+41),DR(NRAD+41),DP(NRAD+41),DQ(NRAD+41), &
                    DPAS,Z,NSTOP,NES,TETS,NP,NUC     
      COMMON  /PS2/ DEXV,DEXE,DCOP,TITRE(30),TEST,TESTE,                 &
                    TESTY,TESTV,NITER,ION,ICUT,IPRAT                    
      COMMON /CHAR/ BAR                                                 
      COMMON /DEUX/ DVN(NRAD+41),DVF(NRAD+41),D(NRAD+41),DC(NRAD+41), &
                    DGC(NRAD+41,30),DPC(NRAD+41,30) 
!                                                                       
!                                                                       
!--------------------------------------------------------------------   
!               INPUT GENERATING                                        
!                                                                       
!    1. WRITE NUMBER OF ATOMORBITALS TO INPUTFILE                       
!    2. WRITE PRINCIPAL QUANTUMNUMBER, QUANTUMNUMBER KAPPA AND          
!             NUMBER OF ELECTRONS FOR EACH ORBITAL IN ONE LINE          
!             OF INPUTFILE                                              
!    3. REPEAT POINT 1 AND 2 FOR THE FOLLOWING ATOMS                    
!                                                                       
!--------------------------------------------------------------------   
!                                                                       
!      CALL XUFLOW(0)                                                   
      CALL GTFNAM(DEFFN,ERRFN)
      CALL ERRFLG(ERRFN,'Error in CORE')
      OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
 8003 CONTINUE
         READ (1,*,END=8001,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL
         OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
      GOTO 8003
 8001 CONTINUE
      CLOSE (1)
!      call getarg(2,fname)
!      if(fname.eq.'      ') call getarg(1,fname)
!      OPEN(1,FILE=fname,STATUS='OLD',ERR=8000)
! 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
!      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,
!     *ERR=8002)
!      rewind iunit
!      GOTO 8003
! 8000 WRITE(*,*) ' ERROR IN OPENING CORE.DEF !!!!'
!      STOP 'CORE.DEF'
! 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
!      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS,
!     *'  FORM:',FORM
!      STOP 'OPEN FAILED'
! 8001 CONTINUE
!      CLOSE (1)
      OLD=.FALSE.                                                       
      READ(20,1511) NATT,cform,IREL                                           
!      if(cform.eq.'NEW ') then
       assign 2021 to iform1  
!      else
!       assign 2020 to iform1  
!      end if
      RELA=.TRUE.                                                       
      IF(IREL.EQ.'NREL') RELA=.FALSE.                                   
      JATOM=0                                                           
      READ(8,2032) ISCF                                                      
!     READ FROM TAPE8=VSP TO MOVE THE FILEPOINTER                       
!fb
!     set of insv-label
!     force-calculation only if 
!     non-sph. potential exists, i.e., insv=0
      read(19,'(//)',iostat=insv)
!fe
!                                                                       
!.....READ TAPE22=SCFDATA UNTIL BOTTOM AND DEFINE ISCF                  
!                                                                       
! 6    READ(22,700,END=7)  LMARGN,RMARGN                                 
!      IF(LMARGN.NE.0)GO TO 11                                           
!      WRITE(21,4043)RMARGN                                              
!      GO TO 6                                                           
! 11   ISCF=LMARGN                                                       
!      WRITE(21,700)LMARGN,RMARGN                                        
!      GOTO 6                                                            
! 7    CONTINUE                                                          
!      WRITE(21,710)                                                     
!                                                                       
      WRITE(9,1960) ISCF                                                
      WRITE(9,*) '    NORM   OF CLM(R) = CLM*R*R/SQRT(4.*PI)  '         
      WRITE(9,*)                                                        
!                                                                       
!                                                                       
!                                                                       
  1   NSTOP=1                                                           
      CALL INSLD(RELA,npxx,cform)                                                  
      IF(NORB.EQ.0) GOTO 100                                            
 2    ITER=1                                                            
      DO 3 I=1,NP                                                       
      DO 3 J=1,NORB                                                     
      DGC(I,J)=0.                                                       
 3    DPC(I,J)=0.                                                       
      WRITE(6,1001) BAR                                                 
      WRITE (6,1002)                                                    
      N=-(ION+1)                                                        
 4    DO 5 I=1,NP                                                       
 5    D(I)=0.                                                           
!     fix H-case
      if(z.lt. 1.1) go to 4000
      TETS=TEST                                                         
      YMAX=0.                                                           
      VMAX=0.                                                           
      EMAX=0.                                                           
! RESOLUTION DE L EQUATION DE DIRAC POUR CHAQUE ORBITALE                
! **********************************************************************
      DO 25  J=1,NORB                                                   
      DE=DEN(J)                                                         
 16   CALL RESLD (NQN(J),NQL(J),NK(J),IMAX,DEN(J),DFL(J),DQ1(J),J,RELA) 
      IF (NSTOP.EQ.0) GO TO 18                                          
      IF (NSTOP.NE.362.OR.ITER.GE.10.OR.TETS.GT.TEST) GO TO 17          
      TETS=TESTV                                                        
      GO TO 16                                                          
 17   WRITE (6,1010) NSTOP,NQN(J),TITRE(J)                              
      GO TO 1                                                           
 18   VAL=ABS((DEN(J)-DE)/DE)                                           
      IF(VAL.GT.EMAX) EMAX=VAL                                          
      NMAX(J)=IMAX                                                      
      DO 24  I=1,NP                                                     
      VAL=DGC(I,J)-DP(I)                                                
      IF (ABS(DP(I)).GT.1.) VAL=VAL/DP(I)                               
      IF (ABS(VAL).LT.ABS(YMAX)) GO TO 21                               
      YMAX=VAL                                                          
      Y=DP(I)                                                           
      YN=DGC(I,J)                                                       
 21   VAL=DPC(I,J)-DQ(I)                                                
      IF (ABS(DQ(I)).GT.1.) VAL=VAL/DQ(I)                               
      IF (ABS(VAL).LT.ABS(YMAX)) GO TO 23                               
      YMAX=VAL                                                          
      Y=DQ(I)                                                           
      YN=DPC(I,J)                                                       
 23   DGC(I,J)=DP(I)                                                    
      DPC(I,J)=DQ(I)                                                    
 24   D(I)=D(I)+NEL(J)*(DP(I)*DP(I)+DQ(I)*DQ(I))                        
 25   CONTINUE                                                          
      CALL POTSL (DC,D,DP,DR,DPAS,DEXV,Z,NP,ION,ICUT)                   
      IF (NUC.LE.0) GO TO 31                                            
      DO 29 I=1,NUC                                                     
 29   DC(I)=DC(I)+Z/DR(I)+Z*((DR(I)/DR(NUC))**2-3.)/(DR(NUC)+DR(NUC))   
 31   DO 33 I=1,NP                                                      
      DVAL=ABS(DC(I)-DV(I))                                             
      IF ((DR(I)*DC(I)).LE.N) DVAL=-DVAL/DC(I)                          
      IF (DVAL.LE.VMAX) GO TO 33                                        
      VMAX=DVAL                                                         
      J=I                                                               
 33   CONTINUE                                                          
      WRITE(6,1003) ITER,VMAX,DR(J),DV(J),DC(J),EMAX,YMAX,YN,Y          
      IF(TETS.LE.TEST.AND.EMAX.LE.TESTE.AND.VMAX.LE.TESTV.AND.YMAX.LE.TE &
      STY) GO TO 61                                                     
      ITER=ITER+1                                                       
      IF(ITER.LE.NITER) GO TO 35                                        
      WRITE(6,1004) NITER                                               
      NSTOP=2                                                           
      GO TO 61                                                          
! POTENTIEL POUR L ITERATION SUIVANTE                                   
! **********************************************************************
 35   IF(ITER.EQ.2) GO TO 41                                            
      IF (IPRAT) 41,45,41                                               
 41   DVAL=1.-DCOP                                                      
      DO 43 I=1,NP                                                      
      DVN(I)=DV(I)                                                      
      DVF(I)=DC(I)                                                      
 43   DV(I)=DVAL*DV(I)+DCOP*DC(I)                                       
      IF (.NOT.OLD) GO TO 44                                            
      READ(10,END=431) DV                                               
431   OLD=.FALSE.                                                       
      GO TO 4                                                           
44    REWIND 10                                                         
      WRITE(10) DV                                                      
      GO TO 4                                                           
 45   DO 55  I=1,NP                                                     
      DVAL=DALP(DVN(I),DVF(I),DV(I),DC(I))                              
      DVN(I)=DV(I)                                                      
      DVF(I)=DC(I)                                                      
 55   DV(I)=DVAL*DV(I)+(1.-DVAL)*DC(I)                                  
      IF (.NOT.OLD) GO TO 56                                            
      READ(10,END=551) DV                                               
551   OLD=.FALSE.                                                       
      GO TO 4                                                           
56    REWIND 10                                                         
      WRITE(10) DV                                                      
      GO TO 4                                                           
 61   WRITE (6,1001) BAR                                                
!     READ (5,3001) I                                                   
 3001 FORMAT (I3)                                                       
      WRITE(6,3004) (DR(KK),D(KK),KK=1,NP)                              
 3004 FORMAT(//' RADII AND CHARGE DENS. ALTERNATELY'/(6(D15.8)))        
!     IF (I.LE.0) GO TO 63                                              
 4000 JATOM=JATOM+1                                                     
      WRITE(21,720) JATOM,BAR,NORB                                      
      WRITE(9,1990) JATOM,BAR                                           
      WRITE(9,2001) 1                                                   
      WRITE(9,2002) ( NQN(JSTATE),TITRE(JSTATE), JSTATE=1,NORB )        
      WRITE(9,1980)                                                     
      WRITE(9,2011) 0,0                                                 
      NP40=NP-npxx                                                        
      WRITE(9,iform1) ( D(KK), KK=1,NP40)                                 
      WRITE(9,2030)                                                     
      WRITE(9,2031)                                                     
      WRITE (6,1005)                                                    
!fb
      if (force.and.insv.eq.0) call fcore(jatom,np40,dr,d)
!fe
!
! VALEURS MOYENNES DE R                                                 
! **********************************************************************
      DO 67 I=1,NP                                                      
      DVF(I)=DC(I)                                                      
 67   DQ(I)=0.                                                          
      DVAL=0.                                                           
      DO 91 I=1,NORB                                                    
      IM=NMAX(I)                                                        
      DVAL=DVAL+NEL(I)*DEN(I)                                           
      DO 73  J=1,IM                                                     
 73   DC(J)=DGC(J,I)*DGC(J,I)+DPC(J,I)*DPC(J,I)                         
      L=5                                                               
      IF (IABS(NK(I)).EQ.1) L=L-1                                       
      DO 88 J=1,L                                                       
      DP(J)=DFL(I)+DFL(I)                                               
      IF (J-2) 74,75,76                                                 
 74   N=4                                                               
      GO TO 88                                                          
 75   N=2                                                               
      GO TO 88                                                          
 76   IF (J-4) 77,78,79                                                 
 77   N=1                                                               
      GO TO 88                                                          
 78   N=-1                                                              
      GO TO 88                                                          
 79   N=-3                                                              
 88   CALL SOMM (DR,DC,DQ,DPAS,DP(J),N,IM)                              
!.... CONVERT TO RYD                                                    
      DENXX=DEN(I)*2.
!cccccc
!  shift the eigenvalues back back by 1 Ry
!      denxx=denxx+1.0d0
!cccccc
      if(nel(i).eq.0) denxx=0.d0                    
      WRITE (RMARGN,'(a2)') TITRE(i)
      IF (RMARGN(2:2).EQ.'*') RMARGN(2:2)=RMARGN(1:1)
      WRITE(21,730) NQN(I),RMARGN,JATOM,NQN(i),TITRE(I),DENXX   
 91   WRITE (6,1006)NQN(I),TITRE(I),DEN(I),(DP(J),J=1,L)                
! ENERGIE TOTALE EN MOYENNE SPHERIQUE                                   
! **********************************************************************
      DC(1)=1                                                           
      DO 95 I=1,NP                                                      
 95   DP(I)=D(I)/DR(I)                                                  
      IF (NUC.LE.0) GO TO 99                                            
      DO 97 I=1,NUC                                                     
 97   DP(I)=D(I)*(3.-DR(I)*DR(I)/(DR(NUC)*DR(NUC)))/(DR(NUC)+DR(NUC))   
      DC(1)=4                                                           
 99   CALL SOMM (DR,DP,DQ,DPAS,DC(1),0,NP)                              
      DO 105  I=1,NP                                                    
      DP(I)=D(I)*DVF(I)                                                 
 105  D(I)=D(I)*((D(I)*DR(I))**(1./3.))                                 
      DC(2)=3                                                           
      DC(3)=1                                                           
      IF (NUC.NE.0) DC(3)=4                                             
      CALL SOMM (DR,DP,DQ,DPAS,DC(3),0,NP)                              
      CALL SOMM (DR,D,DQ,DPAS,DC(2),-1,NP)                              
      DC(2)=-3.*DC(2)/(105.27578**(1./3.))                              
      DC(1)=-Z*DC(1)                                                    
      DC(4)=DVAL-DC(3)                                                  
      DVAL=DVAL+(DC(1)-DC(3)+(DEXE-DEXV)*DC(2))/2.                      
      DC(3)=(DC(3)-DC(1)-DEXV*DC(2))/2.                                 
      DC(2)=DC(2)*DEXE/2.                                               
      WRITE (6,1007) DVAL,DC(4),DC(3),DC(2),DC(1)                       
      IF(NORB.EQ.1) GO TO 205                                           
      WRITE (6,1001) BAR                                                
      WRITE (6,1008)                                                    
! INTEGRALES DE RECOUVREMENT                                            
! **********************************************************************
      DO 115  I=2,NORB                                                  
      K=I-1                                                             
      DO 115  J=1,K                                                     
      IF(NQL(I).NE.NQL(J).OR.NK(I).NE.NK(J)) GO TO 115                  
      IM=NMAX(J)                                                        
      IF(NMAX(I).LT.IM) IM=NMAX(I)                                      
      DO 111  L=1,IM                                                    
      DQ(L)=DPC(L,I)*DPC(L,J)                                           
 111  DC(L)=DGC(L,I)*DGC(L,J)                                           
      DVAL=DFL(I)+DFL(J)                                                
      CALL SOMM (DR,DC,DQ,DPAS,DVAL,0,IM)                               
      WRITE (6,1009) NQN(I),TITRE(I),NQN(J),TITRE(J),DVAL               
 115  CONTINUE                                                          
!205  CALL CDSLD (TITRE,BAR)                                            
 205  CONTINUE                                                        
  100 NATT=NATT-1                                                      
      IF(NATT.EQ.0) THEN
         CALL ERRCLR(ERRFN)
         STOP ' CORE  END'
      ENDIF                                  
      GO TO 1                                                           
!
!        error handling
!
  910 INFO = 1
!
!        'core.def' couldn't be opened
!
      WRITE (ERRMSG,9000) FNAME
      CALL OUTERR('CORE',ERRMSG)
      GOTO 999
  920 INFO = 2
!
!        file FNAME couldn't be opened
!
      WRITE (ERRMSG,9010) IUNIT
      CALL OUTERR('CORE',ERRMSG)
      WRITE (ERRMSG,9020) FNAME
      CALL OUTERR('CORE',ERRMSG)
      WRITE (ERRMSG,9030) STATUS, FORM
      CALL OUTERR('CORE',ERRMSG)
      GOTO 999
  960 INFO = 7
!
!        Error reading file 'core.def'
!
      WRITE (ERRMSG,9040) FNAME
      CALL OUTERR('CORE',ERRMSG)
      GOTO 999
  999 STOP 'CORE - Error'
!                                                                       
!                                                                       
 700  FORMAT(I3,A77)                                                    
 710  FORMAT(//,5X,'CORE STATES CONTRIBUTION TO TOTAL ENERGY',           &
             '   CALCULATED BY CORE',/)                                 
 720  FORMAT(/7X,I2,'.ATOM ',5X,A10,11X,I2,' CORE STATES')              
 730  FORMAT(':',I1,A2,i2.2,': ',i1,a2,9X,F20.6,' Ry') 
 1000 FORMAT(L7)                                                        
 1001 FORMAT (1H1,40X,A10/)                                             
 1002 FORMAT (' ITER',4X,'DVMAX',10X,'R',14X,'VN-1',13X,'VN',10X,'DEMAX' &
      ,6X,'DPMAX',9X,'PN-1',13X,'PN')                                   
 1003 FORMAT (I5,1PE11.2,3(1PE16.6),2(1PE11.2),2(1PE16.6))              
 1004 FORMAT (' NOMBRE D ITERATIONS SUPERIEUR A',I4)                    
 1005 FORMAT (12X,'ENERGIE',12X,'(R4)',14X,'(R2)',14X,'(R)',15X,'(R-1)', &
      13X,'(R-3)'/)                                                     
 1006 FORMAT (I3,A2,6(1PE18.7))                                         
 1007 FORMAT (1H0,5X,'ET=',1PE14.7,5X,'EC=',1PE14.7,5X,'EE=',1PE14.7,5X, &
      'EX=',1PE14.7,5X,'EN=',1PE14.7)                                   
 1008 FORMAT (1H0,47X,'INTEGRALES DE RECOUVREMENT'/)                    
 1009 FORMAT (34X,I1,A2,I3,A2,F19.7)                                    
 1010 FORMAT ('  NSTOP=',I4,'  POUR L ORBITALE',I3,A2)                  
 1511 FORMAT(/,28X,I2,1x,a4,/,13X,A4,/)                                       
 1960 FORMAT(3X,'SPHERICAL CHARGE DENSITY OF CORESTATES',5X,             &
             'GENERATED BY',I3,' ITERATION' )                           
 1980 FORMAT(3X)                                                        
 1990 FORMAT(3X,'ATOMNUMBER  =',I2,5X,A10)                              
 2001 FORMAT(3X,'NUMBER OF LM=',I2)                                     
 2002 FORMAT(3X,'CORE STATES =', 25(I2,A2) )                            
 2011 FORMAT(3X,'CLM(R) FOR L=',I2,3X,'M=',I2/)                         
 2020 FORMAT(3X,5E13.7)                                                 
 2021 FORMAT(3X,4E19.12)                                                 
 2030 FORMAT(///)                                                       
 2031 FORMAT(/)                                                         
 2032 FORMAT(50X,I2,//) 
 4043 FORMAT(3X,A77)    
 9000 FORMAT('can''t open definition file ',A40)
 9010 FORMAT('can''t open unit: ',I2)
 9020 FORMAT('       filename: ',A50)
 9030 FORMAT('         status: ',A,'  form: ',A)
 9040 FORMAT('Error reading file: ',A47)
      END                                                               

