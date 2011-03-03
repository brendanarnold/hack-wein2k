      PROGRAM LSTART                                                    
!     PROGRAM LSTART(INPUT,OUTPUT,TAPE5=INPUT,TAPE6=OUTPUT,TAPE7,POTDN, 
!    * POTUP,TAPE10=POTUP,TAPE11=POTDN,tape13=sigma)                          
! HARTREE FOCK DIRAC SLATER          J P DESCLAUX     CEA PARIS 1969    
! MODIFIE JUILLET 1970                                                  
!   HEDIN LUNDQUIST PARAMETERS INCLUDET                                 
! **********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      character*1 aplot(30,2)                                             
      CHARACTER*4 IREL,cform                                       
      CHARACTER*80 NAME,fname
      CHARACTER*11 STATUS,FORM					
!
      COMMON DEN(30,2),DQ1(30),DFL(30),NQN(30),NQL(30),NK(30),NMAX(30),  &
      ZEL(30,2),NORB                                                    
      COMMON/DIRA/DV(NPT,2),DR(NPT),DP(NPT),DQ(NPT),DPAS,Z,NSTOP,NES,    &
      TETS,NP,NUC                                                       
      COMMON/PS2/DEXV,DEXE,DCOP,TITRE(30),BAR(10),TEST,TESTE,TESTY,TESTV &
      ,NITER,ION,ICUT,IPRAT                                             
      COMMON/DEUX/DVN(NPT,2),DVF(NPT,2),D(NPT,2),DC(NPT,2),DGC(NPT,30,2) &
      ,DPC(NPT,30,2)
      REAL*8         :: dplot(npt,2),dplot2(npt,2)
      DIMENSION DHL(NPT,2),DEXC(NPT,2)                                  
      DIMENSION VLS(2)
      dimension dcsp(npt,2),ddcsp(npt,2)
        DIMENSION dsh(2,30,npt)
      LOGICAL OLD,RELA,lcore(30)      
      data wats/4HWats/                   
      pi=acos(-1.d0)
!
 1000 FORMAT(L7,L7)                                                     
 1001 FORMAT (/,20X,10A4/)                                            
 1002 FORMAT (' ITER',4X,'DVMAX',10X,'R',14X,'VN-1',13X,'VN',10X,'DEMAX' &
      ,6X,'DPMAX',9X,'PN-1',13X,'PN')                                   
 1003 FORMAT (I5,1PE11.2,3(1PE16.6),2(1PE11.2),2(1PE16.6),               &
      /,30(i3,a2,2F15.5,/))                                                   
 1004 FORMAT (' NUMBER OF ITERATIONS GREATER THAN',I4)                  
 1005 FORMAT(9X,'OCCUPANCY    ',                                         &
              'ENERGY(RYD)',9X,'(R4)',14X,'(R2)',14X,'(R)',15X,'(R-1)',  &
      13X,'(R-3)'/)                                                     
 1006 FORMAT (I3,A2,F10.3,6(1PE18.7))                                   
 1007 FORMAT (/,1X,'TOTAL ENERGY (RYD):',F24.6,/,5X,'SUM OF EI:',     &
      1PE14.7,5X,'NUC:',1PE14.7,5X,'COUL:',1PE14.7,/,5X,'V-XC SPIN 1:',   &
      1PE14.7,5X,'E-XC SPIN 1:',1PE14.7,/,5X,'V-XC SPIN 2:',1PE14.7,       &
      5X,'E-XC SPIN 2:',1PE14.7)                                        
 1008 FORMAT (/,47X,'ORTHOGONALITY INTEGRALS'/)                       
 1009 FORMAT (34X,I1,A2,I3,A2,F19.7)                                    
 1010 FORMAT ('  NSTOP=',I4,'  FOR THE ORBITAL',I3,A2)                  
 1011 FORMAT (1P,5E14.7)                                                
 1511 FORMAT(/,28X,I2,1x,a4,/,13X,A4,/)                                 
!                                                                       
!      CALL XUFLOW(0)                                                   
!                                                                       
      write(*,*) ' SELECT XCPOT:'
      write(*,*) ' recommended: 13: GGA (Perdew-Burke-Ernzerhof 96)'
      write(*,*) '               5: LSDA'
      write(*,*) '              11: GGA (Wu-Cohen 2006)'
      read(*,*) iex
 7999 write(*,*) ' SELECT ENERGY to separate core and valence states:'
      write(*,*) ' recommended: -6.0 Ry (check how much core charge', &
      ' leaks out of MT-sphere)'
      read(*,*) ecore
      if(ecore.gt.-2.0d0.or.ecore.lt.-10.0) goto 7999
!
      iarg=iargc()
      call getarg(iarg,fname)
!      call getarg(2,fname)
!      if(fname.eq.'      ') call getarg(1,fname)
      OPEN(1,FILE=fname,STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,NAME,STATUS,FORM
      OPEN(IUNIT,FILE=NAME,STATUS=STATUS,FORM=FORM,ERR=8002)
!      rewind iunit
      GOTO 8003
 8000 WRITE(*,*) ' ERROR IN OPENING LSTART.DEF !!!!'
      STOP 'LSTART.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',NAME,'  STATUS: ',STATUS,'  FORM:' &
       ,FORM
      STOP 'OPEN FAILED'
 8001 CONTINUE
      CLOSE (1)
      READ(20,1511) NATT,cform,IREL                                     
!      if(cform.eq.'NEW ') then
        assign 2021 to iform1
!      else
!        assign 2020 to iform1
!      end if                                    
 2020 FORMAT(3X,5E13.7)                                                 
 2021 FORMAT(3X,4E19.12)                                                
      RELA=.TRUE.                                                       
      IF(IREL.EQ.'NREL') RELA=.FALSE.                                   
      OLD=.FALSE.                                                       
      JATOM=0  
      valel=0.d0
!.....start writing some input files
      write(18,8102)
      write(19,8103)
      write(15,8100) 
      write(15,8101) 
 8100 format('WFFIL        (WFPRI, SUPWF)')
 8101 format('  7.00       10    4 (R-MT*K-MAX; MAX L IN WF, V-NMT' )
 8102 format('BROYD  0.0   YES  (BROYD/PRATT, extra charge', &
      ' (+1 for additional e), norm)',/, &
      '0.40            mixing FACTOR for BROYD/PRATT scheme',/, &
      '0.10  1.00      PW and CLM-scaling factors')
 8103 format('PRATT  0.0   YES   (BROYD/PRATT, extra charge', &
      ' (+1 for additional e), norm)',/, '0.10      FACTOR')
!
      WRITE(7,*) '  '                                                         
      WRITE(11,*) '    TOTAL SPHERICAL POTENTIAL IN MT SPHERES'         
      WRITE(11,*) '    NORM  OF V(0,0) = V(0,0)*R/SQRT(4.*PI)'          
      WRITE(11,*)                                                       
      WRITE(12,*) '    TOTAL SPHERICAL POTENTIAL IN MT SPHERES'         
      WRITE(12,*) '    NORM  OF V(0,0) = V(0,0)*R/SQRT(4.*PI)'          
      WRITE(12,*)                                                       
  1   NSTOP=1                                                           
      JATOM=JATOM+1   
!                                                  
      CALL INSLD(IEX,JSPIN,RELA,aplot,valel,mult,np0,ecore,&
                RWatson,EWatson) 
!
 2    ITER=1                                                            
      DO 3 I=1,NP                                                       
      DO 3 J=1,NORB                                                     
      DO 3 ISPIN=1,JSPIN                                                
      DGC(I,J,ISPIN)=0.                                                 
 3    DPC(I,J,ISPIN)=0.                                                 
      WRITE(6,1001) BAR
      IF (BAR(2).eq.Wats) then
         write(6,*) 'Watson radius=',RWatson
         write(6,*) 'Watson charge=', EWatson       
         write(6,*)
      endif
                   
      WRITE (6,1002)                                                    
      N=-(ION+1)                                                        
4     dplot(1:np,1:jspin)=0.d0
      dplot2(1:np,1:jspin)=0.d0
      d(1:np,1:jspin)=0.d0
      TETS=TEST                                                         
      YMAX=0.                                                           
      VMAX=0.                                                           
      EMAX=0.

      do ispin=1,jspin 
      do i=1,np         
         DV(I,ispin)=DV(i,ispin)!+ watsonw(dr(i),RWatson,EWatson)
      enddo    
      enddo

! RESOLUTION DE L EQUATION DE DIRAC POUR CHAQUE ORBITALE                
! **********************************************************************
      DO 25  J=1,NORB                                                   
      DO 25  ISPIN=1,JSPIN                                              
      DE=DEN(J,ISPIN)
 16   CALL RESLD (NQN(J),NQL(J),NK(J),IMAX,DEN(J,ISPIN),DFL(J),DQ1(J)    &
      ,J,ISPIN,RELA)  
      IF (NSTOP.EQ.0) GO TO 18                                          
      IF (NSTOP.NE.362.OR.ITER.GE.10.OR.TETS.GT.TEST) then
         write(6,*) 'nstop,iter,tets,test',nstop,iter,tets,test
         goto 17
      endif
      TETS=TESTV                                                        
      GO TO 16                                                          
 17   WRITE (6,1010) NSTOP,NQN(J),TITRE(J)                              
      GO TO 1                                                           
 18   VAL=ABS((DEN(J,ISPIN)-DE)/DE)                                     
      IF(VAL.GT.EMAX) EMAX=VAL                                          
      NMAX(J)=IMAX                                                      
      DO 24  I=1,NP                                                     
      VAL=DGC(I,J,ISPIN)-DP(I)                                          
      IF (ABS(DP(I)).GT.1.) VAL=VAL/DP(I)                               
      IF (ABS(VAL).LT.ABS(YMAX)) GO TO 21                               
      YMAX=VAL                                                          
      Y=DP(I)                                                           
      YN=DGC(I,J,ISPIN)                                                 
 21   VAL=DPC(I,J,ISPIN)-DQ(I)                                          
      IF (ABS(DQ(I)).GT.1.) VAL=VAL/DQ(I)                               
      IF (ABS(VAL).LT.ABS(YMAX)) GO TO 23                               
      YMAX=VAL                                                          
      Y=DQ(I)                                                           
      YN=DPC(I,J,ISPIN)                                                 
 23   DGC(I,J,ISPIN)=DP(I)                                              
      DPC(I,J,ISPIN)=DQ(I) 
      if (aplot(j,ispin).EQ.'P') THEN
         dplot(I,ISPIN)=Dplot(I,ISPIN)+ZEL(J,ISPIN)* &
              (DP(I)*DP(I)+DQ(I)*DQ(I)) 
         dplot2(i,ispin)=DP(I)
      ENDIF
      DSH(ISPIN,J,I)=ZEL(J,ISPIN)*(DP(I)*DP(I)+DQ(I)*DQ(I))
 24   D(I,ISPIN)=D(I,ISPIN)+ZEL(J,ISPIN)*(DP(I)*DP(I)+DQ(I)*DQ(I))      
 25   CONTINUE                                         
      CALL POTSL (DC,D,DP,DQ,DR,DPAS,IEX,DEXV,Z,NP,ION,ICUT,JSPIN)      
      IF (NUC.LE.0) GO TO 31                                            
      DO 29 I=1,NUC                                                     
      DC(I,1)=DC(I,1)+Z/DR(I)+Z*((DR(I)/DR(NUC))**2-3.)/(DR(NUC)+        &
      DR(NUC))                                                          
 29   DC(I,2)=DC(I,2)+Z/DR(I)+Z*((DR(I)/DR(NUC))**2-3.)/(DR(NUC)+        &
      DR(NUC))                                                          
 31   DO 33 I=1,NP                                                      
      DO 33 ISPIN=1,JSPIN 
      DC(I,ISPIN)=DC(I,ISPIN)+watsonw(dr(i),RWatson,EWatson)
      DVAL=ABS(DC(I,ISPIN)-DV(I,ISPIN))
      IF ((DR(I)*DC(I,ISPIN)).LE.N) DVAL=-DVAL/DC(I,ISPIN)              
      IF (DVAL.LE.VMAX) GO TO 33                                        
      VMAX=DVAL                                                         
      J=I                                                               
      KSPIN=ISPIN                                                       
 33   CONTINUE                                                          
      WRITE(6,1003)ITER,VMAX,DR(J),DV(J,KSPIN),DC(J,KSPIN),EMAX,YMAX,YN  &
      ,Y,(NQN(I1),TITRE(I1),DEN(I1,1)/5.,DEN(I1,2)/5.,I1=1,NORB)                          
      IF(TETS.LE.TEST.AND.EMAX.LE.TESTE.AND.YMAX.LE.TESTY                &
      ) GO TO 61                                                        
      ITER=ITER+1                                                       
      write (6,*) iter,niter
      IF(ITER.LE.NITER) GO TO 35                                        
      WRITE(6,1004) NITER                                               
      GO TO 61                                                          
! POTENTIEL POUR L ITERATION SUIVANTE                                   
! **********************************************************************
 
 35   IF(ITER.EQ.2) GO TO 41                                            
      IF (IPRAT) 41,45,41                                               
 41   DVAL=1.-DCOP                                                      
      DO 43 ISPIN=1,JSPIN                                              
        DO 43 I=1,NP                                                      
          DVN(I,ISPIN)=DV(I,ISPIN)                                          
          DVF(I,ISPIN)=DC(I,ISPIN)                                          
 43   DV(I,ISPIN)=DVAL*DV(I,ISPIN)+DCOP*DC(I,ISPIN)
!      write (6,*) 'nach 43,old',old
      IF (.NOT.OLD) GO TO 44                                            
      READ(10,END=431) DV                                               
431   OLD=.FALSE.                                                       
      GO TO 4                                                           
44    REWIND 10                                                         
      WRITE(10) DV                                                      
      GO TO 4                                                           
 45   DO 55  I=1,NP                 
      DO 55 ISPIN=1,JSPIN                                               
      DVAL=DALP(DVN(I,ISPIN),DVF(I,ISPIN),DV(I,ISPIN),DC(I,ISPIN))      
      DVN(I,ISPIN)=DV(I,ISPIN)                                          
      DVF(I,ISPIN)=DC(I,ISPIN)                                            
 55   DV(I,ISPIN)=DVAL*DV(I,ISPIN)+(1.-DVAL)*DC(I,ISPIN)
      IF (.NOT.OLD) GO TO 56                                            
      READ(10,END=551) DV                                               
551   OLD=.FALSE.                                                       
      GO TO 4                                                           
56    REWIND 10                                                         
      WRITE(10) DV                                                      
      GO TO 4                                                           
 61   WRITE (6,1001) BAR                                                
      write(13,3001) name,np,dr(1),dpas                                 
      write(13,3002) (dplot(KK,1)+dplot(KK,2),kk=1,np)
      DO kk=1,np
!         WRITE(94,*) dr(kk),dplot(KK,1)+dplot(KK,2),(dplot(KK,1)+dplot(KK,2))/(dr(kk)*dr(kk))
         WRITE(94,*) dr(kk),dplot2(KK,2)/dr(kk),dplot2(KK,2)
      ENDDO
 3001 format(a80/i4,2E12.6) 
 3002 FORMAT(5D16.9)                                                    
!     WRITE(7,3003) BAR                                                 
!3003 FORMAT(10A4)                                                      
!     WRITE(7,3002) (DR(KK),D(KK,1),D(KK,2),KK=1,NP)                    
! 3003 FORMAT(4D20.12)
 3003 FORMAT(3X,4E19.12)                                                 
 3004 FORMAT(I5)
      write(7,3004) NP
      write(7,3003) (d(kk,1)+d(kk,2),kk=1,np)
 63   CONTINUE                                                          
         REWIND 10                                                      
      DO 4545 ISPIN=1,JSPIN                                             
      RVVV=DV(NP0,ISPIN)*DR(NP0)*2.                                 
 4545 WRITE(6,4547) ISPIN,DR(NP0),RVVV,DV(NP0,ISPIN)*2.             
 4547 FORMAT(' FOR SPIN',I3,'   R-MT=',F10.5,'   R*V=',F10.5,            &
       '  V=',F10.5)                                                    
!      READ(5,*) SHIFT                                                   
       shift=0.d0
!      WRITE(6,1991) SHIFT                                               
 1991 FORMAT(20X,' ARBITRARY SHIFT OF POTENTIAL: ',F10.5)               
      SHIFT=SHIFT/2.                                                    
!                                                                       
!bk941018: output of density for dstart
      IT=80                                                             
      DO 4040 ISPIN=1,JSPIN                                             
          IT=IT+1                                                           
          WRITE(IT,1990)JATOM ,NP,z                                           
          WRITE(IT,2001) 1                                               
          WRITE(IT,2012) 0,0                                             

          WRITE(IT,iform1) (D(KK,ISPIN),KK=1,NP)                     
          WRITE(IT,2031)                                                 
          WRITE(IT,2030) 
 4040     continue
          IT=IT+1                                                           
          WRITE(IT,1990)JATOM ,NP,z                                           
          WRITE(IT,2001) 1                                               
          WRITE(IT,2012) 0,0                                             

          WRITE(IT,iform1) (D(KK,1)+D(KK,2),KK=1,NP)                     
          WRITE(IT,2031)                                                 
          WRITE(IT,2030) 
 2012 FORMAT(3X,'RLM(R) FOR L=',I2,3X,'M=',I2/)                         
! jmz addtion for shell by shell densities
!        itt=90+jatom
!        do jshell=1,norb
!        write(itt,911) nint(z),np,dr(1),dpas
!        write(itt,3002) (dsh(1,jshell,kk)+dsh(2,jshell,kk),kk=1,np)
!        end do
!911     FORMAT(1X,2I8,2E16.8)
!        do kk=1,np
!        write(itt,3005) dr(1)*exp((kk-1)*dpas),(dsh(1,jshell,kk)+
!     *                   dsh(2,jshell,kk),jshell=1,norb)
!        end do
! 3005   format(f12.8,10d13.7)
!                                                                       
      DO 20 ISPIN=1,JSPIN                                               
         DO 20 KK=1,NP                                                  
 20      DV(KK,ISPIN)=(DV(KK,ISPIN)+SHIFT)*DR(KK)*2.                    
!    CONVERT POTENTIAL TO R*V (RYD), SHIFT ARBITRARILY TO               
!       REASONABLE ZERO                                                 
      IT=10                                                             
      DO 4546 ISPIN=1,JSPIN                                             
      IT=IT+1                                                           
         WRITE(IT,1990)JATOM                                            
         WRITE(IT,2001) 1                                               
         WRITE(IT,2011) 0,0                                             
         WRITE(IT,iform1) (DV(KK,ISPIN),KK=1,NP0)                     
         WRITE(IT,2031)                                                 
 4546    WRITE(IT,2030)                                                 
 1990 FORMAT(3X,'ATOMNUMBER =',I3,5X,'RADIAL MESH POINTS: ',i5,f9.3)          
 2001 FORMAT(3X,'NUMBER OF LM',I3//)                                   
 2011 FORMAT(3X,'VLM(R) FOR L',I3,3X,'M=',I2/)                         
 2030 FORMAT(///)                                                       
 2031 FORMAT(/)                                                         
 8888    FORMAT(3X,5E13.7)                                              
      WRITE(6,1001) BAR                                                 
      WRITE (6,1005)                                                    
! VALEURS MOYENNES DE R                                                 
! **********************************************************************
      DO 67 ISPIN=1,JSPIN                                               
      DO 67 I=1,NP                                                      
      DVF(I,ISPIN)=DC(I,ISPIN)                                          
 67   DQ(I)=0.                                                          
      DVAL=0.                                                           
      DO 91 I=1,NORB                                                    
      DO 91 ISPIN=1,JSPIN                                               
      IM=NMAX(I)                                                        
      DVAL=DVAL+ZEL(I,ISPIN)*DEN(I,ISPIN)                               
      DO 73  J=1,IM                                                     
 73   DC(J,ISPIN)=DGC(J,I,ISPIN)*DGC(J,I,ISPIN)+DPC(J,I,ISPIN)*DPC(J,I,  &
      ISPIN)                                                            
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
 88   CALL SOMM (DR,DC(1,ISPIN),DQ,DPAS,DP(J),N,IM)                     
 91   WRITE(6,1006)NQN(I),TITRE(I),ZEL(I,ISPIN),DEN(I,ISPIN)*2.,         &
       (DP(J),J=1,L)                                                    
!
      call write1(norb,nqn,nk,den,zel,valel,mult,lcore,ecore)
! ENERGIE TOTALE EN MOYENNE SPHERIQUE                                   
! **********************************************************************
      DC(1,1)=1.                                                        
      DO 95 I=1,NP                                                      
      DP(I)=0.                                                          
      DO 95 ISPIN=1,JSPIN                                               
 95   DP(I)=DP(I)+D(I,ISPIN)/DR(I)                                      
      IF (NUC.LE.0) GO TO 99                                            
      DO 97 I=1,NUC                                                     
 97   DP(I)=(D(I,1)+D(I,2))*(3.-DR(I)*DR(I)/(DR(NUC)*DR(NUC)))/(DR(NUC)  &
      +DR(NUC))                                                         
      DC(1,1)=4                                                         
 99   CALL SOMM (DR,DP,DQ,DPAS,DC(1,1),0,NP)                            
      WRITE(6,*)                                                        
      DO 100 ISPIN=1,JSPIN                                              
      CTEST=1.                                                          
      CALL SOMM (DR,D(1,ISPIN),DQ,DPAS,CTEST,0,NP)                      
      WRITE(6,*) 'TOTAL CHARGE FOR SPIN ',ISPIN,':               ',CTEST   
      CTEST=1.                                                          
      CALL SOMM (DR,D(1,ISPIN),DQ,DPAS,CTEST,0,NP0)                      
  100 WRITE(6,*) 'TOTAL CHARGE FOR SPIN ',ISPIN,' INSIDE SPHERE: ',CTEST       
      DO 9100 ISPIN=1,JSPIN                                              
      CTEST=1.                                                          
      CALL SOMM (DR,Dplot(1,ISPIN),DQ,DPAS,CTEST,0,NP)                      
      WRITE(6,*) 'TOTAL CHARGE in sigma FOR SPIN ',ISPIN, &
                 ':              ',CTEST          
      CTEST=1.                                                          
      CALL SOMM (DR,Dplot(1,ISPIN),DQ,DPAS,CTEST,0,NP0)                   
 9100 WRITE(6,*) 'TOTAL CHARGE in sigma FOR SPIN ',ISPIN, &
                 'INSIDE SPHERE: ',CTEST          
!
      DO 9151 j=1,NP                                                      
 9151 DP(j)=0.d0                                                          
      do 9120 i=1,norb
        if(lcore(i)) then
        DO 9150 j=1,NP                                                      
 9150   DP(j)=dp(j)+dsh(1,i,j)+dsh(2,i,j)                                 
        endif
 9120 continue
      CTEST=1.d0                                                          
      CALL SOMM (DR,Dp,DQ,DPAS,CTEST,0,NP)                      
      WRITE(6,'('' TOTAL CORE-CHARGE:                '',F12.6)') CTEST
      CTEST1=1.                                                          
      CALL SOMM (DR,Dp,DQ,DPAS,CTEST1,0,NP0)                   
      WRITE(6,'('' TOTAL CORE-CHARGE INSIDE SPHERE:  '',F12.6)') CTEST1
      WRITE(6,'('' TOTAL CORE-CHARGE OUTSIDE SPHERE: '',F12.6)') CTEST-CTEST1
      if(ctest-ctest1.gt.0.01d0) then
       write(6,'(/,''WARNING:'',F12.6,'' CORE electrons leak out of MT-sphere !!!!'')') ctest-ctest1
        
       print*, 'WARNING:',ctest-ctest1,' CORE electrons leak out of MT-sphere !!!!'
      endif  
!
! calculate exchange correlation energy
!
      if(iex.gt.10) then
!
! with gradients 
!
! derivatives dCLM/dR and ddCLM/ddR
!
        call drho(np,d,dcsp,ddcsp,dr)
!
!      write(6,*) 'd',(d(ir,1),ir=371,381)
!      write(6,*) 'dcsp',(dcsp(ir,1),ir=371,381)
!      write(6,*) 'ddcsp',(ddcsp(ir,1),ir=371,381)
!
       endif
!
!
      if(iex.lt.10) then
        DO 105  I=1,NP                                                    
          DP(I)=0.                                                          
          SIGMA=0.                                                          
          DO 101 ISPIN=1,JSPIN                                              
            DP(I)=DP(I)+DVF(I,ISPIN)                                          
 101      SIGMA=SIGMA+D(I,ISPIN)                                            
          FPRHO=SIGMA/DR(I)/DR(I)                                           
          IF(FPRHO.LT.1.E-18) THEN                                          
            DEXC(I,1)=0.0d0                                                     
            DHL(I,1)=0.d0                                                       
            DEXC(I,2)=0.0d0                                                     
            DHL(I,2)=0.d0                                                       
            GOTO 105                                                          
          END IF                                                            
          RS=(3.d0/FPRHO)**(1.d0/3.d0) 
          DO 102 ISPIN=1,JSPIN                                              
            FPRHOS=D(I,ISPIN)/DR(I)/DR(I)                                     
            DEXC(I,ISPIN)=0.d0                                                  
            DHL(I,ISPIN)=0.d0                                                   
            DHL1=0.d0                                                           
            IF(FPRHOS.LT.1.E-18) GOTO 102                                     
            XI=FPRHOS/FPRHO                                                   
            IF(IEX.EQ.5) THEN
               ZET=(D(I,1)-D(I,2))/SIGMA
               CALL CORLSD(RS,ZET,ECLSD,VLS(1),VLS(2),ECRS,ECZET,ALFC)
               DHL1=(VXCSP(FPRHO,FPRHOS,iex)+VLS(ISPIN))
               DEXC(I,ISPIN)=D(I,ISPIN)*(ECLSD &
                  -0.75D0*(0.75d0/PI/pi*2.0d0*FPRHOS)**(1.D0/3.D0))
	    ELSE
               DEXC(I,ISPIN)=EPSXC(RS,XI)*D(I,ISPIN) 
               DHL1=VXCSP(FPRHO,FPRHOS,iex)
            ENDIF
            DHL(I,ISPIN)=DHL1*D(I,ISPIN)                                      
  102     DP(I)=DP(I)-DHL1                                                  
          DP(I)=DP(I)*SIGMA/2.d0
! IF(I.LT.30) WRITE(6,*) DHL(I,1),DHL(I,2),DEXC(I,1),DEXC(I,2)     
  105   CONTINUE                                                          
      else
!
        DO 605  I=1,NP                                                    
          DP(I)=0.                                                          
          SIGMA=0.                                                          
          qrhou=d(i,1)
          qrhod=d(i,2)
          fu=qrhou/4.d0/pi
          fd=qrhod/4.d0/pi
          fu=fu/dr(i)/dr(i)
          fd=fd/dr(i)/dr(i)
          ft=fu+fd
!
          qrhou=dcsp(i,1)
          qrhod=dcsp(i,2)
          dndru=qrhou/4.d0/pi
          dndrd=qrhod/4.d0/pi
          dndrt=dndru+dndrd
!       
          qrhou=ddcsp(i,1)
          qrhod=ddcsp(i,2)
          d2ndru=qrhou/4.d0/pi
          d2ndrd=qrhod/4.d0/pi
          d2ndrt=d2ndru+d2ndrd
!
          dndtt=0.
          dndpt=0.
          d2ndtt=0.
          d2ndpt=0.
          dndrtt=0.
          dndrpt=0.
          dndtpt=0.
!
          dndtu=0.
          dndpu=0.
          d2ndtu=0.
          d2ndpu=0.
          dndrtu=0.
          dndrpu=0.
          dndtpu=0.
!
          dndtd=0.
          dndpd=0.
          d2ndtd=0.
          d2ndpd=0.
          dndrtd=0.
          dndrpd=0.
          dndtpd=0.
          qsin=0.
          qcos=0.
          DO 601 ISPIN=1,JSPIN                                              
            DP(I)=DP(I)+DVF(I,ISPIN)                                          
 601      SIGMA=SIGMA+D(I,ISPIN)                                            
          FPRHO=SIGMA/DR(I)/DR(I)                                           
          IF(FPRHO.LT.1.E-18) THEN                                          
            DEXC(I,1)=0.0                                                     
            DHL(I,1)=0.                                                       
            DEXC(I,2)=0.0                                                     
            DHL(I,2)=0.                                                       
            GOTO 605                                                          
          END IF                                                            
          DEXC(i,1)=0.                                                  
          DEXC(i,2)=0.                                                  
          DHL(I,1)=0.                                                   
          DHL(I,2)=0.                                                   
          DHL1=0.                                                           
          IF(fu.LT.1.E-18) GOTO 602                                     
          IF(fd.LT.1.E-18) GOTO 602                                     
          call grans(dr(i),qsin,qcos,dndrt,dndtt,dndpt, &
                   d2ndrt,d2ndtt,d2ndpt,dndrtt,dndrpt,dndtpt, &
                   gx,gy,gz,gmag,g2,ggx,ggy,ggz,gdggt)
          call grans(dr(i),qsin,qcos,dndru,dndtu,dndpu, &
                   d2ndru,d2ndtu,d2ndpu,dndrtu,dndrpu,dndtpu, &
                   gxu,gyu,gzu,gmagu,g2u,ggxu,ggyu,ggzu,gdggu)
          call grans(dr(i),qsin,qcos,dndrd,dndtd,dndpd, &
                   d2ndrd,d2ndtd,d2ndpd,dndrtd,dndrpd,dndtpd, &
                   gxd,gyd,gzd,gmagd,g2d,ggxd,ggyd,ggzd,gdggd)
!
          CALL VXCLM2(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
                    FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD, &
                    vxcup,vxcdn,dexc(i,1),dexc(i,2),GDGGT,iex)
        dexc(I,1)=dexc(i,1)/2.0d0*D(I,1)                                      
        dexc(I,2)=dexc(i,2)/2.0d0*D(I,2)                                      
!
!        vxcup=VXCSP(FPRHO,fu*4.0d0*pi)                                          
!        vxcdn=VXCSP(FPRHO,fd*4.0d0*pi)                                          
!
        DHL(I,1)=vxcup/2.0d0*D(I,1)                                      
        DP(I)=DP(I)-vxcup/2.0d0                                                  
        DHL(I,2)=vxcdn/2.0d0*D(I,2)                                      
        DP(I)=DP(I)-vxcdn/2.0d0                                                   
  602   CONTINUE                                                  
        DP(I)=DP(I)*SIGMA/2.                                              
! IF(I.LT.30) WRITE(6,*) DHL(I,1),DHL(I,2),DEXC(I,1),DEXC(I,2)     
  605   CONTINUE                                                          
      endif
      DC(2,1)=3                                                         
      DC(3,1)=1                                                         
      IF (NUC.NE.0) DC(3,1)=4                                           
      CALL SOMM (DR,DP,DQ,DPAS,DC(3,1),0,NP)                            
      DC(1,1)=-Z*DC(1,1)                                                
      DC(4,1)=DVAL-DC(3,1)                                              
      DC(8,1)=0.0                                                       
      DO 107 ISPIN=1,JSPIN                                              
      DC(5,ISPIN)=1.                                                    
      DC(6,ISPIN)=1.                                                    
      CALL SOMM (DR,DHL(1,ISPIN),DQ,DPAS,DC(5,ISPIN),0,NP)              
      CALL SOMM (DR,DEXC(1,ISPIN),DQ,DPAS,DC(6,ISPIN),0,NP)             
  107 DC(8,1)=DC(8,1)+DC(6,ISPIN)-DC(5,ISPIN)                           
      DC(8,1)=DC(8,1)+DVAL+0.5*DC(1,1)-0.5*DC(3,1)                      
      DVAL=DVAL*2.                                                      
       DC(8,1)=DC(8,1)*2.                                               
      DC(1,1)=DC(1,1)*2.                                                
      DC(3,1)=DC(3,1)*2.                                                
      DC(6,1)=DC(6,1)*2.                                                
      DC(5,1)=DC(5,1)*2.                                                
      DC(6,2)=DC(6,2)*2.                                                
      DC(5,2)=DC(5,2)*2.                                                
!     IF(IEX.EQ.2) WRITE(6,*)'HEDIN-LUNDQUIST, J.DE.PHYSIQUE'           
!     IF(IEX.EQ.3) WRITE(6,*)'HEDIN-LUNDQUIST, J.PHYSICS'               
      WRITE(6,1007)DC(8,1),DVAL,DC(1,1),DC(3,1),(DC(5,ISPIN),DC(6,ISPIN) &
      ,ISPIN=1,JSPIN)                                                   
      IF(NORB.EQ.1) GO TO 205                                           
      WRITE (6,1001) BAR                                                
      WRITE (6,1008)                                                    
! INTEGRALES DE RECOUVREMENT                                            
! **********************************************************************
      DO 115  I=2,NORB                                                  
      DO 115 ISPIN=1,JSPIN                                              
      K=I-1                                                             
      DO 115  J=1,K                                                     
      IF(NQL(I).NE.NQL(J).OR.NK(I).NE.NK(J)) GO TO 115                  
      IM=NMAX(J)                                                        
      IF(NMAX(I).LT.IM) IM=NMAX(I)                                      
      DO 111  L=1,IM                                                    
      DQ(L)=DPC(L,I,ISPIN)*DPC(L,J,ISPIN)                               
 111  DC(L,1)=DGC(L,I,ISPIN)*DGC(L,J,ISPIN)                             
      DVAL=DFL(I)+DFL(J)                                                
      CALL SOMM (DR,DC,DQ,DPAS,DVAL,0,IM)                               
      WRITE (6,1009) NQN(I),TITRE(I),NQN(J),TITRE(J),DVAL               
 115  CONTINUE                                                          
 205  CONTINUE                                                          
!     CALL CDSLD (TITRE,BAR)                                            
      GO TO 1                                                           
      END                                                               











