      SUBROUTINE INSLD(IEX,JSPIN,RELA,aplot,valel,mult,np0,ecore,&
                        RWatson, EWatson)
! LECTURE DES DONNEES ET POTENTIEL DE DEPART                            
! LA DERNIERE CARTE DONNEE DOIT AVOIR 4 ASTERIX                         
! **********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
!     IMPLICIT REAL (D)                                                 
      LOGICAL RELA 
      character*20 config
      character*1 aplot(30,2)                                     
      COMMON DEN(30,2),DQ1(30),DFL(30),NQN(30),NQL(30),NK(30),NMAX(30),  &
      ZEL(30,2),NORB                                                    
      COMMON/DIRA/DV(NPT,2),DR(NPT),DP(NPT),DQ(NPT),DPAS,Z,NSTOP,NES,    &
      TETS,NP,NUC                                                       
      COMMON/PS2/DEXV,DEXE,DCOP,TITRE(30),BAR(10),TEST,TESTE,TESTY,TESTV &
      ,NITER,ION,ICUT,IPRAT                                             
      DIMENSION TTIRE(9)                                                
      DATA TTIRE/4HS   ,4HP*  ,4HP   ,4HD*  ,4HD   ,4HF*  ,4HF   ,4HG*   &
      ,4HG   /,BAR1/4HRHFS/,ENDA/4H****/,wats/4HWats/,PRAT/4HPRAT/                                
 1000 FORMAT (9A4)                                                      
 1001 FORMAT (12I3)                                                     
 1002 FORMAT (6E12.8)                                                   
 1003 FORMAT (E12.8,I1,2I2)                                             
 1004 FORMAT (8F9.4)                                                    
 2000 FORMAT(////,20X,10A4,//,5X,'NUMBER OF ITERATIONS',I4,//,5X,         &
       'PRECISON OF ENERGIES',                                          &
       1PE9.2,//,18X,'WAVEFUNCTION',1PE9.2,//,18X                        &
      ,'POTENTIAL',1PE9.2,/)                                            
 2001 FORMAT (' NORB=',I3,' STOP')                                      
 2002 FORMAT (' ERROR IN INPUT ',E15.8,I1,I2,F5.2)                      
 2003 FORMAT (' ERROR IN NUMBER OF ELECTRONS')                          
 2004 FORMAT (5X,'ORBITAL ',5X,'OCCUPATION',5X,'TRIAL ENERGIES'/)       
 2005 FORMAT (7X,I1,A2,F16.3,1PE23.7,5x,a1)                                   
 2006 FORMAT (' ERROR IN POTENTIAL')                                    
 2007 FORMAT (' CONFIGURATION FOLKLORIQUE')                             
 2008 FORMAT (' INTEGRATION WITH ',I3,' POINTS STARTING AT',             &
       F7.4,'/',I2,' AND INCREMENT ',F7.4,/)                            
 2009 FORMAT (//,' A RELATIVE ENERGY PRESSION OF',      &
       1PE9.2,/                                                          &
       ' IS USED IN ',I3,' ATTEMPTS',//,                                 &
       ' IEX=',I3,5X,'ICUT=',I3,5X,'IPRAT=',I3,/)                      
 2010 FORMAT(' NUMBER OF R MESH POINTS (',I4,') + NPT00 (',i3,') .GT.',       &
       ' NPT (',I3,')  STOP')                                            
 2011 FORMAT (/,40X,10A4,//,5X,'POTENTIEL DE DEPART MULTIPLIE PAR R'/, &
      10(2X,F9.4))                                                      
 2012 FORMAT ('  DEXV=',1PE14.7,'     DEXE=',1PE14.7,/)                 
 2013 FORMAT (' ERREUR POUR IDEP')                                      
 2014 FORMAT (/,'  THE PRATT SCHEME IS USED'/)                        
 2015 FORMAT (/,'  DCOP=',1PE14.7,/)                                  
 2016 FORMAT (/,30X,'NOYAU DE DIMENSIONS FINIES'/)                    
 2017 FORMAT (' ERROR IN ATOMIC MASS')                                  
 2018 FORMAT ('  CAS SUIVANT')                                          
      IF (NSTOP.EQ.0) GO TO 2                                           
  1   CONTINUE                                                          
      JSPIN=2                                                           
      DSAL=DVC+DVC                                                      
      IZ1=0                                                             
      ION1=0                                                            
      NUC1=-1                                                           
      IF (NSTOP-1) 598,3,598                                            
 2    READ(5,1000) BAR(1)                                               
 3    DPAS=0.05                                                         
      DR1=0.01                                                          
      NES=55                                                            
      NITER=350                                                          

      TESTE=5.E-07
                                                      
      TESTY=1.E-06                                                      
      TESTV=1.E-06                                                      

      TEST=1.E-08                                                       
                                                      
      NP=NPT                                                            
      NSTOP=30                                                          
      DEXV=0.666666667                                                  
      DEXE=1.                                                           
      DCOP=0.3                                                          

      READ(5,1000) (BAR(I),I=1,9)                                       
        
!                                                                       
! BAR TITRE DU CAS TRAITE   36 CARACTERES                               
! **********************************************************************
      IF (BAR(1).EQ.ENDA) THEN
         write(15,8100) ecore-3.d0
 8100    format('K-VECTORS FROM UNIT:4  ',f5.1,'       2.0      emin/emax', &
                ' window')                   
         write(16,*) '0'
         write(17,8101) 
 8101    format('TOT             (TOT,FOR,QTL,EFG,FERMI)')
         write(17,8181) ecore-3.d0,valel,0.5,0.05
 8181    format(2f10.1,2f5.2,16x,'EMIN, NE, ESEPERMIN, ESEPER0',/, &
            'TETRA    0.000',10x, &
            '(GAUSS,ROOT,TEMP,TETRA,ALL      eval)')
         STOP 'LSTART ENDS'                           
      END IF

      IPRAT=0                                                           
      RWatson=10000.0d0
      EWatson=0.0d0
!      IF (BAR(2).eq.4HWats) read(5,*) RWatson, EWatson       
!      IF (BAR(2).eq.4HPRAT) IPRAT=1
      IF (BAR(2).eq.Wats) read(5,*) RWatson, EWatson       
      IF (BAR(2).eq.PRAT) IPRAT=1

      BAR(10)=BAR1    
      norb=10
!      iex=5  
      DVC=137.0359895                                                      
      IF(.NOT.RELA) DVC=1.E30                                           
!                                                                       
!     ALL OTHER PARAMETERS ARE SET BY DEFAULT                           
      ION=0                                                             
      IDEP=0                                                            
      ICUT=1                                                            
      IF (IEX.LT.10) IEX=5                                                     
      J=50                                                              
      K=0                                                               
      L=1                                                               
      NUC=0 
!                                                            
!     READ STRUCT FILE                                                  
! ******************************************************                      
      READ(20,1012) IATNR,POSX,POSY,POSZ,MULT                           
 1012 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,3X,/,15X,I2)              
 1011 FORMAT(4X,I4)                                                     
 1050 FORMAT(15X,I5,5X,F10.9,5X,F10.5,5X,F5.1)                         
      DO 110 MU=1,MULT-1                                                
 110  READ(20,1011) IATNR                                               
      READ(20,1050) NP0,DR1,RMTT,ZZ                                       
      IZ=ZZ+0.0000001                                                   
      READ(20,*)                                                        
      READ(20,*)                                                        
      READ(20,*)                                                        
      DPAS=LOG(RMTT/DR1)/(np0-1)                                          
!                                                                       
! IZ NUMERO ATOMIQUE    ION=IZ-NOMBRE D ELECTRONS                       
! NORB NOMBRE D ORBITALES IDEP DOIT ETRE EGAL A 0 OU 1                  
! IDEP=0 POTENTIEL DE DEPART=POTENTIEL THOMAS FERMI                     
! IDEP=1 POTENTIEL DE DEPART EST LU SUR CARTES                          
! SI ICUT EST NUL ON CORRIGE LE POTENTIEL EN -(ION+1)/R                 
! SI IPRAT EST NUL LE PROCEDE DE PRATT EST UTILISE                      
! SI IEX EST NUL ON UTILISE LE POTENTIEL D ECHANGE DE SLATER NON MODIFIE
!  IEX=0   KSG-EXCHANGE                                                 
!  IEX=1   X-ALPHA                                                      
!  IEX=2   HL (JOURNAL DE PHYSIQUE)                                     
!  IEX=3   HL (J.PHYS, OLD, FREEMAN-VERSION)                            
! I NCMBRE DE POINTS POUR L INTEGRATION  NPT PAR DEFAUT                 
! J NCMBRE D ESSAIS POUR AJUSTER L ENERGIE  15 PAR DEFAUT               
! K NOMBRE D ITERATIONS  50 PAR DEFAUT                                  
! L=0 OPTION STANDARD POUR LE BLOC DES POINTS ET LES PRECISIONS         
! NOYAU DE DIMENSIONS FINIES SI NUC POSITIF                             
! **********************************************************************
      IF(NORB.LE.NSTOP) GO TO 4                                         
      WRITE (6,2001) NORB                                               
      GO TO 598                                                         
 4    IF (np0.LE.0) GO TO 6                                               
      IF (np0+npt00.LE.NP) GO TO 5                                           
      WRITE (6,2010) np0,npt00,NP                                              
      GO TO 598                                                         
 5    NP=np0+npt00                                                           
      np=2*(np/2)+1                                                       
 6    IF (J.GT.0) NES=J                                                 
      IF (K.GT.0) NITER=K                                               
      IF (IEX.NE.1) GO TO 7                                             
      READ(5,1002) DEXV,DEXE                                            
!                                                                       
! DEXV COEFFICIENT POUR LE POTENTIEL D ECHANGE  DEXV=1. POUR SLATER     
! DEXE COEFFICIENT POUR L ENERGIE D ECHANGE                             
! DEXV DOIT ETRE EGAL A 2.*DEXE/3. POUR SATISFAIRE LE VIRIEL            
! **********************************************************************
 7    IF (L.EQ.0) GO TO 8                                               
!     READ(5,*) DR1,DPAS                                                
      DR1=DR1*ZZ                                                        
!                                                                       
! DPAS PAS EXPONENTIEL   DR1 DEFINIT LE PREMIER POINT =DR1/IZ           
! TEST PRECISION SUR L ENERGIE DANS RESLS                               
! TESTE CRITERE DE SELF CONSISTENCE POUR LES ENERGIES MONDELECTRONIQUES 
! TESTY CRITERE DE SELF CONSISTENCE POUR LES FONCTIONS D ONDE           
! TESTV CRITERE DE SELF CONSISTENCE POUR LE POTENTIEL                   
! **********************************************************************
 8    IF (IPRAT.EQ.0) GO TO 9                                           
      READ(5,*) DCOP                                                 
!                                                                       
! VI(N+1)=(1.-DCOP)*VI(N)+DCOP*VF(N)                                    
! **********************************************************************
 9    Z=zZ                                                              
!!$ 9    Z=IZ                                                              
      IF (NUC.LE.0) GO TO 16                                            
      READ(5,1002) DVAL                                                 
!                                                                       
! DVAL MASSE ATOMIQUE SI NUC POSITIF                                    
! **********************************************************************
      DVAL=Z*(DVAL**(1./3.))*2.2677D-05/EXP(4.*DPAS)                    
      IF (DVAL-DR1) 11,11,12                                            
 11   DR1=DVAL                                                          
      NUC=5                                                             
      GO TO 16                                                          
 12   DVAL=DVAL*EXP(4.*DPAS)                                            
      DO 13 I=6,NP                                                      
      D1=DR1*EXP((I-1)*DPAS)                                            
      IF (D1.GE.DVAL) GO TO 14                                          
 13   CONTINUE                                                          
      WRITE (6,2017)                                                    
      GO TO 598                                                         
 14   NUC=I                                                             
      DR1=DR1*DVAL/D1                                                   
 16   WRITE (6,2000) BAR,NITER,TESTE,TESTY,TESTV                        
      WRITE (6,2008) NP,DR1,IZ,DPAS                                     
!      WRITE (6,2012) DEXV,DEXE                                          
      ZK=0                                                              
      DVAL=Z*Z/(DVC*DVC)                                                
      IF (NUC.LE.0) GO TO 17                                            
      WRITE (6,2016)                                                    
!                                                                       
! DONNEES POUR LES ORBITALES                                            
! **********************************************************************
 17   WRITE (6,2004)                                                    
      READ(5,2318) config
      call parse(config,norb0,iex,zk,z,nuc,dval,TITRE, &
      TTIRE,aplot)
       iexout=iex
!      if(iex.eq.14) then
!          IEXout=14 
!      else if(iex.eq.13) then
!          IEXout=13
!      else if(iex.eq.25) then
!          IEXout=25 
!      else if(iex.eq.17) then
!          IEXout=17
!      else if(iex.eq.3) then
!          IEXout=1
!      else if(iex.eq.1) then
!          iex=3
!          iexout=1
!      else
!          iex=5
!          iexout=5
!      end if 
 2318 FORMAT (a20)                                          
      DO 24 I=norb0+1,NORB                                                    
      DO 24 ISPIN=1,JSPIN                                               
      READ(5,1070) NQN(I),NK(I),ZEL(I,ISPIN),aplot(i,ispin)     
! DEN ENERGIE DE L ORBITALE EN UNITE ATOMIQUE ET NEGATIVE               
! NQN NOMBRE QUANTIQUE PRINCIPAL   NK NOMBRE QUANTIQUE KAPPA            
! ZEL OCCUPATION DE L ORBITALE                                          
 1070 FORMAT(I1,1X,i2,1x,f5.3,a1)                                   
! **********************************************************************
      ZK=ZK+ZEL(I,ISPIN)                                                
!     IF (DEN(I,ISPIN)) 19,18,18                                        
 18   DEN(I,ISPIN)=-Z*Z/(4.*NQN(I)*NQN(I))                              
 19   NQL(I)=IABS(NK(I))                                                
      IF (NK(I).LT.0) NQL(I)=NQL(I)-1                                   
      IF (NUC.GT.0) GO TO 21                                            
      DFL(I)=NK(I)*NK(I)                                                
      DFL(I)=SQRT(DFL(I)-DVAL)                                          
      GO TO 22                                                          
 21   DFL(I)=IABS(NK(I))                                                
 22   L=2*IABS(NK(I))                                                   
      NEL=ZEL(I,ISPIN)                                                  
      IF (NQL(I).LE.NQN(I).AND.NEL.LE.L.AND.NQN(I).GT.0.AND.NQL(I).LE.4) GO TO 23                                                      
      WRITE (6,2002) DEN(I,ISPIN),NQN(I),NQL(I),J,ZEL(I,ISPIN)          
      GO TO 598                                                         
 23   J=NQL(I)+IABS(NK(I))                                              
      TITRE(I)=TTIRE(J)                                                 
 24   WRITE (6,2005) NQN(I),TITRE(I),ZEL(I,ISPIN),DEN(I,ISPIN), &
       aplot(i,ispin)          
      K=ZK                                                              
      IF (K.EQ.(IZ-ION)) GO TO 25                                       
      WRITE(6,*) 'NUMBER OF ELECTRONS NE IZ ',ZK,IZ                     
      WRITE(*,*) 'NUMBER OF ELECTRONS NE IZ ',ZK,IZ                     
!     GO TO 598                  
 25   continue
      rewind 14
      write(14,8104) iexout
 8104 format('TOT',I5,'     (5...CA-LDA, 13...PBE-GGA, 11...WC-GGA)',/, &
              'NR2V    ')
!                      ,'IFFT      (R2V)',/,' 2.0d0      ifftfactor')
      WRITE (6,2009) TEST,NES,Iex,ICUT,IPRAT                           
      IF (IPRAT) 27,26,27                                               
 26   WRITE (6,2014)                                                    
      GO TO 28                                                          
 27   WRITE (6,2015) DCOP                                               
 28   IF (NUC.NE.NUC1) GO TO 29                                         
      IF (IZ.EQ.IZ1.AND.ION.EQ.ION1) GO TO 101                          
      IF(IZ.EQ.IZ1) GO TO 35                                            
 29   DR(1)=DR1/Z                                                       
      DO 31  I=2,NP                                                     
      DR(I)=DR(1)*EXP((I-1)*DPAS)                                       
 31   if(dr(i).gt.RMAX0) goto 32
      write(6,*) i,np 
      write(6,*) ' WARNING !!!! For good atomic total energies you ', &
      'should change the radial mesh (reduce NRAD or increase' &
      ,' R0), or increase PARAMETERS NPT and NPT00 (ignore this warning,'&
      ,' if you only need starting densities for scf)'
 32   np=i-1
      np=2*(np/2)+1                                                       
      write(6,*) np,' RADIAL MESH POINTS REACH UP TO',dr(np),' A.U.'
!      if(np.eq.npt) 
      write(6,*)
!                                                                       
! POTENTIEL DE DEPART                                                   
! **********************************************************************
 35   VAL=-ION-1                                                        

      IF (IDEP.EQ.1) GO TO 58                                           
      IF (IDEP.EQ.0) GO TO 45                                           
      WRITE (6,2013)                                                    
      GO TO 598                                                         
 45   IF(IZ.EQ.IZ1.AND.ION.GT.ION1.AND.NUC.EQ.NUC1) GO TO 60            
      DO 47  I=1,NP                                                      
      R=DR(I)                                                           
      DV(I,1)=FPOT(R,Z,VAL)
  47  DV(I,2)=DV(I,1)                                                   
      IF (NUC.LE.0) GO TO 101                                           
      DO 48 I=1,NUC                                                     
 48   DV(I,1)=DV(I,1)+Z/DR(I)+Z*((DR(I)/DR(NUC))**2-3.)/(DR(NUC)+        &
      DR(NUC))                                                          
      GO TO 101                                                         
 58   READ(5,1004) (DV(I,1),I=1,NP)                                     
      READ(5,1004) (DV(I,2),I=1,NP)                                     
!                                                                       
! LECTURE DU POTENTIEL DE DEPART (EN U.A. ET NEGATIF)  SI IDEP=1        
! **********************************************************************
      WRITE (6,2011) BAR,(DV(I,1),I=1,NP)                               
      WRITE (6,2011) BAR,(DV(I,2),I=1,NP)                               
      DVAL=-Z/DV(1,1)                                                   
      IF (NUC.GT.0) DVAL=1.                                             
      DO 59  I=1,NP                                                     
      DV(I,1)=DV(I,1)*DVAL/DR(I)                                        
 59   DV(I,2)=DV(I,2)*DVAL/DR(I)                                        
 60   IF (ICUT.NE.0) GO TO 67                                           
      DO 65  I=1,NP                                                     
      IF ((DR(I)*DV(I,1)).GT.VAL) DV(I,1)=VAL/DR(I)                     
 65   CONTINUE                                                          
 67   VAL=Z+DV(1,1)*DR(1)                                               
      IF (NUC.GT.0) VAL=Z+DV(NUC,1)*DR(NUC)                             
      IF (ABS(VAL).LT.0.1) GO TO 101                                    
      WRITE (6,2006)                                                    
      GO TO 598                                                         
 101  IF(NORB.EQ.1) GO TO 107                                           
      DO 105  I=2,NORB                                                  
      K=I-1                                                             
      DO 105  J=1,K                                                     
      IF (NQN(I).NE.NQN(J).OR.NK(I).NE.NK(J)) GO TO 105                 
      WRITE (6,2007)                                                    
      GO TO 598                                                         
 105  CONTINUE                                                          
 107  IZ1=IZ                                                            
      ION1=ION                                                          
      NUC1=NUC                                                          
      DO 119 I=1,NORB                                                   
      NMAX(I)=NP                                                        
      L=1                                                               
      J=NQN(I)-NQL(I)                                                   
      IF((J-2*(J/2)).EQ.0) L=-L                                         
      DQ1(I)=L*NK(I)/IABS(NK(I))                                        
      IF (NUC) 111,119,111                                              
 111  IF (NK(I)) 112,119,119                                            
 112  DQ1(I)=DQ1(I)*(NK(I)-DFL(I))*DVC/Z                                
 119  CONTINUE                                                          
      RETURN                                                            
 598  WRITE (6,2018)                                                    
 599  READ(5,1000) BAR(1)                                               
      IF (BAR(1)-ENDA) 599,601,599                                      
 601  NSTOP=1                                                           
      GO TO 1                                                           
      END                                                               












