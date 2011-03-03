      SUBROUTINE INSLD(RELA,npxx,cform,shift)                                            
!                                                                       
! LECTURE DES DONNEES ET POTENTIEL DE DEPART                            
! LA DERNIERE CARTE DONNEE DOIT AVOIR 4 ASTERIX                         
! **********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      CHARACTER *10 BAR                                                 
      CHARACTER *4  BARONE,ENDA,cform                                         
      LOGICAL RELA                                                      
      real*8 nel
      COMMON DEN(30),DQ1(30),DFL(30),NQN(30),NQL(30),NK(30),NMAX(30),NEL &
      (30),NORB                                                         
      COMMON/DIRA/DV(NRAD+41),DR(NRAD+41),DP(NRAD+41),DQ(NRAD+41),DPAS, &
      Z,NSTOP,NES,TETS,NP,NUC                                                   
      COMMON/PS2/DEXV,DEXE,DCOP,TITRE(30),TEST,TESTE,TESTY,TESTV         &
      ,NITER,ION,ICUT,IPRAT                                             
      COMMON /CHAR/ BAR                                                 
      DIMENSION TTIRE(9)                                                
      DATA TTIRE/4HS   ,4HP*  ,4HP   ,4HD*  ,4HD   ,4HF*  ,4HF   ,4HG*   &
      ,4HG   /                                                          
!    2,BAR1/4HRHFS/,ENDA/4H****/                                        
      DATA ENDA/'****'/
      integer new
      data new /1/
      save new                                                 
 1000 FORMAT (9A4)                                                      
 1001 FORMAT (12I3)                                                     
 1002 FORMAT (6E12.8)                                                   
 1003 FORMAT (E12.8,I1,2I2)                                             
 1004 FORMAT (5E14.7)                                                   
 1990 FORMAT(1H1,40X,A10,//,5X,'NOMBRE D ITERATIONS',I4,//,5X, &
      'PRECISIO N SUR LES ENERGIES',1PE9.2,//,23X,'FONCTIONS D ONDE',1PE9.2,//,23X &
      ,'POTENTIEL',1PE9.2,/)                                            
 2001 FORMAT (' NORB=',I3,'TROP GRAND')                                 
 2002 FORMAT (' ERREUR DANS LA CARTE ',E15.8,I1,2I2,f5.2)                    
 2003 FORMAT (' ERREUR POUR LE NOMBRE D ELECTRONS')                     
 2004 FORMAT (5X,'ORBITALE',5X,'OCCUPATION',5X,'ENERGIE DE DEPART'/)    
 2005 FORMAT (7X,I1,A2,f16.2,1PE23.7)                                     
 2006 FORMAT (' ERREUR POUR LE POTENTIEL')                              
 2007 FORMAT (' CONFIGURATION FOLKLORIQUE')                             
 2008 FORMAT (' L INTEGRATION EST FAITE EN ',I3,' POINTS LE PREMIER EST EGAL A ',F7.4,'/',I2,' ET LE PAS A ',F7.4,/)                      
 2009 FORMAT (' DANS LE SOUS PROGRAMME RESLD LA PRECISION RELATIVE A OBTENIR SUR L ENERGIE EST ' &
      ,1PE9.2,' ET LE NOMBRE D ESSAIS ',I3,//,  &
      'IDEP=',I3,5X,'ICUT=',I3,5X,'IPRAT=',I3,/)                         
 2010 FORMAT ('  NP=',I3,'  TROP GRAND')                                
 2011 FORMAT (1H1,40X,A10,//,5X,'POTENTIEL DE DEPART '/,                 &
      10(2X,E9.3))                                                      
 2012 FORMAT ('  DEXV=',1PE14.7,'     DEXE=',1PE14.7,/)                 
 2013 FORMAT (' ERREUR POUR IDEP')                                      
 2014 FORMAT (1H0,'  LE PROCEDE DE PRATT EST UTILISE'/)                 
 2015 FORMAT (1H0,'  DCOP=',1PE14.7,/)                                  
 2016 FORMAT (1H0,30X,'NOYAU DE DIMENSIONS FINIES'/)                    
 2017 FORMAT (' ERREUR POUR LA MASSE ATOMIQUE')                         
 2018 FORMAT ('  CAS SUIVANT')                                          
 2020 FORMAT(3X,5E13.7)                                                 
 2021 FORMAT(3X,4E19.12)                                                 
!      if(cform.eq.'NEW ') then
       assign 2021 to iform1  
!      else
!       assign 2020 to iform1  
!      end if
      IF (NSTOP.EQ.0) GO TO 2                                           
 1    DVC=137.0359895d0                                                     
      IF(.NOT.RELA) DVC=1.d30                                           
      DSAL=DVC+DVC                                                      
      IZ1=0                                                             
      ION1=0                                                            
      NUC1=-1                                                           
      IF (NSTOP-1) 598,3,598                                            
 2    CONTINUE                                                          
 3    DPAS=0.05d0                                                         
      DR1=0.01d0                                                          
      NES=15                                                            
      NITER=50                                                          
      TESTE=5.d-06                                                      
      TESTY=1.d-05                                                      
      TESTV=1.d-05                                                      
      TEST=1.d-07                                                       
      NP=381                                                            
      NSTOP=30                                                          
      DEXV=2.d0/3.d0                                                  
      DEXE=1.d0                                                           
      DCOP=0.3d0                                                          
!     READ (5,1000) (BAR(I),I=1,9)                                      
!                                                                       
! BAR TITRE DU CAS TRAITE   36 CARACTERES                               
! **********************************************************************
!     IF (BAR(1).EQ.ENDA) CALL EXIT                                     
!     BAR(10)=BAR1                                                      
!     READ (20,1012) BAR,NP,DR1,DPAS,ZZ                                 
!1012 FORMAT(/A10,I5,5X,2F10.9,5X,F5.2)                                 
!                                                                       
!.....READ FILE POTE                                                    
      READ(20,1012) IATNR,POSX,POSY,POSZ,MULT                           
         DO 110 MU=1,MULT-1                                             
            READ(20,1011) IATNR                                         
 110     CONTINUE                                                       
      READ(20,1050) BAR,NP,DR1,RMTT,ZZ                                  
      DPAS=LOG(RMTT/DR1) / (NP-1)                                       
      READ(20,*)                                                        
      READ(20,*)                                                        
      READ(20,*)                                                        
 1011 FORMAT(4X,I4)                                                     
 1012 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,7X,I2,8X,I2)     
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
!                                                                       
      IZ=ZZ                                                             
      IDEP=1
      shift=0.d0  
      if(new.eq.0) goto 87                                                   
      READ(5,*,err=88) NORB,shift
      goto 89
 88   continue
      new=0
      rewind (5)
 87   READ(5,*) NORB     
 89   continue
      write(6,*) 'NORB, SHIFT ',norb,shift
      ICUT=1                                                            
      IPRAT=0                                                           
      IEX=0                                                             
      I=NP                                                              
      J=30                                                              
      L=0                                                               
      NUC=0                                                             
      K=1                                                               
!                                                                       
! IZ NUMERO ATOMIQUE    ION=IZ-NOMBRE D ELECTRONS                       
! NORB NOMBRE D ORBITALES IDEP DOIT ETRE EGAL A 0 OU 1                  
! IDEP=0 POTENTIEL DE DEPART=POTENTIEL THOMAS FERMI                     
! IDEP=1 POTENTIEL DE DEPART EST LU SUR CARTES                          
! SI ICUT EST NUL ON CORRIGE LE POTENTIEL EN -(ION+1)/R                 
! SI IPRAT EST NUL LE PROCEDE DE PRATT EST UTILISE                      
! SI IEX EST NUL ON UTILISE LE POTENTIEL D ECHANGE DE SLATER NON MODIFIE
! I NCMBRE DE POINTS POUR L INTEGRATION  NRAD PAR DEFAUT                 
! J NCMBRE D ESSAIS POUR AJUSTER L ENERGIE  15 PAR DEFAUT               
! K NOMBRE D ITERATIONS  50 PAR DEFAUT                                  
! L=0 OPTION STANDARD POUR LE BLOC DES POINTS ET LES PRECISIONS         
! NOYAU DE DIMENSIONS FINIES SI NUC POSITIF                             
! **********************************************************************
      IF(NORB.LE.NSTOP) GO TO 4                                         
      WRITE (6,2001) NORB                                               
      GO TO 598                                                         
 4    IF (I.LE.0) GO TO 6                                               
!      IF (I.LE.NP) GO TO 5                                              
!      WRITE (6,2010) I                                                  
!      GO TO 598                                                         
 5    NP=I                                                              
 6    IF (J.GT.0) NES=J                                                 
      IF (K.GT.0) NITER=K                                               
      IF (IEX.EQ.0) GO TO 7                                             
      READ(5,1002) DEXV,DEXE                                           
!                                                                       
! DEXV COEFFICIENT POUR LE POTENTIEL D ECHANGE  DEXV=1. POUR SLATER     
! DEXE COEFFICIENT POUR L ENERGIE D ECHANGE                             
! DEXV DOIT ETRE EGAL A 2.*DEXE/3. POUR SATISFAIRE LE VIRIEL            
! **********************************************************************
 7    IF (L.EQ.0) GO TO 8                                               
      READ(5,1002) DPAS,DR1                                            
      DR1=DR1*IZ                                                        
!                                                                       
! DPAS PAS EXPONENTIEL   DR1 DEFINIT LE PREMIER POINT =DR1/IZ           
! TEST PRECISION SUR L ENERGIE DANS RESLS                               
! TESTE CRITERE DE SELF CONSISTENCE POUR LES ENERGIES MONDELECTRONIQUES 
! TESTY CRITERE DE SELF CONSISTENCE POUR LES FONCTIONS D ONDE           
! TESTV CRITERE DE SELF CONSISTENCE POUR LE POTENTIEL                   
! **********************************************************************
 8    IF (IPRAT.EQ.0) GO TO 9                                           
      READ(5,1002) DCOP                                                
!                                                                       
! VI(N+1)=(1.-DCOP)*VI(N)+DCOP*VF(N)                                    
! **********************************************************************
 9    Z=ZZ
!     z=iz                                                              
      IF (NUC.LE.0) GO TO 16                                            
      READ(5,1002) DVAL                                                
!                                                                       
! DVAL MASSE ATOMIQUE SI NUC POSITIF                                    
! **********************************************************************
      DVAL=Z*(DVAL**(1.d0/3.d0))*2.2677D-05/EXP(4.d0*DPAS)                    
      IF (DVAL-DR1) 11,11,12                                            
 11   DR1=DVAL                                                          
      NUC=5                                                             
      GO TO 16                                                          
 12   DVAL=DVAL*EXP(4.d0*DPAS)                                            
      DO 13 I=6,NP                                                      
      D1=DR1*EXP((I-1)*DPAS)                                            
      IF (D1.GE.DVAL) GO TO 14                                          
 13   CONTINUE                                                          
      WRITE (6,2017)                                                    
      GO TO 598                                                         
 14   NUC=I                                                             
      DR1=DR1*DVAL/D1                                                   
 16   WRITE (6,1990) BAR,NITER,TESTE,TESTY,TESTV                        
      WRITE (6,2008) NP,DR1,IZ,DPAS                                     
      WRITE (6,2009) TEST,NES,IDEP,ICUT,IPRAT                           
!     WRITE (6,2012) DEXV,DEXE                                          
      K=0                                                               
      DVAL=Z*Z/(DVC*DVC)                                                
      IF (NUC.LE.0) GO TO 17                                            
      WRITE (6,2016)                                                    
!                                                                       
! DONNEES POUR LES ORBITALES                                            
! **********************************************************************
 17   WRITE (6,2004)                                                    
      IF(NORB.EQ.0) GOTO 58                                             
      DO 24 I=1,NORB                                                    
      READ(5,*) NQN(I),NK(I),NEL(I)                                    
!                                                                       
! DEN ENERGIE DE L ORBITALE EN UNITE ATOMIQUE ET NEGATIVE               
! NQN NOMBRE QUANTIQUE PRINCIPAL   NK NOMBRE QUANTIQUE KAPPA            
! NEL OCCUPATION DE L ORBITALE                                          
! **********************************************************************
      K=K+NEL(I)                                                        
!     IF (DEN(I)) 19,18,18                                              
 18   DEN(I)=-Z*Z/(4.d0*NQN(I)*NQN(I))                                    
 19   NQL(I)=IABS(NK(I))                                                
      IF (NK(I).LT.0) NQL(I)=NQL(I)-1                                   
      IF (NUC.GT.0) GO TO 21                                            
      DFL(I)=NK(I)*NK(I)                                                
      DFL(I)=SQRT(DFL(I)-DVAL)                                          
      GO TO 22                                                          
 21   DFL(I)=IABS(NK(I))                                                
 22   L=2*IABS(NK(I))                                                   
      IF (NQL(I).LE.NQN(I).AND.NEL(I).LE.L.AND.NQN(I).GT.0.AND.NQL(I).LE.4) GO TO 23                                                      
      WRITE (6,2002) DEN(I),NQN(I),NQL(I),J,NEL(I)                      
      GO TO 598                                                         
 23   J=NQL(I)+IABS(NK(I))                                              
      TITRE(I)=TTIRE(J)                                                 
 24   WRITE (6,2005) NQN(I),TITRE(I),NEL(I),DEN(I)                      
      ION=IZ-K                                                          
      IF (K.EQ.(IZ-ION)) GO TO 25                                       
      WRITE(6,2003)                                                     
      GO TO 598                                                         
 25   IF (IPRAT) 27,26,27                                               
 26   WRITE (6,2014)                                                    
      GO TO 28                                                          
 27   WRITE (6,2015) DCOP                                               
 28   IF (NUC.NE.NUC1) GO TO 29                                         
      IF (IZ.EQ.IZ1.AND.ION.EQ.ION1) GO TO 101                          
      IF(IZ.EQ.IZ1) GO TO 35                                            
 29   DR(1)=DR1                                                         
      DO 31  I=2,NRAD+41                                                    
 31   DR(I)=DR(1)*EXP((I-1)*DPAS)                                       
!                                                                       
! POTENTIEL DE DEPART                                                   
! **********************************************************************
 35   VAL=-ION-1                                                        
      IF (IDEP.EQ.1) GO TO 58                                           
      IF (IDEP.EQ.0) GO TO 45                                           
      WRITE (6,2013)                                                    
      GO TO 598                                                         
 45   IF(IZ.EQ.IZ1.AND.ION.GT.ION1.AND.NUC.EQ.NUC1) GO TO 101           
      DO 47  I=1,NP                                                      
      R=DR(I)                                                           
 47   DV(I)=FPOT(R,Z,VAL)                                               
      IF (NUC.LE.0) GO TO 101                                           
      DO 48 I=1,NUC                                                     
 48   DV(I)=DV(I)+Z/DR(I)+Z*((DR(I)/DR(NUC))**2-3.d0)/(DR(NUC)+DR(NUC))   
      GO TO 101                                                         
 58   READ(8,1980)                                                      
      READ(8,2000) IDUMMY                                               
      READ(8,2031)                                                      
      READ(8,iform1) ( DV(I), I=1,NP )                                    
      READ(8,2033)                                                      
      READ(8,2031)                                                      
      NPXX=40                                                           
      I=2*(np/2)+1  
      if(i.gt.np) npxx=41                                                     
      DVR=DV(NP)
      DO 10 I=1,NP+NPXX                                                 
      IF(I.GT.NP)DV(I)=DVR                                           
  10  DV(I)= DV(I)/2.d0/DR(I)                                             
      NP=NP+NPXX                                                        
!cccccccc
!  constant shift of the potential by 1 Ry, the resulting eigenvalues
!  should also be shifted back.
      do 911 i=1,np
 911      dv(i)=dv(i)-shift/2.d0
!cccccccc
!                                                                       
! LECTURE DU POTENTIEL DE DEPART (EN U.A. ET NEGATIF)  SI IDEP=1        
! **********************************************************************
      WRITE (6,2011) BAR,(DV(I),I=1,NP)                                 
      IF(NORB.EQ.0) RETURN                                              
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
 599  READ(5,1000) BARONE                                              
      IF (BARONE.NE.ENDA) GO TO 599                                     
 601  NSTOP=1                                                           
      GO TO 1                                                           
 1980 FORMAT(3X)                                                        
 2000 FORMAT(16X,I2,//)                                                 
 2031 FORMAT(/)                                                         
 2033 FORMAT(///)                                                       
      END                                                               
