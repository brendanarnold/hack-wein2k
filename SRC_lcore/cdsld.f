      SUBROUTINE CDSLD (TITRE,BAR)                                      
!                                                                       
! TITRE IDENTIFICATION DES FONCTIONS D ONDE   S,P*,P,........           
! BAR TITRE DU CAS TRAITE  10 MEMOIRES  TITRE(10)= RHFS                 
! **********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      real*8 nel
      COMMON DEN(30),DQ1(30),DFL(30),NQN(30),NQL(30),NK(30),NMAX(30),NEL &
      (30),NORB                                                         
      COMMON/DIRA/DV(NRAD+41),DR(NRAD+41),DP(NRAD+41),DQ(NRAD+41), &
                  DPAS,Z,NSTOP,NES,TETS,NP,NUC 
      COMMON/DEUX/DVN(NRAD+41),DVF(NRAD+41),D(NRAD+41),DC(NRAD+41), &
      DGC(NRAD+41,30),DPC(NRAD+41,30)
      DIMENSION TITRE(*),BAR(10)                                        
 1000 FORMAT (8I3)                                                      
 2000 FORMAT (1H1,40X,10A1,/)                                           
 2001 FORMAT (1H0,'    ERREUR POUR LA CARTE   ',8I3)                    
 2002 FORMAT(24X,'(P1P2+Q1Q2)R**',I2,' POUR ',I1,A2,I3,A2,5X,'=',1PE14.7 &
      ,/)                                                               
 2003 FORMAT(24X,'(P1Q2+Q1P2)R**',I2,' POUR ',I1,A2,I3,A2,5X,'=',1PE14.7 &
      ,/)                                                               
 2004 FORMAT (9X,'R',14X,3(I1,A2,'G.C.',I11,A2,'P.C.',10X))             
 2005 FORMAT (7(1PE17.7))                                               
 2006 FORMAT (8F9.4)                                                    
 2007 FORMAT (10A4)                                                     
 2008 FORMAT (/,30X,'INTEGRALES MAGNETIQUES DIRECTES ET D ECHANGE'//)   
 2009 FORMAT (20X,'FM',I2,' (',I1,A2,',',I1,A2,') =',1PE14.7)           
 2010 FORMAT (' GM',I2,' (',I1,A2,',',I1,A2,')',5X,'(-1)=',1PE14.7,5X,'(0)=', &
      1PE14.7,5X,'(+1)=',1PE14.7)                                  
 2011 FORMAT (/,20X,'INTEGRALES MAGNETIQUES RK=INTEGRALE DE P1(1)*Q2(1)* UK(1,2)*P3(2)*Q4(2)'//)                                           
 2012 FORMAT (20X,'RM',I2,' (',I1,A2,',',I1,A2,',',I1,A2,',',I1,A2,') =' &
      ,1PE14.7)                                                         
      READ(5,1000) IRM,INS,NPUN,NMFG,NMRK                              
!                                                                       
! VALEURS MOYENNES DE R**J  SI IRM NON NUL                              
! TABULATION DES FONCTIONS D ONDE SI INS NON NUL                        
! LE POTENTIEL MULTIPLIE PAR R EST PERFORE SI NPUN NON NUL              
! **********************************************************************
      IF (IRM.EQ.0) GO TO 21                                            
      WRITE (6,2000) BAR                                                
 1    READ(5,1000) J,L,N1,L1,J1,N2,L2,J2                               
      IF (L.EQ.0) GO TO 21                                              
!                                                                       
! VALEUR MOYENNE DE (P1*P2+Q1*Q2)*R**J  SI L POSITIF                    
! VALEUR MOYENNE DE (P1*Q2+P2*Q1)*R**J  SI L NEGATIF                    
! **********************************************************************
      IF (N1.GT.0) GO TO 2                                              
      IF (((N1+1)*(N1+2)).NE.0) GO TO 6                                 
      I1=1                                                              
      I2=1                                                              
      GO TO 7                                                           
 2    I1=0                                                              
      I2=0                                                              
      DO 5 I=1,NORB                                                     
      IF (NQN(I).EQ.N1.AND.NQL(I).EQ.L1.AND.(J1-1).EQ.(-NK(I)/IABS(NK(I) &
      ))) I1=I                                                          
      IF (NQN(I).EQ.N2.AND.NQL(I).EQ.L2.AND.(J2-1).EQ.(-NK(I)/IABS(NK(I) &
      ))) I2=I                                                          
 5    CONTINUE                                                          
      IF (I1.NE.0.AND.I2.NE.0) GO TO 7                                  
 6    WRITE (6,2001) J,L,N1,L1,J1,N2,L2,J2                              
      GO TO 1                                                           
 7    DVAL=DFL(I1)+DFL(I2)                                              
      IF ((DVAL+J).GT.-1.d0) GO TO 8                                      
      IF (N1) 16,16,6                                                   
 8    IM=NMAX(I1)                                                       
      IF (NMAX(I2).LT.IM) IM=NMAX(I2)                                   
      IF (L.LT.0) GO TO 11                                              
      DO 9 I=1,IM                                                       
      DV(I)=DGC(I,I1)*DGC(I,I2)                                         
 9    DQ(I)=DPC(I,I1)*DPC(I,I2)                                         
      GO TO 14                                                          
 11   DO 13 I=1,IM                                                      
      DV(I)=DGC(I,I1)*DPC(I,I2)                                         
 13   DQ(I)=DGC(I,I2)*DPC(I,I1)                                         
 14   CALL SOMM (DR,DV,DQ,DPAS,DVAL,J,IM)                               
      IF (L.LT.0) GO TO 15                                              
      WRITE (6,2002) J,NQN(I1),TITRE(I1),NQN(I2),TITRE(I2),DVAL         
      GO TO 16                                                          
 15   WRITE (6,2003) J,NQN(I1),TITRE(I1),NQN(I2),TITRE(I2),DVAL         
 16   IF (N1+1) 19,17,1                                                 
 17   I1=I1+1                                                           
      I2=I1                                                             
      IF (I1-NORB) 7,7,1                                                
 19   I2=I2+1                                                           
      IF (I2-NORB) 7,7,17                                               
 21   IF (INS.EQ.0) GO TO 31                                            
      DO 29 I=1,NORB,3                                                  
      J=I+2                                                             
      IF (J.GT.NORB) J=NORB                                             
      IM=0                                                              
      DO 22 L=I,J                                                       
      IF (NMAX(L).GT.IM) IM=NMAX(L)                                     
 22   CONTINUE                                                          
      DO 25 K=1,IM                                                      
      IF (((K-1)*(K-48*(K/48))).NE.0) GO TO 25                          
      WRITE (6,2000) BAR                                                
      WRITE (6,2004) (NQN(L),TITRE(L),NQN(L),TITRE(L),L=I,J)            
 25   WRITE (6,2005) DR(K),(DGC(K,L),DPC(K,L),L=I,J)                    
 29   CONTINUE                                                          
 31   IF (NPUN.EQ.0) GO TO 99                                           
      DO 35 I=1,NP                                                      
 35   DP(I)=DVF(I)*DR(I)                                                
 99   DO 101 I=1,NP                                                     
 101  D(I)=0.d0                                                           
      NAG=1                                                             
      IF (NMFG.EQ.0) GO TO 199                                          
      WRITE (6,2000) BAR                                                
      WRITE (6,2008)                                                    
 103  READ(5,1000) I1,I2,N1                                            
      IF (I1.LE.0) GO TO 199                                            
      IF (I2.GT.0) GO TO 105                                            
      IF (((I2+1)*(I2+2)).NE.0) GO TO 104                               
      IF (N1.LE.0) N1=1                                                 
      I1=N1                                                             
      N1=I2                                                             
      I2=I1                                                             
      GO TO 107                                                         
 104  WRITE (6,2001) I1,I2,N1                                           
      GO TO 103                                                         
 105  IF (I1.GT.NORB.OR.I2.GT.NORB) GO TO 104                           
      N1=1                                                              
 107  J1=2*IABS(NK(I1))-1                                               
      J2=2*IABS(NK(I2))-1                                               
      KMA=MIN0(J1,J2)                                                   
      NM=NMAX(I2)                                                       
      DO 121 J=1,KMA,2                                                  
      CALL YKDIR (I1,I1,J,NAG)                                          
      DO 111 I=1,NM                                                     
 111  DP(I)=DQ(I)*DGC(I,I2)*DPC(I,I2)                                   
      DVAL=J+1                                                          
      CALL SOMM (DR,D,DP,DPAS,DVAL,-1,NM)                               
 121  WRITE (6,2009) J,NQN(I1),TITRE(I1),NQN(I2),TITRE(I2),DVAL         
      IF (I1.EQ.I2) GO TO 161                                           
      J1=(IABS(1-2*NK(I1))-1)/2                                         
      J2=(IABS(1-2*NK(I2))-1)/2                                         
      KMA=MAX0(NQL(I1)+J2,NQL(I2)+J1)                                   
      J1=IABS(NQL(I2)-J1)                                               
      J2=IABS(NQL(I1)-J2)                                               
      KMI=MIN0(J1,J2)                                                   
      J1=KMI+NQL(I1)+NQL(I2)                                            
      J1=J1-2*(J1/2)                                                    
      IF (J1.EQ.0) KMI=KMI+1                                            
      NM=MIN0(NMAX(I1),NMAX(I2))                                        
      DO 141 J=KMI,KMA,2                                                
      CALL YKDIR (I1,I2,J,NAG)                                          
      DO 131 I=1,NM                                                     
      DP(I)=DQ(I)*DGC(I,I1)*DPC(I,I2)                                   
 131  DC(I)=DQ(I)*DGC(I,I2)*DPC(I,I1)                                   
      DVAL=J+1                                                          
      DVALP=DVAL                                                        
      DVALM=DVAL                                                        
      CALL SOMM (DR,D,DP,DPAS,DVALP,-1,NM)                              
      CALL SOMM (DR,D,DC,DPAS,DVAL,-1,NM)                               
      CALL YKDIR (I2,I1,J,NAG)                                          
      DO 135 I=1,NM                                                     
 135  DP(I)=DQ(I)*DGC(I,I2)*DPC(I,I1)                                   
      CALL SOMM (DR,D,DP,DPAS,DVALM,-1,NM)                              
 141  WRITE (6,2010) J,NQN(I1),TITRE(I1),NQN(I2),TITRE(I2),DVALM,DVAL,DVALP   
 161  IF (N1+1) 165,163,103                                             
 163  I1=I1+1                                                           
      I2=I1                                                             
      IF (I1-NORB) 107,107,103                                          
 165  I2=I2+1                                                           
      IF (I2-NORB) 107,107,163                                          
 199  IF (NMRK.EQ.0) GO TO 999                                          
      WRITE (6,2000) BAR                                                
      WRITE (6,2011)                                                    
 201  READ(5,1000) I1,I2,I3,I4,K                                       
      IF (I1.LE.0) GO TO 999                                            
      IF (I1.LE.NORB.AND.I2.GT.0.AND.I2.LE.NORB.AND.I3.GT.0.AND.I3.LE.NORB &
      .AND.I4.GT.0.AND.I4.LE.NORB.AND.K.GE.0) GO TO 205               
      WRITE (6,2001) I1,I2,I3,I4,K                                      
      GO TO 201                                                         
 205  CALL YKDIR (I1,I2,K,NAG)                                          
      DO 207 I=1,NP                                                     
 207  DP(I)=DQ(I)*DGC(I,I3)*DPC(I,I4)                                   
      DVAL=K+1                                                          
      CALL SOMM (DR,D,DP,DPAS,DVAL,-1,NP)                               
      WRITE (6,2012) K,NQN(I1),TITRE(I1),NQN(I2),TITRE(I2),NQN(I3),TITRE &
      (I3),NQN(I4),TITRE(I4),DVAL                                       
      GO TO 201                                                         
 999  RETURN                                                            
      END                                                               
