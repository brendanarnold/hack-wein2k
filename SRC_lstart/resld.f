      SUBROUTINE RESLD (NQN,NQL,NK,IMAX,DE,DFL,DQ1,JC,ISPIN,RELA)       
!                                                                       
! RESOLUTION DE L EQUATION DE DIRAC                                     
! NQN NOMBRE QUANTIQUE PRINCIPAL NQL NOMBRE QUANTIQUE ORBITAL           
! NK NOMBRE QUANTIQUE KAPPA   IMAX DERNIER POINT DE TABULATION DE LA    
! FONCTION D ONDE   DE ENERGIE   DFL PUISSANCE DU PREMIER TERME DU      
! DEVELOPPEMENT LIMITE   DQ1 PENTE A L ORIGINE DE DP OU DQ              
! **********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
!     IMPLICIT REAL (D)                                                 
      LOGICAL RELA                                                      
      COMMON/DIRA/DV(NPT,2),DR(NPT),DP(NPT),DQ(NPT),DPAS,Z,NSTOP,NES,    &
      TEST,NP,NUC                                                       
!                                                                       
! DV POTENTIEL EN U.A. ET NEGATIF  DR BLOC DES POINTS  DE TABULATION    
! DP GRANDE COMPOSANTE   DQ PETITE COMPOSANTE   DPAS PAS EXPONENTIEL    
! Z NUMERO ATOMIQUE   NSTOP CONTROLE L INTEGRATION NUMERIQUE            
! NES NOMBRE D ESSAIS POUR AJUSTER L ENERGIE                            
! TEST PRECISION A OBTENIR SUR L ENERGIE   NP NOMBRE MAXIMUM DE POINTS  
! NOYAU DE DIMENSIONS FINIES SI NUC NON NUL                             
! **********************************************************************
      COMMON/PS1/DEP(5),DEQ(5),DB,DVC,DSAL,DK,DM                        
!                                                                       
! DEP,DEQ DERIVEES DE OP ET DQ   DB=ENERGIE/DVC    DVC VITESSE DE LA    
! LUMIERE EN U.A.   DSAL=2.*DVC   DK NOMBRE QUANTIQUE KAPPA             
! DM=PAS EXPONENTIEL/720., DKOEF=1./720.                                
! **********************************************************************
      COMMON/TROIS/ DPNO(4,30),DQNO(4,30)                               
      DATA DKOEF/.1388888888888888D-2/                                  
      NSTOP=0                                                           
      DVC=137.0359895                                                      
      IF(.NOT.RELA) DVC=1.E30                                           
      DSAL=DVC+DVC                                                      
      IMM=0                                                             
      IES=0                                                             
      DK=NK                                                             
      LLL=(NQL*(NQL+1))/2                                               
      ND=0                                                              
      NOEUD=NQN-NQL                                                     
      IF (LLL.NE.0) GO TO 11                                            
      ELIM=-Z*Z/(1.5*NQN*NQN)                                           
      GO TO 19                                                          
 11   ELIM=DV(1,ISPIN)+LLL/(DR(1)*DR(1))                                
      DO 15 I=2,NP                                                      
      VAL=DV(I,ISPIN)+LLL/(DR(I)*DR(I))                                 
      IF (VAL.LE.ELIM) ELIM=VAL                                         
 15   CONTINUE                                                          
      IF(ELIM) 19,17,17                                                 
 17   NSTOP=17                                                          
! 2*V+L*(L+1)/R**2 EST PARTOUT POSITIF                                  
! **********************************************************************
      RETURN                                                            
 19   IF(DE.LE.ELIM) DE=ELIM*0.5D0                                      
 21   IF (IMM.EQ.1) GO TO 35                                            
      DO 25 I=7,NP,2                                                    
      IMAT=NP+1-I                                                       
      IF ((DV(IMAT,ISPIN)+LLL/(DR(IMAT)*DR(IMAT))-DE).LE.0.) GO TO 26   
 25   CONTINUE                                                          
 26   IF (IMAT.GT.5) GO TO 35                                           
      DE=DE*0.5D0                                                       
      IF(DE.LT.-TEST.AND.ND.LE.NOEUD) GO TO 21                          
 28   NSTOP=28                                                          
! 2*V+L*(L+1)/R**2-2*E EST PARTOUT POSITIF                              
! **********************************************************************
      RETURN                                                            
! VALEURS INITIALES POUR L INTEGRATION VERS L EXTERIEUR                 
! **********************************************************************
 35   DB=DE/DVC                                                         
      CALL INOUH (DP,DQ,DR,DQ1,DFL,DV(1,ISPIN),Z,TEST,NUC,NSTOP,JC)     
      IF (NSTOP) 399,47,399                                             
!     NSTOP=45                                                          
! LE DEVELOPPEMENT A L ORIGINE NE CONVERGE PAS                          
! **********************************************************************
 47   ND=1                                                              
      DO 51 I=1,5                                                       
      DVAL=DR(I)**DFL                                                   
      IF (I.EQ.1) GO TO 50                                              
      IF (DP(I-1).EQ.0.) GO TO 50                                       
      IF ((DP(I)/DP(I-1)).GT.0.) GO TO 50                               
      ND=ND+1                                                           
 50   DP(I)=DP(I)*DVAL                                                  
      DQ(I)=DQ(I)*DVAL                                                  
      DEP(I)=DEP(I)*DVAL                                                
 51   DEQ(I)=DEQ(I)*DVAL                                                
      K=-1+2*(NOEUD-2*(NOEUD/2))                                        
      IF ((DP(1)*K).GT.0.) GO TO 54                                     
 53   NSTOP=53                                                          
! ERREUR POUR LE DEVELOPPEMENT A L ORIGINE                              
! **********************************************************************
      RETURN                                                            
 54   IF ((K*NK*DQ(1)).LT.0.) GO TO 53                                  
      DM=DPAS*DKOEF                                                     
! INTEGRATION VERS L EXTERIEUR                                          
! **********************************************************************
      DO 195 I=6,IMAT                                                   
      DP(I)=DP(I-1)                                                     
      DQ(I)=DQ(I-1)                                                     
      CALL INTH (DP(I),DQ(I),DV(I,ISPIN),DR(I))                         
      IF (DP(I-1).EQ.0.) GO TO 195                                      
      IF ((DP(I)/DP(I-1)).GT.0.) GO TO 195                              
      ND=ND+1                                                           
      IF(ND.GT.NOEUD) GO TO 209                                         
 195  CONTINUE                                                          
      IF (ND.EQ.NOEUD) GO TO 240                                        
      DE=0.8*DE                                                         
      IF(DE.LT.-TEST) GO TO 21                                          
 206  NSTOP=206                                                         
! LE NOMBRE DE NOEUDS EST TROP PETIT                                    
! **********************************************************************
      RETURN                                                            
 209  DE=1.2*DE                                                         
      IF(DE.GT.ELIM) GO TO 21                                           
 210  NSTOP=210                                                         
! LE NOMBRE DE NOEUDS EST TROP GRAND                                    
! **********************************************************************
      RETURN                                                            
! VALEURS INITIALES POUR L INTEGRATION VERS L INTERIEUR                 
! **********************************************************************
 240  DQM=DQ(IMAT)                                                      
      DPM=DP(IMAT)                                                      
      IF (IMM.EQ.1) GO TO 258                                           
      DO 255 I=1,NP,2                                                   
      IMAX=NP+1-I                                                       
      IF (((DV(IMAX,ISPIN)-DE)*DR(IMAX)*DR(IMAX)).LE.300.) GO TO 258    
 255  CONTINUE                                                          
 258  DD=SQRT(-DE*(2.+DB/DVC))                                          
      DPQ=-DD/(DSAL+DB)                                                 
      DM=-DM                                                            
      DO 277 I=1,5                                                      
      J=IMAX+1-I                                                        
      DP(J)=EXP(-DD*DR(J))                                              
      DEP(I)=-DD*DP(J)*DR(J)                                            
      DQ(J)=DPQ*DP(J)                                                   
 277  DEQ(I)=DPQ*DEP(I)                                                 
      M=IMAX-5                                                          
! INTEGRATION VERS L INTERIEUR                                          
!***********************************************************************
      DO 301 I=IMAT,M                                                   
      J=M+IMAT-I                                                        
      DP(J)=DP(J+1)                                                     
      DQ(J)=DQ(J+1)                                                     
 301  CALL INTH (DP(J),DQ(J),DV(J,ISPIN),DR(J))                         
! RACCORDEMENT DE LA GRANDE COMPOSANTE                                  
! **********************************************************************
      DVAL=DPM/DP(IMAT)                                                 
      IF (DVAL.GT.0.) GO TO 313                                         
      NSTOP=312                                                         
! ERREUR SUR LE SIGNE DE LA GRANDE COMPOSANTE                           
! **********************************************************************
      RETURN                                                            
 313  DO 315 I=IMAT,IMAX                                                
      DP(I)=DP(I)*DVAL                                                  
 315  DQ(I)=DQ(I)*DVAL                                                  
! CALCUL DE LA NORME                                                    
! **********************************************************************
      DSUM=3.*DR(1)*(DP(1)**2+DQ(1)**2)/(DPAS*(DFL+DFL+1.))             
      DO 333 I=3,IMAX,2                                                 
 333  DSUM=DSUM+DR(I)*(DP(I)**2+DQ(I)**2)+4.*DR(I-1)*(DP(I-1)**2+DQ(I-1) &
      **2)+DR(I-2)*(DP(I-2)**2+DQ(I-2)**2)                              
      DSUM=DPAS*(DSUM+DR(IMAT)*(DQM*DQM-DQ(IMAT)*DQ(IMAT)))              &
      *0.3333333333333333D0                                             
! MODIFICATION DE L ENERGIE                                             
! **********************************************************************
      DBE=DP(IMAT)*(DQM-DQ(IMAT))*DVC/DSUM                              
      IMM=0                                                             
      VAL=ABS(DBE/DE)                                                   
      IF (VAL.LE.TEST) GO TO 365                                        
 340  DVAL=DE+DBE                                                       
      IF(DVAL.LT.0.) GO TO 360                                          
      DBE=DBE*0.5D0                                                     
      VAL=VAL*0.5D0                                                     
      IF (VAL.GT.TEST) GO TO 340                                        
 345  NSTOP=345                                                         
! ENERGIE NULLE                                                         
! **********************************************************************
      RETURN                                                            
 360  DE=DVAL                                                           
      IF (VAL.LE.0.1) IMM=1                                             
      IES=IES+1                                                         
      IF(IES.LE.NES) GO TO 21                                           
 362  NSTOP=362                                                         
! NOMBRE D ITERATIONS TROP GRAND                                        
! **********************************************************************
      RETURN                                                            
 365  DSUM=SQRT(DSUM)                                                   
      DQ1=DQ1/DSUM                                                      
      DO 367 I=1,IMAX                                                   
      DP(I)=DP(I)/DSUM                                                  
 367  DQ(I)=DQ(I)/DSUM                                                  
      DO 368 I=1,4                                                      
      DPNO(I,JC)=DPNO(I,JC)/DSUM                                        
 368  DQNO(I,JC)=DQNO(I,JC)/DSUM                                        
      IF(IMAX.EQ.NP) GO TO 398                                          
      J=IMAX+1                                                          
      DO 371 I=J,NP                                                     
      DP(I)=0.                                                          
 371  DQ(I)=0.                                                          
 398  NSTOP=0                                                           
 399  RETURN                                                            
      END                                                               

