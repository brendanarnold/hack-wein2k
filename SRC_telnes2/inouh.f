!BOP
! !ROUTINE: Inouh
! !INTERFACE:
      SUBROUTINE INOUH (DQ1,DFL,TEST,JC)        
! !USES:
      use initialstate1
! !INPUT/OUTPUT PARAMETERS:
!   dq1      :   pente à l'origine de dp ou dq
!   dfl      :   puissance du premier terme du développement limité
!   test     :   test de précision
!   jc       :
! !DESCRIPTION:
! VALEURS INITIALES POUR L INTEGRATION VERS L EXTERIEUR                 
! !REVISION HISTORY:
!     Taken from SRC_lcore.
!     Updated November 2004 (Kevin Jorissen)
!EOP


! DP GRANDE COMPOSANTE    DQ PETITE COMPOSANTE    DR BLOC DES POINTS    
! DV POTENTIEL AU PREMIER POINT               
! Z NUMERO ATOMIQUE
! NOYAU DE DIMENSIONS FINIES SI NUC NON NUL                             
! NSTOP CONTROLE LA CONVERGENCE DU DEVELOPPEMENT LIMITE                 
! **********************************************************************
	  implicit none
!  INPUT VARIABLES
      integer jc
	  real*8 dq1, dfl,test
!  COMMON VARIABLES
      real*8 DEP(5),DEQ(5),DD,DVC,DSAL,DK,DM,DPNO(4,30),DQNO(4,30)
      COMMON/PS1/DEP,DEQ,DD,DVC,DSAL,DK,DM                        
!                                                                       
! DEP,DEQ DERIVEES DE DP ET DQ   DD=ENERGIE/DVC    DVC VITESSE DE LA    
! LUMIERE EN U.A.   DSAL=2.*DVC   DK NOMBRE QUANTIQUE KAPPA             
! DM=PAS EXPONENTIEL/720.                                               
! **********************************************************************
      COMMON/TROIS/ DPNO,DQNO      
!  LOCAL VARIABLES
      integer i,j,m
	  real*8 dpr,dqr,dsum,dbe,deva3,deva2,deva1,dval

      DP(1:10)=dble(0)
      DQ(1:10)=dble(0)
      DVAL=Z/DVC                                                        
      DEVA1=-DVAL
      DEVA2=DV(1)/DVC+DVAL/DR(1)-DD                                        
      DEVA3=dble(0)
      IF (DK) 9,9,11                                                    
 9    DBE=(DK-DFL)/DVAL                                                 
      GO TO 15                                                          
 11   DBE=DVAL/(DK+DFL)                                                 
 15   DQ(10)=DQ1                                                        
      DP(10)=DBE*DQ1                                                    
      GO TO 39
 21   DVAL=DV(1)+Z*(3.-DR(1)*DR(1)/(DR(NUC)*DR(NUC)))/(DR(NUC)+DR(NUC))    
      DEVA1=0.                                                          
      DEVA2=(DVAL-3.*Z/(DR(NUC)+DR(NUC)))/DVC-DD                        
      DEVA3=Z/(DR(NUC)*DR(NUC)*DR(NUC)*DSAL)                            
      IF (DK) 33,33,35                                                  
 33   DP(10)=DQ1                                                        
      GO TO 39                                                          
 35   DQ(10)=DQ1                                                        
 39   do I=1,5                                                       
        DP(I)=DP(10)                                                      
        DQ(I)=DQ(10)                                                      
        DEP(I)=DP(I)*DFL                                                  
        DEQ(I)=DQ(I)*DFL                                                  
      enddo
      M=1                                                               
 41   DM=M+DFL                                                          
      DSUM=DM*DM-DK*DK+DEVA1*DEVA1                                      
      DQR=(DSAL-DEVA2)*DQ(M+9)-DEVA3*DQ(M+7)                            
      DPR=DEVA2*DP(M+9)+DEVA3*DP(M+7)                                   
      DVAL=((DM-DK)*DQR-DEVA1*DPR)/DSUM                                 
      DSUM=((DM+DK)*DPR+DEVA1*DQR)/DSUM                                 
      J=-1                                                              
      do I=1,5                                                       
        DPR=DR(I)**M                                                      
        DQR=DSUM*DPR                                                      
        DPR=DVAL*DPR                                                      
        IF (M.EQ.1) GO TO 43                                              
        IF (ABS(DPR/DP(I)).LE.TEST.AND.ABS(DQR/DQ(I)).LE.TEST) J=1        
 43     DP(I)=DP(I)+DPR                                                   
        DQ(I)=DQ(I)+DQR                                                   
        DEP(I)=DEP(I)+DPR*DM                                              
        DEQ(I)=DEQ(I)+DQR*DM
      enddo
      IF (J.EQ.1) GO TO 99                                              
      DP(M+10)=DVAL                                                     
      DQ(M+10)=DSUM                                                     
      M=M+1                                                             
      IF (M.LE.20) GO TO 41                                             
      NSTOP=45                                                          
 99   do I=1,4                                                      
        DPNO(I,JC)=DP(I+9)                                                
        DQNO(I,JC)=DQ(I+9)
      enddo
 999  RETURN                                                            
      END                                                               
