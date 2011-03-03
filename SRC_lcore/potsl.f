      SUBROUTINE POTSL (DV,D,DP,DR,DPAS,DEXV,Z,NP,ION,ICUT)             
!                                                                       
! POTENTIEL INTEGRATION PAR UNE METHODE A 4 POINTS                      
! DV POTENTIEL   D DENSITE   DP BLOC DE TRAVAIL  DR POINTS DE TABULATION
! DPAS PAS EXPONENTIEL   DEXV COEFFICIENT MULTIPLICATIF POUR L ECHANGE  
! Z NUMERO ATOMIQUE   NP NOMBRE DE POINTS   ION=Z-NOMBRE D ELECTRONS    
! SI ICUT EST NUL ON CORRIGE EVENTUELLEMENT LE POTENTIEL EN -(ION+1)/R  
! **********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DV(*),D(*),DP(*),DR(*)                                  
      DAS=DPAS/24.d0                                                      
      DO 1 I=1,NP                                                       
 1    DV(I)=D(I)*DR(I)                                                  
      DLO=EXP(DPAS)                                                     
      DLO2=DLO*DLO                                                      
      DP(2)=DR(1)*(D(2)-D(1)*DLO2)/(12.d0*(DLO-1.d0))                       
      DP(1)=DV(1)/3.d0-DP(2)/DLO2                                         
      DP(2)=DV(2)/3.d0-DP(2)*DLO2                                         
      J=NP-1                                                            
      DO 3 I=3,J                                                        
 3    DP(I)=DP(I-1)+DAS*(13.d0*(DV(I)+DV(I-1))-(DV(I-2)+DV(I+1)))         
      DP(NP)=DP(J)                                                      
      DV(J)=DP(J)                                                       
      DV(NP)=DP(J)                                                      
      DO 5 I=3,J                                                        
      K=NP+1-I                                                          
 5    DV(K)=DV(K+1)/DLO+DAS*(13.d0*(DP(K+1)/DLO+DP(K))-(DP(K+2)/DLO2+DP(K- &
      1)*DLO))                                                          
      DV(1)=DV(3)/DLO2+DPAS*(DP(1)+4.d0*DP(2)/DLO+DP(3)/DLO2)/3.d0          
      DLO=-(ION+1)                                                      
      DO 7 I=1,NP                                                       
      DV(I)=DV(I)-(Z+3.d0*DEXV*((DR(I)*D(I)/105.27578d0)**(1.d0/3.d0)))       
      IF (ICUT.NE.0) GO TO 7                                            
      IF (DV(I).GT.DLO) DV(I)=DLO                                       
 7    DV(I)=DV(I)/DR(I)                                                 
      RETURN                                                            
      END                                                               
