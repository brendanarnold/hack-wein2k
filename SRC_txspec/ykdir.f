      SUBROUTINE YKDIR (IA,IB,NK1,NAG)                                  
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      COMMON DEN(30),DQ1(30),DFL(30),NQN(30),NQL(30), &
      NK(30),NMAX(30),NEL(30),NORB                                                         
      COMMON/DIRA/DV(NRAD+41),DR(NRAD+41),DP(NRAD+41), &
                  DQ(NRAD+41),DPAS,Z,NSTOP,NES,TETS,NP,NUC
      COMMON/DEUX/DVN(NRAD+41),DVF(NRAD+41),D(NRAD+41), &
                  DC(NRAD+41),DGC(NRAD+41,30),DPC(NRAD+41,30)
      COMMON/TROIS/ DPNO(4,30),DQNO(4,30)                               
      DIMENSION DPN1(4)                                                 
      DPAH=EXP(DPAS)                                                    
      DPYK=DPAS/24.                                                     
      ID=MIN0(NMAX(IA)+2,NMAX(IB)+2,NP)                                 
      IDM1=ID-1                                                         
      IF (NAG.NE.0) GO TO 4                                             
      DO 2 I=1,ID                                                       
 2    DQ(I)=DR(I)*(DGC(I,IA)*DGC(I,IB)+DPC(I,IA)*DPC(I,IB))             
      DO 3 I=1,4                                                        
      DPN1(I)=0.                                                        
      DO 3 J=1,I                                                        
 3    DPN1(I)=DPN1(I)+DPNO(J,IA)*DPNO(I+1-J,IB)+DQNO(J,IA)*DQNO(I+1-J,IB &
      )                                                                 
      GO TO 7                                                           
 4    DO 5 I=1,ID                                                       
 5    DQ(I)=DR(I)*DGC(I,IA)*DPC(I,IB)                                   
      DO 6 I=1,4                                                        
      DPN1(I)=0.                                                        
      DO 6 J=1,I                                                        
 6    DPN1(I)=DPN1(I)+DPNO(J,IA)*DQNO(I+1-J,IB)                         
 7    DI=DFL(IA)+DFL(IB)+NK1                                            
      DP(1)=0.                                                          
      DP(2)=0.                                                          
      DO 8 I=1,4                                                        
      DI=DI+1.                                                          
      DP(1)=DP(1)+(DR(1)**DI)*DPN1(I)/DI                                
 8    DP(2)=DP(2)+(DR(2)**DI)*DPN1(I)/DI                                
      DM=DPAH**(-NK1)                                                   
      DIM2=-DPYK*DM*DM                                                  
      DIM1=13.*DPYK*DM                                                  
      DI=13.*DPYK                                                       
      DIP1=-DPYK/DM                                                     
      DO 9 I=3,IDM1                                                     
 9    DP(I)=DP(I-1)*DM+DIM2*DQ(I-2)+DIP1*DQ(I+1)+DIM1*DQ(I-1)+DI*DQ(I)  
      DQ(ID-2)=DP(ID-2)                                                 
      DO 11 I=IDM1,NP                                                   
 11   DQ(I)=DQ(I-1)*DM                                                  
      I=NK1+NK1+1                                                       
      DM=DM/DPAH                                                        
      DIM2=I*DIM2/(DPAH*DPAH)                                           
      DIM1=I*DIM1/DPAH                                                  
      DI=I*DI                                                           
      DIP1=I*DIP1*DPAH                                                  
      I=ID-3                                                            
 14   DQ(I)=DQ(I+1)*DM+DIM2*DP(I+2)+DIP1*DP(I-1)+DIM1*DP(I+1)+DI*DP(I)  
      I=I-1                                                             
      IF (I-1) 15,15,14                                                 
 15   DQ(1)=DQ(3)*DM*DM+8.*((DI*DP(1)+4.*DIM1*DP(2))/13.-DIM2*DP(3))    
      RETURN                                                            
      END                                                               
