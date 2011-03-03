      SUBROUTINE STERN(G,IZ,TAU,IORD)
!                                                                       
!.... GENERATE STAR OF G                                                
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
      INTEGER        STG,IND(noper)                                    
      COMPLEX*16        TAUP,IMAG                                       
      COMMON /FOUR/  NST,STG(3,noper),TAUP(noper)                    
      INTEGER        G(3)
      DIMENSION       IZ(3,3,noper),TAU(3,noper)     
      DATA           TPI/6.2831853071796d0/, IMAG/(0.D0,1.D0)/          
!---------------------------------------------------------------------  
!                                         
      tpi= 2.d0 * acos(-1.d0)                                           
      NST=0                                                             
      DO 1 I=1,IORD                                                     
         TK=0.                                                          
         DO 2 J=1,3                                                     
            TK=TK+TAU(J,I)*G(J)*TPI                                     
            K=0                                                         
            DO 3 L=1,3                                                  
               K=IZ(J,L,I)*G(L)+K                                       
 3          CONTINUE                                                    
            STG(J,I)=K                                                  
 2       CONTINUE                                                       
         IF(NST.EQ.0) GOTO 7                                            
         DO 4 M=1,NST                                                   
            DO 5 J=1,3                                                  
               IF(STG(J,M).NE.STG(J,I)) GOTO 4                          
 5          CONTINUE                                                    
            IND(M)=IND(M)+1                                             
            TAUP(M)=TAUP(M) + EXP(IMAG*TK)                              
            GOTO 1                                                      
 4       CONTINUE                                                       
 7       NST=NST+1                                                      
         DO 6 J=1,3                                                     
            STG(J,NST)=STG(J,I)                                         
 6       CONTINUE                                                       
         IND(NST)=1                                                     
         TAUP(NST)=EXP(IMAG*TK)                                         
 1    CONTINUE                                                          
      DO 10 I=1,NST                                                     
  10  TAUP(I)=TAUP(I)/IND(I)                                            
      RETURN                                                            
      END                                                               
