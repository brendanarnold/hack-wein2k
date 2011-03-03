      SUBROUTINE STERN                                                  
!                                                                       
!.... GENERATE STAR OF G                                                
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16          TAUP,IMAG                                     
      COMMON /FOUR/ G(3),NST,STG(3,48),TAUP(48)                         
      COMMON /SYM2/ IORD,IZ(3,3,48),ITAU,TAU(3,48)                      
      INTEGER G,STG,IND(48)                                             
!      DATA TPI/6.2831853071796d0/                                         
!      DATA IMAG        /(0.D0,1.D0)/                                    
      tpi=8.D0*DATAN(1.D0)
      NST=0                                                             
      DO 1 I=1,IORD                                                     
           TK=0.d0                                                 
      DO 2 J=1,3                                                        
               TK=TK+TAU(J,I)*G(J)*TPI                                  
      K=0                                                               
      DO 3 L=1,3                                                        
   3  K=IZ(J,L,I)*G(L)+K                                                
   2  STG(J,I)=K                                                        
      IF(NST.EQ.0) GOTO 7                                               
      DO 4 M=1,NST                                                      
      DO 5 J=1,3                                                        
      IF(STG(J,M).NE.STG(J,I)) GOTO 4                                   
   5  CONTINUE                                                          
      IND(M)=IND(M)+1                                                   
!      TAUP(M)=TAUP(M)+EXP(TK*IMAG)                                            
      TAUP(M)=TAUP(M)+cmplx(cos(TK),sin(TK))
      GOTO 1                                                            
   4  CONTINUE                                                          
   7  NST=NST+1                                                         
      DO 6 J=1,3                                                        
   6  STG(J,NST)=STG(J,I)                                               
!       TAUP(NST)=EXP(TK*IMAG)                                                 
      TAUP(NST)=cmplx(cos(tk),sin(tk))
      IND(NST)=1                                                        
   1  CONTINUE                                                          
      DO 10 I=1,NST                                                     
  10  TAUP(I)=TAUP(I)/IND(I)                                            
      RETURN                                                            
      END                                                               
