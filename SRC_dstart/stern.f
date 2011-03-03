      SUBROUTINE STERN                                                  
!                                                                       
!.... GENERATE STAR OF G                                                
!                                                                       
      use struct
      IMPLICIT REAL*8 (A-H,O-Z)
!      include 'param.inc'
!p    PARAMETER (NSYM= 48)                                              
      INTEGER        G,STG,IND(NSYM)                                    
      COMPLEX*16        TAUP,IMAG , TMP                                      
      COMMON /FOUR/  G(3),NST,STG(3,NSYM),TAUP(NSYM)                    
!      COMMON /SYM2/  TAU(3,NSYM),IORD,IZ(3,3,NSYM)                      
      DATA           TPI/6.2831853071796d0/, IMAG/(0.D0,1.D0)/          
!---------------------------------------------------------------------  
!                                         
      tpi= 2.d0 * acos(-1.d0)                                           
      NST=0                                                             
      DO 1 I=1,IORD                                                     
         TK=0.                                                          
         DO 2 J=1,3                                                     
            TK=TK+TAU(J,I)*G(J)                                    
            K=0                                                         
            DO 3 L=1,3                                                  
               K=IZ(J,L,I)*G(L)+K                                       
 3          CONTINUE                                                    
            STG(J,I)=K                                                  
 2       CONTINUE            
         TK=TK*TPI
         TMP=DCMPLX(COS(TK),SIN(TK))                  
         IF(NST.EQ.0) GOTO 7                                            
         DO 4 M=1,NST                                                   
            DO 5 J=1,3                                                  
               IF(STG(J,M).NE.STG(J,I)) GOTO 4                          
 5          CONTINUE                                                    
            IND(M)=IND(M)+1                                             
            TAUP(M)=TAUP(M) + TMP                             
            GOTO 1                                                      
 4       CONTINUE                                                       
 7       NST=NST+1                                                      
         DO 6 J=1,3                                                     
            STG(J,NST)=STG(J,I)                                         
 6       CONTINUE                                                       
         IND(NST)=1                                                     
         TAUP(NST)=TMP                                       
 1    CONTINUE                                                          
      DO 10 I=1,NST                                                     
  10  TAUP(I)=TAUP(I)/IND(I)                                            
      RETURN                                                            
      END                                                               
