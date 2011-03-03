SUBROUTINE STERN(G,NST,STG,TAUP)
  !.... GENERATE STAR OF G 
  USE param
  use defs
  USE sym2
  IMPLICIT NONE
  !
  INTEGER I,J,K,L,M
  REAL*8  TK,TPI
  INTEGER        G(3),NST,STG(3,NSYM),IND(NSYM)
  COMPLEX*16     TAUP(NSYM)
!---------------------------------------------------------------------  
!                                         
  tpi=two*pi
  NST=0                                                             
  DO 1 I=1,IORD                                                     
     TK=0.                                                          
     DO 2 J=1,3                                                     
        TK=TK+TAU(J,I)*G(J)*TPI                                     
        K=0                                                         
        DO L=1,3                                                  
           K=IZ(J,L,I)*G(L)+K                                       
        ENDDO
        STG(J,I)=K                                                  
2    CONTINUE                                                       
     IF(NST.EQ.0) GOTO 7                                            
     DO 4 M=1,NST                                                   
        DO 5 J=1,3                                                  
           IF(STG(J,M).NE.STG(J,I)) GOTO 4                          
 5      CONTINUE                                                    
        IND(M)=IND(M)+1                                             
!        TAUP(M)=TAUP(M) + EXP(IMAG*TK)                              
!_REAL         TAUP(M) = TAUP(M) + COS(TK)
!_COMPLEX         TAUP(M) = TAUP(M) + DCMPLX(COS(TK),SIN(TK))
        GOTO 1                                                      
 4   CONTINUE                                                       
7    NST=NST+1                                                      
     DO J=1,3                                                     
        STG(J,NST)=STG(J,I)                                         
     ENDDO
     IND(NST)=1                                                     
!     TAUP(NST)=EXP(IMAG*TK)                                         
!_REAL         TAUP(NST) = COS(TK)
!_COMPLEX         TAUP(NST) = DCMPLX(COS(TK),SIN(TK))
1 CONTINUE                                                          
  DO I=1,NST                                                     
     TAUP(I)=TAUP(I)/IND(I)                                            
  ENDDO
  RETURN                                                            
END
