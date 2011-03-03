      SUBROUTINE GENER(NSA,NSB,NSC,BR1)                           
      use reallocate
      use atpos
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
      DIMENSION BR1(3,3)                                                
      nnpos=3000                                    
      allocate ( ATP(3,nnpos))                                    
      NPOS=0                                                            
      JA=-NSA-1                                                         
   10 JA=JA+1                                                           
      IF(JA.GT.NSA) GOTO 1                                              
      JB=-NSB-1                                                         
   11 JB=JB+1                                                           
      IF(JB.GT.NSB) GOTO 10                                             
      JC=-NSC-1                                                         
   12 JC=JC+1                                                           
      IF(JC.GT.NSC) GOTO 11                                             
      NPOS=NPOS+1                                                       
      IF(NPOS.GT.nnpos) then
        nnpos=nnpos*2
        call doreallocate(atp,3,nnpos)
      endif                                 
      DO 2 I=1,3                                                        
   2  ATP(I,NPOS)=BR1(1,I)*JA+BR1(2,I)*JB+BR1(3,I)*JC                   
      GOTO 12                                                           
    1 CONTINUE                                                          
              WRITE(6,20) NPOS                                          
  20  FORMAT('0NUMBER OF CELL ORIGINS:',I5)                             
      RETURN                                                            
      END                                                               
