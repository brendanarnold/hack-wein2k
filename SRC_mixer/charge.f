      SUBROUTINE CHARGE (RX,DELTA,FUNCT,JSTART,JEND,SUM)                
!                                                                       
!-------------------------------------------------------------------    
!---  CHARGE CALCULATES THE INTEGRAL OF A NUMERIC GIVEN FUNCTION  --    
!---  BY THE SIMPSON ALGORITHM                                    --    
!---  THE FUNCTION VALUES HAVE TO BE STORED IN THE  ARRAY  FUNCT  --    
!---  AND CORRESPOND TO THE ARGUMENT VALUES IN RX. THE ALGORITHM  --    
!---  ASSUMES THAT THE ARGUMENTS IN RX  FOLLOW  AN   EXPONENTIAL  --    
!---  EXPRESSION RX(J)=R0*EXP(DELTA*(J-1)) , J=JSTART,JEND        --    
!---                                                              --    
!---                                 P.I.  SORANTIN               --    
!-------------------------------------------------------------------    
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      DIMENSION RX(NRAD),FUNCT(NRAD)                                    
!                                                                       
      SUM   =0.0                                                        
      NPOINT=JEND-JSTART+1                                              
      JS    =JSTART+1-MOD(NPOINT,2)                                     
      JF    =JEND-2                                                     
!                                                                       
      DO 1 J=JS,JF,2                                                    
  1      SUM=SUM+FUNCT(J)*RX(J)+4.*FUNCT(J+1)*RX(J+1)+FUNCT(J+2)*RX(J+2)
      SUM   =SUM*DELTA/3.0                                              
      SUM   =SUM + (1.-MOD(NPOINT,2))*( FUNCT(JSTART)+FUNCT(JS) ) /2.0   &
                                     *( RX(JS)-RX(JSTART)       )       
!                                                                       
      RETURN                                                            
      END                                                               
