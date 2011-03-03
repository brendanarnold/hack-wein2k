      SUBROUTINE rel_CHARGE2(RX,S,FPIV,FUNCT,DELTA,JEND,SUM,zz)
      use param
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      DIMENSION RX(NRAD),FUNCT(NRAD)                                    
!                                                                       
      SUM=0.0D0                                                         
      NPOINT=JEND                                                       
      JS    =2-MOD(NPOINT,2)                                            
      JF    =JEND-2                                                     
!                                                                       
      DO J=JS,JF,2                                                    
         IF(rx(j).LT.zz*0.00005325d0) cycle
         SUM=SUM+FUNCT(J)*RX(J)+4.*FUNCT(J+1)*RX(J+1)+FUNCT(J+2)*RX(J+2)   
      ENDDO
!                                                                       
      SUM=SUM*DELTA/3.0D0                                               
      SUM=SUM+(1-MOD(NPOINT,2))*( 0.0d0+FUNCT(JS) )/ 2.0D0            &
                               *( RX(JS)-RX(1)       )                  
      SUM=SUM*FPIV*S                                                    
!                                                                       
      RETURN                                                            
      END                                                               
