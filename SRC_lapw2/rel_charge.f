      SUBROUTINE rel_charge(RX,S,FPIV,FUNCT,DELTA,JEND,SUM,zz)                 
      use param
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      DIMENSION RX(NRAD),FUNCT(NRAD)                                    
      real*8   :: funct_rel(nrad)
!                                                                       
      SUM=0.0D0                                                         
      NPOINT=JEND                                                       
      JS    =2-MOD(NPOINT,2)                                            
      JF    =JEND-2                                                     
!                                                                       
      DO j=1,jend
         IF(rx(j).LT.zz*0.00005325d0) THEN 
            funct_rel(j)=funct(j)*2.d0/(zz*0.00005325d0)*rx(j)
         ELSE
            funct_rel(j)=funct(j)
         ENDIF
      ENDDO
      DO J=JS,JF,2                                                    
         SUM=SUM+FUNCT_rel(J)*RX(J)+4.*FUNCT_rel(J+1)*RX(J+1)+FUNCT_rel(J+2)*RX(J+2)   
      ENDDO
!                                                                       
      SUM=SUM*DELTA/3.0D0                                               
      SUM=SUM+(1-MOD(NPOINT,2))*( FUNCT_rel(1)+FUNCT_rel(JS) )/ 2.0D0            &
                               *( RX(JS)-RX(1)       )                  
      SUM=SUM*FPIV*S                                                    
!                                                                       
      RETURN                                                            
      END                                                               
