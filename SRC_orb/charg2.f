      SUBROUTINE CHARG2(RX,DELTA,FUNCT,JSTART,JEND,SUM)                 
      IMPLICIT REAL*8 (A-H,O-Z)
!                    ... DIMENSION OF YKA, QQ,..                        
      INCLUDE 'param.inc'
!
      DIMENSION RX(NRAD),FUNCT(NRAD)                                    
!     for even points 3/8 formula at beginning
      SUM=0.0D0                                                         
      NPOINT=JEND-JSTART+1                                              
      JS    =JSTART+3-3*MOD(NPOINT,2)                                     
      JF    =JEND-2                                                     
!                                                                       
      DO 1 J=JS,JF,2                                                    
 1    SUM=SUM+FUNCT(J)*RX(J)+4.*FUNCT(J+1)*RX(J+1)+FUNCT(J+2)*RX(J+2)   
!                                                                       
      SUM=SUM*DELTA/3.0D0                                               
      if(npoint.eq.2) then
      SUM=SUM+(1-MOD(NPOINT,2))*( FUNCT(JSTART)+FUNCT(JStart+1) )/ 2.0D0       &
                               *( RX(JStart+1)-RX(JSTART)       )             
      else if(mod(npoint,2).eq.0) then
      sum=sum+(FUNCT(jstart)*RX(jstart)+3.*FUNCT(jstart+1)*RX(jstart+1)+ &
       3.*FUNCT(jstart+2)*RX(jstart+2)+ &
       FUNCT(jstart+3)*RX(jstart+3))*delta*3./8.
      end if
!                                                                       
      RETURN                                                            
      END                                                               
