      SUBROUTINE CHARG3(RX,DELTA,FUNCT,JSTART,JEND,SUM)                 
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
!
      DIMENSION RX(NRAD),FUNCT(NRAD)                                    
!     for even points 3/8 formula at end
      SUM=0.0D0                                                         
      NPOINT=JEND-JSTART+1                                              
      JS    =JSTART                                   
      JF    =JEND-2-3+3*MOD(NPOINT,2)                                                       
!                                                                       
      DO 1 J=JS,JF,2                                                    
 1    SUM=SUM+FUNCT(J)*RX(J)+4.*FUNCT(J+1)*RX(J+1)+FUNCT(J+2)*RX(J+2)   
!                                                                       
      SUM=SUM*DELTA/3.0D0                                               
      if(npoint.eq.2) then
      SUM=SUM+(1-MOD(NPOINT,2))*( FUNCT(Jend)+FUNCT(Jend-1) )/ 2.0D0       &
                               *( RX(Jend)-RX(Jend-1)       )             
      else if(mod(npoint,2).eq.0) then
      sum=sum+(FUNCT(jend)*RX(jend)+3.*FUNCT(jend-1)*RX(jend-1)+ &
       3.*FUNCT(jend-2)*RX(jend-2)+ &
       FUNCT(jend-3)*RX(jend-3))*delta*3./8.
      end if
!                                                                       
      RETURN                                                            
      END                                                               
