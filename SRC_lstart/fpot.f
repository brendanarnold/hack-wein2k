      FUNCTION FPOT (R,Z,WA)
!                                                                       
! POTENTIEL THOMAS FERMI AU POINT R  Z NUMERO ATOMIQUE                  
! WA NOMBRE D ELECTRONS-Z-1                                             
! **********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      WC=SQRT((R*(Z+WA)**(1./3.))/0.8853)                               
      WD=WC*(0.60112*WC+1.81061)+1.                                     
      WE=WC*(WC*(WC*(WC*(0.04793*WC+0.21465)+0.77112)+1.39515)+1.81061)+ &
      1.                                                                
      WC=(Z+WA)*(WD/WE)**2-WA                                           

      FPOT=-WC/R
      
      RETURN                                                            
      END                                                               
