      FUNCTION FPOT (R,Z,WA)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 FPOT
!                                                                       
! POTENTIEL THOMAS FERMI AU POINT R  Z NUMERO ATOMIQUE                  
! WA NOMBRE D ELECTRONS-Z-1                                             
! **********************************************************************
      WC=SQRT((R*(Z+WA)**(1.d0/3.d0))/0.8853d0)                               
      WD=WC*(0.60112d0*WC+1.81061d0)+1.d0                                     
      WE=WC*(WC*(WC*(WC*(0.04793d0*WC+0.21465d0)+0.77112d0)+1.39515d0)+1.81061d0)+ &
      1.d0                                                                
      WC=(Z+WA)*(WD/WE)**2-WA                                           
      FPOT=-WC/R                                                        
      RETURN                                                            
      END                                                               
