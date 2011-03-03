!BOP
! !ROUTINE: FPot
! !INTERFACE:
      real*8 function fpot(R,Z,WA)
! !INPUT/OUTPUT PARAMETERS:
!   r   :  radial point
!   z   :  atomic number
!   wa  :  nombre d'electrons-z-1
! !DESCRIPTION:
!   Gives the Thomas Fermi potential at point R
! !REVISION HISTORY:
!     Taken from SRC_lcore.
!     Heavily modified November 2004 (Kevin Jorissen)
!EOP
      implicit none
!   INPUT
      real*8 r,z,wa
!   LOCALS
      real*8 wc,wd,we
                                                                       
! POTENTIEL THOMAS FERMI AU POINT R  Z NUMERO ATOMIQUE                  
! WA NOMBRE D ELECTRONS-Z-1                                             
! **********************************************************************
      WC=DSQRT((R*(Z+WA)**(dble(1)/dble(3)))/dble(0.8853))                               
      WD=WC*(0.60112*WC+1.81061)+1.                                     
      WE=WC*(WC*(WC*(WC*(0.04793*WC+0.21465)+0.77112)+1.39515)+1.81061)+1.                                                                
      WC=(Z+WA)*(WD/WE)**2-WA                                           
      FPOT=-WC/R                                                        
      RETURN                                                            
      END                                                               
