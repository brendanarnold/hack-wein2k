!BOP
! !ROUTINE: ArbDos
! !INTERFACE:
      SUBROUTINE ARBDOS (NFLOW,NFUN1,NFIRST,NLAST,EMIN,EMAX,NEMIN,NEMAX  &
       ,DET,ebs,fc,nst,RNUMB,DENSTY)                                                            
! !USES:
       use tetra_params  
       use program_control,only : verbosity
! !INPUT/OUTPUT PARAMETERS:
!     nflow  : always zero
!     nfun1  : number of dos cases
!     nfirst : index of first dos case (always = 1 in this program)
!     nlast  : index of last dos case
!     emin   : lower end of energy range
!     emax   : upper end of energy range
!     nemin  : index of lowest energy point
!     nemax  : index of highest energy point
!     det    : step size of energy mesh
!     ebs    :
!     fc     :
!     nst    : number of k-points
!     rnumb  :
!     densty :
! !DESCRIPTION:
!  DENSITY OF STATES PROGRAM                                        *  
!  CALCULATES THE PROJECTED DENSITIES OF STATES, CORRESPONDING      *  
!  TO THE BAND ENERGIES AND ANGULAR MOMENTUM WEIGHTS, STORED ON     *  
!  THE FILE 8                                                       *  
!
!   Note Kevin Jorissen : I allowed some shit to remain in the code for 
!   compatibility reasons.
! !REVISION HISTORY:
!     NEW TETRA.-SUBR. (P.BLOECHL, OKA, O. JEPSEN, M. ALOUANI              
!     Taken from SRC\_tetra.
!     Updated November 2004 (Kevin Jorissen)
!     Cleanup February 2005 (Kevin Jorissen)
!EOP

!  IN/OUTPUT  !KJ
      implicit none  ! added KJ
      integer nflow,nfun1,nfirst,nlast,nst,nemin,nemax !KJ
      real*8 emin,emax,det,ebs(nst,*),fc(mg,nst,*),rnumb(nlast,mg),densty(nlast,mg)  !KJ

!  LOCALS   !KJ
       real*8,allocatable :: snumb(:,:)
       COMMON /EME/ EMIN1, EMAX1, EFACTR, ESTEP, NO, NFUN, NU            
	   integer i,l,nu1,nu,no,nfun         !KJ declaration added (variables existed already)
	   real*8 estep, efactr, emin1, emax1 !KJ declaration added (variables existed already)
       allocate(snumb(nlast,mg))  ! added KJ


      snumb(:,:)=dble(0)
      ESTEP=DET                                                         
      EFACTR=dble(1)/DET                                                   
      NFUN=NFUN1                                                        
      EMIN1=EMIN                                                        
      EMAX1=EMAX                                                        
      NU=NFIRST                                                         
      NO=NLAST                                                          

      if (verbosity.ge.2) WRITE (6,10) EMIN,EMAX,EFACTR,ESTEP,NEMIN,NEMAX,NU,NO             
      CALL ADOS (NFLOW,NFUN1,NEMIN,NEMAX,ebs,fc,nst, RNUMB,DENSTY,SNUMB,nlast)
      NU1=NU+1                                                                                                                            
      
      do L=1,NFUN1                                                   
        do I=NU1,NO                                      
          SNUMB(I,L+NFLOW)=SNUMB(I,L+NFLOW)+SNUMB(I-1,L+NFLOW)              
          RNUMB(I,L+NFLOW)=(RNUMB(I,L+NFLOW)+SNUMB(I,L+NFLOW))*2.           
        enddo
      enddo                  

      RETURN                                                            

   10 FORMAT (1H ,' EMIN=',F10.5,' EMAX=',F10.5,' EFACTR=',F16.8, &
      ' ESTEP =',F10.5/' ENERGY BAND',I5,' THROUGH',I5,' ENERGY CHANNEL:',I5,'   TO',I5)                                                            
      END                                                               

