      SUBROUTINE ARBDOS (NFLOW,NFUN1,NFIRST,NLAST,EMIN,EMAX,NEMIN,NEMAX  &
       ,DET,ebs,fc,nst,RNUMB,DENSTY)                                                            
!  NEW TETRA.-SUBR. (P.BLOECHL, OKA, O. JEPSEN, M. ALOUANI              
!********************************************************************** 
!*                                                                   *  
!*  DENSITY OF STATES PROGRAM                                        *  
!*  CALCULATES THE PROJECTED DENSITIES OF STATES, CORRESPONDING      *  
!*  TO THE BAND ENERGIES AND ANGULAR MOMENTUM WEIGHTS, STORED ON     *  
!*  THE FILE 8                                                       *  
!*********************************************************************  
      INCLUDE 'param.inc'
      real*4 EBS(Nst,*), FC(MG,Nst,*)
      real*4 RNUMB(nlast,mg), DENSTY(nlast,mg)
      real*4,allocatable :: SNUMB(:,:)                                        
      COMMON /EME/ EMIN1, EMAX1, EFACTR, ESTEP, NO, NFUN, NU            
!                         
      allocate (SNUMB(nlast,MG))                           
      snumb=0.d0
      ESTEP=DET                                                         
      EFACTR=1.E0/DET                                                   
      NFUN=NFUN1                                                        
      EMIN1=EMIN                                                        
      EMAX1=EMAX                                                        
      NU=NFIRST                                                         
      NO=NLAST                                                          
      WRITE (6,10) EMIN,EMAX,EFACTR,ESTEP,NEMIN,NEMAX,NU,NO             
   10 FORMAT (1H ,' EMIN=',F10.5,' EMAX=',F10.5,' EFACTR=',F16.8, &
      ' ESTEP =',F10.5/' ENERGY BAND',I5,' THROUGH',I5,' ENERGY CHANNEL:', &
      I5,'   TO',I5)                                                            
      CALL ADOS (NFLOW,NFUN1,NEMIN,NEMAX,ebs,fc,nst, RNUMB,DENSTY,SNUMB,nlast)
!                                                                       
      NU1=NU+1                                                          
      DO 30 L=1,NFUN1                                                   
      DO 20 I=NU1,NO                                                    
      SNUMB(I,L+NFLOW)=SNUMB(I,L+NFLOW)+SNUMB(I-1,L+NFLOW)              
   20 RNUMB(I,L+NFLOW)=(RNUMB(I,L+NFLOW)+SNUMB(I,L+NFLOW))*2.           
   30 CONTINUE                                                          
!      DO 40 L=1,NFUN1                                                  
!      DO 40 I=NU,NO                                                    
!      TDOS(I)=TDOS(I)+DENSTY(I,L+NFLOW)                                
!   40 TNOS(I)=TNOS(I)+RNUMB(I,L+NFLOW)                                 
      RETURN                                                            
!   50 WRITE (6,60) Nlast                                                
!   60 FORMAT (' ARRAYS TOO SMALL.   Nlast=',I7)                         
!      STOP                                                              
      END                                                               
