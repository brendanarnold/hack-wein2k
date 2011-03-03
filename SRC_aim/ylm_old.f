      SUBROUTINE YLM(cth,sth,cfi,sfi,LOMAX,Y)                        
      implicit none
!      IMPLICIT REAL*8 (A-H,O-Z)

      COMPLEX*16 Y(*)

      real*8 cth,sth,cfi,sfi,p,pi,tpi,fpi,sqpi,sqtpi,zero,snull,one
      real*8 rf,yr,csr,cyp,sqfpi,c2l,tcth,cmfi,smfi,cd,sgnm,c1l,yi,cn

      integer i,l,m,l1,m1,l2,idwn,lm,lm2,lomax

      INCLUDE 'param.inc'
      DIMENSION P(lmax2+4,lmax2+4)                                      
      COMMON /CTES/ pi,tpi,fpi,sqpi,sqtpi
      DATA ZERO/0.00/,SNULL/1.E-6/,ONE/1.00/         
!.....CALCULATE SINES AND COSINES OF THE POLAR ANGLES OF THE VECTOR V.  
      sqfpi=2.0*sqpi
!      FPI=4.*ACOS(-1.D0)
!      XY=V(1)**2+V(2)**2                                                
!      XYZ=XY+V(3)**2                                                    
!      IF (XY.GT.SNULL)  GO TO 6                                         
!      CTH=ONE                                                           
!      IF (V(3).LT.ZERO)  CTH=-CTH                                       
!      STH=ZERO                                                          
!      CFI=ONE                                                           
!      SFI=ZERO                                                          
!      GO TO 7                                                           
! 6    XY=SQRT(XY)                                                       
!      XYZ=SQRT(XYZ)                                                     
!      CTH=V(3)/XYZ                                                      
!      STH=XY/XYZ                                                        
!      CFI=V(1)/XY                                                       
!      SFI=V(2)/XY                                                       
 7    RF=ONE/SQFPI                                                  
      YR=RF                                                             
      Y(1)=YR                                                           
      I=1                                                               
      P(1,1)=ONE                                                        
      P(2,1)=CTH                                                        
      C2L=CTH                                                           
      TCTH=CTH+CTH                                                      
      L1=2                                                              
      L=1                                                               
 5    M=1                                                               
      I=I+L                                                             
      IDWN=I+2                                                          
      M1=2                                                              
      L2=L                                                              
      L=L1                                                              
      CMFI=ONE                                                          
      SMFI=ZERO                                                         
      L1=L+1                                                            
      LM=L2                                                             
      LM2=L                                                             
      CD=ONE                                                            
      C2L=C2L+TCTH                                                      
      SGNM=ONE                                                          
!.....RECURSE UPWARD IN L.                                              
 42   P(L1,M)=(C2L*P(L,M)-LM*P(L2,M))/LM2                               
 44   C1L=(LM+1)*CTH                                                    
      P(L,M1)=ZERO                                                      
!.....RECURSE UPWARD IN M.                                              
      IF (ABS(STH).LT.SNULL)  GO TO 3                                   
      P(L,M1)=(C1L*P(L,M)-LM2*P(L1,M))/STH                              
 3    I=I+1                                                             
      IDWN=IDWN-1                                                       
    9 CSR=SQRT((2*L-ONE)/(FPI*CD))                                      
      CYP=SGNM*CSR*P(L,M)                                               
      YR=CYP*CMFI                                                       
      YI=CYP*SMFI                                                       
      Y(I)=DCMPLX(YR,YI)
      IF(IDWN.NE.I) Y(IDWN)=SGNM*DCMPLX(YR,-YI)
!      write(6,*) 'ylm ',i,yr,yi,y(i)
!      if(idwn.ne.i) write(6,*) 'ylm(idwn) ',idwn,yr,yi,y(idwn)
      CN=CMFI                                                           
      CMFI=CN*CFI-SFI*SMFI                                              
      SMFI=SFI*CN+SMFI*CFI                                              
      M=M1                                                              
      M1=M+1                                                            
      LM2=LM2-1                                                         
      LM=LM+1                                                           
      CD=CD*LM*LM2                                                      
      SGNM=-SGNM                                                        
      IF (M-L)  42,3,12                                                 
 12   IF (L.LE.LOMAX)  GOTO 5                                           
!      write (6,*) 'ylm ylm= ',
!     $   (y(jj),jj=1,(lomax)*(lomax+1))
      RETURN                                                            
      END                                                               


