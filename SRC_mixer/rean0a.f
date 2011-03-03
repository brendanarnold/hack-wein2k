      SUBROUTINE REAN0a(NKK,KZZ,NAT,IFFT1,IFFT2,IFFT3,U,mult,rmt,v,pos,inversion) 
!                                                                       
!     SETS UP IN THE ARRAY U THE STEP FUNCTION FOR ALL K VECTORS        
!     USED IN REAN1 (REANALYSING THE FOURIER COEFFICIENTS)              
!                                                                       
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
!
      COMPLEX*16         U,IMAG,P,phs                                   
!                                                                       
!
!        ALAT(1:3)   - lattice constants (a,b,c)
!        ALPHA(1:3)  - angles of the unit-cell's unit-vectors
!        MULT(i)     - number of equivalent atoms for inequiv. atom i
!        PIA(1:3)    - reciprocal lattice constants (2*pi/a, 2*pi/b,
!                      2*pi/c)
!        POS(1:3,i)  - position vector of atom i (including equivalent
!                      atoms) in the unit-cell
!        RMT(i)      - muffin tin radius of atom i
!        V(i)        - relative muffin tin spherevolume for atom i
!        VI          - volume of the direct lattice unit-cell
!
      logical inversion
      INTEGER            MULT(NAT)
      DOUBLE PRECISION   VI
      DOUBLE PRECISION   ALAT(3), ALPHA(3), PIA(3), POS(3,Nat*48)
      DOUBLE PRECISION   RMT(NAT), V(NAT)
      COMMON  /STRUK/    ALAT, ALPHA, PIA, VI
      SAVE    /STRUK/
      COMMON /GENER/  BR1(3,3)                                          
      DIMENSION U(-IFFT1/2:IFFT1/2,-IFFT2/2:IFFT2/2,-IFFT3/2:IFFT3/2)         
      DIMENSION SK(3),SSK(3),KZZ(3,NKK)                 
!      complex*16,allocatable :: PHS(:)           
!                                                                       
      DATA ZILCH/1.D-10/,IMAG/(0.D0,1.D0)/                  
!--------------------------------------------------------------------   
!                                                                       
!      allocate ( PHS(NAT))           
      PI=ACOS(-1.D0)
      PI2=PI*2.D0
!.....DEFINE STEPFUNCTION U FOR K=0                                     
         U0=VI                                                         
         DO 5 JATOM=1,NAT                                               
            U0=U0-V(JATOM)*MULT(JATOM)                                  
 5       CONTINUE                                                       
      if(inversion)then
         DO 5000 I1=-IFFT1/2,IFFT1/2                                        
         DO 5000 I2=-IFFT2/2,IFFT2/2                                        
         DO 5000 I3=-IFFT3/2,IFFT3/2                                        
               SK(1)=I1                                                 
               SK(2)=I2                                                 
               SK(3)=I3                                                 
      SSK(1)=SK(1)*BR1(1,1)+SK(2)*BR1(1,2)+SK(3)*BR1(1,3)               
      SSK(2)=SK(1)*BR1(2,1)+SK(2)*BR1(2,2)+SK(3)*BR1(2,3)               
      SSK(3)=SK(1)*BR1(3,1)+SK(2)*BR1(3,2)+SK(3)*BR1(3,3)               
      S=SQRT(SSK(1)**2+SSK(2)**2+SSK(3)**2)                             
      IF(S.LT.ZILCH) THEN                                               
!                                                                       
!.....DEFINE STEPFUNCTION U FOR K=0                                     
         U(I1,I2,I3)=U0                                                 
!                                                                       
!.....DEFINE STEPFUNCTION FOR K#0                                       
      ELSE                                                              
   10 U(I1,I2,I3)=0.d0                                                    
      INDEX=0    
      utmp=0.d0                                                       
      DO 35 JATOM=1,NAT                                                 
         PHSR=0.d0                                            
         DO 40 MU=1,MULT(JATOM)                                         
            INDEX=INDEX+1                                               
!            ARG1=2.D0*PI*POS(1,INDEX)*SK(1)                             
!            ARG2=2.D0*PI*POS(2,INDEX)*SK(2)                             
!            ARG3=2.D0*PI*POS(3,INDEX)*SK(3)                             
            ARG1=POS(1,INDEX)*SK(1)                             
            ARG2=POS(2,INDEX)*SK(2)                             
            ARG3=POS(3,INDEX)*SK(3)  
            ARG=(ARG1+ARG2+ARG3)*PI2
!            P=cmplx(cos(arg),sin(arg)
!            P=EXP( IMAG*(ARG1+ARG2+ARG3) )                              
            PHSR=PHSR + cos(arg)                                   
 40      CONTINUE                                                       
                     Q=S*RMT(JATOM)
                     BESR=(Q*COS(Q)-SIN(Q))/(Q*Q*Q)
                     UTMP=UTMP+V(JATOM)*BESR*PHSR
 35   CONTINUE                                                          
!                                                                       
!  21  DO 20 JATOM=1,NAT                                                 
!         Q=S*RMT(JATOM)                                                   
!         BESR=(Q*COS(Q)-SIN(Q))/(Q*Q*Q)                                  
!         U(I1,I2,I3)=U(I1,I2,I3)+V(JATOM)*PHS(JATOM)*BESR               
! 20   CONTINUE                                                          
      U(I1,I2,I3)=3.0D0*UTMP                                            
      ENDIF                                                             
 5000 CONTINUE
                                                          
      else

         DO 5001 I1=-IFFT1/2,IFFT1/2                                        
         DO 5001 I2=-IFFT2/2,IFFT2/2                                        
         DO 5001 I3=-IFFT3/2,IFFT3/2                                        
               SK(1)=I1                                                 
               SK(2)=I2                                                 
               SK(3)=I3                                                 
      SSK(1)=SK(1)*BR1(1,1)+SK(2)*BR1(1,2)+SK(3)*BR1(1,3)               
      SSK(2)=SK(1)*BR1(2,1)+SK(2)*BR1(2,2)+SK(3)*BR1(2,3)               
      SSK(3)=SK(1)*BR1(3,1)+SK(2)*BR1(3,2)+SK(3)*BR1(3,3)               
      S=SQRT(SSK(1)**2+SSK(2)**2+SSK(3)**2)                             
      IF(S.LT.ZILCH) THEN                                               
!                                                                       
!.....DEFINE STEPFUNCTION U FOR K=0                                     
         U(I1,I2,I3)=U0                                                 
!                                                                       
!.....DEFINE STEPFUNCTION FOR K#0                                       
      ELSE                                                              
      U(I1,I2,I3)=0.d0                                                    
      INDEX=0                                                           
      DO 36 JATOM=1,NAT                                                 
         PHS=(0.,0.)                                             
         DO 41 MU=1,MULT(JATOM)                                         
            INDEX=INDEX+1                                               
!            ARG1=2.D0*PI*POS(1,INDEX)*SK(1)                             
!            ARG2=2.D0*PI*POS(2,INDEX)*SK(2)                             
!            ARG3=2.D0*PI*POS(3,INDEX)*SK(3)                             
            ARG1=POS(1,INDEX)*SK(1)                             
            ARG2=POS(2,INDEX)*SK(2)                             
            ARG3=POS(3,INDEX)*SK(3)  
            ARG=(ARG1+ARG2+ARG3)*PI2
!            P=cmplx(cos(arg),sin(arg)
!            P=EXP( IMAG*(ARG1+ARG2+ARG3) )                              
            PHS=PHS + DCMPLX(cos(arg),sin(arg))                                   
 41      CONTINUE                                                       
         Q=S*RMT(JATOM)
         BESR=(Q*COS(Q)-SIN(Q))/(Q*Q*Q)
         U(I1,I2,I3)=U(I1,I2,I3)+V(JATOM)*PHS*BESR
 36   CONTINUE                                                          
!                                                                       
!      DO 22 JATOM=1,NAT                                                 
!         Q=S*RMT(JATOM)                                                   
!         BESR=(Q*COS(Q)-SIN(Q))/(Q*Q*Q)                                  
!         U(I1,I2,I3)=U(I1,I2,I3)+V(JATOM)*PHS(JATOM)*BESR               
! 22   CONTINUE                                                          
      U(I1,I2,I3)=3.0D0*U(I1,I2,I3)                                     
      ENDIF                                                             
 5001 CONTINUE                                                          
      endif
      RETURN                                                            
      END                                                               
