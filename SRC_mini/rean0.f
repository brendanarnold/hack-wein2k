      SUBROUTINE REAN0 (NKK,KZZ,NAT,IFFT1,IFFT2,IFFT3,U) 
!                                                                       
!     SETS UP IN THE ARRAY U THE STEP FUNCTION FOR ALL K VECTORS        
!     USED IN REAN1 (REANALYSING THE FOURIER COEFFICIENTS)              
!                                                                       
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
!
      COMPLEX*16         U,IMAG,PHS,P                                   
!                                                                       
!
!        ALAT(1:3)   - lattice constants (a,b,c)
!        ALPHA(1:3)  - angles of the unit-cell's unit-vectors
!        IATNR(i)    - atom index of (inequivalent) atom i
!                      also indicates cubic and non-cubic symmetries
!                      IATNR(i) .GT. 0 ... cubic symmetry
!                      IATNR(i) .LT. 0 ... non-cubic symmetry
!        MULT(i)     - number of equivalent atoms for inequiv. atom i
!        PIA(1:3)    - reciprocal lattice constants (2*pi/a, 2*pi/b,
!                      2*pi/c)
!        POS(1:3,i)  - position vector of atom i (including equivalent
!                      atoms) in the unit-cell
!        RMT(i)      - muffin tin radius of atom i
!        V(i)        - relative muffin tin spherevolume for atom i
!        VI          - volume of the direct lattice unit-cell
!
      INTEGER            IATNR(NATO), MULT(NATO)
      DOUBLE PRECISION   VI
      DOUBLE PRECISION   ALAT(3), ALPHA(3), PIA(3), POS(3,NDIF)
      DOUBLE PRECISION   RMT(NATO), V(NATO)
      COMMON  /STRUK/    POS, ALAT, ALPHA, RMT, V, PIA, VI, IATNR, MULT
      SAVE    /STRUK/
      COMMON /GENER/  BR1(3,3)                                          
      DIMENSION U(-IFFT1/2:IFFT1/2,-IFFT2/2:IFFT2/2,-IFFT3/2:IFFT3/2)         
      DIMENSION SK(3),PHS(NATO),SSK(3),KZZ(3,NKK)                 
!                                                                       
      DATA ZILCH/1.D-10/,IMAG/(0.D0,1.D0)/                  
!--------------------------------------------------------------------   
!                                                                       
      PID=DACOS(-1.d0)
!.....DEFINE STEPFUNCTION U FOR K=0                                     
         U0=VI                                                         
         DO 5 JATOM=1,NAT                                               
            U0=U0-V(JATOM)*MULT(JATOM)                                  
 5       CONTINUE                                                       
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
   10 U(I1,I2,I3)=0.                                                    
      INDEX=0                                                           
      DO 35 JATOM=1,NAT                                                 
         PHS(JATOM)=(0.,0.)                                             
         DO 40 MU=1,MULT(JATOM)                                         
            INDEX=INDEX+1                                               
            ARG1=2.D0*PID*POS(1,INDEX)*SK(1)                             
            ARG2=2.D0*PID*POS(2,INDEX)*SK(2)                             
            ARG3=2.D0*PID*POS(3,INDEX)*SK(3)                             
            P=EXP( IMAG*(ARG1+ARG2+ARG3) )                              
            PHS(JATOM)=PHS(JATOM) + P                                   
 40      CONTINUE                                                       
 35   CONTINUE                                                          
!                                                                       
  21  DO 20 JATOM=1,NAT                                                 
         Q=S*RMT(JATOM)                                                   
         BESR=(Q*COS(Q)-SIN(Q))/(Q**3)                                  
         U(I1,I2,I3)=U(I1,I2,I3)+V(JATOM)*PHS(JATOM)*BESR               
 20   CONTINUE                                                          
      U(I1,I2,I3)=3.0D0*U(I1,I2,I3)                                     
      ENDIF                                                             
 5000 CONTINUE                                                          
      RETURN                                                            
      END                                                               
