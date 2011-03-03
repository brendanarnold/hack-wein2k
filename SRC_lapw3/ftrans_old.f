      SUBROUTINE FTRANS (BCOF,LM,NPOINT,R,H,P,NDIM)                     
!CC                                                                     
!CC SPHERICAL BESSEL TRANSFORM                                          
!CC METHOD: EXTENDED FILON INTEGRATION                                  
!CC                                                                     
!CC     SOMMER & ZABOLIIZKY, COMP.PHYS.COMM.,VOL16,383(1979)            
!CC                                                                     
!CC PARAMETER LIST                                                      
!CC BCOF  : OUTPUT ARRAY CONTAINING THE INTEGRATION COEFFICIENTS        
!CC LM    : MAXIMUM L TO BE CALCULATED + 1                              
!CC            LM .GE. 1                                                
!CC NPOINT: NUMBER OF RADIUS POINTS, MUST BE ODD                        
!CC R     : STARTING RADIUS            .GE. 0                           
!CC H     : INCREMENT   .GT. 0                                          
!CC P     : TRANSFORMED VARIABLE       .GE. 0                           
!CC BINT  : AUXILLIARY ARRAY, MINIMUM LENGTH 9*LM                       
!CC B     : AUXILLIARY ARRAY, MINIMUM LENGTH   LM                       
!CC NDIM  : DECLARED FIRST DIMENSION OF ARRAY BCOF                      
!CC                                                                     
!CC APPLICATION                                                         
!CC   1. BUILD BCOF BY CALL TO FTRANS                                   
!CC   2. BESSEL TRANSFORM OF VECTOR F IS CALCULATED AS                  
!CC      SCALAR PRODUCT OF VECTORS F AND BCOF(*,L+1)                    
!CC                                                                     
!CC SUBROUTINES CALLED :  SININT, DFAC.                                 
!CC                                                                     
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION BCOF(NDIM,LM),  BINT(7,3,3),  F(2),  BFAC(2),  LL(2),    &
         FI(3),  B(7),  BDIF(3)                                         
      DATA EPS, EPS1, EPS4 /1.E-3, 1.E-13, 1.E-1/                       
!CC EPS  MAX PRECISION LOST IN BESSEL FUNCTION EVALUATION BY RECURSION  
!CC EPS1 RELATIVE ACCURACY DESIRED FOR BESSEL FN EVALUATION BY SERIES   
!CC EPS4 CRITERION WHEN TO USE STANDARD SIMPSON INTEGRATION, NOT FILON  
!CC        IF H*P .LT. EPS4, SIMPSON IS USED.                           
      L1 = LM-1                                                         
      L2 = LM-2                                                         
      D = H*P                                                           
      DFLM = DFAC(L1)                                                   
      IF (LM.GT.1) DFLM1 = DFAC(L2)                                     
      XEPS = 0.0                                                        
      IF (L1 .GT. 0) XEPS = (EPS*DFLM)**(1./FLOAT(L1))                  
!CC         FOR ARG .LT. XEPS,  USE SERIES EVALUATION FOR BESSEL FN.    
!CC                                                                     
!CC DECIDE WHETHER TO USE FILON OR SIMPSON                              
      IF (ABS(D).LT.EPS4) GOTO 300                                      
!CC                                                                     
!CC                                                                     
!CC                                                                     
!CC USE EXTENDED FILON METHOD                                           
!CC                                                                     
      D2 = D+D                                                          
      DD = D*D                                                          
      XFAC = 1./(P*DD)                                                  
      X = R*P                                                           
      N = (NPOINT+1)/2                                                  
!CC                                                                     
!CC                                                                     
      DO 5 IN=1, N                                                      
      XX = X*X                                                          
      M = MOD(IN,3)+1                                                   
      M1 = MOD(IN-1,3)+1                                                
      M2 = MOD(IN-2,3)+1                                                
!CC                                                                     
!CC DECIDE WHETHER TO USE RECURSION UP OR DOWN                          
      IF (X .EQ. 0.0) GO TO 8                                           
      IF (LM .EQ. 1 .OR. ABS(X) .GE. XEPS) GO TO 100                    
!CC                                                                     
!CC SMALL ARG CASE, USE DOWNWARD RECURSION                              
      IF (ABS(X).GT.EPS1) GOTO 11                                       
!CC ARG TOO SMALL, TREAT AS ZERO                                        
    8 I2=LM*3                                                           
      DO 9 I=1,I2                                                       
    9 BINT(I,1,M)=0.0                                                   
      GOTO 200                                                          
!CC                                                                     
!CC DOWNWARD RECURSION. USE SERIES EXPANSION FOR STARTING VALUES        
   11 DO 10 I=1, 2                                                      
      LL(I) = LM-I                                                      
      F(I) = 1.                                                         
   10 BFAC(I) = 1.                                                      
      DO 15 I2=1, 3                                                     
   15 FI(I2) = 1./(LL(1)+I2)                                            
!CC SERIES EVALUATION                                                   
      IND = 0                                                           
   20 IND = IND+1                                                       
      DO 30 I=1, 2                                                      
      BFAC(I) = -BFAC(I)/(4.*IND*(LL(I)+IND+.5))*XX                     
   30 F(I) = F(I)+BFAC(I)                                               
      DO 35 I2=1, 3                                                     
   35 FI(I2) = FI(I2)+BFAC(1)/(LL(1)+I2+IND+IND)                        
      IF (ABS(BFAC(2)).GT.EPS1) GOTO 20                                 
!CC                                                                     
      XEX = X**L2 * XFAC                                                
      B(LM) = XEX*X*F(1)/DFLM                                           
      B(L1) = XEX  *F(2)/DFLM1                                          
      DO 60 I2=1, 3                                                     
   60 BINT(LM,I2,M) = XEX*X**(I2+1)*FI(I2)/DFLM                         
      BINT(L1,1,M) = (BINT(LM,2,M)+X*B(L1))/L1                          
      BINT(L1,2,M) = L1*BINT(LM,1,M)+X*B(LM)                            
      BINT(L1,3,M) = L2*BINT(LM,2,M)+XX*B(LM)                           
      IF (LM.LE.2) GOTO 200                                             
!CC                                                                     
!CC DOWNWARD RECURSION                                                  
      LI = L1                                                           
      DO 40 L=1, L2                                                     
      LI = LI-1                                                         
      B(LI) = (LI+LI+1)/X*B(LI+1)-B(LI+2)                               
      BINT(LI,1,M) = ((LI+1)*BINT(LI+2,1,M)+(LI+LI+1)                    &
       *B(LI+1))/LI                                                     
      BINT(LI,2,M) = LI*BINT(LI+1,1,M) + X*B(LI+1)                      
   40 BINT(LI,3,M) = (LI-1)*BINT(LI+1,2,M) + XX*B(LI+1)                 
      GOTO 200                                                          
!CC                                                                     
!CC                                                                     
!CC LARGE ARGUMENT CASE, UPWARD RECURSION                               
  100 S = SIN(X)                                                        
      C = COS(X)                                                        
      B(1) = S/X*XFAC                                                   
      B(2) = (B(1)-C*XFAC)/X                                            
      BINT(1,1,M) = SININT(X)*XFAC                                      
      BINT(1,2,M) = (-C+1.)*XFAC                                        
      BINT(1,3,M) =(S-X*C)*XFAC                                         
      IF (LM.EQ.1) GOTO 200                                             
      BINT(2,1,M) = -B(1)+XFAC                                          
      BINT(2,2,M) = BINT(1,1,M)-S*XFAC                                  
      BINT(2,3,M) =(-2.*C-X*S+2.)*XFAC                                  
      IF (LM.EQ.2) GOTO 200                                             
!CC                                                                     
      DO 110 L=3, LM                                                    
      B(L) = (L+L-3)/X*B(L-1)-B(L-2)                                    
      BINT(L,1,M) = ((L-2)*BINT(L-2,1,M) - (L+L-3)*B(L-1))/(L-1)        
      BINT(L,2,M) = (L-1)*BINT(L-1,1,M) - X*B(L-1)                      
  110 BINT(L,3,M) = L*BINT(L-1,2,M) - XX*B(L-1)                         
!CC                                                                     
!CC                                                                     
!CC CONSTRUCT BCOF ARRAY                                                
  200 CONTINUE                                                          
      IF (IN.EQ.1) GOTO 5                                               
      X1 = X-D                                                          
      X2 = X-D2                                                         
      DO 230 L=1, LM                                                    
      DO 210 I2=1, 3                                                    
  210 BDIF(I2) = BINT(L,I2,M)-BINT(L,I2,M1)                             
      BCOF(2*IN-2,L) = -BDIF(3) + 2.*X1*BDIF(2)                          &
       - (X1*X1-DD)*BDIF(1)                                             
      IF (IN.GT.2) GOTO 240                                             
      BCOF(1,L) = .5*BDIF(3) - (X2+1.5*D)*BDIF(2)                        &
       + ((.5*X2+1.5*D)*X2+DD)*BDIF(1)                                  
      GOTO 230                                                          
  240 BCOF(2*IN-3,L) = .5*(BINT(L,3,M)-BINT(L,3,M2))                     &
       -X2*(BINT(L,2,M)-BINT(L,2,M2))                                    &
       -1.5*D*(BINT(L,2,M)-2.*BINT(L,2,M1)+BINT(L,2,M2))                 &
       +(.5*X2*X2+DD)*(BINT(L,1,M)-BINT(L,1,M2))                         &
       +1.5*X2*D*(BINT(L,1,M)-2.*BINT(L,1,M1)+BINT(L,1,M2))             
      IF (IN.EQ.N) BCOF(NPOINT,L) = .5*BDIF(3) - (X-1.5*D)               &
       *BDIF(2) + ((.5*X-1.5*D)*X+DD)*BDIF(1)                           
  230 CONTINUE                                                          
!CC                                                                     
!CC                                                                     
    5 X = X+D2                                                          
      RETURN                                                            
!CC                                                                     
!CC                                                                     
!CC                                                                     
!CC STANDARD SIMPSON METHOD, CASE H*P .LT. EPS4                         
!CC                                                                     
  300 CONTINUE                                                          
      XFAC = H/3.                                                       
      X = P*R                                                           
      IVF = 4                                                           
      DO 340 IN=1, NPOINT                                               
!CC                                                                     
!CC DECIDE WHETHER TO USE RECURSION UP OR DOWN                          
      IF (X .EQ. 0.0) GO TO 301                                         
      IF (LM .EQ. 1 .OR. ABS(X) .GE. XEPS) GO TO 400                    
!CC                                                                     
!CC SMALL ARGUMENT CASE, DOWNWARD RECURSION                             
  301 IVF = 6-IVF                                                       
      IF (ABS(X).GT.EPS1) GOTO 360                                      
!CC ARGUMENT TOO SMALL, TREAT AS ZERO                                   
      BCOF(IN,1) = IVF   *  XFAC                                        
      DO 370 I1=2, LM                                                   
  370 BCOF(IN,I1) = 0.                                                  
      GOTO 340                                                          
!CC                                                                     
  360 XX = X*X                                                          
      DO 310 I=1, 2                                                     
      LL(I) = LM-I                                                      
      F(I) = 1.                                                         
  310 BFAC(I) = 1.                                                      
!CC                                                                     
      IND = 0                                                           
  320 IND = IND+1                                                       
      DO 330 I=1, 2                                                     
      BFAC(I) = -BFAC(I)/(4.*IND*(LL(I)+IND+.5))*XX                     
  330 F(I) = F(I)+BFAC(I)                                               
      IF (ABS(BFAC(2)).GT.EPS1) GOTO 320                                
!CC                                                                     
      XEX = X**L2*XFAC                                                  
      BCOF(IN,LM) = XEX*X*F(1)*IVF/DFLM                                 
      BCOF(IN,L1) = XEX  *F(2)*IVF/DFLM1                                
      IF (LM.LE.2) GOTO 340                                             
!CC                                                                     
      LI = L1                                                           
      DO 345 L=1, L2                                                    
      LI = LI-1                                                         
  345 BCOF(IN,LI) = (LI+LI+1)/X*BCOF(IN,LI+1)-BCOF(IN,LI+2)             
!CC                                                                     
  340 X = X+D                                                           
      GOTO 410                                                          
!CC                                                                     
!CC                                                                     
!CC LARGE ARGUMENT CASE, UPWARD RECURSION                               
  400 N1= IN                                                            
      DO 351 IN=N1, NPOINT                                              
      IVF = 6-IVF                                                       
      S = SIN(X)*XFAC                                                   
      C = COS(X)*XFAC                                                   
      BCOF(IN,1) = S/X*IVF                                              
      IF (LM .EQ. 1) GO TO 351                                          
      BCOF(IN,2) = (BCOF(IN,1)-C*IVF)/X                                 
      IF (LM .EQ. 2) GO TO 351                                          
      DO 350 L=3, LM                                                    
  350 BCOF(IN,L) = (L+L-3)/X*BCOF(IN,L-1)-BCOF(IN,L-2)                  
  351 X = X+D                                                           
!CC                                                                     
  410 DO 420 L=1, LM                                                    
      BCOF(1,L) = .5*BCOF(1,L)                                          
  420 BCOF(NPOINT,L) = .5*BCOF(NPOINT,L)                                
      RETURN                                                            
      END                                                               
