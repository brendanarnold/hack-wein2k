      FUNCTION AINTER(R,P,RS)                                           
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(4),P(4)                                               
      N=4                                                               
      PS=0.                                                             
      DO 1 J=1,N                                                        
      TERM=1.                                                           
      DENOM=1.                                                          
      DO 2 I=1,N                                                        
      IF(I.EQ.J) GOTO 2                                                 
      DENOM=DENOM*(R(J)-R(I))                                           
      TERM=TERM*(RS-R(I))                                               
  2   CONTINUE                                                          
  1   PS=PS+TERM*P(J)/DENOM                                             
       AINTER=PS                                                        
      RETURN                                                            
      END                                                               
      SUBROUTINE ATPAR(NAT)                                             
!-----X----X----X----X----X----X----X----X----X----X----X----X----X----X
!-----X LINEAR EXPANDED ORBITAL APW ROUTINE--------               X----X
!-----X SET UP THE RADIAL INTEGRAL PARAMETERS.                    X----X
!-----X SPIN-ORBITLESS RELATIVISTIC EQUATIONS USED.               X----X
!-----X L IS 7  (NT)                                              X----X
!-----X REQUIRED SUBROUTINES ARE OUTWIN AND RINT13.               X----X
!-----X SET UP FOR 4 ATOMS, AND CALL 4 TIMES                      X----X
!-----X                                      P.BLAHA              X----X
!-----X----X----X----X----X----X----X----X----X----X----X----X----X----X
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      CHARACTER*10 ANAME                                                
      COMMON /STRUK/ POS(ndif,3),AA,BB,CC,alpha(3), &
       RMT(nato),V(nato),PIA(3),VOL,IATNR(nato)   &
      ,MULT(nato),ISPLIT(nato)   
      COMMON /POTNLC/ VR(885),RNOT(nato),DX(nato),JRI(nato)             
      COMMON /CUBWK1/ A(885),B(885),AP(885),BP(885),AE(885),BE(885)     
      COMMON /ROTMAT/   ROTLOC(3,3,NATO),ROTIJ(3,3,NDIF),TAUIJ(3,NDIF)  
    1 FORMAT(A10,I5,5X,2E10.5,I5,F5.2)                                  
    2 FORMAT(5E14.7)                                                    
   12 FORMAT(8I2)                                                       
   13 FORMAT(/,' AT.NR.:',I2,5X,A10,5X,'POSITION:.',3F8.3,10X,         &
      'MULTIPLICITY: ',I5)                                              
!                                                                       
!     INDEX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE    
!     NOTEQUIVALENT ATOMS                                               
      INDEX=0                                                           
      DO 50 JATOM = 1,NAT                                               
         INDEX=INDEX+1                                                  
         READ(20,1030) IATNR(JATOM),( POS(INDEX,j),J=1,3 ),MULT(JATOM),  &
                       ISPLIT(JATOM)                                    
         IF (MULT(JATOM).EQ.0) THEN                                     
            WRITE(66,1040) JATOM,INDEX,MULT(JATOM)                      
            STOP ' LAPW2: MULT EQ 0'                                    
         ENDIF                                                          
            DO 55 M=1,MULT(JATOM)-1                                     
               INDEX=INDEX+1                                            
               READ(20,1031) IATNR(JATOM),( POS(INDEX,j),J=1,3)         
 55         CONTINUE                                                    
         READ(20,1050) ANAME,JRI(JATOM),Rnot(JATOM),RMT(JATOM)     
         DX(JATOM)=LOG(RMT(JATOM)/Rnot(JATOM)) / (JRI(JATOM)-1)           
         READ(20,1051) ((ROTLOC(I,J,JATOM),I=1,3),J=1,3)                
 1030 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
 1031 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1040 FORMAT(///,3X,'ERROR IN LAPW2 : MULT(JATOM)=0 ...',                &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)           
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
 1051    FORMAT(20X,3F10.8)                                             
 50   CONTINUE                                                          
!                                                                       
   51 CONTINUE                                                          
!                                                                       
      RETURN                                                            
      END                                                               

!     .....................................................GBASS........
      SUBROUTINE DETER(GMAX,RBAS,kxmax,kymax,kzmax)                                       
!     **                                                              **
!     **  CALCULATE RECIPROCAL LATTICE VECTORS FROM REAL SPACE        **
!     **  LATTICE VECTORS OR VICE VERSA                               **
!     **                                                              **
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION RBAS(3,3),GBAS(3,3)                                     
      GBAS(1,1)=RBAS(2,2)*RBAS(3,3)-RBAS(3,2)*RBAS(2,3)                 
      GBAS(2,1)=RBAS(3,2)*RBAS(1,3)-RBAS(1,2)*RBAS(3,3)                 
      GBAS(3,1)=RBAS(1,2)*RBAS(2,3)-RBAS(2,2)*RBAS(1,3)                 
      GBAS(1,2)=RBAS(2,3)*RBAS(3,1)-RBAS(3,3)*RBAS(2,1)                 
      GBAS(2,2)=RBAS(3,3)*RBAS(1,1)-RBAS(1,3)*RBAS(3,1)                 
      GBAS(3,2)=RBAS(1,3)*RBAS(2,1)-RBAS(2,3)*RBAS(1,1)                 
      GBAS(1,3)=RBAS(2,1)*RBAS(3,2)-RBAS(3,1)*RBAS(2,2)                 
      GBAS(2,3)=RBAS(3,1)*RBAS(1,2)-RBAS(1,1)*RBAS(3,2)                 
      GBAS(3,3)=RBAS(1,1)*RBAS(2,2)-RBAS(2,1)*RBAS(1,2)                 
      DET=0.D0                                                          
      DO 100 I=1,3                                                      
      DET=DET+GBAS(I,1)*RBAS(I,1)                                       
100   CONTINUE                                                          
      kxmax=DABS(GBAS(1,1)/DET*GMAX)
      kymax=DABS(GBAS(1,2)/DET*GMAX)
      kzmax=DABS(GBAS(1,3)/DET*GMAX)
      xtest=DABS(GBAS(2,1)/DET*GMAX)
      if (xtest.gt.kxmax) kxmax= xtest
      xtest=DABS(GBAS(3,1)/DET*GMAX)
      if (xtest.gt.kxmax) kxmax= xtest
      ytest=DABS(GBAS(2,2)/DET*GMAX)
      if (ytest.gt.kymax) kymax= ytest
      ytest=DABS(GBAS(3,2)/DET*GMAX)
      if (ytest.gt.kymax) kymax= ytest
      ztest=DABS(GBAS(2,3)/DET*GMAX)
      if (ztest.gt.kzmax) kzmax= ztest
      ztest=DABS(GBAS(3,3)/DET*GMAX)
      if (ztest.gt.kzmax) kzmax= ztest
      RETURN                                                            
      END                                                               
      FUNCTION DFAC(L)                                                  
!CC CALCULATES DOUBLE FACTORIALS OF ARGUMENT 2*L+1                      
      IMPLICIT REAL*8 (A-H,O-Z)
      DFAC = 1.                                                         
      IF (L.EQ.0) RETURN                                                
      X = 1.                                                            
      DO 10 I=1, L                                                      
      X = X+2.                                                          
   10 DFAC = DFAC*X                                                     
      RETURN                                                            
      END                                                               
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
      SUBROUTINE GBASS(RBAS,GBAS)                                       
!     **                                                              **
!     **  CALCULATE RECIPROCAL LATTICE VECTORS FROM REAL SPACE        **
!     **  LATTICE VECTORS OR VICE VERSA (without 2 Pi factor !!!      **
!     **                                                              **
      INTEGER i,j
      DOUBLE PRECISION PI, DET
      DOUBLE PRECISION RBAS(3,3),GBAS(3,3) ,help(3,3)                                 
      PI=4.D0*DATAN(1.D0)                                               
      GBAS(1,1)=RBAS(2,2)*RBAS(3,3)-RBAS(3,2)*RBAS(2,3)                 
      GBAS(2,1)=RBAS(3,2)*RBAS(1,3)-RBAS(1,2)*RBAS(3,3)                 
      GBAS(3,1)=RBAS(1,2)*RBAS(2,3)-RBAS(2,2)*RBAS(1,3)                 
      GBAS(1,2)=RBAS(2,3)*RBAS(3,1)-RBAS(3,3)*RBAS(2,1)                 
      GBAS(2,2)=RBAS(3,3)*RBAS(1,1)-RBAS(1,3)*RBAS(3,1)                 
      GBAS(3,2)=RBAS(1,3)*RBAS(2,1)-RBAS(2,3)*RBAS(1,1)                 
      GBAS(1,3)=RBAS(2,1)*RBAS(3,2)-RBAS(3,1)*RBAS(2,2)                 
      GBAS(2,3)=RBAS(3,1)*RBAS(1,2)-RBAS(1,1)*RBAS(3,2)                 
      GBAS(3,3)=RBAS(1,1)*RBAS(2,2)-RBAS(2,1)*RBAS(1,2)                 
      DET=0.D0                                                          
      DO 100 I=1,3                                                      
      DET=DET+GBAS(I,1)*RBAS(I,1)                                       
100   CONTINUE                                                          
      DO 110 I=1,3                                                      
      DO 110 J=1,3     
 110     help(i,j)=gbas(i,j)/det                                                 
!      GBAS(I,J)=GBAS(I,J)*2.D0*PI/DET                                   
      DO 111 I=1,3                                                      
      DO 111 J=1,3     
      GBAS(I,J)=help(j,i)                                   
111   CONTINUE                                                          
      RETURN                                                            
      END                                                               
      LOGICAL FUNCTION KDELTA(K)                                        
!                                                                       
!.... TEST, IF K IS IN STAR OF G  (GENERATED IN STERN)                  
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /FOUR/ G(3),NST,STG(3,48)                                  
      DIMENSION K(3)                                                    
      INTEGER G,STG                                                     
      DO 1 I=1,NST                                                      
      DO 2 J=1,3                                                        
      IF(STG(J,I).NE.K(J)) GOTO 1                                       
   2  CONTINUE                                                          
      KDELTA=.TRUE.                                                     
      RETURN                                                            
  1   CONTINUE                                                          
      KDELTA=.FALSE.                                                    
      RETURN                                                            
      END                                                               
!$hp9000_800 intrinsics on
      PROGRAM LAPW3                                                     
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
!                                                                       
!                                                                       
      CHARACTER*80     FNAME                                      
      CHARACTER*11     STATUS,FORM                                      
      CHARACTER*4 LATTIC,ISYM,IREL,cform                                
      CHARACTER*10 TITLE(8)                                             
      CHARACTER*10     ANAME                                     
      LOGICAL REL                                                       
      COMMON /CHAR/   TITLE,LATTIC                   
      COMMON /COM/ XWT,EMIN,EF,ELECN,REL,MODUS,NBAND,NAT                
      COMMON /STRUK/ POS(ndif,3),AA,BB,CC,alpha(3), &
       RMT(nato),V(nato),PIA(3),VOL,IATNR(nato)   &
      ,MULT(nato),ISPLIT(nato)   
      COMMON /GENER/ GX1,GY1,GZ1,GX2,GY2,GZ2,GX3,GY3,GZ3,BR2(3,3)       
      COMMON /CUBWK1/ DUM(999,ncom)                                     
      COMMON /POTNLC/ DUMMY2(885),R0(nato),DX(nato),JRI(nato)           
      DATA PI/3.1415926535898d0/                                        
!
      call getarg(2,fname)
      if(fname.eq.'      ') call getarg(1,fname)
      OPEN(1,FILE=fname,STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM, &
      ERR=8002)
      GOTO 8003
 8000 WRITE(*,*) ' ERROR IN OPENING LAPW3.DEF !!!!'
      STOP 'LAPW2.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, &
      '  FORM:',FORM
      STOP 'OPEN FAILED'
 8001 CONTINUE
      CLOSE (1)
!
      DO 9575 IND=1,nato                                                
      RMT(IND)=0.0                                                      
      V(IND)=0.0                                                        
 9575 MULT(IND)=0                                                       
      SQRT3=SQRT(3.)                                                    
      READ(20,1000) TITLE                                               
      READ(20,1001) LATTIC,NAT,cform,IREL                               
      WRITE(66,1004) TITLE,LATTIC,IREL                                  
      WRITE(66,1006) NAT                                                
        assign 2021 to iform1
        assign 2071 to iform2
 2021 FORMAT(3X,4E19.12)                                                
 2071 FORMAT(3X,3I5,2E19.12)                                            
!                                                                       
      REL=.FALSE.                                                       
      IF(IREL.EQ.'RELA') REL=.TRUE.                                     
      READ(20,1007) AA,BB,CC,alpha                                            
      WRITE(66,1008) AA,BB,CC,alpha                                           
!
      call latgen(nat)
! 
      CALL ATPAR(NAT)                                                   
!
        CALL FOURIR                                              
      STOP                                                              
 1000 FORMAT(8A10)                                                      
 1001 FORMAT(A4,24X,I2,1x,a4,/,13X,A4,18X,A4)                           
 1004 FORMAT(8A10,//,10X,A4,' LATTIC ASSUMED',//,10X,2A4,            &
      '-CALCULATION',//)                                                
 1006 FORMAT(10X,'NUMBER OF ATOMS:',I3)                                 
 1007 FORMAT(6F10.7)                                                    
 1008 FORMAT('  LATTICE CONSTANTS  ',6F10.5)                        
      END                                                               
      SUBROUTINE LATGEN (NAT)                                           
!                                                                       
!     LATGEN GENERATES TWO BRAVAIS MATRICES, DEFINES THE VOLUME OF      
!     THE UNIT CELL AND CALLS ROTDEF                                    
!     BR1(3,3)  : TRANSFORMS INTEGER RECIPROCAL LATTICE VECTORS AS      
!                 GIVEN IN THE VECTORLIST OF LAPW1  ( GENERATED IN      
!                 COORS, TRANSFORMED IN BASISO, AND WRITTEN OUT IN      
!                 WFTAPE) INTO CARTESIAN SYSTEM                         
!     BR2(3,3) :  TRANSFORMS A RECIPROCAL LATTICE VECTOR OF A SPE-      
!                 CIAL COORDINATE SYSTEM ( IN UNITS OF 2 PI / A )       
!                 TO CARTESIAN SYSTEM                                   
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
      LOGICAL          ORTHO           
      CHARACTER*4      LATTIC
      CHARACTER*80     TITLE    
!                           
      COMMON /ORTH/   ORTHO                
      COMMON /STRUK/  POS(3,NDIF),AA,BB,CC,ALPHA(3),RMT(NATO),V(NATO), &
                      PIA(3),VOL,IATNR(NATO),MULT(NATO),ISPLIT(NATO)           
      COMMON /CHAR/   TITLE,LATTIC                    
      COMMON /GENER/  BR1(3,3),BR2(3,3)                                 
!---------------------------------------------------------------------  
!                      
      PI=ACOS(-1.0D0)                                                   
      SQRT3=SQRT(3.D0)
      ALPHA(1)=ALPHA(1)*PI/180.0D0                                             
      ALPHA(2)=ALPHA(2)*PI/180.0D0                                             
      ALPHA(3)=ALPHA(3)*PI/180.0D0                                             
      PIA(1)=2.D0*PI/AA                                                 
      PIA(2)=2.D0*PI/BB                                                 
      PIA(3)=2.D0*PI/CC                                                 
      IF(LATTIC(1:1).EQ.'H') GOTO 10                                    
      IF(LATTIC(1:1).EQ.'S') GOTO 20                                    
      IF(LATTIC(1:1).EQ.'P') GOTO 20                                    
      IF(LATTIC(1:1).EQ.'F') GOTO 30                                    
      IF(LATTIC(1:1).EQ.'B') GOTO 40                                    
      IF(LATTIC(1:1).EQ.'C') GOTO 50                                    
      IF(LATTIC(1:1).EQ.'R') GOTO 60                                    
!
!        Error: wrong lattice, stop execution
!
         GOTO 900
!                                                                       
!.....HEXAGONAL LATTICE                                                 
 10   CONTINUE                                                          
      BR1(1,1)=2.D0/SQRT3*PIA(1)                                        
      BR1(1,2)=1.D0/SQRT3*PIA(1)                                        
      BR1(1,3)=0.0D0                                                    
      BR1(2,1)=0.0D0                                                    
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=0.0D0                                                    
      BR1(3,1)=0.0D0                                                    
      BR1(3,2)=0.0D0                                                    
      BR1(3,3)=PIA(3)                                                   
!                                                                       
      BR2(1,1)=2.D0/SQRT3*PIA(1)                                        
      BR2(1,2)=1.D0/SQRT3*PIA(1)                                        
      BR2(1,3)=0.0D0                                                    
      BR2(2,1)=0.0D0                                                    
      BR2(2,2)=PIA(2)                                                   
      BR2(2,3)=0.0D0                                                    
      BR2(3,1)=0.0D0                                                    
      BR2(3,2)=0.0D0                                                    
      BR2(3,3)=PIA(3)                                                   
!                                                                       
      RVFAC=2.D0/SQRT(3.D0)                                             
      ORTHO=.FALSE.                                             
      GOTO 100                                                          
!                                                                       
!.....RHOMBOHEDRAL CASE                                                    
 60   BR1(1,1)=1.D0/SQRT(3.D0)*PIA(1)                                          
      BR1(1,2)=1.D0/SQRT(3.D0)*PIA(1)                                          
      BR1(1,3)=-2.d0/sqrt(3.d0)*PIA(1)                                         
      BR1(2,1)=-1.0d0*PIA(2)                                                  
      BR1(2,2)=1.0d0*PIA(2)                                                    
      BR1(2,3)=0.0d0*PIA(2)                                                    
      BR1(3,1)=1.0d0*PIA(3)                                                    
      BR1(3,2)=1.0d0*PIA(3)                                                    
      BR1(3,3)=1.0d0*PIA(3)                                                    
!
      BR2(1,1)=1.D0/SQRT(3.D0)*PIA(1)                                          
      BR2(1,2)=1.D0/SQRT(3.D0)*PIA(1)                                          
      BR2(1,3)=-2.d0/sqrt(3.d0)*PIA(1)                                         
      BR2(2,1)=-1.0d0*PIA(2)                                                   
      BR2(2,2)=1.0d0*PIA(2)                                                    
      BR2(2,3)=0.0d0*PIA(2)                                                    
      BR2(3,1)=1.0d0*PIA(3)                                                    
      BR2(3,2)=1.0d0*PIA(3)                                                    
      BR2(3,3)=1.0d0*PIA(3)                                                    
      RVFAC=6.D0/SQRT(3.D0)
      ORTHO=.FALSE.                                             
      GOTO 100                                                          
!                                                                       
!.....PRIMITIVE LATTICE                                                 
 20   IF(ABS(ALPHA(3)-PI/2.0D0).GT.0.0001) then
      write(*,*) '  gamma not equal 90'
      GOTO 21     
      endif
      IF(ABS(ALPHA(2)-PI/2.0D0).GT.0.0001) then
      write(*,*) '  beta not equal 90'
      GOTO 22     
      endif
      BR1(1,1)=PIA(1)                                                   
      BR1(1,2)=0.0D0                                                    
      BR1(1,3)=0.0D0                                                    
      BR1(2,1)=0.0D0                                                    
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=0.0D0                                                    
      BR1(3,1)=0.0D0                                                    
      BR1(3,2)=0.0D0                                                    
      BR1(3,3)=PIA(3)                                                   
!                                                                       
      BR2(1,1)=PIA(1)                                                   
      BR2(1,2)=0.0D0                                                    
      BR2(1,3)=0.0D0                                                    
      BR2(2,1)=0.0D0                                                    
      BR2(2,2)=PIA(2)                                                   
      BR2(2,3)=0.0D0                                                    
      BR2(3,1)=0.0D0                                                    
      BR2(3,2)=0.0D0                                                    
      BR2(3,3)=PIA(3)                                                   
!                                                                       
      RVFAC=1.D0                                                        
      ORTHO=.TRUE.                                             
      GOTO 100                                                          
!                                                                       
!.....FC LATTICE                                                        
 30   CONTINUE                                                          
      BR1(1,1)=PIA(1)                                                   
      BR1(1,2)=0.0D0                                                    
      BR1(1,3)=0.0D0                                                    
      BR1(2,1)=0.0D0                                                    
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=0.0D0                                                    
      BR1(3,2)=0.0D0                                                    
      BR1(3,1)=0.0D0                                                    
      BR1(3,3)=PIA(3)                                                   
!                    
!     definitions according to column, rows convention for BR2
!
      BR2(1,1)=-PIA(1)                                                  
      BR2(1,2)= PIA(1)                                                  
      BR2(1,3)= PIA(1)                                                  
      BR2(2,1)= PIA(2)                                                  
      BR2(2,2)=-PIA(2)                                                  
      BR2(2,3)= PIA(2)                                                  
      BR2(3,1)= PIA(3)                                                  
      BR2(3,2)= PIA(3)                                                  
      BR2(3,3)=-PIA(3)                                                  
!                                                                       
      RVFAC=4.D0                                                        
      ORTHO=.TRUE.                                             
      GOTO 100                                                          
!                                                                       
!.....BC LATTICE                                                        
 40   CONTINUE                                                          
      BR1(1,1)=PIA(1)                                                   
      BR1(1,2)=0.0D0                                                    
      BR1(1,3)=0.0D0                                                    
      BR1(2,1)=0.0D0                                                    
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=0.0D0                                                    
      BR1(3,1)=0.0D0                                                    
      BR1(3,2)=0.0D0                                                    
      BR1(3,3)=PIA(3)                                                   
!                                                                       
      BR2(1,1)= 0.0D0                                                    
      BR2(1,2)= PIA(1)                                                  
      BR2(1,3)= PIA(1)                                                  
      BR2(2,1)= PIA(2)                                                  
      BR2(2,2)= 0.0D0                                                    
      BR2(2,3)= PIA(2)                                                  
      BR2(3,1)= PIA(3)                                                  
      BR2(3,2)= PIA(3)                                                  
      BR2(3,3)= 0.0D0
!
      RVFAC=2.D0
      ORTHO=.TRUE.                                             
      GOTO 100                                                    
!             
 50   CONTINUE                                                          
      IF(LATTIC(2:3).EQ.'XZ') GOTO 51                                    
      IF(LATTIC(2:3).EQ.'YZ') GOTO 52                                    
!.....CXY LATTICE                                                          
      BR1(1,1)=PIA(1)                                                   
      BR1(1,2)=0.0D0                                                    
      BR1(1,3)=0.0D0                                                    
      BR1(2,1)=0.0D0                                                    
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=0.0D0                                                    
      BR1(3,1)=0.0D0                                                    
      BR1(3,2)=0.0D0                                                    
      BR1(3,3)=PIA(3)                                                   
!                                                                       
      BR2(1,1)= PIA(1)                                                    
      BR2(1,2)= PIA(1)                                                  
      BR2(1,3)= 0.0D0                                                  
      BR2(2,1)=-PIA(2)                                                  
      BR2(2,2)= PIA(2)                                                    
      BR2(2,3)= 0.0D0                                                 
      BR2(3,1)= 0.0D0                                                  
      BR2(3,2)= 0.0D0                                                 
      BR2(3,3)= PIA(3)                                                    
!                                                                       
      RVFAC=2.D0                                                        
      ORTHO=.TRUE.                                             
      GOTO 100                                                          
!                                                                       
!.....CXZ CASE (CXZ LATTICE BUILD UP)                                     
 51   CONTINUE                                     
!.....CXZ ORTHOROMBIC CASE 
      IF(ABS(ALPHA(3)-PI/2.0D0).LT.0.0001) then
         BR1(1,1)=PIA(1)                                                   
         BR1(1,2)=0.0D0                                                    
         BR1(1,3)=0.0D0                                                    
         BR1(2,1)=0.0D0                                                    
         BR1(2,2)=PIA(2)                                                   
         BR1(2,3)=0.0D0                                                    
         BR1(3,1)=0.0D0                                                    
         BR1(3,2)=0.0D0                                                    
         BR1(3,3)=PIA(3)                                                   
!                                                                       
         BR2(1,1)= PIA(1)                                                   
         BR2(1,2)= 0.0                                                   
         BR2(1,3)= PIA(1)                                                      
         BR2(2,1)= 0.0                                                      
         BR2(2,2)= PIA(2)                                                     
         BR2(2,3)= 0.0                                                     
         BR2(3,1)=-PIA(3)                                                     
         BR2(3,2)= 0.0                                                     
         BR2(3,3)= PIA(3)                                                     
!                                                                       
         RVFAC=2.0                                                         
         ORTHO=.TRUE.                                             
         GOTO 100                                                          
      ELSE
!.....CXZ MONOCLINIC CASE 
         write(*,*) '  gamma not equal 90'
         SINAB=SIN(ALPHA(3))
         COSAB=COS(ALPHA(3))
!                                                                       
         BR1(1,1)= PIA(1)/SINAB 
         BR1(1,2)= -PIA(1)*COSAB/SINAB
         BR1(1,3)= PIA(1)/SINAB 
         BR1(2,1)= 0.0                                                      
         BR1(2,2)= PIA(2)                                                     
         BR1(2,3)= 0.0                                                     
         BR1(3,1)=-PIA(3)                                                     
         BR1(3,2)= 0.0                                                     
         BR1(3,3)= PIA(3)                                                     
!                                                                       
         BR2(1,1)= PIA(1)/SINAB 
         BR2(1,2)= -PIA(1)*COSAB/SINAB
         BR2(1,3)= PIA(1)/SINAB 
         BR2(2,1)= 0.0                                                      
         BR2(2,2)= PIA(2)                                                     
         BR2(2,3)= 0.0                                                     
         BR2(3,1)=-PIA(3)                                                     
         BR2(3,2)= 0.0                                                     
         BR2(3,3)= PIA(3)                                                     
!                                                                       
         RVFAC=2.0/SINAB                                                   
         ORTHO=.FALSE.                                             
         GOTO 100                                                          
      ENDIF
!                                                                       
!.....CYZ CASE (CYZ LATTICE BUILD UP)                                     
 52   CONTINUE                                     
      BR1(1,1)=PIA(1)                                                   
      BR1(1,2)=0.0D0                                                    
      BR1(1,3)=0.0D0                                                    
      BR1(2,1)=0.0D0                                                    
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=0.0D0                                                    
      BR1(3,1)=0.0D0                                                    
      BR1(3,2)=0.0D0                                                    
      BR1(3,3)=PIA(3)                                                   
!                                                                       
      BR2(1,1)= PIA(1)                                                      
      BR2(1,2)= 0.0                                                   
      BR2(1,3)= 0.0                                                      
      BR2(2,1)= 0.0                                                      
      BR2(2,2)= PIA(2)                                                     
      BR2(2,3)= PIA(2)                                                     
      BR2(3,1)= 0.0                                                     
      BR2(3,2)=-PIA(3)                                                     
      BR2(3,3)= PIA(3)                                                     
!                                                                       
      RVFAC=2.0                                                         
      ORTHO=.TRUE.                                             
      GOTO 100                                                          
!                                                                       
!.....M CASE (MONOCLINIC LATTICE BUILD UP)
 21   CONTINUE
!
!     IF(ABS(ALPHA(2)-PI/2.).GT.0.0001) GOTO 61
!     GAMMA NE 90                                     
      SINAB=SIN(ALPHA(3))
      COSAB=COS(ALPHA(3))
      BR1(1,1)= PIA(1)/SINAB                                                      
      BR1(1,2)= -PIA(2)*COSAB/SINAB                                                     
      BR1(1,3)= 0.0                                                   
      BR1(2,1)= 0.0                                                      
      BR1(2,2)= PIA(2)                                                     
      BR1(2,3)= 0.0                                                     
      BR1(3,1)= 0.0                                                     
      BR1(3,2)= 0.0                                                     
      BR1(3,3)= PIA(3)                                                     
!
      BR2(1,1)= PIA(1)/SINAB                                                      
      BR2(1,2)= -PIA(2)*COSAB/SINAB                                                     
      BR2(1,3)= 0.0                                                   
      BR2(2,1)= 0.0                                                      
      BR2(2,2)= PIA(2)                                                     
      BR2(2,3)= 0.0                                                     
      BR2(3,1)= 0.0                                                     
      BR2(3,2)= 0.0                                                     
      BR2(3,3)= PIA(3)                                                     
!
      RVFAC=1.d0/SINAB                                                         
      ORTHO=.FALSE.                                             
      GOTO 100                                                          
!                                                                       
!.....M CASE (MONOCLINIC LATTICE BUILD UP)                                     
!     BETA NE 90                                     
 22   CONTINUE
      SINAC=SIN(ALPHA(2))
      COSAC=COS(ALPHA(2))
      BR2(1,1)= PIA(1)/SINAC                                                      
      BR2(1,2)= 0.0                                                   
      BR2(1,3)= -PIA(3)*COSAC/SINAC                                                     
      BR2(2,1)= 0.0                                                      
      BR2(2,2)= PIA(2)                                                     
      BR2(2,3)= 0.0                                                     
      BR2(3,1)= 0.0                                                     
      BR2(3,2)= 0.0                                                     
      BR2(3,3)= PIA(3)                                                     
!
      BR1(1,1)= PIA(1)/SINAC                                                      
      BR1(1,2)= 0.0                                                   
      BR1(1,3)= -PIA(3)*COSAC/SINAC                                                     
      BR1(2,1)= 0.0                                                      
      BR1(2,2)= PIA(2)                                                     
      BR1(2,3)= 0.0                                                     
      BR1(3,1)= 0.0                                                     
      BR1(3,2)= 0.0                                                     
      BR1(3,3)= PIA(3)                                                     
!
      RVFAC=1.d0/SINAC                                                         
      ORTHO=.FALSE.                                             
      GOTO 100                                                          
!                                                                       
!                                                                       
!.....DEFINE VOLUME OF UNIT CELL                                        
 100  CONTINUE                                                          
      VOL=AA*BB*CC/RVFAC                                                
!                                                                       
!.....DEFINE ROTATION MATRICES IN NONSYMMORPHIC CASE                    
!      CALL ROTDEF (NAT)                                            
!      this is done later in strfac
!
      RETURN
!
!        Error messages
!
!  900 CALL OUTERR('LATGEN','wrong lattice.')
  900 WRITE(*,*) 'LATGEN','wrong lattice.'
      STOP 'LATGEN - Error'
!
!        End of 'LATGEN'
!
      END                                                               
      SUBROUTINE LOCDEF(RBAS,GBAS,ROTLOC)
!
!    Redefines local rotation matix from struct-file with
!    unitary transformation  u(-1) * S * u  for non-orthogonal lattice
!
      INTEGER j
      DOUBLE PRECISION rbas(3,3),gbas(3,3),ROTLOC(3,3)
      DOUBLE PRECISION B(3,3)
!      write(6,111) (rotloc(1,j),j=1,3)
!      write(6,111) (rotloc(2,j),j=1,3)
!      write(6,111) (rotloc(3,j),j=1,3)
             CALL MATMM(B,GBAS,ROTLOC)
!      write(6,111) (b(1,j),j=1,3)
!      write(6,111) (b(2,j),j=1,3)
!      write(6,111) (b(3,j),j=1,3)
             CALL MATMM(ROTLOC,B,RBAS)
!      write(6,111) (rotloc(1,j),j=1,3)
!      write(6,111) (rotloc(2,j),j=1,3)
!      write(6,111) (rotloc(3,j),j=1,3)
111   format(3f10.4)
      RETURN
      END
      SUBROUTINE  MATMM (C, A,B)                                       
!     (MULTIPLY 3 X 3 MATRICES)                                        
!     FORM C = A B, WHERE C MAY OVERLAP WITH EITHER A OR B, OR BOTH,   
!     SINCE THE PRODUCT IS DEVELOPED IN A TEMPORARY MATRIX.            
!     (06-JUL-80)                                                      
      REAL*8           A(3,3),      AB(3,3),     B(3,3),      C(3,3)   
      AB(1,1) = A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
      AB(1,2) = A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)
      AB(1,3) = A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3)
      AB(2,1) = A(2,1)*B(1,1)+A(2,2)*B(2,1)+A(2,3)*B(3,1)
      AB(2,2) = A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2)
      AB(2,3) = A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3)
      AB(3,1) = A(3,1)*B(1,1)+A(3,2)*B(2,1)+A(3,3)*B(3,1)
      AB(3,2) = A(3,1)*B(1,2)+A(3,2)*B(2,2)+A(3,3)*B(3,2)
      AB(3,3) = A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3)
      C(1,1) = AB(1,1)                                                 
      C(2,1) = AB(2,1)                                                 
      C(3,1) = AB(3,1)                                                 
      C(1,2) = AB(1,2)                                                 
      C(2,2) = AB(2,2)                                                 
      C(3,2) = AB(3,2)                                                 
      C(1,3) = AB(1,3)                                                 
      C(2,3) = AB(2,3)                                                 
      C(3,3) = AB(3,3)                  
      RETURN                                                           
      END                                                              
      SUBROUTINE RECPR(NKK,INDMAX,INST,KXMAX,KYMAX,KZMAX,GMAX)                                 
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      LOGICAL KDELTA,ORTHO
!                                                                       
      COMMON /ORTH/   ORTHO                                            
      COMMON /GENER/ BR1(3,3),BR2(3,3)                                  
      COMMON /STRUK/ POS(ndif,3),AA,BB,CC,alpha(3), &
       RMT(nato),V(nato),PIA(3),VOL,IATNR(nato)   &
      ,MULT(nato),ISPLIT(nato)   
      COMMON /GITT/ ABSK(NKD+1),NK,KZZ(3,NKD+1)                         
      COMMON /FOUR/ IM(3),NST,ISTM(3,48),TAUP(48)                       
      COMMON /SYM2/ IORD,IZ(3,3,48),ITAU,TAU(3,48)                      
      COMMON /INVERS/ TAUK(NFK),KREC(3,NFK)                            
      DIMENSION AM(3),M(3),INST(NKD)                                    
      save l,test,nk1,ll
!C    L=16                                                              
      L=32                                                              
!      L=08                                                             
      TEST=1.E-08                                                       
      NK=nkd                                                          
      LL=2*L+1                                                          
!                                                                       
!..... READ IN SYMMETRY OPERATIONS AND NONPRIMITIVE TRANSLATIONS        
!                                                                       
      READ(20,100) IORD                                                 
      ITAU=0                                                            
      DO 15 I=1,IORD                                                    
      READ(20,101) ((IZ(I1,I2,I),I1=1,3),TAU(I2,I),I2=1,3)              
      DO 16 IND=1,3                                                     
   16 IF(TAU(IND,I).NE.0.0) GOTO 7777                                   
      GOTO 15                                                           
 7777 ITAU=ITAU+1                                                       
   15 CONTINUE                                                          
!                                                                       
      ENTRY RECPR2(NKK,INDMAX,INST,KXMAX,KYMAX,KZMAX,GMAX)                 
      NK1=NK+1                                                          
      J=1                                                               
      NK=0                                                              
      LL1=KXMAX+1                                                       
      LL2=KYMAX+1                                                       
      LL3=KZMAX+1                                                       
      ILL1=-ll1+2                                                       
      ILL2=-ll2+2                                                       
      ILL3=-ll3+2                                                       
!      
      DO 1 I1=ILL1,LL1                                                 
      M(1)=I1-1                                                         
      DO 1 I2=ILL2,LL2                                                   
      M(2)=I2-1                                                         
      DO 1 I3=ILL3,LL3                                                
      M(3)=I3-1                                                         
!                                                                       
      ABSM=0.0d0                                                          
      DO 22 IND1=1,3                                                    
      AM(IND1)=0.0                                                      
      DO 23 IND2=1,3                                                    
  23  AM(IND1)=AM(IND1)+M(IND2)*BR2(IND1,IND2)                          
      ABSM=ABSM+AM(IND1)**2                                             
      AHELP=AM(IND1)/PIA(IND1)                                          
      IM(IND1)=AHELP+SIGN(0.01d0,AHELP)                                 
      IF(ORTHO) GOTO 22                                 
      IM(1)=M(1)                                                        
      IM(2)=M(2)                                                        
      IM(3)=M(3)
   22 CONTINUE                                                          
!                                                                       
       ABSM=SQRT(ABSM)                                                  
      if(absm.gt.gmax) goto 1
!
      IF(NK.EQ.0) GOTO 7                                                
      DO 2 J=1,NK                                                       
      ABST=ABSM-ABSK(J)                                                 
      IF(ABS(ABST).LT.TEST)GOTO 5                                       
      IF(ABST.LT.0.0) GOTO 4                                            
   2  CONTINUE                                                          
! .... NEW VECTOR GT. OLD MAX                                           
      J=NK                                                              
      IF(J.GE.NK1) GOTO 1                                               
      J=J+1                                                             
      GOTO 7                                                            
  5   CALL STERN                                                        
! .... NEW VECTOR EQU. OLD ONE                                          
      DO 3 I=J,NK                                                       
      IF(KDELTA(KZZ(1,I))) GOTO 1                                       
  3   CONTINUE                                                          
  4   NJ=NK-J+1                                                         
! .... NEW VECTOR .NE. OLD ONE                                          
      DO 6 JJ=1,NJ                                                      
      NJJ=NK-JJ+1                                                       
      ABSK(NJJ+1)=ABSK(NJJ)                                             
      DO 6 N=1,3                                                        
  6   KZZ(N,NJJ+1)=KZZ(N,NJJ)                                           
   7  ABSK(J)=ABSM                                                      
! .... PUT NEW VECTOR IN LIST                                           
      DO 8 N=1,3                                                        
   8  KZZ(N,J)=IM(N)                                                    
!      KZZ(3,J)=IABS(KZZ(3,J))                                          
      NK=NK+1                                                           
      IF(NK.GE.NK1) NK=NK1-1                                            
   1  CONTINUE                                                          
!
      IND=0                                                             
      DO 30 I=1,NKK
      DO 31 J=1,3                                                       
   31 IM(J)=KZZ(J,I)                                                    
      CALL STERN                                                        
      INST(I)=NST                                                       
      IF(IND+NST.GT.NFK) THEN                                          
      WRITE(66,102) IND,I                                               
 102  FORMAT('0 MORE THAN NFK K REQUESTED, TRUNCATED TO',I5,            &
      '   NKK+1=',I5,//)                                                
      GOTO 33                                                           
      END IF                                                            
      DO 32 JJ=1,NST                                                    
      IND=IND+1                                                         
      TAUK(IND)=TAUP(JJ)                                                
      DO 32 JJJ=1,3                                                     
  32  KREC(JJJ,IND)=ISTM(JJJ,JJ)                                        
   30 CONTINUE                                                          
!
 33   NKK=I-1                                                           
      INDMAX=IND                                                        
      RETURN                                                            
  100 FORMAT(I4)                                                        
  101 FORMAT(3(3I2,F5.2/))                                              
      END                                                               
      SUBROUTINE ROTATE (VECTOR,ROTMAT,ROTVEC)                          
!                                                                       
!     ROTATE PERFORMS A ROTATION OF THE VECTOR FROM THE GENERAL         
!     CARTESIAN COORDINATION SYSTEM INTO THE  LOCAL ONE  OF THE         
!     JATOM-TH SPHERE.                                                  
!     THIS SUBROUTINE IS ONLY REQUIRED FOR NONSYMMORPHIC CASES.         
!                                                                       
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      DIMENSION VECTOR(3),ROTVEC(3),ROTMAT(3,3)                         
!---------------------------------------------------------------------- 
!                                                                       
      DO 10 JCOORD=1,3                                                  
         DOTPRO=0.0D00                                                  
         DO 20 J=1,3                                                    
            DOTPRO=DOTPRO + VECTOR(J)*ROTMAT(JCOORD,J)                  
 20      CONTINUE                                                       
         ROTVEC(JCOORD)=DOTPRO                                          
 10   CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE ROTDEF (NAT,br2,ortho)                                
!                                                                       
!     ROTDEF GENERATES THE ROTATION-MATRICES ROTLOC(3,3,JATOM) FOR      
!     NONSYMMORPHIC STRUCTURES. THESE MATRICES PERFORM ROTATIONS        
!     FROM THE GENERAL COORDINATION-SYSTEM INTO LOCAL SYSTEMS WITH      
!     SPECIFIED POINTSYMMETRY.                                          
!     THE MATRICES  ROTIJ(3,3,INDEX)  TRANSFORM A POSITION IN   AN      
!     EQUIVALENT ATOM TO IT'S CORRESPONDING POSITION IN A NOTEQUIV      
!     ATOM.                                                             
!                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'param.inc'
!
      logical ortho
      COMMON /STRUK/ POS(ndif,3),AA,BB,CC,alpha(3), &
       RMT(nato),V(nato),PIA(3),VOL,IATNR(nato)   &
      ,MULT(nato),ISPLIT(nato)   
      COMMON /SYM2/ IORD,IZ(3,3,48),ITAU,TAU(3,48)                      
      COMMON /ROTMAT/   ROTLOC(3,3,NATO),ROTIJ(3,3,NDIF),TAUIJ(3,NDIF)  
!                                                                       
      DOUBLE PRECISION   GBAS(3,3), RBAS(3,3)
      dimension br2(3,3)
      DATA    TOLER    /1.D-4/                                          
      DATA    ONE      /1.D00/                                          
!---------------------------------------------------------------------  
!                                                                       
      PI = 4.0D+0*ATAN(1.0D+0)
      DO 15 I = 1,3
         DO 12 J = 1,3
            GBAS(I,J) = BR2(I,J)/2.0D+0/PI
 12      CONTINUE
 15   CONTINUE
      CALL GBASS(GBAS,RBAS)
      WRITE (66,6000) GBAS, RBAS
 6000 FORMAT ('RECIPROCAL AND REAL BRAVAIS MATRIX:',/,2(3(3F10.5,/)/))
!
      INDEX=0                                                           
      NCOUNT=0                                                          
      DO 20 JATOM=1,NAT
         INDEX1=INDEX+1                                                 
         DO 30 M=1,MULT(JATOM)                                          
            INDEX=INDEX+1                                               
            DO 25 I=1,IORD                                              
            X=0.D00                                                     
            DO 26 J=1,3                                                 
 26         X=X+IZ(J,1,I)*POS(INDEX1,j)                                 
            X=X+TAU(1,I) + 1.D00                                        
            X= MOD (X,ONE)                                              
            Y=0.D00                                                     
            DO 27 J=1,3                                                 
 27         Y=Y+IZ(J,2,I)*POS(INDEX1,j)                                 
            Y=Y+TAU(2,I) + 1.D00                                        
            Y= MOD (Y,ONE)                                              
            Z=0.D00                                                     
            DO 28 J=1,3                                                 
 28         Z=Z+IZ(J,3,I)*POS(INDEX1,j)                                 
            Z=Z+TAU(3,I) + 1.D00
            Z= MOD (Z,ONE)                                              
            X1=ABS(X-POS(INDEX,1))                                      
            Y1=ABS(Y-POS(INDEX,2))                                      
            Z1=ABS(Z-POS(INDEX,3))                                      
!          WRITE(*,*) 'JATOM,INDEX1,INDEX,I',JATOM,INDEX1,INDEX,I       
!            WRITE(*,*) X1,Y1,Z1                                        
            IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) THEN        
            NCOUNT=NCOUNT+1                                             
            DO 29 J=1,3                                                 
            TAUIJ(J,INDEX)= TAU(J,I)                                    
            DO 29 L=1,3                                                 
 29         ROTIJ(J,L,INDEX)=IZ(J,L,I)                                  
!
!     Redefine rotation matrix for non-orthogonal case
!
         IF(.not.ortho) CALL LOCDEF(RBAS,GBAS,ROTIJ(1,1,INDEX))
!
         write(66,6002) jatom,m,(rotij(1,i1,index),i1=1,3)           
         write(66,6001) (rotij(2,i1,index),i1=1,3)                      
         write(66,6001) (rotij(3,i1,index),i1=1,3)                      
 6001    format(3f10.5)
 6002    format(/,'LOCAL ROTATION MATRIX for ATOM',i3,'  INDEX:',i3,/ &
        3f10.5)
            GOTO 30                                                     
            END IF                                                      
 25         CONTINUE                                                    
            WRITE(*,*) 'NO SYMMETRY OPERATION FOUND IN ROTDEF'          
            WRITE(*,*) JATOM,INDEX                                      
            WRITE(*,*) (POS(JATOM,i1),I1=1,3)                           
            WRITE(*,*) (POS(INDEX,i1),I1=1,3)                           
            STOP 'ROTDEF'                                               
!                    
 30      CONTINUE            
 20   CONTINUE                                                          
      IF(NCOUNT.NE.INDEX) THEN                                          
         WRITE(66,1000) NCOUNT                                          
         STOP ' ROTDEF: NCOUNT NE INDEX'                                
      ENDIF                                                             

      RETURN                                                            
 1000 FORMAT(///,3X,'ERROR IN ROTDEF: ROTIJ NOT DEFINED FOR ALL ',       &
             'ATOMS OF BASIS',/,20X,'NCOUNT=',I2)                       
      END                                                               
      FUNCTION SININT(X)                                                
!CC SINE INTEGRAL FROM ZERO TO ARGUMENT X                               
      IMPLICIT REAL*8 (A-H,O-Z)
      IF (ABS(X).GE.16.) GOTO 1                                         
      Z = X/16.                                                         
      Y = 4.*Z**2-2.                                                    
      B =      +0.00000 00000 00007                                     
      A = Y*B  -0.00000 00000 00185                                     
      B = Y*A-B+0.00000 00000 04185                                     
      A = Y*B-A-0.00000 00000 84710                                     
      B = Y*A-B+0.00000 00015 22370                                     
      A = Y*B-A-0.00000 00241 00076                                     
      B = Y*A-B+0.00000 03329 88589                                     
      A = Y*B-A-0.00000 39729 08746                                     
      B = Y*A-B+0.00004 04202 71419                                     
      A = Y*B-A-0.00034 54691 55569                                     
      B = Y*A-B+0.00243 62214 04749                                     
      A = Y*B-A-0.01386 74455 89417                                     
      B = Y*A-B+0.06203 36794 32003                                     
      A = Y*B-A-0.21126 37809 76555                                     
      B = Y*A-B+0.53014 88479 1652                                      
      A = Y*B-A-0.96832 22369 8708                                      
      B = Y*A-B+1.38930 87711 7188                                      
      A = Y*B-A-1.92656 50911 5065                                      
      B = Y*A-B+2.77875 63817 4266                                      
      A = Y*B-A-4.06398 08449 119                                       
      A = Y*A-B+8.10585 29553 612                                       
      SININT = Z*.5*(A-B)                                               
      RETURN                                                            
!                                                                       
    1 Z = 16./X                                                         
      G = Z**2                                                          
      Y = 4.*G-2.                                                       
      B =      +0.00000 00000 00002                                     
      A = Y*B  -0.00000 00000 00014                                     
      B = Y*A-B+0.00000 00000 00107                                     
      A = Y*B-A-0.00000 00000 00964                                     
      B = Y*A-B+0.00000 00000 10308                                     
      A = Y*B-A-0.00000 00001 36096                                     
      B = Y*A-B+0.00000 00023 56196                                     
      A = Y*B-A-0.00000 00586 70317                                     
      B = Y*A-B+0.00000 24537 55677                                     
      A = Y*B-A-0.00023 37560 41393                                     
      A = Y*A-B+0.12452 74580 57854                                     
      F = Z*.5*(A-B)                                                    
      B =      +0.00000 00000 00002                                     
      A = Y*B  -0.00000 00000 00012                                     
      B = Y*A-B+0.00000 00000 00087                                     
      A = Y*B-A-0.00000 00000 00717                                     
      B = Y*A-B+0.00000 00000 06875                                     
      A = Y*B-A-0.00000 00000 79604                                     
      B = Y*A-B+0.00000 00011 69202                                     
      A = Y*B-A-0.00000 00234 68225                                     
      B = Y*A-B+0.00000 07249 95950                                     
      A = Y*B-A-0.00004 26441 82622                                     
      A = Y*A-B+0.00772 57121 93407                                     
      G = G*.5*(A-B)                                                    
      B = 1.57079 63267 9489                                            
      IF (X.LT.0.) B = -B                                               
      SININT = B-F*COS(X)-G*SIN(X)                                      
      RETURN                                                            
      END                                                               
      SUBROUTINE STERN                                                  
!                                                                       
!.... GENERATE STAR OF G                                                
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /FOUR/ G(3),NST,STG(3,48),TAUP(48)                         
      COMMON /SYM2/ IORD,IZ(3,3,48),ITAU,TAU(3,48)                      
      INTEGER G,STG,IND(48)                                             
      DATA TPI/6.2831853071796/                                         
      NST=0                                                             
      DO 1 I=1,IORD                                                     
                  TK=0.                                                 
      DO 2 J=1,3                                                        
               TK=TK+TAU(J,I)*G(J)*TPI                                  
      K=0                                                               
      DO 3 L=1,3                                                        
   3  K=IZ(J,L,I)*G(L)+K                                                
   2  STG(J,I)=K                                                        
      IF(NST.EQ.0) GOTO 7                                               
      DO 4 M=1,NST                                                      
      DO 5 J=1,3                                                        
      IF(STG(J,M).NE.STG(J,I)) GOTO 4                                   
   5  CONTINUE                                                          
      IND(M)=IND(M)+1                                                   
      TAUP(M)=TAUP(M)+COS(TK)                                           
      GOTO 1                                                            
   4  CONTINUE                                                          
   7  NST=NST+1                                                         
      DO 6 J=1,3                                                        
   6  STG(J,NST)=STG(J,I)                                               
       TAUP(NST)=COS(TK)                                                
      IND(NST)=1                                                        
   1  CONTINUE                                                          
      DO 10 I=1,NST                                                     
  10  TAUP(I)=TAUP(I)/IND(I)                                            
      RETURN                                                            
      END                                                               
      SUBROUTINE SUML(LMMAX,TCC,YL,F,FA,LM,FFX)                         
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 YL(49),IMAG,IMAG1                                      
      DIMENSION TCC(*),LM(2,*),FFX(*)                          
        REAL*8 MINU                                                     
      DATA IMAG /(0.0,1.0)/                                             
      F=0.0                                                             
      FA=0.0                                                            
      DO 10 ILM=1,LMMAX                                                 
      IMAG1=(1.,0.)                                                     
      LL=IABS(LM(1,ILM))                                                
      MM=LM(2,ILM)                                                      
      MINU=1.                                                           
      IF(LM(1,ILM).LT.0) THEN                                           
      IMAG1=-IMAG                                                       
      MINU=-1.                                                          
      END IF                                                            
      IF(MOD(MM,2).EQ.1) THEN                                           
      IMAG1=-IMAG1                                                      
      MINU=-MINU                                                        
      END IF                                                            
      LY1=LL*(LL+1)+MM+1                                                
      IF(MM.EQ.0) GOTO 11                                               
      LY2=LL*(LL+1)-MM+1                                                
      IF(MOD(LL,2).NE.0) GOTO 12                                        
! ...SYMMETRIC PART, MM NE.0                                            
      FFX(ILM)=TCC(ILM)*IMAG1*(YL(LY1)+YL(LY2)*MINU)*0.707107           
       F=F+FFX(ILM)                                                     
      GOTO 10                                                           
  11  IF(MOD(LL,2).NE.0) GOTO 13                                        
! ...SYMMETRIC PART,  MM EQ.0                                           
      FFX(ILM)=TCC(ILM)*IMAG1*YL(LY1)                                   
      F=F+FFX(ILM)                                                      
      GOTO 10                                                           
! ...ANTISYMMETRIC PART,  MM NE.0                                       
  12  FFX(ILM)=TCC(ILM)*IMAG1*(YL(LY1)+YL(LY2)*MINU)*0.707107           
      FA=FA+FFX(ILM)                                                    
      GOTO 10                                                           
! ...ANTISYMMETRIC PART, MM EQ.0                                        
  13  FFX(ILM)=TCC(ILM)*IMAG1*YL(LY1)                                   
      FA=FA+FFX(ILM)                                                    
!                                                                       
  10  CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE SUMLM(LMMAX,TCC,YL,F,FA,JATOM,FFX)                     
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 YL(49),ANG,IMAGM                                       
      DIMENSION TCC(*),FFX(*)                                     
      DATA IMAGM/(0.0,-1.0)/                                            
      FA=0.0                                                            
!                                                                       
!.....WRITTEN FOR CUBIC STRUCTUR                                        
!.....LM:  00  40  44  60  66                                           
!                                                                       
      CFO=SQRT(7./12.)                                                  
      CFF=SQRT(5./12.)                                                  
      CSO=SQRT(2.)/4.                                                   
      CSF=-SQRT(14.)/4.                                                 
      F=0.0
      fa=0.0                                                            
      F=F+TCC(1)*YL(1)                                                  
      FFX(1)=F                                                          
      ANG=(YL(25)+YL(17))/SQRT(2.)                                      
      F1=CFO*TCC(2)+CFF*TCC(3)                                          
      F2=CFO*YL(21)+CFF*ANG                                             
      F=F+F1*F2                                                         
      FFX(2)=F1*F2                                                      
      ANG=(YL(47)+YL(39))/SQRT(2.)                                      
      F1=CSO*TCC(4)+CSF*TCC(5)                                          
      F2=CSO*YL(43)+CSF*ANG                                             
      F=F+F1*F2                                                         
      FFX(3)=F1*F2                                                      
      FFX(4)=0.                                                         
      FFX(5)=0.                                                         
      FFX(6)=0.                                                         
      FFX(7)=0.                                                         
      IF(LMMAX.EQ.5) RETURN                                             
!ccc      IF(JATOM.EQ.1) RETURN                                         
!                                                                       
!.....ANFL-STRUCTUR                                                     
!.....AT.NR.2  LM=32                                                    
!                                                                       
      ANG=IMAGM*(YL(15)-YL(11))*0.707107                                
      FA=TCC(6)*ANG                                                     
      FFX(6)=FA                                                         
!......FA IS THE ANTISYMMETRIC PART OF THE STRUCTUR-AMPLIDUDE           
!......SEE KURKU-... ISR.J.CHEM.                                        
      RETURN                                                            
      END                                                               
      SUBROUTINE YLM(V,LMAX,Y)
      INTEGER            LMAX
      DOUBLE PRECISION   V(3)
      COMPLEX*16         Y(*)
!
!     ..................................................................
! 1.     PROGRAM UNIT 'YLM'
!           Calculates spherical harmonics
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           The spherical harmonics (Condon and Shortley convention)
!             Y(0,0),Y(1,-1),Y(1,0),Y(1,1),Y(2,-2) ... Y(LMAX,LMAX)
!           for vector V (given in Cartesian coordinates)
!           are calculated. In the Condon Shortley convention the
!           spherical harmonics are defined as
!                             +------+ 
!                        m    |   1     m              im(Phi)
!           Y(l,m) = (-1)  -+ | -----  P (cos(Theta)) e
!                            \| 2(Pi)   l
!                  m      
!           where P (cos(Theta)) is the normalized Associated Legendre
!                  l
!           function. Thus,
!                                          m      *
!                            Y(l,-m) = (-1) Y(l,m)
!           
!
! 3.     USAGE
!           DOUBLE PRECISION V(3), Y(5*5)
!           V(1) = ...
!           V(2) = ...
!           V(3) = ...
!           CALL YLM(V,4,Y)
!
!        ARGUMENT-DESCRIPTION
!           V      - DOUBLE PRECISION vector, dimension 3        (input)
!                    Must be given in Cartesian coordinates.
!                    Conversion of V to polar coordinates gives the
!                    angles Theta and Phi necessary for the calculation
!                    of the spherical harmonics.
!           LMAX   - INTEGER value                               (input)
!                    upper bound of L for which spherical harmonics
!                    will be calculated
!                    constraint:
!                       LMAX .GE. 0 (not checked)
!           Y      - COMPLEX*16 array, dimension (LMAX+1)**2    (output)
!                    contains the calculated spherical harmonics
!                    Y(1)                   for L .EQ. 0 (M = 0)
!                    Y(2), ..., Y(4)        for L .EQ. 1 (M = -1, 0, 1)
!                    ...
!                    Y(LMAX*LMAX+1), ..., Y((LMAX+1)*(LMAX+1))
!                                           for L .EQ. LMAX
!                                              (M = -L,...,L)
!                    constraint:
!                       Dimension of Y .GE. (LMAX+1)**2 (not checked)
!
!        USED SUBROUTINES (DIRECTLY CALLED)
!           none
!
!        INDIRECTLY CALLED SUBROUTINES
!           none
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           none
!
!        INPUT/OUTPUT (READ/WRITE)
!           none
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           Type COMPLEX*16 is used which does not conform to the
!           FORTRAN 77 standard.
!           Also the non-standard type conversion function DCMPLX()
!           is used which combines two double precision values into
!           one double complex value.
!
! 4.     REMARKS
!           none
!
! 5.     METHOD
!           The basic algorithm used to calculate the spherical
!           harmonics for vector V is as follows:
!
!           Y(0,0)
!           Y(1,0)
!           Y(1,1)
!           Y(1,-1) = -Y(1,1)
!           DO L = 2, LMAX
!              Y(L,L)   = f(Y(L-1,L-1)) ... Formula 1
!              Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
!              DO M = L-2, 0, -1
!                 Y(L,M) = f(Y(L-1,M),Y(L-2,M)) ... Formula 2
!                 Y(L,-M)= (-1)**M*Y(L,M)
!              ENDDO
!           ENDDO
!
!           In the following the necessary recursion formulas and
!           starting values are given:
!
!        Start:
!                        +------+
!                        |   1     
!           Y(0,0) =  -+ | -----  
!                       \| 4(Pi)  
!
!                                   +------+
!                                   |   3     
!           Y(1,0) =  cos(Theta) -+ | -----  
!                                  \| 4(Pi)  
!
!                                     +------+
!                                     |   3    i(Phi)
!           Y(1,1) =  - sin(Theta) -+ | ----- e
!                                    \| 8(Pi)  
!
!        Formula 1:
!
!           Y(l,l) =
!                           +--------+
!                           | (2l+1)   i(Phi)
!            -sin(Theta) -+ | ------  e       Y(l-1,l-1)
!                          \|   2l  
!
!        Formula 2:
!                                  +---------------+  
!                                  |  (2l-1)(2l+1)   
!           Y(l,m) = cos(Theta) -+ | -------------- Y(l-1,m)  -
!                                 \|   (l-m)(l+m)       
!
!                                    +--------------------+  
!                                    |(l-1+m)(l-1-m)(2l+1)
!                              -  -+ |-------------------- Y(l-2,m)
!                                   \|  (2l-3)(l-m)(l+m)                 
!
!        Formula 3: (not used in the algorithm because of the division
!                    by sin(Theta) which may be zero)
!
!                                    +--------------+  
!                      cos(Theta)    |  4(m+1)(m+1)   -i(Phi)
!           Y(l,m) = - ---------- -+ | ------------  e       Y(l,m+1) -
!                      sin(Theta)   \| (l+m+1)(l-m)       
!
!                                    +--------------+  
!                                    |(l-m-1)(l+m+2)  -2i(Phi)
!                              -  -+ |-------------- e        Y(l,m+2)
!                                   \| (l-m)(l+m+1)                         
!                                  
! 6.     DATE
!           26. April 1994                                   Version 1.2
!
!        INSTITUT FUER TECHNISCHE ELEKTROCHEMIE            --  TU VIENNA
!     ..................................................................
!
!
      INTEGER            I2L, I4L2, INDEX, INDEX2, L, M, MSIGN
      DOUBLE PRECISION   A, B, C, AB, ABC, ABMAX, ABCMAX
      DOUBLE PRECISION   D4LL1C, D2L13, PI
      DOUBLE PRECISION   COSTH, SINTH, COSPH, SINPH
      DOUBLE PRECISION   TEMP1, TEMP2, TEMP3
      DOUBLE PRECISION   YLLR, YLL1R, YL1L1R, YLMR
      DOUBLE PRECISION   YLLI, YLL1I, YL1L1I, YLMI
!
      PI = (4.0D+0)*ATAN(1.0D+0)
!
!        Y(0,0)
!
      YLLR = 1.0D+0/SQRT(4.0D+0*PI)
      YLLI = 0.0D+0
      Y(1) = DCMPLX(YLLR,YLLI)
!
!        continue only if spherical harmonics for (L .GT. 0) are desired
!
      IF (LMAX .LE. 0) GOTO 999
!
!        calculate sin(Phi), cos(Phi), sin(Theta), cos(Theta)
!        Theta, Phi ... polar angles of vector V
!
      ABMAX  = MAX(ABS(V(1)),ABS(V(2)))
      IF (ABMAX .GT. 0.0D+0) THEN
         A = V(1)/ABMAX
         B = V(2)/ABMAX
         AB = SQRT(A*A+B*B)
         COSPH = A/AB
         SINPH = B/AB
      ELSE
         COSPH = 1.0D+0
         SINPH = 0.0D+0
      ENDIF
      ABCMAX = MAX(ABMAX,ABS(V(3)))
      IF (ABCMAX .GT. 0.0D+0) THEN
         A = V(1)/ABCMAX
         B = V(2)/ABCMAX
         C = V(3)/ABCMAX
         AB = A*A + B*B
         ABC = SQRT(AB + C*C)
         COSTH = C/ABC
         SINTH = SQRT(AB)/ABC
      ELSE
         COSTH = 1.0D+0
         SINTH = 0.0D+0
      ENDIF
!
!        Y(1,0)
!
      Y(3) = DCMPLX(SQRT(3.0D+0)*YLLR*COSTH,0.0D+0)
!
!        Y(1,1) ( = -DCONJG(Y(1,-1)))
!
      TEMP1 = -SQRT(1.5D+0)*YLLR*SINTH
      Y(4) = DCMPLX(TEMP1*COSPH,TEMP1*SINPH)
      Y(2) = -DCONJG(Y(4))
!
      DO 20 L = 2, LMAX
         INDEX  = L*L+1
         INDEX2 = INDEX + 2*L
         MSIGN  = 1 - 2*MOD(L,2)
!
!        YLL = Y(L,L) = f(Y(L-1,L-1)) ... Formula 1
!
         YL1L1R = DBLE(Y(INDEX-1))
         YL1L1I = DIMAG(Y(INDEX-1))
         TEMP1 = -SQRT(DBLE(2*L+1)/DBLE(2*L))*SINTH
         YLLR = TEMP1*(COSPH*YL1L1R - SINPH*YL1L1I)
         YLLI = TEMP1*(COSPH*YL1L1I + SINPH*YL1L1R)
         Y(INDEX2) = DCMPLX(YLLR,YLLI)
         Y(INDEX)  = MSIGN*DCONJG(Y(INDEX2))
         INDEX2 = INDEX2 - 1
         INDEX  = INDEX  + 1
!
!        YLL1 = Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
!               (the coefficient for Y(L-2,L-1) in Formula 2 is zero)
!
         TEMP2 = SQRT(DBLE(2*L+1))*COSTH
         YLL1R = TEMP2*YL1L1R
         YLL1I = TEMP2*YL1L1I
         Y(INDEX2) = DCMPLX(YLL1R,YLL1I)
         Y(INDEX)  = -MSIGN*DCONJG(Y(INDEX2))
         INDEX2 = INDEX2 - 1
         INDEX  = INDEX  + 1
!
         I4L2 = INDEX2 - 4*L + 2
         I2L  = INDEX2 - 2*L
         D4LL1C = COSTH*SQRT(DBLE(4*L*L-1))
         D2L13  = -SQRT(DBLE(2*L+1)/DBLE(2*L-3))
!
         DO 10 M = L - 2, 0, -1
!
!        YLM = Y(L,M) = f(Y(L-2,M),Y(L-1,M)) ... Formula 2
!
            TEMP1 = 1.0D+0/SQRT(DBLE((L+M)*(L-M)))
            TEMP2 = D4LL1C*TEMP1
            TEMP3 = D2L13*SQRT(DBLE((L+M-1)*(L-M-1)))*TEMP1
            YLMR = TEMP2*DBLE(Y(I2L))  + TEMP3*DBLE(Y(I4L2))
            YLMI = TEMP2*DIMAG(Y(I2L)) + TEMP3*DIMAG(Y(I4L2))
            Y(INDEX2) = DCMPLX(YLMR,YLMI)
            Y(INDEX)  = MSIGN*DCONJG(Y(INDEX2))
!
            MSIGN  = -MSIGN
            INDEX2 = INDEX2 - 1
            INDEX  = INDEX  + 1
            I4L2   = I4L2   - 1
            I2L    = I2L    - 1
   10    CONTINUE
   20 CONTINUE
!
  999 RETURN
!
!        End of 'YLM'
!
      END
