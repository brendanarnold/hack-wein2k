      SUBROUTINE gpoint(spt,weight,npt1,LUSE)
!
!************************************************
!
!     GENERATE POINTS ON UNIT SPHERE
!
!************************************************
!
      IMPLICIT REAL*8(A-H,O-Z)
!
      INCLUDE 'param.inc'
!
      dimension spt(3,nspd),weight(nspd)
      DIMENSION PM(3),POS0(3),XX(LMAX1D)
!
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, FOUR=4.D0)
!      REAL*16 SUM
      REAL*8 SUM
!
      PI=FOUR*ATAN(ONE)
      FPI=FOUR*PI
      SFP=SQRT(FPI)
      LBIG=min(LUSE+1,LMAX1D)
      IF(LUSE+1 .gt. LMAX1D)then
        CALL OUTERR('LAPW0','Too large an L value for compiled LMAX1D')
        CALL OUTERR('LAPW0','Increase LMAX1D in param.inc or reduced harmonics')
        STOP 'LAPW0 - Error. Check file lapw0.error.'
      ENDIF  
!
      NPT1=0
      CALL GRULE(LBIG,XX,WEIGHT)
      DO 115 I=(LBIG+3)/2,LBIG
      J=I-(LBIG+1)/2
      WEIGHT(I)=WEIGHT(J)
      XX(I)=-XX(J)
  115 CONTINUE
      DO 110 J=1,LBIG
      Z=XX(J)
      RXY=SQRT(1.D0-Z*Z)
      DO 110 I=0,2*LBIG-2
      NPT1=NPT1+1
      PHI=PI*DBLE(2*I)/DBLE(2*LBIG-1)
      INDEX=LBIG*I+J
      WEIGHT(INDEX)=WEIGHT(J)
      SPT(1,INDEX)=COS(PHI)*RXY
      SPT(2,INDEX)=SIN(PHI)*RXY
      SPT(3,INDEX)=Z
!     write(6,*) index,(spt(j1,index),j1=1,3),weight(index)
 110  continue
!
      SUM=0.D0
      DO J=1,NPT1
        SUM=SUM+weight(j)
      enddo
      sum=4.D0*PI/sum
      do J=1,NPT1
        weight(j)=weight(j)*sum
      enddo
      write(6,*)'Gauss-Legendre grid of ',LBIG,'x',2*LBIG-1
      return
      end
      SUBROUTINE GRULE(N,X,W)
      IMPLICIT REAL*8 (A-H,O-Z)
!
!
!     DETERMINES THE (N+1)/2 NONNEGATIVE POINTS X(I) AND
!     THE CORRESPONDING WEIGHTS W(I) OF THE N-POINT
!     GAUSS-LEGENDRE INTEGRATION RULE, NORMALIZED TO THE
!     INTERVAL \-1,1\. THE X(I) APPEAR IN DESCENDING ORDER.
!
!     THIS ROUTINE IS FROM 'METHODS OF NUMERICAL INTEGRATION',
!     P.J. DAVIS AND P. RABINOWITZ, PAGE 369.
!
!
      DIMENSION X(N),W(N)
!      DATA ZERO,ONE,TWO,THREE,FOUR /0.0D0,1.0D0,2.0D0,3.0D0,4.0D0/
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, FOUR=4.D0)

      PI=FOUR*ATAN(ONE)
      M=(N+1)/2
      E1=N*(N+1)
      DO 1 I=1,M
      T=(4*I-1)*PI/(4*N+2)
      X0=(1.D0-(1.D0-1.D0/N)/(8.D0*N*N))*COS(T)
!--->    ITERATE ON THE VALUE  (M.W. JAN. 1982)
      DO 2 IT=1,3
      PKM1=ONE
      PK=X0
      DO 3 K=2,N
      T1=X0*PK
      PKP1=T1-PKM1-(T1-PKM1)/K+T1
      PKM1=PK
      PK=PKP1
    3 CONTINUE
      DEN=ONE-X0*X0
      D1=N*(PKM1-X0*PK)
      DPN=D1/DEN
      D2PN=(2.D0*X0*DPN-E1*PK)/DEN
      D3PN=(4.D0*X0*D2PN+(2.D0-E1)*DPN)/DEN
      D4PN=(6.D0*X0*D3PN+(6.D0-E1)*D2PN)/DEN
      U=PK/DPN
      V=D2PN/DPN
      H=-U*(1.D0+.5D0*U*(V+U*(V*V-U*D3PN/(3.D0*DPN))))
      P=PK+H*(DPN+.5D0*H*(D2PN+H/3.D0*(D3PN+.25D0*H*D4PN)))
      DP=DPN+H*(D2PN+.5D0*H*(D3PN+H*D4PN/3.D0))
      H=H-P/DP
      X0=X0+H
    2 CONTINUE
      X(I)=X0
      FX=D1- &
       H*E1*(PK+.5D0*H*(DPN+H/3.D0*(D2PN+.25D0*H*(D3PN+.2D0*H*D4PN))))
      W(I)=TWO*(ONE-X(I)*X(I))/(FX*FX)
    1 CONTINUE
      IF(M+M.GT.N) X(M)=ZERO
      RETURN
      END

