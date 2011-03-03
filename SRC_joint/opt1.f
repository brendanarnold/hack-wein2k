!........1.........2.........3.........4.........5.........6.........7
!234567890123456789012345678901234567890123456789012345678901234567890
!............................................................NOCULC....
!
      SUBROUTINE OPT1 (IB)
      use felder
      IMPLICIT REAL*8 (A-H,O-Z)
!  THE PARTIAL DS AND NS FROM ONE BAND AND ONE MICROZONE IS CALC.
      COMMON /EME/ EEF,EMIN, EMAX, EFA, DET, NU, NFUN, NO
      COMMON /EMICRO/ D(4), F(MG0,4), V,D1(4)
      COMMON /INTNEW/ VL,EPS,BEGIN,STEP
      real*8, allocatable ::  x(:),y(:)
      DIMENSION E0(4),H0(4),P0(4)
      allocate (x(NO),y(NO))
      DO NF=1,NFUN
      P0=F(NF,:)
      H0=D
      E0=D1
      VL=V
      EPS=EEF
      BEGIN=EMIN
      STEP=DET
      DO I=1,NO
       x(i)=BEGIN+(i-1)*DET
      ENDDO
	y=0
      CALL OPT0(E0,H0,P0,NO,X,Y)
      DO 80 I=1,NO
   80 DENSTY(IB,I,NF)=DENSTY(IB,I,NF)+Y(I)

	
      ENDDO
      END

      SUBROUTINE OPT0(E0,H0,P0,NUMB,X,Y)
!   Igor Mazin, written in 1984
!   This program computes the
!   following integrals over a microtetrahedron:(OM0=H0-E0)
!   Integral(P*Delta(OMG-OM)*(Theta(EPS-H)-Theta(EPS-E))*d^3k)
!   i.e. for optics, E0 are the energies of the lower band,
!   H0 are the energies of the upper band, EPS is the Fermi energy,
!   P0 are the matrix elements squared for a given transition at 
!   the four vertices of the tetrahedron, BEGIN is the minimal frequency,
!   STEP is the frequency step, NUMB is the total number of frequencies,
!   and the whole list of frequencies is stored in X.
!   The integration result is in the array Y
!
!   Optionally the program can compute the second derivative of the 
!   same quanitity with respect to the Fermi energy, D^2(Y)/D^EPS,
!   which e.g. controls termoreflectance, result returned in Z. 
!   For this, uncomment c (not C).
!
!
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION E0(4),E(4),A0(3),OM0(4),OM(4),P0(4),H(4),H0(4),P(4),PI0(3),B0(3),X(NUMB),Y(NUMB)
      INTEGER      IND(4)
      COMMON /INTNEW/ V,EPS,BEGIN,STEP,PP,DOM,PIMEAN,AA,A0,BB,B0,PI0,G,G1,G2
             common /tsting/iib,i
!  V-VOLUME OF MKTET (INPUT),BEG=X(1),STEP=X(I+1)-I),OTHER-INTERNAL VARS
! ****************************************************
      DO 11 I=1,4
11      OM0(I)=H0(I)-E0(I)
!  ORDERING OF OM0 (OM - ORDERED)
      DO 10 I=1,4
10    IND(I)=1
      DO 1 I=1,3
      DO 1 J=I+1,4
      IF(OM0(I).GT.OM0(J))GOTO 2
      IND(I)=IND(I)+1
      GOTO 1
2     IND(J)=IND(J)+1
1     CONTINUE
      DO 20 I=1,4
      JJ=IND(I)
      OM(JJ)=OM0((I))
       E(JJ)=E0((I))
       H(JJ)=H0((I))
20     P(JJ)=P0((I))
      MN=MAX(INT((OM(4)-BEGIN)/STEP)+1,1)
      MX=MIN(INT((OM(1)-BEGIN)/STEP)+1,NUMB)
      DO 100 I=mn,mx
      OMG=X(I)
      IF(OMG.GT.OM(4)) GOTO 4
!      TYPE*,'CASE EMPT'
      Y(I)=0.
      GOTO 100
4     IF(OMG.GT.OM(3)) GOTO 5
!      TYPE*,'CASE A'
      DOM=OMG-OM(4)
      DOM1=OM(1)-OM(4)
      DOM2=OM(2)-OM(4)
      DOM3=OM(3)-OM(4)
      AA=(EPS-E(4))/DOM
      BB=(EPS-H(4))/DOM
      A0(1)=(E(1)-E(4))/DOM1
      B0(1)=(H(1)-H(4))/DOM1
      A0(2)=(E(2)-E(4))/DOM2
      B0(2)=(H(2)-H(4))/DOM2
      A0(3)=(E(3)-E(4))/DOM3
      B0(3)=(H(3)-H(4))/DOM3
      PI0(1)=(P(1)-P(4))/DOM1
      PI0(2)=(P(2)-P(4))/DOM2
      PI0(3)=(P(3)-P(4))/DOM3
      F=DOM/DOM1/DOM2/DOM3
      PP=P(4)
      PIMEAN=(PI0(1)+PI0(2)+PI0(3))*.333333333
      CALL INTA
      TT=G
      CALL INTB
      Y(I)=-6*V*F*(G-TT)
      GOTO 100
5     IF(OMG.GT.OM(2)) GOTO 6
!      TYPE*,'CASE B'
      OMG3=OMG-OM(3)
      OM23=OM(2)-OM(3)
      EMDL=E(3)+(E(2)-E(3))*OMG3/OM23
      HMDL=H(3)+(H(2)-H(3))*OMG3/OM23
      PMDL=P(3)+(P(2)-P(3))*OMG3/OM23
      DOM=OMG-OM(4)
      DOM1=OM(1)-OM(4)
      DOM2=OM(2)-OM(4)
      DOM3=DOM
      AA=(EPS-E(4))/DOM
      A0(1)=(E(1)-E(4))/DOM1
      A0(2)=(E(2)-E(4))/DOM2
      A0(3)=(EMDL-E(4))/DOM3
      BB=(EPS-H(4))/DOM
      B0(1)=(H(1)-H(4))/DOM1
      B0(2)=(H(2)-H(4))/DOM2
      B0(3)=(HMDL-H(4))/DOM3
      PI0(1)=(P(1)-P(4))/DOM1
      PI0(2)=(P(2)-P(4))/DOM2
      PI0(3)=(PMDL-P(4))/DOM3
      F=1./DOM1/DOM2 *(OM(2)-OMG)/OM23
      PP=P(4)
      PIMEAN=(PI0(1)+PI0(2)+PI0(3))*.333333333
      CALL INTA
      TT=G
      CALL INTB
      TMPRR=-F*(G-TT)
      DOM=OM(1)-OMG
      DOM2=OM(1)-OM(3)
      DOM3=DOM1
      DOM1=DOM
      AA=- (E(1)-EPS)/DOM
      A0(1)=- (E(1)-EMDL)/DOM1
      A0(2)=- (E(1)-E(3))/DOM2
      A0(3)=- (E(1)-E(4))/DOM3
      BB=- (H(1)-EPS)/DOM
      B0(1)=- (H(1)-HMDL)/DOM1
      B0(2)=- (H(1)-H(3))/DOM2
      B0(3)=- (H(1)-H(4))/DOM3
      PI0(1)=- (P(1)-PMDL)/DOM1
      PI0(2)=- (P(1)-P(3))/DOM2
      PI0(3)=- (P(1)-P(4))/DOM3
      F=1./DOM2/DOM3  *(OMG-OM(3))/OM23
      PP=P(1)
      PIMEAN=(PI0(1)+PI0(2)+PI0(3))*.333333333
      CALL INTA
      TT=G
      CALL INTB
      Y(I)=6*V*(TMPRR-F*(G-TT))
      GOTO 100
6     IF(OMG.GT.OM(1)) GOTO 7
!      TYPE*,'CASE C'
      DOM=OM(1)-OMG
      DOM1=OM(1)-OM(2)
      DOM2=OM(1)-OM(3)
      DOM3=OM(1)-OM(4)
      AA=- (E(1)-EPS)/DOM
      A0(1)=- (E(1)-E(2))/DOM1
      A0(2)=- (E(1)-E(3))/DOM2
      A0(3)=- (E(1)-E(4))/DOM3
      BB=- (H(1)-EPS)/DOM
      B0(1)=- (H(1)-H(2))/DOM1
      B0(2)=- (H(1)-H(3))/DOM2
      B0(3)=- (H(1)-H(4))/DOM3
      PI0(1)=- (P(1)-P(2))/DOM1
      PI0(2)=- (P(1)-P(3))/DOM2
      PI0(3)=- (P(1)-P(4))/DOM3
      F=DOM/DOM1/DOM2/DOM3
      PP=P(1)
      PIMEAN=(PI0(1)+PI0(2)+PI0(3))*.333333333
      CALL INTA
      TT=G
      CALL INTB
      Y(I)=-6*V*F*(G-TT)
      GOTO 100
7      Y(I)=0.
!      TYPE*,'CASE FU'
100      CONTINUE
      END
      SUBROUTINE INTA
!
!
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION A0(3),A(3),PI0(3),PI(3),B0(3),B(3)
      INTEGER IND(3)
      COMMON /INTNEW/ V,EPS,BEG,STEP,PP,DOM,PIMEAN,AA,A0,BB,B0,PI0,G,G1,G2
             common /tsting/iib,indx
      DO 10 I=1,3
10    IND(I)=1
      DO 1 I=1,2
      DO 1 J=I+1,3
      IF (A0(I).GT. A0(J))GOTO 2
      IND(I)=IND(I)+1
      GOTO 1
2     IND(J)=IND(J)+1
1     CONTINUE
      DO 20 I=1,3
      JJ=IND(I)
      PI(JJ)=PI0((I))
20     A(JJ)=A0((I))
!      TYPE*,'INTA'
!      TYPE*,A,AA
      IF(AA .GT. A(3)) GOTO 4
!      TYPE*,'CASE -'
       G=0.
       G2=0.
      RETURN
4     IF(AA .GT. A(2)) GOTO 5
      DAA= AA- A(3)
      DAM = A(2)- A(3)
      DAL = A(1)- A(3)
      DPI=DAA*((PI(2)-PI(3))/DAM+(PI(1)-PI(3))/DAL)
!      TYPE*,'CASE 1'
      AR=1./DAM/DAL
      G=.5*DAA*DAA*AR*DOM*(PP+DOM*(PI(3)+DPI*.33333))
      RETURN
5     IF(AA .GT. A(1)) GOTO 6
!      TYPE*,'CASE  2'
      DAA= A(1)-AA
      DAM = A(1)- A(2)
      DAS = A(1)- A(3)
      DPI=DAA*((PI(2)-PI(1))/DAM+(PI(3)-PI(1))/DAS)
      AR=1./DAM/DAS
      G=.5*DOM*(PP+DOM*PIMEAN-DAA*DAA*AR*(PP+DOM*(PI(1)+DPI*.33333)))
      RETURN
6      G=.5*(PP+DOM*PIMEAN)*DOM
!      TYPE*,'CASE 3'
      RETURN
      END
      SUBROUTINE INTB
!
      IMPLICIT double precision (A-H,O-Z)
!
      DIMENSION B0(3),B(3),PI0(3),PI(3),A0(3),A(3)
      INTEGER      IND(3)
      COMMON /INTNEW/ V,EPS,BEG,STEP,PP,DOM,PIMEAN,AA,A0,BB,B0,PI0,G,G1,G2
             common /tsting/iib,indx
      DO 10 I=1,3
10    IND(I)=1
      DO 1 I=1,2
      DO 1 J=I+1,3
      IF (B0(I).GT. B0(J))GOTO 2
      IND(I)=IND(I)+1
      GOTO 1
2     IND(J)=IND(J)+1
1     CONTINUE
      DO 20 I=1,3
      JJ=IND(I)
      PI(JJ)=PI0((I))
20     B(JJ)=B0((I))
!      TYPE*,'INTB'
      IF(BB .GT. B(3)) GOTO 4
       G=0.
!      TYPE*,'CASE -1'
      RETURN
4     IF(BB .GT. B(2)) GOTO 5
      DBB= BB- B(3)
      DBM = B(2)- B(3)
      DBL = B(1)- B(3)
      DPI=DBB*((PI(2)-PI(3))/DBM+(PI(1)-PI(3))/DBL)
!      TYPE*,'CASE 1'
      BR=1./DBM/DBL
      G=.5*DBB*DBB*BR*DOM*(PP+DOM*(PI(3)+DPI*.33333))
      RETURN
5     IF(BB .GT. B(1)) GOTO 6
!      TYPE*,'CASE  2'
      DBB= B(1)-BB
      DBM = B(1)- B(2)
      DBS = B(1)- B(3)
      DPI=DBB*((PI(2)-PI(1))/DBM+(PI(3)-PI(1))/DBS)
      BR=1./DBM/DBS
      G=.5*DOM*(PP+DOM*PIMEAN-DBB*DBB*BR*(PP+DOM*(PI(1)+DPI*.33333)))
      RETURN
6      G=.5*(PP+DOM*PIMEAN)*DOM
!      TYPE*,'CASE 3'
      RETURN
      END
