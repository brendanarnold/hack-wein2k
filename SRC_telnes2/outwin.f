!BOP
! !ROUTINE: Outwin
! !INTERFACE:
      SUBROUTINE OUTWIN(REL,V,RNOT,DH,JRI,EH,FL,VAL,SLO,Nodes,Z)
! !USES:
      use work_atpar
      use dimension_constants,only : nrad
! !INPUT/OUTPUT PARAMETERS:
!  Rydberg Einheiten 
!  Input:
!    REL   switch for relativistic / nonrelativistic calculations
!    EH    Energie in Hartree 
!    FL    Drehimpuls 
!    Z     Kernladung
!    V     rad.sym. Potential in Hartree
!    RNOT  erster radialer Netzpunkt
!    DH    log. Schrittweite
!    JRI   Anzahl radialer Netzpunkte 
!  Output:
!    VAL,SLO:  Wellenfunktion und Steigung am Kugelrand
!    Nodes:    Anzahl Knoten 
! !DESCRIPTION:
!         Integration der skalarrel. Schroedingergleichung
! !REVISION HISTORY:
!     Taken from SRC\_lapw2 (?)
!     Updated November 2004 (Kevin Jorissen)
!EOP
      IMPLICIT none
      logical rel
! IN/OUT
      real*8,intent(in)   :: v(nrad),rnot,dh,eh,fl,z
	  real*8,intent(out)  :: val,slo
	  integer,intent(in)  :: jri
	  integer,intent(out) :: Nodes
! LOCAL VARIABLES
      real*8 e,rnet(nrad),zz,c,fllp1,r83sq,r1,r2,r3,h83,g0,s,sf,f0,aa,r,drdi,d(2,3),dg1,dg2,dg3,b22,b11,det,y,x,u,phi,df1,df2,df3
	  integer iiij,k

!     Hartree in Ryd
      E=EH*dble(2)

      do iiij=1,JRI
         RNET(iiij)=RNOT*(dexp(DH*(iiij-dble(1))))
      enddo

      Nodes = 0
      ZZ = Z + Z
      C = dble(274.074)
      if(.not.rel) C=dble(1.D+10)
      FLLP1 = FL*(FL + dble(1))
      R83SQ = dble(64)/dble(9)
      R1 = dble(1)/dble(9)
      R2 = dble(-5)*R1
      R3 = dble(19)*R1
      H83 = dble(8)/dble(3)
 

      G0 = dble(1)
      IF (Z .LT. dble(0.9)) THEN
        S = FL+dble(1)
        SF = FL
        F0 = FL/C
      ELSE
        AA = ZZ/C
        S = DSQRT(FLLP1 + dble(1) - AA*AA)
        SF = S
        F0 = G0*(S - dble(1))/AA
      ENDIF
      do K = 1,3
        R = RNET(K)
        DRDI = DH*R
        A(K) = (R**S)*G0
        B(K) = (R**SF)*F0
        D(1,K) = DRDI*A(K)*S/R
        D(2,K) = DRDI*B(K)*SF/R
      enddo


      DG1 = D(1,1)
      DG2 = D(1,2)
      DG3 = D(1,3)
      DF1 = D(2,1)
      DF2 = D(2,2)
      DF3 = D(2,3)
      do K = 4, JRI
        R = RNET(K)
        DRDI = DH*R

!       Faktor zwei vor V wegen Hartree-Rydberg !
        PHI = (E - dble(2)*V(K)/R)*DRDI/C
        U = DRDI*C + PHI
        X = -DRDI/R
        Y = -FLLP1*X*X/U + PHI
        DET = R83SQ - X*X + U*Y
        B11 = A(K-1)*H83 + R1*DG1 + R2*DG2 + R3*DG3
        B22 = B(K-1)*H83 + R1*DF1 + R2*DF2 + R3*DF3
        A(K) = (B11*(H83-X) + B22*U)/DET
        B(K) = (B22*(H83+X) - B11*Y)/DET
        IF (A(K)*A(K-1) .LT. dble(0)) Nodes = Nodes + 1
        DG1 = DG2
        DG2 = DG3
        DG3 = U*B(K) - X*A(K)
        DF1 = DF2
        DF2 = DF3
        DF3 = X*B(K) - Y*A(K)
      enddo
     
      do iiij=1,JRI
         B(iiij)=B(iiij)*c/dble(2)
      enddo

      VAL = A(JRI)/RNET(JRI)
      SLO = DG3/(DH*RNET(JRI))
      SLO = (SLO-VAL)/RNET(JRI) 
      RETURN
      END

