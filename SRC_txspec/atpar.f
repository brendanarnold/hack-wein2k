      SUBROUTINE ATPAR(NATOM,NT,NAT,REL,JEMIN,JEMAX,a1,b1,ene)
      use atomgrid
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'param.inc'
!
!        Arguments
!
      INTEGER            LNSMAX, NAT, NT
      LOGICAL            REL
!
!     ..................................................................
!
!        ATPAR calculates the solutions u(l) of the radial Schroedinger
!        Equation, the energy derivatives ue(l), the radial non muffin-
!        tin integrals uvu, uevu, and the Gaunt-coefficients.
!        Spin-orbitless relativistic equations are used.
!
!     ..................................................................
!
!        Common blocks
!
!
!        DP(l,i)      - derivative of the radial-function ul(r,E) with
!                       respect to radius r
!                       (evaluated at r = muffin tin radius of atom i)
!        DPE(l,i)     - derivative of the energy derivative ul_dot(r,E)
!                       of the radial-function ul(r,E) with respect to
!                       radius r
!                       (evaluated at r = muffin tin radius of atom i)
!        GFAC(:,i)    - gaunt-factors for atom i
!        INS          - flag if nonspherical correction (warpin) should
!                       take place
!                       if (INS .EQ. 0) then do warpin otherwise don't
!        LM(:,j,i)    - L,M combination for lattice harmonics j for
!                       atom i
!                       LM(1,j,i) ... L for j-th gaunt factor
!                       LM(2,j,i) ... M for j-th gaunt factor
!        LMMAX(i)     - number of L,M combinations for the
!                       nonspherical contribution of atom i
!        LQIND(i)     - number of gaunt factors for atom i
!        LQNS(:,j,i)  - l,m,LM,l',m' combination corresponding to
!                       the gaunt-factor GFAC(j,i) for atom i
!                       LQNS(1,1:LQIND(i),i) ... l'+1 for atom i
!                       LQNS(2,1:LQIND(i),i) ... L
!                       LQNS(3,1:LQIND(i),i) ... l+1 for atom i
!                       LQNS(4,1:LQIND(i),i) ... index of LM for atom i
!                       LQNS(5,1:LQIND(i),i) ... m' for atom i
!                       LQNS(6,1:LQIND(i),i) ... m for atom i
!        E(l,i)       - expansion energy El for atom i
!        P(l,i)       - radial-function ul(r,E)
!                       (evaluated at r = muffin tin radius of atom i)
!        PE(l,i)      - derivative of the radial-function ul(r,E) with
!                       respect to energy E
!                       (evaluated at r = muffin tin radius of atom i)
!        PEI(l,i)     - norm of ul_dot(r,E) integrated over the muffin
!                       tin sphere
!
      INTEGER            INS
!      real*8,allocatable ::  DP(:,:), DPE(:,:), E(:,:) 
!      real*8,allocatable ::  P(:,:), PE(:,:), PEI(:,:)
!
!        DH(i)   - step size of the logarithmic radial mesh for atom i
!        JRJ(i)  - number of radial (logarithmic) mesh points for atom i
!        RO(i)   - first radial mesh point for atom i
!        VR(j,i) - spherical part (l=0, m=0) of the total potential r*V
!                  at mesh point j for atom i 
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
!        VI          - inverse volume of the direct lattice unit-cell
!
      DOUBLE PRECISION   VI
      DOUBLE PRECISION   ALAT(3), ALPHA(3), PIA(3) 
      COMMON  /STRUK/    ALAT, ALPHA, PIA, VI
      SAVE    /STRUK/
!
!        A, A1, A1LO, AE, AE1,
!        AE1LO, AP, B, B1, B1LO,  ..... various radial functions
!        BE, BE1, BE1LO, BP, VLM
!
!        Common block WORK is used as a work-space and to exchange data
!        between the following routines:
!
!           ATPAR, OUTWIN
!
!        This common block will be reused with other data elements
!        in the routines:
!
!           HAMILT, HNS, PRTKPT, SECLR4
!
!        These routines are not conflicting with the above routines -
!        no fancy overlaying of data takes place, just a reuse of
!        the common block name (and therefore its storage) with
!        completely different data.
!
      DOUBLE PRECISION   A(NRAD), A1(NRAD,jemax)
      DOUBLE PRECISION   AE(NRAD)
      DOUBLE PRECISION   AP(NRAD), B(NRAD), B1(NRAD,jemax)
      DOUBLE PRECISION   BE(NRAD)
      DOUBLE PRECISION   BP(NRAD)
      COMMON  /WORK/     A,B,AP,BP,AE,BE,VLM            
      SAVE    /WORK/

      DOUBLE PRECISION ENE(jemax),EF
      COMMON  /ENE/      EF
      SAVE    /ENE/
!
!        Local Scalars
!
      INTEGER            I, ICOUNT, IDUMMY, INDEX, J, JATOM, JC, JCOL
      INTEGER            JRI, JROW, L, L0, L0BEFO, LALT, LATOM, LL
      INTEGER            LLBEFO, LLL, LLMM, LM1, LM2, LMBEFO, LMX, LP
      INTEGER            LPBEFO, LQX, M, MINU, MM, MP, NLR, NODEL
      INTEGER            NODES, NODEU
      DOUBLE PRECISION   AINT, CFEIN1, CFEIN2, CROSS, DE, DELE, DELEI
      DOUBLE PRECISION   DUV, DUVB, DX, E1, EI, FL, OVLP, RNOT
      DOUBLE PRECISION   TRX, TRY, UV, UVB
      COMPLEX*16         IMAG, IMAG1
      CHARACTER*4        EMAIN
      CHARACTER*67       ERRMSG
!fb
!     variables for [yu91]-force formalism
      logical     force
!fe
!
!        External Subroutines
!
      EXTERNAL           OUTERR,OUTWIN,RINT13
!
!        Intrinsic Functions
!
      INTRINSIC          ABS, DBLE, MOD, SQRT
!
!        Data statements
!
      DATA DELE /2.0D-3/
      DATA IMAG /(0.0D+0,1.0D+0)/
!
!fb
      force=.false.
!      if(allocated(dp) ) then
!      deallocate (  DP, DPE, E) 
!      deallocate (  P, PE, PEI)
!      endif
!      allocate (  DP(LMAX,NAT), DPE(LMAX,NAT), E(LMAX,NAT)) 
!      allocate (  P(LMAX,NAT), PE(LMAX,NAT), PEI(LMAX,NAT))
!fe
!        initialize constants
!        CFEIN1 and CFEIN2 are used in RINT13 to calculate the
!        integrals of radialfunctions and are derivates of the
!        Feinstruktur-constant (in Hartrees)
!
      IF (REL) THEN
         CFEIN1 = 1.0D+0
	 CFEIN2 = 1.d0/137.0359895d0**2
      ELSE
         CFEIN1 = 1.0D+0
         CFEIN2 = 1.0D-22
      ENDIF
!
!        read total spherical potential V(0,0) of TAPE18=VSP
!        norm of V(0,0)=V00(R)*R/SQRT(4.0D+0*PI)
!
      REWIND 18
      READ (18,5070)
      DO 20 JATOM = 1, NAT
         READ (18,5010)
         READ (18,5020) IDUMMY
         READ (18,5060)
         READ (18,5040) (VR(J,JATOM),J=1,JRJ(JATOM))
         READ (18,5060)
         READ (18,5050)
         DO 10 J = 1, JRJ(JATOM)
            VR(J,JATOM) = VR(J,JATOM)/2.0D+0
   10    CONTINUE
   20 CONTINUE
!

   30 CONTINUE
!
!        start loop over atoms in unitcell
!
      LATOM = 0
      INDEX = 0
      NLO = 0
!pi
      JATOM = NATOM
      RNOT = R0(JATOM)
      JRI = JRJ(JATOM)
      DX = DH(JATOM)
        
!pi
!     loop over Energy-Range
!     Energy steps are stored in ENE()

      DO 150 J = JEMIN,JEMAX
         
         FL = NT

!     get Energy in Ry
         EI = (ENE(j)/13.6058+EF)/2.0D+0


!         write(50,*)EI*2


!
!        calculate function at EI
!
            CALL OUTWIN(REL,VR(1,JATOM),RNOT,DX,JRI,EI,FL,UV,DUV,NODES &
                        ,ZZ(JATOM))
            CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM)
            TRX = 1.0D+0/SQRT(OVLP)
!            P(J,JATOM) = TRX*UV
!            DP(J,JATOM) = TRX*DUV
            DO 120 M = 1, JRI
               A(M) = TRX*A(M)
               B(M) = TRX*B(M)
  120       CONTINUE
!
!        insure orthogonalization
!
            CALL RINT13(CFEIN1,CFEIN2,A,B,AE,BE,CROSS,JATOM)
            TRY = -CROSS
!            IF (TRY .LT. -0.05D+0) WRITE (6,6070) L, TRY, OVLP
            DO 130 M = 1, JRI
               AE(M) = AE(M) + TRY*A(M)
               BE(M) = BE(M) + TRY*B(M)
              
  130       CONTINUE
!            PE(J,JATOM) =  TRY*P(J,JATOM)
!            DPE(J,JATOM) = TRY*DP(J,JATOM)
!            CALL RINT13(CFEIN1,CFEIN2,AE,BE,AE,BE,PEI(J,JATOM),JATOM)
            CALL RINT13(CFEIN1,CFEIN2,AE,BE,AE,BE,PEI,JATOM)

!            IF (J .LE. NSLMAX) THEN
               DO 140 M = 1, JRI

                  A1(M,J) = A(M)
                  B1(M,J) = B(M)
  140          CONTINUE
             
!      

!            WRITE (6,6030) L, P(J,JATOM), DP(J,JATOM), PE(J,JATOM),
!     ,                     DPE(J,JATOM), PEI(J,JATOM), NODEL, NODES,
!     ,                     NODEU
            
!           

  150    CONTINUE
!
!     PI end atpar here!
!
         RETURN
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!

!
!
 5010 FORMAT (3X)
 5020 FORMAT (15X,I3,/,/)
 5040 FORMAT (3X,4E19.12)
 5050 FORMAT (/,/,/)
 5060 FORMAT (/)
 5070 FORMAT (/,/)
 6070 FORMAT (10X,'FOR L=',I2,' CORRECTION=',E14.5,' OVERLAP=',E14.5)
!
!        End of 'ATPAR'
!
      END
