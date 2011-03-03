!BOP
! !ROUTINE: AveragedAngularSpectrum
! !INTERFACE:
      SUBROUTINE ATPAR(jATOM,nt)
! !USES:
      use work_atpar
      use energygrid
	  use potential
	  use struct
	  use dimension_constants,only : nrad
      use leftovers,only : cfein1, cfein2
! !INPUT/OUTPUT PARAMETERS:
!     jatom  :  atom position for which radial functions are needed (number
!                as in case.struct)
!     nt     :  l-value of radial function (I think)
! !DESCRIPTION:
!     Calculates APW radial basis functions for the wave function.
!     Input is taken from file (total spherical potential).
! !REVISION HISTORY:
!     Originally taken from SRC\_lapw2.
!     Created November 2004 (Kevin Jorissen)
!EOP
      implicit none
!  INPUT
      integer,intent(in) :: jatom,nt
!  LOCALS
	  integer nlo,jri,index,latom,j,m,nodes,idummy
      real*8 ei,fl,ovlp,duv,uv,dx,trx,rnot,cross,try,pei


	  call make_pot(nat,nrad)
!        read total spherical potential V(0,0) of TAPE18=VSP
!        norm of V(0,0)=V00(R)*R/SQRT(4.0D+0*PI)
      REWIND 18
      READ (18,5070)
      DO LATOM = 1, NAT
         READ (18,5010)
         READ (18,5020) IDUMMY
         READ (18,5060)
         READ (18,5040) (VR(J,LATOM),J=1,JRJ(LATOM))
         READ (18,5060)
         READ (18,5050)
         do J = 1, JRJ(LATOM)
            VR(J,LATOM) = VR(J,LATOM)/dble(2)
         enddo
      enddo

      LATOM = 0
      INDEX = 0
      NLO = 0
      RNOT = RO(JATOM)
      JRI = JRJ(JATOM)
      DX = DH(JATOM)

      do J = JEMIN,JEMAX

         FL = dble(NT)
!     get Energy in Ry
         EI = (ENE(j)/dble(13.6058)+EF)/dble(2)

!        calculate function at EI
         CALL OUTWIN(REL,VR(1,JATOM),RNOT,DX,JRI,EI,FL,UV,DUV,NODES,ZZ(JATOM))
         CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM)
         TRX = dble(1)/DSQRT(OVLP)
         do M = 1, JRI
            A(M) = TRX*A(M)
            B(M) = TRX*B(M)
         enddo

!        ensure orthogonalization
         CALL RINT13(CFEIN1,CFEIN2,A,B,AE,BE,CROSS,JATOM)
         TRY = -CROSS
         do M = 1, JRI
            AE(M) = AE(M) + TRY*A(M)
            BE(M) = BE(M) + TRY*B(M)
         enddo
         CALL RINT13(CFEIN1,CFEIN2,AE,BE,AE,BE,PEI,JATOM)
         do  M = 1, JRI
            A1(M,J) = A(M)
            B1(M,J) = B(M)
         enddo

      enddo

      call destroy_pot        
      RETURN

 5000 FORMAT (1X,I1,2F10.5,A4)
 5010 FORMAT (3X)
 5020 FORMAT (15X,I3,/,/)
 5040 FORMAT (3X,4E19.12)
 5050 FORMAT (/,/,/)
 5060 FORMAT (/)
 5070 FORMAT (/,/)
      END





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
!        E(l,i)       - expansion energy El for atom i
!        P(l,i)       - radial-function ul(r,E)
!                       (evaluated at r = muffin tin radius of atom i)
!        PE(l,i)      - derivative of the radial-function ul(r,E) with
!                       respect to energy E
!                       (evaluated at r = muffin tin radius of atom i)
!        PEI(l,i)     - norm of ul_dot(r,E) integrated over the muffin
!                       tin sphere
!
!
!        DH(i)   - step size of the logarithmic radial mesh for atom i
!        JRJ(i)  - number of radial (logarithmic) mesh points for atom i
!        RO(i)   - first radial mesh point for atom i
!        VR(j,i) - spherical part (l=0, m=0) of the total potential r*V
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
!        A, A1, A1LO, AE, AE1,
!        AE1LO, AP, B, B1, B1LO,  ..... various radial functions
!        BE, BE1, BE1LO, BP, VLM
