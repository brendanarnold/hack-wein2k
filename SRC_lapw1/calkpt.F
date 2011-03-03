      SUBROUTINE CALKPT(nmat_only)
!
#ifdef Parallel
      use matrices, only : S, H, HSROWS, Z, EIGVAL, KZZ
#else
      use matrices, only : HS, HSDIAG, HSROWS, Z, EIGVAL, KZZ
#endif
      use parallel, only : barrier, ictxtall, myid, deschs, close_parallel
      use lapw_timer, only : READ_CPU_TIME, time_coor, time_hamilt,  &
                             time_hns, time_diag,  &
                             START_TIMER, STOP_TIMER
      use lolog, only   : nlo, loor, lapw, ilo
      use orb, only     : BEXT, VORB, NMOD, NSP, NATORB, IATOM, NLORB, LORB
      use atspdt, only  : INS, NL, LMMAX, LQIND, LM, LQNS, DP, DPE, E, P, PE, PEI, VNS1, VNS2, VNS3, GFAC
      use comi, only    : NVAA, NE, NBELW, NAT, NT, LNSMAX
      use coml, only    : SPRSWF, PRNTWF, WFTAPE, REL, ITER, NOHNS
      use comr, only    : RKM, EL, EU, WEIGHT
      use lstapw, only  : RKMT, NV
      use onekpt, only  : SX1, SY1, SZ1
      use rotmat, only  : ROTIJ, ROTLOC
      use struk, only   : POS, ALAT, ALPHA, RMT, V, PIA, VI, IATNR, MULT
      use albl, only    : END_ALBL
      IMPLICIT NONE
      INCLUDE 'param.inc'
!
!     ..................................................................
! 1.     PROGRAM UNIT 'CALKPT'
!           Calculate Eigenvalues and Eigenvectors for one k-point.
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           The Hamilton- and Overlap matrices (HP and SP) are setup and
!           the corresponding general Eigenvalue problem Z*HP = E*SP*Z
!           is solved. Eigenvalues E and corresponding Eigenvectors Z
!           are computed for a specified energy window.
!
! 3.     USAGE
!           INTEGER NUMKPT, INFO
!           CALL INILPW(NUMKPT,INFO)
!           CALL INIKPT(INFO)
!           DO I = 1, NUMKPT
!              CALL SETKPT(I)
!              CALL CALKPT
!              CALL PRTKPT(I,INFO)
!           ENDDO
!
!        ARGUMENT-DESCRIPTION
!           none
!
!        USED SUBROUTINES (DIRECTLY CALLED)
!           COORS  - generate linearized augmented plane wave list
!           CPUTIM - measure used CPU-time
!           HAMILT - setup (spherical) Hamilton and Overlap-matrix
!           HNS    - add nonspherical contribution to Hamilton-matrix
!           SECLR4 - solve generalized eigenvalue problem
!           SECLR5 - solve generalized eigenvalue problem iteratively
!
!        INDIRECTLY CALLED SUBROUTINES
!           BISECT - determine eigenvalues of a symm. tridiag. matrix
!           DASUM  - sum of the magnitudes of vector elements (BLAS)
!           DCOPY  - copy a vector (BLAS)
!           DDOT   - dot product of two vectors (BLAS)
!           DSCAL  - multiply a vector by a scalar (BLAS)
!           DSUPMV - symmetric matrix (packed) times vector
!           DSUPR2 - rank two update of a symm. matrix (packed)
!           DVBES1 - calc. the derivatives of the spher. Bessel function
!           LMSORT - reorder (l' L l lmx m' m) to (l m l' m' lmx)
!           LOHNS  - Hamilton matrix local orbital index range
!           LOPW   - LAPW (K+G)-vector for local orbitals
!           NEBLW1 - compute number of eigenvalues below energy value
!           RDC_12 - blocked routine needed for REBK5
!           RDC_22 - blocked routine needed for REBK5
!           RDC_N2 - blocked routine needed for REBK5
!           REBK5  - backtrans. Eigenvectors to original general problem
!           REDC5  - reduce general symmetric Eigenvalue problem
!           ROTATE - rotate a vector according to a rotation matrix
!           SPHBES - compute spherical Bessel functions
!           TINVIT - determine Eigenvectors of symm. tridiag. matrix
!           TRBAK5 - backtransform Eigenvectors from tridiag form
!           TRED5  - tridiagonalize a real symmetric packed matrix
!           USTPHX - step function U(K)
!           YLM    - calculates spherical harmonics
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           INILPW - preliminalry calculations (set COMMON blocks)
!           INIKPT - setup values (valid for all k-points) not covered
!                    by 'INILPW'
!           SETKPT - set current k-point
!           PRTKPT - print out/save the results obtained by 'CALKPT'
!
!        INPUT/OUTPUT (READ/WRITE)
!           none
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           COMPLEX*16 declaration for GFAC is used
!
! 4.     REMARKS
!           'CALKPT' requires that data in the following COMMON blocks
!           is present:  ATSPDT, COMI,   COMR,   KVEC,
!                        LSTAPW, ONEKPT, ROTMAT, STRUK
!           These precalculations are done in routines INILPW, INIKPT
!           and SETKPT.
!           The resulting Eigenvalues and Eigenvectors are available
!           through COMMON block WORK (see description below in the
!           'Common blocks' section).
!
!           Hamilton- and Overlap-matrices are stored in a COMMON block
!           in order to avoid allocating them on the local stack. This
!           is necessary because they are far too large to fit onto the
!           local stack area.
!
! 5.     METHOD
!           - compute reciprocal lattice vectors for all considered
!             atoms in the crystal structure (for the current k-point)
!           - setup Hamilton matrix (spherical contribution) and
!             Overlap matrix
!           - add nonspherical contribution to the Hamilton matrix
!           - solve the generalized eigenvalue problem
!                HP*Z = EIGVAL*SP*Z       (HP ... Hamilton matrix)
!                                         (SP ... Overlap matrix )
!
!        INSTITUT FUER THEORETISCHE CHEMIE            --  TU VIENNA
!     ..................................................................
!
!        Locals
!
      logical            nmat_only
      INTEGER            NU, NLDAU, LLDAU
      DOUBLE PRECISION   TSEC1, TSEC2
      DOUBLE PRECISION   WTIME(3)
!
!        External Subroutines
!
      EXTERNAL           COORS, CPUTIM, HAMILT, HNS, SECLR4, SECLR5
!
!        setup timing
!
      CALL START_TIMER(time_coor)
!
!
!        set up array warp and generate plane wave list
!
      CALL COORS(RKM,NT,SX1,SY1,SZ1,NAT,nmat_only)
      if(nmat_only) return
      NVAA = NV
!
      EIGVAL(1) = EL
      EIGVAL(2) = EU
!
!        construct secular equation
!
      CALL STOP_TIMER(time_coor)
      CALL START_TIMER(time_hamilt)
!      write(*,*) warp(-1,-1,1)
      CALL HAMILT (NT,NAT,NV,P,PE,DP,DPE,MULT,RMT, &
                   POS,V,VI,PEI,E,INS)
!     Diagonal elements of Hamilt in HSDIAG, of Overlap in HS
      CALL STOP_TIMER(time_hamilt)
!
!        swap matrix diagonal
!
#ifndef Parallel
!_REAL      CALL DSWAP(HSROWS,HS,HSROWS+1,HSDIAG,1)
!_COMPLEX      CALL ZSWAP(HSROWS,HS,HSROWS+1,HSDIAG,1)
#endif
!     Diagonal elements of Hamilt in HS, of Overlap in HSDIAG

      CALL START_TIMER(time_hns)
      if(.not.nohns) &
      CALL HNS(NV,NAT,LNSMAX+1,POS,VI,RMT,MULT, &
               LQIND,LQNS,GFAC,ROTLOC,ROTIJ,VNS1(1,1,1,1), &
               VNS2(1,1,1,1),VNS3(1,1,1,1))

! Add orbital potential to the hamiltonian
#ifndef Parallel
       if(Nmod.ne.0)then
	 CALL  horb(NV,NAT,LNSMAX+1,POS,VI,RMT,MULT, &
                      LQIND,LQNS,PEI,ROTLOC,ROTIJ,HS)
!! Add LDA+U correction to the hamiltonian
!      if(NU.eq.1)then
!      CALL  Hldau(NV,NAT,LNSMAX+1,POS,VI,RMT,MULT, &
!                      LQIND,LQNS,PEI,ROTLOC,ROTIJ)
      endif
#endif
      CALL END_ALBL
      CALL STOP_TIMER(time_hns)
      CALL START_TIMER(time_diag)
!
!        solve secular equation
!
      if(ITER) then
!.....iterative diagonalization
           CALL SECLR5
      else
!.....full diagonalization
           CALL SECLR4
      endif
      CALL CPUTIM(TSEC1)
      CALL STOP_TIMER(time_diag)
      CALL BARRIER
!
      RETURN
!
!        End of 'CALKPT'
!
      END
