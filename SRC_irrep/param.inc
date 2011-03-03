!     Constant parameter definition
!
      INTEGER    NATO, NDIF, NUME, NMAT, LOMAX
      INTEGER    MSTP, MAXDG, MAXIRDG, FLMAX
      INTEGER    MAXIR, NSYM, LMAX2
      DOUBLE PRECISION   PI, TOLDG
!
!     NATO     - max. no. of nonequivalent atoms
!     NDIF     - max. no. of atoms
!     NUME     - max. no. of eigenvalues
!     NMAT     - max. no. of plane waves
!     LOMAX    - max. no. of local orbital. should be 
!                consistent with lapw1 and lapwso
!     MSTP     - max. step to describe k as a fraction
!     MAXDG    - max. no. of degenerate eigenfunctions
!     MAXIRDG  - max. no. of degenerate irr. representations
!     FLMAX    - size of flag (FL) array (should be 4)
!     MAXIR    - max. no. of irreducible representations 
!     NSYM     - max. no. of symmetry operations
!     TOLDG    - min. energy deviation of degenerate states,
!                in units of Rydberg
!
!      PARAMETER (NATO=     18)
!      PARAMETER (NDIF=     64)
!      PARAMETER (NUME=    300)
!      PARAMETER (NMAT=   3500)
      PARAMETER (LOMAX=     3) 
      PARAMETER (NLOAT=     3) 
      PARAMETER (LMAX2=    10) 
      PARAMETER (MSTP=   8000)
      PARAMETER (MAXDG=     6)
      PARAMETER (MAXIRDG=   4)
      PARAMETER (FLMAX=     4)
      PARAMETER (MAXIR=    48)
      PARAMETER (NSYM=     48)
      PARAMETER (TOLDG=0.5E-8)
!      PARAMETER (TOLDG=0.5E-5)
      PARAMETER (PI=3.141592653589793d0)      
