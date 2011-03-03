      SUBROUTINE INILPW(NAT,REL)
      use atomgrid
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'param.inc'
!
!        Arguments
!
!      INTEGER            INFO, NUMKPT
!      CHARACTER*(*)      DEFFN
!
!     ..................................................................
! 1.     PROGRAM UNIT 'INILPW'
!           Setup COMMON block data necessary for K-point processing.
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           Initialize the processing of K-points - setup data necessary
!           for the calculation of Hamilton and Overlap matrices and the
!           calculation of Eigenvalues and Eigenvectors).
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
!           DEFFN  - CHARACTER*(*) string                        (input)
!                    name of the LAPW1 definition file 'lapw1.def'
!           INFO   - INTEGER value                              (output)
!                    returns the error condition:
!                    INFO .EQ. 0 ... no error
!                    INFO .NE. 0 ... currently not supported
!           NUMKPT - INTEGER value                              (output)
!                    returns the number of k-points to process.
!
!        USED SUBROUTINES (DIRECTLY CALLED)
!           ATPAR  - solve radial Schroedinger equation, calc. gaunts
!           CPUTIM - measure used CPU-time
!           GAUNT2 - initialize GAUNT routine
!           LATGEN - setup lattice and, if required, rotation matrices
!           NN     - checks if overlapping atomic spheres exist
!           OUTERR - output a message on the standard error device
!           RDSWAR - initialize warpin calculation
!
!        INDIRECTLY CALLED SUBROUTINES
!           ABC    - calculate local orbital coefficients A,B and C
!           CBCOMB - combine the radial total potential coefficients
!           GAUNT1 - calculate gaunt factors
!           OUTWIN - solve an IVP ODE using a Runge-Kutta method
!           RINT13 - perform radial integrals
!           ROTDEF - find symmetry operations, define rotation matrices
!           SELECT - automatic energy search
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           CALKPT - calculate Eigenvalues/-vectors for one k-point
!           INIKPT - setup values (valid for all k-points) not covered
!                    by 'INILPW' (e.g. warpin)
!           SETKPT - set current k-point
!           PRTKPT - print out/save the results obtained by 'CALKPT'
!
!        INPUT/OUTPUT (READ/WRITE)
!        Input:
!           the following files are read by 'INILPW':
!              - a file containing all files to open
!                (given as a command-line parameter)
!              - structure-file of the crystal
!              - input-file for LAPW1
!        Output:
!           'INILPW' writes to the following files:
!              - output-file for LAPW1
!              - vector-file (energy parameters for the first k-point
!                are written)
!        Error messages (to an error file):
!           - 'INILPW - can''t open definition file ' FNAME
!             where FNAME is the name of the file (given as command-line
!             parameter) which contains all files to be opened.
!           - 'INILPW - can''t open unit: ' IUNIT
!             'INILPW -         filename: ' FNAME
!             'INILPW -           status: ' STATUS '  form: ' FORM
!             where IUNIT is the intended unitnumber for the file,
!                   FNAME is the filename of the file to be opened,
!                   STATUS is the status of the file to be opened,
!                   FORM is the formatting parameter of the file.
!           - 'INILPW - MULT .EQ. 0'
!             The number of equivalent atoms has to be greater 0.
!           - 'INILPW - Maximum number of K-points exceeded.'
!             The number of K-points to be processed must be less
!             than or equal to NKPT.
!           - 'INILPW - Too many atoms (NATO too small)'
!             The number of inequivalent atoms must be less than or
!             equal to NATO.
!           - 'INILPW - Error reading file: ' FNAME
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           - COMPLEX*16 declaration for GFAC is used
!
! 4.     REMARKS
!           none
!
! 5.     METHOD
!           - setup timing
!           - open all required files
!           - initialize gaunt routine
!           - read lattice parameters
!           - read crystal structure
!           - determine output-format parameters
!           - read convergence parameters
!           - setup lattice and, if required, rotation matrices
!           - solve radial Schroedinger equation, calculate gaunts
!             and radial integrals of nonspherical potential
!           - initialize warpin calculations
!           - read first K-point
!           - write energy parameters on vector-file for the
!             first K-point
!           - read the remaining K-points
!           - perform timing
!           - return the number of read K-points to the calling
!             routine
!
! 6.     DATE
!          25. August 1993                                  Version 1.05
!
!        INSTITUT FUER TECHNISCHE ELEKTROCHEMIE            --  TU VIENNA
!     ..................................................................
!
!
!        Common Blocks
!
!

!
!        LATTIC  - lattice type
!                  'P' ...... primitive latice (cubic, tetragonal,
!                             orthorhombic, monoclinic, triclin)
!                  'FC' ..... face centered
!                  'BC' ..... body centered
!                  'HEX' .... hexagonal
!                  'CXY' .... c-base centered (only orthorombic)
!                  'CYZ' .... a-base centered (only orthorombic)
!                  'CXZ' .... b-base centered (only orthorombic)
!                  'R' ...... rhombohedral
!        NAME(i) - name of atom i in the unit-cell (e.g. 'Titanium')
!
      CHARACTER*4        LATTIC
      CHARACTER*10,allocatable ::       NAME(:)
      COMMON  /CHAR/     LATTIC
      SAVE    /CHAR/
!

!
!        LNSMAX - maximum considered azimuthal quantum number l of ul
!                 for non-spherical hamilton matrix contributions
!        NAT    - number of inequivalent atoms
!        NBELW  - number of Eigenvalues below given energy window
!        NE     - number of calculated Eigenvalues
!        NT     - maximum considered azimuthal quantum number + 1
!                 (l+1 of ul) for the spherical contribution
!        NVAA   - number of plane waves
!
!      INTEGER            NVAA, NE, NBELW, NAT, NT, LNSMAX
!      COMMON  /COMI/     NVAA, NE, NBELW, NAT, NT, LNSMAX
!      SAVE    /COMI/
!
!        PRNTWF - if .TRUE. print wave functions to output-file
!        SPRSWF - if .TRUE. suppress wave function calculation
!        WFTAPE - if .TRUE. write wave functions to vector-file
!        REL    - if .TRUE. perform relativistic calculations
!
      LOGICAL            REL

!



!
!        DH(i)   - step size of the logarithmic radial mesh for atom i
!        JRJ(i)  - number of radial (logarithmic) mesh points for atom i
!        RO(i)   - first radial mesh point for atom i
!        VR(j,i) - spherical part (l=0, m=0) of the total potential r*V
!                  at mesh point j for atom i 
!
!
!        ROTIJ(1:3,1:3,i)  - rotate atom i (including equivalent ones)
!                            according to the symmetry. For atoms with
!                            only one equiv. atom this is always the
!                            identity matrix.
!        ROTLOC(1:3,1:3,i) - local rotation matrix for atom i (from
!                            global to local coordinate system)
!
      DOUBLE PRECISION   ROTLOC(3,3)
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

!
!        Local Scalars
!
      INTEGER            KINDEX, ITAPE, IUNIT, INDEX, JATOM, J
      INTEGER            I, I1, I2, NT1, IA, M, IRECL
      INTEGER            ISX, ISY, ISZ, IDV
      DOUBLE PRECISION   PI, E1, E2, RNN
      DOUBLE PRECISION   DTIME1, DTIME2, DTIMG1, DTIMG2
      CHARACTER*4        IREL
      CHARACTER*5        MODUS
      CHARACTER*11       STATUS, FORM
      CHARACTER*67       ERRMSG
      CHARACTER*80       FNAME, TITLE
!
!        External Subroutines
!
      EXTERNAL           OUTERR,  ATPAR
!
!        Intrinsic Functions
!
      INTRINSIC          LOG, EXP, ATAN
!

!
      PI = 4.0D+0*ATAN(1.0D+0)
!
      REWIND 20
!
!        calculate initialization data
!
      NT = LMAX
!


!
      READ (20,5000) TITLE
      WRITE (6,5000) TITLE
      READ (20,5010) LATTIC, NAT, IREL
      nato=nat
      allocate ( name(nato),pos(3,48*nato),jrj(nato),iatnr(nato),mult(nato))
      allocate ( v(nato),rmt(nato),zz(nato),vr(nrad,nato),r0(nato),dh(nato))
      DO 30 I = 1, NATO
         V(I)   = 0.0D+0
         RMT(I) = 0.0D+0
   30 CONTINUE
!
!        read in lattice constants
!
      READ (20,5020) ALAT(1), ALAT(2), ALAT(3), &
                     ALPHA(1), ALPHA(2), ALPHA(3)
      IF (ALPHA(1) .EQ. 0) ALPHA(1) = 90.0D+0
      IF (ALPHA(2) .EQ. 0) ALPHA(2) = 90.0D+0
      IF (ALPHA(3) .EQ. 0) ALPHA(3) = 90.0D+0
      ALPHA(1) = ALPHA(1)*PI/180.0D0
      ALPHA(2) = ALPHA(2)*PI/180.0D0
      ALPHA(3) = ALPHA(3)*PI/180.0D0
      IF (IREL .EQ. 'RELA') REL = .TRUE.
      IF (IREL .EQ. 'NREL') REL = .FALSE.
!
!        read crystal-structure (atompositions, symmetry-parameters,
!                                muffin-tin radius, ...)
!        'INDEX' counts all atoms in the unit cell,
!        'JATOM' counts only the notequivalent atoms
!
      INDEX = 0
      DO 50 JATOM = 1,NAT
         INDEX = INDEX + 1
!
         READ (20,5030) IATNR(JATOM),(POS(J,INDEX),J=1,3),MULT(JATOM)
         IF (MULT(JATOM) .EQ. 0) THEN
!
!           illegal number of equivalent atoms
!
            WRITE (6,6000) JATOM, INDEX, MULT(JATOM)
            GOTO 930
         ENDIF
         DO 40 M = 1, MULT(JATOM) - 1
            INDEX = INDEX + 1
            READ (20,5040) IATNR(JATOM), (POS(J,INDEX),J=1,3)
   40    CONTINUE
         READ (20,5050) NAME(JATOM),JRJ(JATOM),R0(JATOM),RMT(JATOM) &
                       ,ZZ(JATOM)
         DH(JATOM)  = LOG(RMT(JATOM)/R0(JATOM))/(JRJ(JATOM) - 1)
         RMT(JATOM) = R0(JATOM)*EXP(DH(JATOM)*(JRJ(JATOM)-1))
         READ (20,5060) ((ROTLOC(I1,I2),I1=1,3),I2=1,3)
   50 CONTINUE
!

!pi      WRITE (6,6010) LATTIC,IREL
!

      RETURN
      


!
!        error handling
!
  910 INFO = 1
!
!        'lapw1.def' couldn't be opened
!
!      WRITE (ERRMSG,9000) FNAME
!      CALL OUTERR('INILPW',ERRMSG)
      GOTO 999
  920 INFO = 2
!
!        file FNAME couldn't be opened
!
!      WRITE (ERRMSG,9010) IUNIT
!      CALL OUTERR('INILPW',ERRMSG)
!      WRITE (ERRMSG,9020) FNAME
!      CALL OUTERR('INILPW',ERRMSG)
!      WRITE (ERRMSG,9030) STATUS, FORM
!      CALL OUTERR('INILPW',ERRMSG)
      GOTO 999
  930 INFO = 3
!
!        illegal number of equivalent atoms
!
      WRITE (ERRMSG,9050) JATOM, MULT(JATOM), INDEX
      CALL OUTERR('INILPW',ERRMSG)
      CALL OUTERR('INILPW','MULT .EQ. 0')
      GOTO 999
  940 INFO = 4
!      CALL OUTERR('INILPW','Maximum number of K-points exceeded.')
      GOTO 999
  960 INFO = 7
!
!        Error reading file 'lapw1.def'
!
!      WRITE (ERRMSG,9040) FNAME
!      CALL OUTERR('INILPW',ERRMSG)
      GOTO 999
!
  999 RETURN
!
 5000 FORMAT(A80)
 5010 FORMAT(A4,23X,I3,/,13X,A4)
 5020 FORMAT(6F10.7)
 5030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)
 5040 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)
 5050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
 5060 FORMAT(20X,3F10.8)
 5070 FORMAT(A5)
 5080 FORMAT(F10.7,2I5)
 5090 FORMAT(20X,I1)
 5100 FORMAT(A10,4I5,3F5.2,A3)
 6000 FORMAT(///,3X,'ERROR IN EBND  : MULT(JATOM)=0 ...', &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
 6010 FORMAT(1H0,9X,A4,3X,'TYPE LATTICE ASSUMED'/10X,A4, &
              '-CALCULATION',/)
 6020 FORMAT(10X,'R-MT TIMES K-MAX IS',F5.2,/10X,'MAX L IS',I3,5X, &
      'MAX L IN NONSPHERICAL MATRIXELEMENTS:',I3,/ &
             ' NUMBER OF ATOMS IS',I5)
 6030 FORMAT(40X,'CPTIME ATPAR  =',F8.1)
 6040 FORMAT(///10X,'R-MT=',4F10.7/15X,4F10.7/15X,4F10.7)
 6050 FORMAT(10X,' FRACTIONAL VOLUME WITHIN MT=',4F10.7,2(/39X,4F10.7))
 6060 FORMAT(10X,' ONE/UNIT CELL VOLUME=',E16.9)
 6070 FORMAT(/10X,'LATTICE CONSTANTS ARE:',3F9.5)
 9000 FORMAT('can''t open definition file ',A40)
 9010 FORMAT('can''t open unit: ',I2)
 9020 FORMAT('       filename: ',A50)
 9030 FORMAT('         status: ',A,'  form: ',A)
 9040 FORMAT('Error reading file: ',A47)
 9050 FORMAT('MULT(',I3,'=',I3,', INDEX=',I3)
!
!        End of 'INILPW'
!
      END
