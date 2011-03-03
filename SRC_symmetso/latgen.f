      SUBROUTINE LATGEN(NAT,ndif,rotij,tauij,rotloc,pos,mult)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'param.inc'
!
!        Arguments
!
      INTEGER            NAT
!
!     ..................................................................
!
!        LATGEN generates the Bravais matrix BR2(3,3), which transforms
!        a reciprocal lattice vector of a special coordinate system (in
!        units of (2*PI)/ALAT(I)) into Cartesian system
!
!     ..................................................................
!
!        Common blocks
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
      COMMON  /CHAR/     LATTIC
      SAVE    /CHAR/
!
!        BR2(1:3,1:3) - reciprocal Bravais matrix
!
      DOUBLE PRECISION   BR2(3,3)
      COMMON  /GENER/    BR2
      SAVE    /GENER/
!
!        ORTH - .TRUE. for orthogonal lattice
!               .FALSE. otherwise
!
      LOGICAL            ORTHO
      COMMON  /ORTH/     ORTHO
      SAVE    /ORTH/
!
!        ALAT(1:3)   - lattice constants (a,b,c)
!        ALPHA(1:3)  - angles of the unit-cell's unit-vectors
!        MULT(i)     - number of equivalent atoms for inequiv. atom i
!        PIA(1:3)    - reciprocal lattice constants (2*pi/a, 2*pi/b,
!                      2*pi/c)
!        POS(1:3,i)  - position vector of atom i (including equivalent
!                      atoms) in the unit-cell
!        VI          - inverse volume of the direct lattice unit-cell
!
      INTEGER            MULT(*)
      DOUBLE PRECISION   VI
      DOUBLE PRECISION   ALAT(3), ALPHA(3), PIA(3), POS(3,*) 
      COMMON  /STRUK/    ALAT, ALPHA, PIA, VI
      SAVE    /STRUK/
      COMMON /SYM2/   IZ(3,3,noper),tau(3,noper),iord                   
      dimension rotij(3,3,*),tauij(3,*),rotloc(3,3,*)
!
!        Locals
!
      DOUBLE PRECISION   COSAB, COSAC, PI, RVFAC, SINAB, SINAC
!
!        External Subroutines
!
      EXTERNAL            ROTDEF
!
!        Intrinsic Functions
!
      INTRINSIC          ABS, ATAN, COS, SIN, SQRT
!
      ORTHO = .FALSE.
      PI = 4.0D+0*ATAN(1.0D+0)
      PIA(1) = 2.0D+0*PI/ALAT(1)
      PIA(2) = 2.0D+0*PI/ALAT(2)
      PIA(3) = 2.0D+0*PI/ALAT(3)
!
      IF (LATTIC(1:1) .EQ. 'H') THEN
!
!        hexagonal case
!
         BR2(1,1) = 2.0D+0/SQRT(3.0D+0)
         BR2(1,2) = 1.0D+0/SQRT(3.0D+0)
         BR2(1,3) = 0.0D+0
         BR2(2,1) = 0.0D+0
         BR2(2,2) = 1.0D+0
         BR2(2,3) = 0.0D+0
         BR2(3,1) = 0.0D+0
         BR2(3,2) = 0.0D+0
         BR2(3,3) = 1.0D+0
         RVFAC = 2.0D+0/SQRT(3.0D+0)
         ORTHO = .FALSE.
      ELSEIF ((LATTIC(1:1) .EQ. 'S') .OR. (LATTIC(1:1) .EQ. 'P')) THEN
!
!        primitive case (primitive lattice build up)
!
         IF (ABS(ALPHA(3)-PI/2.0D+0) .GT. 0.0001) THEN
            WRITE (*,*) 'gamma not equal 90'
!
!        M case (monoclinic lattice build up)
!        gamma .NE. 90
!
            SINAB = SIN(ALPHA(3))
            COSAB = COS(ALPHA(3))
            BR2(1,1) = PIA(1)/SINAB
            BR2(1,2) = -PIA(2)*COSAB/SINAB
            BR2(1,3) = 0.0D+0
            BR2(2,1) = 0.0D+0
            BR2(2,2) = PIA(2)
            BR2(2,3) = 0.0D+0
            BR2(3,1) = 0.0D+0
            BR2(3,2) = 0.0D+0
            BR2(3,3) = PIA(3)
            RVFAC = 1.0D+0/SINAB
            ORTHO = .FALSE.
!
!        do not normalize Bravais matrix
!
            GOTO 20
         ELSEIF (ABS(ALPHA(2)-PI/2.0D+0) .GT. 0.0001) THEN
            WRITE (*,*) 'beta not equal 90'
!
!        M case (monoclinic lattice build up)
!        beta .NE. 90
!
            SINAC = SIN(ALPHA(2))
            COSAC = COS(ALPHA(2))
            BR2(1,1) = PIA(1)/SINAC
            BR2(1,2) = 0.0D+0
            BR2(1,3) = -PIA(3)*COSAC/SINAC
            BR2(2,1) = 0.0D+0
            BR2(2,2) = PIA(2)
            BR2(2,3) = 0.0D+0
            BR2(3,1) = 0.0D+0
            BR2(3,2) = 0.0D+0
            BR2(3,3) = PIA(3)
            RVFAC = 1.0D+0/SINAC
!
!        do not normalize Bravais matrix
!
            ORTHO = .FALSE.
            GOTO 20
         ELSE
            BR2(1,1) = 1.0D+0
            BR2(1,2) = 0.0D+0
            BR2(1,3) = 0.0D+0
            BR2(2,1) = 0.0D+0
            BR2(2,2) = 1.0D+0
            BR2(2,3) = 0.0D+0
            BR2(3,1) = 0.0D+0
            BR2(3,2) = 0.0D+0
            BR2(3,3) = 1.0D+0
            RVFAC = 1.0D+0
            ORTHO = .TRUE.
         ENDIF
      ELSEIF (LATTIC(1:1) .EQ. 'F') THEN
!
!        fc case (bc lattice build up )
!
         BR2(1,1) = -1.0D+0
         BR2(1,2) =  1.0D+0
         BR2(1,3) =  1.0D+0
         BR2(2,1) =  1.0D+0
         BR2(2,2) = -1.0D+0
         BR2(2,3) =  1.0D+0
         BR2(3,1) =  1.0D+0
         BR2(3,2) =  1.0D+0
         BR2(3,3) = -1.0D+0
         RVFAC = 4.0D+0
         ORTHO = .TRUE.
      ELSEIF (LATTIC(1:1) .EQ. 'B') THEN
!
!        bc case (fc lattice build up)
!
         BR2(1,1) = 0.0D+0
         BR2(1,2) = 1.0D+0
         BR2(1,3) = 1.0D+0
         BR2(2,1) = 1.0D+0
         BR2(2,2) = 0.0D+0
         BR2(2,3) = 1.0D+0
         BR2(3,1) = 1.0D+0
         BR2(3,2) = 1.0D+0
         BR2(3,3) = 0.0D+0
         RVFAC = 2.0D+0
         ORTHO = .TRUE.
      ELSEIF (LATTIC(1:1) .EQ. 'C') THEN
         IF (LATTIC(2:3) .EQ. 'XZ') THEN
!
!        cxz case (cxz lattice build up)
!
!.....CXZ ORTHOROMBIC CASE
            IF(ABS(ALPHA(3)-PI/2.0D0).LT.0.0001) then
               BR2(1,1) =  1.0D+0
               BR2(1,2) =  0.0D+0
               BR2(1,3) =  1.0D+0
               BR2(2,1) =  0.0D+0
               BR2(2,2) =  1.0D+0
               BR2(2,3) =  0.0D+0
               BR2(3,1) = -1.0D+0
               BR2(3,2) =  0.0D+0
               BR2(3,3) =  1.0D+0
               RVFAC = 2.0D+0
               ORTHO = .TRUE.
            ELSE
!.....CXZ MONOCLINIC CASE
               write(*,*) 'gamma not equal 90'
               SINAB=SIN(ALPHA(3))
               COSAB=COS(ALPHA(3))
               BR2(1,1)= 1.0D+0/SINAB
               BR2(1,2)= -COSAB/SINAB
               BR2(1,3)= 1.0D+0/SINAB
               BR2(2,1)= 0.0
               BR2(2,2)= 1.0D+0
               BR2(2,3)= 0.0
               BR2(3,1)=-1.0D+0
               BR2(3,2)= 0.0
               BR2(3,3)= 1.0D+0
               RVFAC=2.0/SINAB
               ORTHO=.FALSE.
            ENDIF
!
         ELSEIF (LATTIC(2:3) .EQ. 'YZ') THEN
!
!        cyz case (cyz lattice build up)
!
            BR2(1,1) =  1.0D+0
            BR2(1,2) =  0.0D+0
            BR2(1,3) =  0.0D+0
            BR2(2,1) =  0.0D+0
            BR2(2,2) =  1.0D+0
            BR2(2,3) =  1.0D+0
            BR2(3,1) =  0.0D+0
            BR2(3,2) = -1.0D+0
            BR2(3,3) =  1.0D+0
            RVFAC = 2.0D+0
            ORTHO = .TRUE.
         ELSE
!
!        cxy case (cxy lattice build up)
!
            BR2(1,1) =  1.0D+0
            BR2(1,2) =  1.0D+0
            BR2(1,3) =  0.0D+0
            BR2(2,1) = -1.0D+0
            BR2(2,2) =  1.0D+0
            BR2(2,3) =  0.0D+0
            BR2(3,1) =  0.0D+0
            BR2(3,2) =  0.0D+0
            BR2(3,3) =  1.0D+0
            RVFAC = 2.0D+0
            ORTHO = .TRUE.
         ENDIF
      ELSEIF (LATTIC(1:1) .EQ. 'R') THEN
!
!        rhombohedral case
!
         BR2(1,1) =  1.0D+0/SQRT(3.0D+0)
         BR2(1,2) =  1.0D+0/SQRT(3.0D+0)
         BR2(1,3) = -2.0D+0/SQRT(3.0D+0)
         BR2(2,1) = -1.0D+0
         BR2(2,2) =  1.0D+0
         BR2(2,3) =  0.0D+0
         BR2(3,1) =  1.0D+0
         BR2(3,2) =  1.0D+0
         BR2(3,3) =  1.0D+0
         RVFAC = 6.0D+0/SQRT(3.0D+0)
         ORTHO = .FALSE.
      ELSE
!
!        Error: wrong lattice, stop execution
!
         GOTO 900
      ENDIF
!
!        normalize Bravais matrix
!        definitions according to column,rows convention of BR2
!
      BR2(1,1) = BR2(1,1)*PIA(1)
      BR2(2,1) = BR2(2,1)*PIA(2)
      BR2(3,1) = BR2(3,1)*PIA(3)
      BR2(1,2) = BR2(1,2)*PIA(1)
      BR2(2,2) = BR2(2,2)*PIA(2)
      BR2(3,2) = BR2(3,2)*PIA(3)
      BR2(1,3) = BR2(1,3)*PIA(1)
      BR2(2,3) = BR2(2,3)*PIA(2)
      BR2(3,3) = BR2(3,3)*PIA(3)
!
   20 CONTINUE
!
!       define rotation matrices if required
!
!      CALL ROTDEF(NAT,BR2,ORTHO,rotij,rotloc,pos,mult)
      CALL ROTDEF(iz,tau,iord,nat,pos,ndif,rotij,tauij,mult,lattic)
!
!        define inverse of cellvolume
!
      VI = RVFAC/ (ALAT(1)*ALAT(2)*ALAT(3))
!
      RETURN
!
  900 continue                               
      write(6,*)' LATGEN error'
      STOP 'LATGEN - Error'
!
!        End of 'LATGEN'
!
      END
