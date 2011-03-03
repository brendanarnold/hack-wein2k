      SUBROUTINE LATGEN(NAT)
!
      use char, only  : LATTIC
      use gener, only : BR2
      use orth, only  : ORTHO
      use struk, only : ALAT, ALPHA, PIA, VI, mult,pos,ndf
      use out, only   : imat,iord,tau, init_out
      use rotmat, only: rotij,tauij
      IMPLICIT NONE
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
!        Locals
!
      integer i,j,j1,j2,index,jatom
      DOUBLE PRECISION   COSAB, COSAC, PI, RVFAC, SINAB, SINAC, WURZEL
      DOUBLE PRECISION   COSBC, SINBC
      DOUBLE PRECISION   GBAS(3,3), RBAS(3,3)
!
!        External Subroutines
!
      EXTERNAL           OUTERR, ROTDEF
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
      IF(ABS(ALPHA(1)-PI/2.0D0).GT.0.0001) THEN
        WRITE(6,*) 'ALPHA NOT EQUAL 90'
        IF(ABS(ALPHA(2)-PI/2.0D0).GT.0.0001) THEN
          WRITE(6,*) 'BETA NOT EQUAL 90'
	  IF(ABS(ALPHA(3)-PI/2.0D0).GT.0.0001) THEN
	    WRITE(6,*) 'GAMMA NOT EQUAL 90'
          ENDIF
          GOTO 10
        ENDIF
        IF(ABS(ALPHA(3)-PI/2.0D0).GT.0.0001) THEN
          WRITE(6,*) 'GAMMA NOT EQUAL 90'
          GOTO 10
        ENDIF
      ENDIF
      IF(ABS(ALPHA(2)-PI/2.0D0).GT.0.0001) THEN
        WRITE(6,*) 'BETA NOT EQUAL 90'
        IF(ABS(ALPHA(3)-PI/2.0D0).GT.0.0001) THEN
          WRITE(6,*) 'GAMMA NOT EQUAL 90'
          GOTO 10
        ENDIF
      ENDIF
!     MONOCLINIC CASE -> GOTO 11
!cc      GOTO 11
!     TRICLINIC CASE
10    SINBC=SIN(ALPHA(1))
      COSAB=COS(ALPHA(3))
      COSAC=COS(ALPHA(2))
      COSBC=COS(ALPHA(1))
      WURZEL=SQRT(SINBC**2-COSAC**2-COSAB**2+2*COSBC*COSAC*COSAB)
      BR2(1,1)= SINBC/WURZEL*PIA(1)
      BR2(1,2)= (-COSAB+COSBC*COSAC)/(SINBC*WURZEL)*PIA(2)
      BR2(1,3)= (COSBC*COSAB-COSAC)/(SINBC*WURZEL)*PIA(3)
      BR2(2,1)= 0.0
      BR2(2,2)= PIA(2)/SINBC
      BR2(2,3)= -PIA(3)*COSBC/SINBC
      BR2(3,1)= 0.0
      BR2(3,2)= 0.0
      BR2(3,3)= PIA(3)
!
      RVFAC= 1.d0/WURZEL
      ORTHO=.FALSE.
!
      GOTO 20
!
!     MONOCLINIC CASE
!
  11  IF (ABS(ALPHA(3)-PI/2.0D+0) .GT. 0.0001) THEN
            WRITE (6,*) 'gamma not equal 90'
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
            WRITE (6,*) 'beta not equal 90'
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
               write(6,*) 'gamma not equal 90'
               SINAB=SIN(ALPHA(3))
               COSAB=COS(ALPHA(3))
               BR2(1,1)= pia(1)*1.0D+0/SINAB
               BR2(1,2)= -pia(2)*COSAB/SINAB
               BR2(1,3)= pia(1)*1.0D+0/SINAB
               BR2(2,1)= 0.0
               BR2(2,2)= pia(2)*1.0D+0
               BR2(2,3)= 0.0
               BR2(3,1)=-pia(3)*1.0D+0
               BR2(3,2)= 0.0
               BR2(3,3)= pia(3)*1.0D+0
               RVFAC=2.0/SINAB
               ORTHO=.FALSE.
!              do not renormalize BR2
               goto 20
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
!      CALL ROTDEF(NAT,BR2,ORTHO,lattic,mult)

 5000 FORMAT (I4)
 5010 FORMAT (3 (3I2,F11.8,/))
      READ (20,5000) IORD
      CALL init_out(IORD,0,0,0)
      DO  J = 1, IORD
         READ (20,5010) ((IMAT(J1,J2,J),J1=1,3),TAU(J2,J),J2=1,3)
      enddo
!
      DO  I = 1,3
         DO  J = 1,3
            GBAS(I,J) = BR2(I,J)/2.0D+0/PI
         enddo
      enddo
      CALL GBASS(GBAS,RBAS)
      WRITE (6,6000) GBAS, RBAS
 6000 FORMAT (3F10.5)
      CALL ROTDEF(imat,tau,iord,nat,pos,ndf,rotij,tauij,mult,lattic)
!
!     Redefine rotation matrix for non-orthogonal case
!     for  mon.CXZ type rotation skipped (caution!!!)
      IF(.NOT. ORTHO.and.lattic(1:3).ne.'CXZ')  then
       INDEX = 0
       DO JATOM = 1, NAT
         DO J = 1, MULT(JATOM)
            INDEX = INDEX + 1
            CALL LOCDEF(RBAS,GBAS,ROTIJ(1,1,INDEX))
         enddo
       enddo
      endif
!
!        define inverse of cellvolume
!
      VI = RVFAC/ (ALAT(1)*ALAT(2)*ALAT(3))
!
      RETURN
!
!        Error messages
!
  900 CALL OUTERR('LATGEN','wrong lattice.')
      STOP 'LATGEN - Error'
!
!        End of 'LATGEN'
!
      END

