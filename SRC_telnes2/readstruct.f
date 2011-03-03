!BOP
! !ROUTINE: ReadStruct
! !INTERFACE:
      SUBROUTINE ReadStruct
! !USES:
      use input, only : natom
      use dimension_constants
	  use struct
	  use rotation_matrices
	  use leftovers,only : cfein1, cfein2
      use crash
! !DESCRIPTION:
!    Reads structural data from case.struct.  Calls latgen, rotdef and readrotij.
! !REVISION HISTORY:
!     Ripped from wien inilpw routine.
!     Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none

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

!        Local Scalars
      INTEGER            INDEX,JATOM,J,I,I1,I2,M,mul
      CHARACTER*4        IREL


      REWIND 20


      write(6,'(/,a)') 'Output from subroutine inilpw :'

      READ (20,*)
      READ (20,5010) LATTIC, NAT, IREL
      
      READ (20,5020) ALAT(1), ALAT(2), ALAT(3), &
                     ALPHA(1), ALPHA(2), ALPHA(3)
      ALPHA(1) = ALPHA(1)*3.14159265359D0/180.0D0
      ALPHA(2) = ALPHA(2)*3.14159265359D0/180.0D0
      ALPHA(3) = ALPHA(3)*3.14159265359D0/180.0D0
      IF (IREL .EQ. 'RELA') REL = .TRUE.
      IF (IREL .EQ. 'NREL') REL = .FALSE.

!     CFEIN1 and CFEIN2 are used in RINT13 and RINT14 to calculate the
!     integrals of radialfunctions and are derivates of the Feinstruktur-constant (in Hartrees).
      IF (REL) THEN
         CFEIN1 = dble(1)
         CFEIN2 = dble(4)*dble(1.331258D-5)
      ELSE
         CFEIN1 = dble(1)
         CFEIN2 = dble(0)
      ENDIF      

! now we first get the total number of atoms
      index=0
	  do jatom=1,nat
	    index=index+1
	    read(20,5031) mul
		if(jatom.eq.natom) j=mul
	    do m=1,mul-1
	      index=index+1
	      read(20,*)
	    enddo
	    do m=1,4;read(20,*);enddo
	  enddo
	  rewind(20)
      do m=1,4;read(20,*);enddo
      nats=index
! now we allocate arrays
	  call make_struct(nat,nats)
	  call make_rot(j,nat)  ! used to be (nats,nat)

!   now we read all data properly
      INDEX = 0
      DO JATOM = 1,NAT
         INDEX = INDEX + 1

         READ (20,5030) IATNR(JATOM),(POS(J,INDEX),J=1,3),MULT(JATOM),ISPLIT(JATOM)
         IF (MULT(JATOM) .EQ. 0) THEN
            WRITE (6,6000) JATOM, INDEX, MULT(JATOM)
            GOTO 930
         ENDIF
          DO M = 1, MULT(JATOM) - 1
            INDEX = INDEX + 1
            READ (20,5040) IATNR(JATOM), (POS(J,INDEX),J=1,3)
         enddo
         READ (20,5050) NAME(JATOM),JRJ(JATOM),RO(JATOM),RMT(JATOM) &
                       ,ZZ(JATOM)
         DH(JATOM)  = DLOG(RMT(JATOM)/RO(JATOM))/(JRJ(JATOM) - 1)
         RMT(JATOM) = RO(JATOM)*DEXP(DH(JATOM)*(JRJ(JATOM)-1))
         READ (20,5060) ((ROTLOC(I1,I2,JATOM),I1=1,3),I2=1,3)
      enddo

!   now we generate the Bravais matrix :
      call errflg(errfn,'Error while making the Bravais matrix in latgen.')
      call latgen
!   we calculate the ROTIJ rotation matrices :
      call errflg(errfn,'Error while making the ROTIJ matrices in rotdef.')
	  call rotdef
!   if necessary, read/write the rotation matrices for equivalent atoms
      call errflg(errfn,'Error at the processing of the ROTIJ Matrices in ReadROTIJ')
      call ReadROTIJ




      RETURN
 5010 FORMAT(A4,23X,I3,/,13X,A4)
 5020 FORMAT(6F10.7)
 5030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)
 5031 FORMAT(/,15x,I2)
 5040 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)
 5050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
 5060 FORMAT(20X,3F10.8)
 6000 FORMAT(///,3X,'ERROR IN EBND  : MULT(JATOM)=0 ...', &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)

  930 WRITE (ERRMSG,9050) JATOM, MULT(JATOM), INDEX
      CALL OUTERR('INILPW',ERRMSG)
      CALL OUTERR('INILPW','MULT .EQ. 0')
  999 RETURN
 9050 FORMAT('MULT(',I3,'=',I3,', INDEX=',I3)

      END
















