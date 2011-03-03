!BOP
! !ROUTINE: ReadROTIJ
! !INTERFACE:
      SUBROUTINE ReadROTIJ
! !USES:
      use input, only : natom, neqatdeb, neqat
	  use rotation_matrices, only : rotij
	  use program_control, only : mrot, wrot, verbosity
	  use struct, only : mult
! !DESCRIPTION:
!     Read the ROTIJ Matrices from the file case.rotij
! !REVISION HISTORY:
!     Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none
!   LOCALS :
      INTEGER JR, JC, i


      write(6,'(/,a)') 'Output from subroutine ReadROTIJ:'

!     Sanity check :
!     If we get any funny input for the ATOMS key, we just calculate all atoms.
      if(NEqAt.lt.1.or.NEqAt.gt.mult(natom)) then
	    if(NEqAt.ne.0) &   ! This is the automatic value - let's not blame that on the user
	    write(6,'(a)') ':INFO - resetting the atom summation because of tricky input.'
	    NEqAT=mult(natom)
	  endif
      if(NEqAtDeb.lt.1.or.NeqAtDeb.gt.NEqAt) then
	    write(6,'(a)') ':INFO - resetting the atom summation because of tricky input.'
	    NEqAtDeb=1
	  endif
!     This check could not be done in lectureinnes because it does not know mult yet.


      if(.not.mrot) then
!     read matrices from file :

      do i=1,NEqAt
        READ(31,*)   ! this line is for indexing.
        DO JR=1, 3
          READ(31, 101) ROTIJ(1,JR,i), ROTIJ(2,JR,i), ROTIJ(3,JR,i)
        ENDDO
      enddo

      endif

	  if(wrot.or.verbosity.ge.1) then
	  write(6,'(a)') 'ROTIJ rotation matrices: '
      do i=NEqAtDeb,NEqAt
		 write(6,'(a,i4)') 'for equivalent position ',i
		 do jr=1,3
           write(6,101) (ROTIJ(JC,JR,i), JC=1, 3)
		 enddo
      enddo
      endif


      if (wrot) then
	    rewind(31)
		do i=1,mult(natom)            !  we write some that were possibly not used, but this way, case.rotlm is useful regardless of the ATOM input.
		  write(31,*) 'for equivalent position',i
          do jr=1,3
		    write(31,101) (ROTIJ(JC,JR,i), JC=1, 3)
	      enddo
		enddo
      endif
      
 910  CONTINUE 
 100  FORMAT(24X,I3,8X,I2,7X,I4)
 101  FORMAT(3F10.5)
      
      RETURN
      END



