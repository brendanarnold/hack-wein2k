!BOP
! !ROUTINE: Inverse
! !INTERFACE:
      SUBROUTINE Inverse(M, Mout)
! !USES:
	  use program_control,only : verbosity
! !INPUT/OUTPUT PARAMETERS:
!   M    :   3*3 input matrix
!   Mout :   the inverse of M
! !DESCRIPTION:
!   Calculates the inverse of a 3*3 real matrix.
! !REVISION HISTORY:
!     Updated November 2004 (Kevin Jorissen)
!EOP
	  implicit none
      DOUBLE PRECISION M(3,3), Mout(3,3), Det, CoMat(3,3)
      INTEGER I, J

      CALL Determinant(M, Det)
      IF (Det.NE.dble(0)) THEN
         CALL CoMatrix(M, CoMat)
         CALL Transpose(CoMat, Mout)
         Mout=Mout / Det
      ELSE
	     write(6,'(a)') 'Output of subroutine Inverse :'
         WRITE(6,'(a)') 'Determinant null !!'
      ENDIF
         
!     to verify the inverse:
	  if (verbosity.ge.2) then
      CALL ProductMatMat(M, Mout, CoMat)
        WRITE(6,'(a)') 'Identite in subroutine Inverse ? :'
        DO I=1, 3
           WRITE(6,'(3(f8.3,x))') (CoMat(I, J), J=1, 3)
        ENDDO
	  endif ! verbosity
         
      RETURN
      END
      
      
