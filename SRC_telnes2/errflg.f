!BOP
! !ROUTINE: ErrFlg
! !INTERFACE:
      SUBROUTINE ERRFLG(FNAME,MSG)
! !INPUT/OUTPUT PARAMETERS:
!           FNAME  - CHARACTER*(*) string                        (input)
!                    The name of the file acting as error-flag.
!           MSG    - CHARACTER*(*) string                        (input)
!                    The (error) message which should be written to the
!                    errorflag-file.
! !DESCRIPTION:
!
! 1.     PURPOSE
!           Because there is no standard (or even semi-standard) way to
!           generate exit codes in FORTRAN 77, this routine writes a
!           non-empty file to the current subdirectory as an indication
!           that some serious error has occured. Other programs can then
!           check the contents of this file to determine whether an
!           error has occured. The errorflag-file is left opened when
!           returning from this routine to enable writing other
!           errormessages to it.
!
! 1.     USAGE
!           CALL ERRFLG('lapw2.error','Error in OUTWIN')
!!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           ERRCLR - clears the contents of a file
!
!        INPUT/OUTPUT (READ/WRITE)
!           A message given as argument MSG is written to a the file
!           given as argument FNAME. File FNAME is created if not
!           existing and otherwise overwritten.
!
!        MACHINEDEPENDENT PROGRAMPARTS
!           none
!
! 3.     REMARKS
!           The best way to use this routine is to call ERRFLG at the
!           start of a program writing some message to the errorflag-
!           file and ERRCLR before a successful exit of the program. By
!           checking the contents of the errorflag-file it is possible
!           to determine if the program was successfully completed.
!
!           This method has the advantage of working even if some
!           runtime-error occurs which is not taken care of in the
!           program.
!
! 4.     METHOD
!           - open errorflag-file
!           - write some message to this errorflag-file
!           - exit leaving the file opened

! !REVISION HISTORY:
!     24. August 1993  -  Version 1.02  -  INSTITUT FUER TECHNISCHE ELEKTROCHEMIE  --  TU VIENNA
!     Modified November 2004 (Kevin Jorissen)
!EOP




      CHARACTER*(*),intent(in) ::      FNAME, MSG

      OPEN (99,FILE=FNAME,ERR=900,STATUS='UNKNOWN')
!      WRITE (99,9000) 'TELNES2 ERROR'
      WRITE (99,9000) MSG
      CLOSE (99)
      OPEN (99,FILE=FNAME,ERR=900,STATUS='UNKNOWN')

      RETURN
!
!        Errors
!
  900 STOP 'ERRFLG - couldn''t open errorflag-file.'
!
 9000 FORMAT (A)
!
!        End of 'ERRFLG'
!
      END
