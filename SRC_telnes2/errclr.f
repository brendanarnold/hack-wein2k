!BOP
! !ROUTINE: ErrClr
! !INTERFACE:
      SUBROUTINE ERRCLR(FNAME)
! !INPUT/OUTPUT PARAMETERS:
!           FNAME  - CHARACTER*(*) string                        (input)
!                    The name of the file acting as error-flag.
! !DESCRIPTION:
! 1.     PURPOSE
!           Clears the contents of a file indicating an error condition,
!           remove the error flagging. ERRCLR works in conjunction with
!           the routine ERRFLG.
!
! 2.     USAGE
!           CALL ERRCLR('lapw1.error')
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           ERRFLG - create an errorflag-file (notify error conditions)
!
!        INPUT/OUTPUT (READ/WRITE)
!           The contents of the file given as argument FNAME is deleted.
!           File FNAME is created if not existing and otherwise
!           overwritten.
!
!        MACHINENDEPENDENT PROGRAMPARTS
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
!           - close the errorflag-file (it was left opened by 'ERRFLG')
!           - open the errorflag-file
!           - clear the contents of the file by writing an end-of-file
!             marker at the beginning
!           - close the errorflag-file
!           - exit leaving the file opened

! !REVISION HISTORY:
!     24. August 1993  -  Version 1.01  -  INSTITUT FUER TECHNISCHE ELEKTROCHEMIE  --  TU VIENNA
!     Modified November 2004 (Kevin Jorissen)
!EOP



      CHARACTER*(*),intent(in) ::      FNAME

      CLOSE (99)
      OPEN (99,FILE=FNAME,STATUS='REPLACE')

      RETURN
      END
