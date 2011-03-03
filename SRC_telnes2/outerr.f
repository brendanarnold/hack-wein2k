!BOP
! !ROUTINE: Outerr
! !INTERFACE:
      SUBROUTINE OUTERR(SRNAME,ERRMSG)
! !INPUT/OUTPUT PARAMETERS:
!   SRNAME    The name of the routine which called 'OUTERR'.
!   ERRMSG    The (error) message to be printed.
! !DESCRIPTION:
! 1.     PURPOSE
!           OUTERR is a gemeral purpose message output routine for
!           FORTRAN 77 subroutines. A message given as an argument
!           is printed on standard error and execution of the calling
!           routine resumes.
!
! 2.     USAGE
!           SUBROUTINE CALLER
!           ...
!           EXTERNAL OUTERR
!           ...
!           CHARACTER*67 ERRMSG
!           WRITE (ERRMSG,FMT='(''error number:'',I3)') 10
!           CALL OUTERR('CALLER',ERRMSG)
!           ...
!           END
!
!        INPUT/OUTPUT (READ/WRITE)
!           If OUTERR is called by
!              CALL OUTERR('FOO','couldn''t open file.')
!           the resulting (output) error message will look like:
!
!           'FOO' - couldn't open file.
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           The actual value of the standard error device (STDERR)
!           has to be setup for different brand's of computers.
!           For instance: STDERR = 7 on HP9000 (800/700) under HPUX
!                                = 0 on IBM RS6000 under AIX
!                                = 0 on Siemens/Fujitsu S100 under UXP/M
!           
! 3.     REMARKS
!           In the program LAPW1 (part of the WIEN93 software package),
!           error messages are written to a special error-file, accessed
!           as logical filenumber 99.
!
! 4.     METHOD
!           Write error message on output device 'standard error'.

! !REVISION HISTORY:
!     24. August 1993  -  Version 1.21  -  INSTITUT FUER TECHNISCHE ELEKTROCHEMIE  --  TU VIENNA
!     Modified November 2004 (Kevin Jorissen)
!EOP

!   INPUT :
      CHARACTER*(*),intent(in) ::      SRNAME,ERRMSG
!   LOCALS :
      INTEGER,parameter ::            STDERR=99

      WRITE (STDERR,9010) SRNAME, ERRMSG

      RETURN

 9010 FORMAT (a,' - ',a)
      END
