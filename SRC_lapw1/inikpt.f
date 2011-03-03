      SUBROUTINE INIKPT(INFO)
!
      use lapw_timer, only: START_TIMER, STOP_TIMER, time_setwar, &
                            READ_CPU_TIME
      use atspdt, only  : INS
      IMPLICIT NONE
      INCLUDE 'param.inc'
!
!        Arguments
!
      INTEGER            INFO
!
!     ..................................................................
! 1.     PROGRAM UNIT 'INIKPT'
!           Setup data values (valid for all k-points) not covered
!           by 'INILPW' (e.g. warpin).
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           Perform calculations valid for all k-points, therefore
!           initializing the node processors in the parallel version
!           of 'LAPW1'.
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
!           INFO   - INTEGER value                              (output)
!                    returns the error condition:
!                    INFO .EQ. 0 ... no error
!                    INFO .NE. 0 ... currently not supported
!
!        USED SUBROUTINES (DIRECTLY CALLED)
!           CPUTIM - measure used CPU-time
!           SETWAR - setup interstitial warpin factors
!
!        INDIRECTLY CALLED SUBROUTINES
!           STERN  - generate the star of a reciprocal lattice vector
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           INILPW - preliminalry calculations (set COMMON blocks)
!           CALKPT - calculate Eigenvalues/-vectors for one k-point
!           SETKPT - set current k-point
!           PRTKPT - print out/save the results obtained by 'CALKPT'
!
!        INPUT/OUTPUT (READ/WRITE)
!           none
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           COMPLEX*16 declaration for GFAC is used.
!
! 4.     REMARKS
!           The distinction between subroutines 'SETKPT' and 'INIKPT'
!           originates from the parallel version of 'LAPW1'. Although
!           the calculations done in INIKPT are the same for all
!           k-points it is more appropriate to calculate INIKPT on all
!           nodes, because otherwise huge amounts of data would have
!           to be communicated to every node.
!           Therefore the original intention of INIKPT is to setup
!           a node processor with calculations valid for all k-points.
!           INIKPT is performed only once per node.
!           The separation from INILPW would not be necessary for the
!           sequential version.
!           
!
! 5.     METHOD
!           If interstitial warpin should be applied (INS .EQ. 0) then
!           call a routine which sets up warpin factors. Additionally
!           the used CPU-time during executing 'INIKPT' is estimated.
!
!        INSTITUT FUER THEORETISCHE CHEMIE            --  TU VIENNA
!     ..................................................................
!
!
!        External subroutines
!
      EXTERNAL           SETWAR
!
!        setup timing
!
      CALL START_TIMER(time_setwar)
!
      IF (INS .EQ. 0) CALL SETWAR
      INFO = 0
!
!        timing of 'INIKPT'
!
      CALL STOP_TIMER(time_setwar)
!      TIMINI = TIMINI + READ_CPU_TIME(time_setwar)
!      TIMGES = TIMGES + READ_CPU_TIME(time_setwar)
!
      RETURN
!
!        End of 'INIKPT'
!
      END
