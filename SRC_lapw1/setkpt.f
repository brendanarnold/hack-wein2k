      SUBROUTINE SETKPT(KPOINT)
!
      use kpts, only    : SX, SY, SZ
      use onekpt, only  : SX1, SY1, SZ1
      IMPLICIT NONE
      INCLUDE 'param.inc'
!
!        Arguments
!
      INTEGER            KPOINT
!
!     ..................................................................
!
! 1.     PROGRAM UNIT 'SETKPT'
!           Set the current K-point.
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           The the current K-point (identified by the argument KPOINT)
!           is put into the associated data structure.
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
!           KPOINT - INTEGER value                               (input)
!                    the index of the K-point to process
!                    constraint: (1 .LE. KPOINT .LE. MAXKPT)
!                                this constraint won't be checked
!
!        USED SUBROUTINES (DIRECTLY CALLED)
!           CPUTIM - measure used CPU-time
!
!        INDIRECTLY CALLED SUBROUTINES
!           none
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           CALKPT - calculate Eigenvalues/-vectors for one k-point
!           INILPW - preliminalry calculations (set COMMON blocks)
!           INIKPT - setup values (valid for all k-points) not covered
!                    by 'INILPW' (e.g. warpin)
!           PRTKPT - print out/save the results obtained by 'CALKPT'
!
!        INPUT/OUTPUT (READ/WRITE)
!           none
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           none
!
! 4.     REMARKS
!           This routine originates from the parallel version of LAPW1.
!           In the parallel version of 'SETKPT' the current K-point is
!           received from the master and stored in the associated
!           COMMON block. To maintain the same program structure as in
!           the parallel version, 'SETKPT' is introduced in the
!           sequential version also.
!
! 5.     METHOD
!           - get the current K-point and store it in the associated
!             COMMON-block
!
!        INSTITUT FUER THEORETISCHE CHEMIE            --  TU VIENNA
!     ..................................................................
!
!        Local Scalars
!
!
!        setup timing of 'SETKPT'
!
!
      SX1 = SX(KPOINT)
      SY1 = SY(KPOINT)
      SZ1 = SZ(KPOINT)
!
!        measure cpu-time used by 'SETKPT'
!
      RETURN
!
!        End of 'SETKPT'
!
      END
