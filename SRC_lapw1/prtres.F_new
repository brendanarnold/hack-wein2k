      SUBROUTINE PRTRES(KPOINT,INFO)
!
!      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      use lapw_timer, only : READ_WALL_TIME, READ_CPU_TIME, &
                       time_setwar, time_g1, time_hamilt, &
                       time_hns, time_diag, time_pkp
      use cut1, only  : CUT
      use parallel, only: MYID

      IMPLICIT NONE
      INCLUDE 'param.inc'
!
!        Arguments
!
      INTEGER            KPOINT, INFO
      Double Precision r1, r2, r3
!
!     ..................................................................
!
! 1.     PROGRAM UNIT 'PRTRES'
!           Process (print) overall results (taking all K-points into
!           account).
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           This routine is intended to perform calculations where the
!           computed results of all K-points have to be taken into
!           account.
!           Currently only the number of K-points and overall timing
!           results are printed out.
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
!           KPOINT - INTEGER value                               (input)
!                    the index of the K-point to process
!                    constraint: (1 .LE. KPOINT .LE. MAXKPT)
!                                this constraint wont be checked
!
!        USED SUBROUTINES (DIRECTLY CALLED)
!           none
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
!           SETKPT - set current k-point
!
!        INPUT/OUTPUT (READ/WRITE)
!           'PRTKPT' writes to the following files:
!              - output-file
!              - SCF-file (history-file)
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           none
!
! 4.     REMARKS
!           No timing is performed in this routine, because there is
!           very little done here. The consumed CPU-time can be
!           neglected.
!
! 5.     METHOD
!           - write the number of processed K-points to the output- and
!             SCF-file
!           - write the overall used CPU-time to the output-file
!
!        INSTITUT FUER THEORETISCHE CHEMIE            --  TU VIENNA
!     ..................................................................
!
!        ARGUMENT-DESCRIPTION
!
!        External Subroutines
!
      EXTERNAL           CPUTIM
!
      IF(myid.eq.0) WRITE(21,*) '      NUMBER OF K-POINTS:',KPOINT
      IF(myid.eq.0) WRITE(6,*)  '      NUMBER OF K-POINTS:',KPOINT
      if(cut) then
      IF(myid.eq.0) WRITE(21,99) 
 99   format(':WARN :      WARNING: RKmax reduced due to NMATMAX')
      IF(myid.eq.0) WRITE(6,*)  '      WARNING: RKmax reduced due to NMATMAX'
      endif
!
!        timing results
!
!!$      IF(myid.eq.0) WRITE(6,6000) READ_CPU_TIME(time_setwar) +  &
!!$                      READ_CPU_TIME(time_g1) +  &
!!$                      READ_CPU_TIME(time_hamilt) +  &
!!$                      READ_CPU_TIME(time_hns) +  &
!!$                      READ_CPU_TIME(time_diag) +  &
!!$                      READ_CPU_TIME(time_pkp),  &
!!$                    READ_CPU_TIME(time_setwar) + &
!!$                      READ_CPU_TIME(time_g1), &
!!$                    READ_CPU_TIME(time_hamilt) +  &
!!$                      READ_CPU_TIME(time_hns) +  &
!!$                      READ_CPU_TIME(time_diag) +  &
!!$                      READ_CPU_TIME(time_pkp) 
      r1= READ_CPU_TIME(time_setwar) +  READ_CPU_TIME(time_g1) + &
           & READ_CPU_TIME(time_hamilt) +  READ_CPU_TIME(time_hns) + &
           & READ_CPU_TIME(time_diag) + READ_CPU_TIME(time_pkp)
      r2 = READ_CPU_TIME(time_setwar) + READ_CPU_TIME(time_g1)
      r3 = READ_CPU_TIME(time_hamilt) + READ_CPU_TIME(time_hns) + &
           & READ_CPU_TIME(time_diag) +  READ_CPU_TIME(time_pkp)  
      If(myid.Eq.0) Write(6,6000) r1, r2, r3
!!$      IF(myid.eq.0) WRITE(6,6010) READ_WALL_TIME(time_setwar) +  &
!!$                      READ_WALL_TIME(time_g1) +  &
!!$                      READ_WALL_TIME(time_hamilt) +  &
!!$                      READ_WALL_TIME(time_hns) +  &
!!$                      READ_WALL_TIME(time_diag) +  &
!!$                      READ_WALL_TIME(time_pkp),  &
!!$                    READ_WALL_TIME(time_setwar) + &
!!$                      READ_WALL_TIME(time_g1), &
!!$                    READ_WALL_TIME(time_hamilt) +  &
!!$                      READ_WALL_TIME(time_hns) +  &
!!$                      READ_WALL_TIME(time_diag) +  &
!!$                      READ_WALL_TIME(time_pkp) 

 r1 = READ_WALL_TIME(time_setwar) +   READ_WALL_TIME(time_g1) + &
      & READ_WALL_TIME(time_hamilt) +   READ_WALL_TIME(time_hns) +   &
      & READ_WALL_TIME(time_diag) +  READ_WALL_TIME(time_pkp) 
 r2 = READ_WALL_TIME(time_setwar) +  READ_WALL_TIME(time_g1)
 r3 =  READ_WALL_TIME(time_hamilt) +  READ_WALL_TIME(time_hns) + &
      & READ_WALL_TIME(time_diag) +  READ_WALL_TIME(time_pkp)  
 If(myid.Eq.0) Write(6,6010)r1,r2,r3
      INFO = 0
!
      RETURN
!
 6010 FORMAT(3X,   '> SUM OF WALL CLOCK TIMES: ',F8.1,' (INIT = ',F8.1, &
            & ' + K-POINTS = ',F8.1,')')
 6000 FORMAT(3X,'===> TOTAL CPU       TIME: ',F8.1, &
            & ' (INIT = ',F8.1, &
            & ' + K-POINTS = ',F8.1, ')' )
!
!        End of 'PRTRES'
!
      END
