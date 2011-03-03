      PROGRAM LAPW1
      use parallel, only : init_parallel, close_parallel, ictxtall, myid
      use lapw_timer, only : init_all_timer, start_timer, time_total, &
                             stop_timer, read_wall_time, read_cpu_time
      use cut1, only  : CUT

!
      IMPLICIT NONE
!
      CHARACTER*5        MNAM
      PARAMETER          (MNAM = 'LAPW1')
!
!        Local Scalars
!
      logical            nmat_only
      INTEGER            NUMWK, INFO, KPOINT
      CHARACTER*67       ERRMSG
      CHARACTER*80       DEFFN, ERRFN
      DOUBLE PRECISION   WTIME(3)
!
!        External Subroutines
!
      EXTERNAL           CALKPT, ERRCLR, ERRFLG, GTFNAM, INIKPT
      EXTERNAL           INILPW, OUTERR, PRTKPT, PRTRES, SETKPT

!      EXTERNAL           SLBOOT, SLTIMER, SLCOMBINE, SLENABLE
!      EXTERNAL           DCPUTIME00
!      DOUBLE PRECISION   DCPUTIME00
!
      cut=.false.
!
!      write(*,*) "In lapw1"
      CALL GTFNAM(DEFFN,ERRFN)
!      CALL INIT_PARALLEL
      CALL INIT_ALL_TIMER
      CALL START_TIMER(time_total)
      CALL ERRFLG(ERRFN,'Error in LAPW1')

!
      CALL INILPW(NUMWK,DEFFN,INFO,nmat_only)
      IF (INFO .NE. 0) THEN
         CALL OUTERR(MNAM,'INILPW aborted unsuccessfully.')
         GOTO 999
      ENDIF
      CALL INIKPT(INFO)
      IF (INFO .NE. 0) THEN
         CALL OUTERR(MNAM,'INIKPT aborted unsuccessfully.')
         GOTO 999
      ENDIF
      DO 10 KPOINT = 1, NUMWK
         CALL SETKPT(KPOINT)
         CALL CALKPT(nmat_only)
         if(nmat_only) then
          if (myid == 0) CALL ERRCLR(ERRFN)
          stop 'NMAT_ONLY'
         endif
         CALL PRTKPT(KPOINT,INFO)
         IF (INFO .NE. 0) THEN
!
!        couldn't process results of work-piece
!
            WRITE(ERRMSG,1000) KPOINT
            CALL OUTERR(MNAM,ERRMSG)
            GOTO 999
         ENDIF
   10 CONTINUE
!
!        print summary
!
      CALL PRTRES(NUMWK,INFO)
      IF (INFO .NE. 0) THEN
         CALL OUTERR(MNAM,'couldn''t process overall results.')
         GOTO 999
      ENDIF
      CALL STOP_TIMER(time_total)
      WRITE(6,*) '     Maximum WALL clock time: ', &
                     READ_WALL_TIME(time_total)
      WRITE(6,*) '     Maximum CPU time:        ', &
                     READ_CPU_TIME(time_total)
      if (myid == 0) CALL ERRCLR(ERRFN)
      CALL CLOSE_PARALLEL
!
!      CALL ERRCLR(ERRFN)
      STOP ' LAPW1 END'
!
  999 STOP 'LAPW1 - Error'
!
 1000 FORMAT ('couldn''t process results of work-piece ',I5,'.')
!
!        End of 'LAPW1'
!
      END
