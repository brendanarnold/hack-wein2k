! This file uses LAPACK routines 
!
      SUBROUTINE SECLR4
!
#ifdef Parallel
      use matrices, only : H, S, HSROWS, Z, EIGVAL
#else 
      use matrices, only : HS, HSDIAG, HSROWS, Z, EIGVAL
#endif
      use parallel, only : DESCHS, DESCZ, NPROW, NPCOL, BARRIER,  &
                           BLOCKSIZE, ICTXTALL
      use lapw_timer, only : READ_CPU_TIME, READ_WALL_TIME,  &
                             START_TIMER, STOP_TIMER, &
                             time_cholesky, time_transform, time_sep,  &
                             time_backtransform
      use lolog, only : nlo, loor, lapw, ilo
      use comi, only    : NVAA, NE, NBELW, NAT, NT, LNSMAX
      use coml, only    : SPRSWF, PRNTWF, WFTAPE, REL, ITER, NOHNS
      use comr, only    : RKM, EL, EU, WEIGHT
      IMPLICIT NONE
      INCLUDE 'param.inc'
!
!     ..................................................................
!
!        interface between matrix set up in LAPW1 and LAPACK (BLAS) matrix
!        subroutines; DOUBLE PRECISION version
!
!        Solve the secular equations H*Z = EIGVAL*S*Z
!        H and S are stored in the array HS,
!        H is stored below the diagonal, whereas S is stored above
!        the diagonal. The diagonal of H is stored in the 
!        one dimensional array HSDIAG.
!        NVAA is the dimension of the secular matrices
!        there is a limit on NE because of the dimensions used
!        NS is the total space available for the eigenvectors Z
!        NE is returned as the number of Eigenvalues between EL1 and EU1
!        except NE is limited to (NS/NVAA)
!        EL1 and EU1 are input as EIGVAL(1) and EIGVAL(2).
!        IF EIGVAL(2) .LT. EIGVAL(1) THEN no bounds are applied
!
!                                              D.D.Koelling   April,1974
!
!     ..................................................................
!
!
!
!        Locals
!
      INTEGER            IER1, IER2, NE1, NMAT1, NSPC
      INTEGER            I, BLOCK1, BLOCK2
      DOUBLE PRECISION   EL1, EU1, Q
      DOUBLE PRECISION   Flop
      DOUBLE PRECISION   WTIME(4)
      INTEGER            STDERR
      PARAMETER          (STDERR = 99)
!
!        External Functions
!
!      INTEGER            NEBLW1
!      EXTERNAL           NEBLW1
!
!        Local Parameters
!
      INTEGER            NBMIN, NBVAL, NXVAL
      PARAMETER          (NBMIN = 1, NBVAL = 2, NXVAL = 1)
!
!        Locals
!
      CHARACTER*67       ERRMSG
      CHARACTER          JOBZ, RANG, UPLO
      INTEGER            IL, INFO, IU, J
      INTEGER, allocatable :: IFAIL(:), IWORK(:)
      INTEGER            LWORK, LIWORK, LRWORK
      INTEGER            NEV
      INTEGER, allocatable :: ICLUSTR(:)
      DOUBLE PRECISION, allocatable :: GAP(:)
      DOUBLE PRECISION   SCAL
      INTEGER            NN, NP0, MQ0, NQ0, CLUSTERSIZE
      INTEGER            IERR
      DOUBLE PRECISION   ABSTOL, ORFAC
!_REAL      DOUBLE PRECISION, allocatable ::   WORK(:)
!_COMPLEX      DOUBLE COMPLEX, allocatable ::  WORK(:)
!_COMPLEX      DOUBLE PRECISION, allocatable ::   RWORK(:)

      DOUBLE PRECISION   CP(5)

      logical out,out1
      data out/.false./
      data out1/.false./

!
!        External Subroutines
!
!_REAL      EXTERNAL           OUTERR, DPOTRF, DSCGST, DSYEVX2
!_COMPLEX      EXTERNAL           OUTERR, ZPOTRF, ZHCGST, ZHEEVX2
#ifdef Parallel
      INTEGER NUMROC, ICEIL
      EXTERNAL NUMROC, ICEIL
#endif
!
!        Intrinsic Functions
!
      INTRINSIC          MIN
!
!
      read(5,'(L1)',err=111,end=111) out1
      if(out1)  out=out1
      print*, 'no out-error'
 111  continue
!      out=.true.
!      out=.false.

      CALL BARRIER
#ifdef Parallel
      NN = MAX(HSROWS, BLOCKSIZE, 2)
      NP0 = NUMROC(NN,BLOCKSIZE, 0, 0, NPROW)
      MQ0 = NUMROC(MAX(NUME, BLOCKSIZE, 2), BLOCKSIZE, 0, 0, NPCOL)
      NQ0 = NUMROC( NN, BLOCKSIZE, 0, 0, NPCOL )
      CLUSTERSIZE = 200
! LRWORK from Documentation
!      LRWORK = 4*HSROWS + MAX( 5*NN, NP0*MQ0) + &
!              ICEIL(NUME,NPROW*NPCOL)*NN + (CLUSTERSIZE-1) * HSROWS
! LRWORK from PZHEEVX
!_COMPLEX      LRWORK = 4*HSROWS + &
!_COMPLEX              MAX( ICEIL( NUME, NPROW*NPCOL )*NN + MAX( 5*NN, NP0*MQ0 ), &
!_COMPLEX                   MAX( ( BLOCKSIZE*( BLOCKSIZE-1 ) ) / 2, &
!_COMPLEX                        ( MQ0+NP0 )*BLOCKSIZE ) + BLOCKSIZE*BLOCKSIZE) + &
!_COMPLEX              (CLUSTERSIZE-1) * HSROWS

!_COMPLEX      LWORK = HSROWS + MAX(BLOCKSIZE*(NP0+1),3) + CLUSTERSIZE * HSROWS
!_COMPLEX      LWORK = HSROWS + (NP0 + NQ0 + BLOCKSIZE) * BLOCKSIZE + &
!_COMPLEX              CLUSTERSIZE * HSROWS
!_REAL      LWORK = 5*HSROWS + MAX(5*NN, NP0*MQ0 + 2*BLOCKSIZE*BLOCKSIZE) +  &
!_REAL              ICEIL(NUME, NPROW*NPCOL)*NN + CLUSTERSIZE * HSROWS
      LIWORK = 6*MAX(HSROWS, NPROW*NPCOL+1,4)
      allocate( ICLUSTR(2*NPROW*NPCOL),GAP(NPROW*NPCOL) )
#else
      LWORK = HSROWS*(BLOCKSIZE+3)
      LIWORK = 5*HSROWS
      LRWORK = 7*HSROWS
#endif
!_COMPLEX      allocate( RWORK(LRWORK) )
      allocate( IFAIL(HSROWS), IWORK(LIWORK), WORK(LWORK), STAT=IERR )
      if(IERR.ne.0) THEN 
        write(6,*) 'not enough memory: ', IERR
        stop 'not enough memory' 
      end if

      IF (EIGVAL(1) .GE. EIGVAL(2)) THEN
         EIGVAL(1) = -1.0D+70
         EIGVAL(2) =  1.0D+70
      ENDIF
      EL1 = EIGVAL(1)
      EU1 = EIGVAL(2)

!	
!      conjugate matrix
!
#ifndef Parallel
!_COMPLEX       DO J = 1, NVAA+NLO
!_COMPLEX         DO I = 1, NVAA+NLO
!_COMPLEX	   HS(I,J) = DCONJG(HS(I,J))
!_COMPLEX	 END DO
!_COMPLEX	 HSDIAG(J) = DCONJG(HSDIAG(J))
!_COMPLEX       END DO
#endif
!
!        form a Cholesky factorization of SP (packed overlap-matrix S).
!
!
!        swap Diagonal of Overlap matrix into matrix
!
!print
#ifndef Parallel
if(out) then
write(87,*) HSrows,nvaa,nlo
write(88,1234) HS
write(87,1234) HSdiag
 1234  format(4e23.15)
endif
#endif
      
#ifndef Parallel
!_REAL      CALL DSWAP(HSROWS,HS,HSROWS+1,HSDIAG,1)
!_COMPLEX      CALL ZSWAP(HSROWS,HS,HSROWS+1,HSDIAG,1)
#endif
!
!      write (6,*)'matrix size',NVAA+NLO
      CALL START_TIMER(time_cholesky)
#ifdef Parallel
      if(deschs(2).gt.0) then
!_REAL      CALL PDPOTRF('U',NVAA+NLO,S,1,1,DESCHS,INFO)
!_COMPLEX      CALL PZPOTRF('U',NVAA+NLO,S,1,1,DESCHS,INFO)
      end if
#else
!_REAL      CALL DPOTRF('U',NVAA+NLO,HS,HSROWS,INFO)
!_COMPLEX      CALL ZPOTRF('U',NVAA+NLO,HS,HSROWS,INFO)
#endif
      CALL STOP_TIMER(time_cholesky)
      IF(INFO .NE. 0) THEN
         write(6,*) 'Cholesky INFO = ', INFO
         write(stderr,*) 'Cholesky INFO = ', INFO
         CALL OUTERR('SECLR4','POTRF (Scalapack/LAPACK) failed.')
         STOP 'SECLR4 - Error'
      ENDIF

!
!        swap Diagonal of Hamilton matrix into matrix
!      
#ifndef Parallel
!_REAL      CALL DSWAP(HSROWS,HS,HSROWS+1,HSDIAG,1)
!_COMPLEX      CALL ZSWAP(HSROWS,HS,HSROWS+1,HSDIAG,1)
#endif
      CALL BARRIER
!
!        transform the problem to a standard eigenvalue problem.
!
      CALL START_TIMER(time_transform)
#ifdef Parallel
      if(deschs(2).gt.0) then
!_REAL      CALL PDSYGST(1,'U',NVAA+NLO,H,1,1,DESCHS,S,1,1,DESCHS,SCAL,INFO)
!_COMPLEX      CALL PZHEGST(1,'U',NVAA+NLO,H,1,1,DESCHS,S,1,1,DESCHS,SCAL,INFO)
      end if
#else
!_REAL      CALL DSCGST(NVAA+NLO,HS,HSROWS,HSDIAG,WORK,INFO)
!_COMPLEX      CALL ZHCGST(NVAA+NLO,HS,HSROWS,HSDIAG,WORK,INFO)
#endif

      CALL STOP_TIMER(time_transform)
      IF(INFO .NE. 0) THEN
         write(6,*) 'GST INFO = ', INFO
         write(stderr,*) 'GST INFO = ', INFO
         CALL OUTERR('SECLR4','SCGST (Scalapack/LAPACK) failed.')
         STOP 'SECLR4 - Error'
      ENDIF
      CALL BARRIER
!
!        compute the eigenvalues within energy-range EL1, EU1
!
      CALL START_TIMER(time_sep)
      IL = 0
      IU = 0
      ABSTOL = 0.0D+0
      ORFAC = 1.0D-3
      ORFAC = 1.0D-2
!_REAL      BLOCK1 = BLOCKSIZE
!_REAL      BLOCK2 = 12      
!_REAL      IF(HSROWS.LT.2500) BLOCK2=8      
!_REAL      IF(HSROWS.LT.1500) BLOCK2=6      
!_REAL      IF(HSROWS.LT.1200) BLOCK2=4      
!_REAL      IF(HSROWS.LT.1000) BLOCK2=1      
!_REAL      BLOCK2 = 1
#ifdef Parallel
            JOBZ='V'
	    RANG='V'
	    UPLO='U'
      if(deschs(2).gt.0) then
!C_REAL      CALL PDSYEVX_1_7(JOBZ,RANG,UPLO,NVAA+NLO,H,1,1,
!C_REAL      CALL PDSYEVX2(JOBZ,RANG,UPLO,NVAA+NLO,H,1,1,
!_REAL      CALL PDSYEVX16(JOBZ,RANG,UPLO,NVAA+NLO,H,1,1, &
!_REAL                  DESCHS,EL1,EU1,IL,IU, &
!_REAL                  ABSTOL,NE,NEV,EIGVAL,ORFAC,Z,1,1,DESCZ, &
!_REAL                  WORK,LWORK,IWORK,LIWORK,IFAIL,ICLUSTR,GAP,INFO)
!_COMPLEX      CALL PZHEEVX16(JOBZ,RANG,UPLO,NVAA+NLO,H,1,1, &
!_COMPLEX                  DESCHS,EL1,EU1,IL,IU, &
!_COMPLEX                  ABSTOL,NE,NEV,EIGVAL,ORFAC,Z,1,1,DESCZ, &
!_COMPLEX                  WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,IFAIL, &
!_COMPLEX                  ICLUSTR,GAP,INFO)
      end if
      IF (INFO .NE. 0) THEN
         write(6,*) 'SEP INFO = ', INFO
	 if(info.eq.2) then
 99        format(':WARN :      WARNING: Not all eigenvectors are orthogonal') 
!           WRITE(21,99)
           WRITE(6,99)
!	   write(6,*) 'Numbers of eigenvectors that failed to converge =', IFAIL
	   write(6,*) 'Eigenvalue clusters in eigensolver =', ICLUSTR
	   write(6,*) 'Gap in eigensolver =', Gap
           write(6,*) '  Increase CLUSTERSIZE in SECLR4'
	 else
           write(stderr,*) 'SEP INFO = ', INFO
           CALL OUTERR('SECLR4','SYEVX (Scalapack/LAPACK) failed.')
           IF(INFO.EQ.4) THEN
             write(6,*) 'NEV, NE = ', NEV, NE
             write(6,*) 'ICLUSTR = ', ICLUSTR
             write(6,*) 'WORK(1) = ', WORK(1:2), LWORK
!_COMPLEX           write(6,*) 'RWORK(1) = ', RWORK(1:2), LRWORK
             write(6,*) 'IWORK(1) = ', IWORK(1:2), LIWORK
           END IF
           STOP 'SECLR4 - Error'
	 end if
      ENDIF
#else
!_REAL      CALL DSYXEV4('V','V','L',NVAA+NLO,HS,HSROWS,EL1,EU1,IL,IU, &
!_REAL                  ABSTOL,NE,EIGVAL,Z,HSROWS,WORK,LWORK,IWORK,IFAIL, &
!_REAL                  BLOCK1,BLOCK2,0.0D0,NUME,NBELW,INFO)
!_COMPLEX!      CALL ZHEEVX('V','V','L',NVAA+NLO,HS,HSROWS,EL1,EU1,IL,IU,ABSTOL,
!_COMPLEX!     ,            NE,EIGVAL,Z,HSROWS,WORK,LWORK,RWORK,IWORK,IFAIL,
!_COMPLEX!     ,            INFO)
!_COMPLEX      CALL ZHEEVX2('V','V','L',NVAA+NLO,HS,HSROWS,EL1,EU1,IL,IU,ABSTOL, &
!_COMPLEX                  NE,EIGVAL,Z,HSROWS,WORK,LWORK,RWORK,IWORK,IFAIL, &
!_COMPLEX                  BLOCKSIZE, NUME, NBELW, INFO)
!     compare: Complex routines need additional parameter RWORK (after WORK)
      IF (INFO .NE. 0) THEN
         write(6,*) 'SEP INFO = ', INFO
         write(stderr,*) 'SEP INFO = ', INFO
         CALL OUTERR('SECLR4','SYEVX (Scalapack/LAPACK) failed.')
         STOP 'SECLR4 - Error'
      ENDIF
      write(6,*) NE, ' Eigenvalues computed '
      IF (NE .GT. NUME) THEN
         CALL OUTERR('SECLR4','More than NUME eigenvalues.')
         WRITE (ERRMSG,9000) NUME,NE
 9000    FORMAT ('Increase NUME from',i4,' to',i4,'   OR')
         CALL OUTERR('SECLR4',ERRMSG)
         CALL OUTERR('SECLR4','reduce E2 in *.klist.(see *.output1)')
         WRITE(6,9001) (eigval(j),j=1,ne)
 9001    FORMAT (5f15.6)
         STOP 'SECLR4 - Error'
      ENDIF
#endif
      CALL STOP_TIMER(time_sep)
      CALL BARRIER
!
!        swap Diagonal of Overlap matrix into matrix
!      
#ifndef Parallel
!_REAL      CALL DSWAP(HSROWS,HS,HSROWS+1,HSDIAG,1)
!_COMPLEX      CALL ZSWAP(HSROWS,HS,HSROWS+1,HSDIAG,1)
!
!	
!      conjugate overlap matrix
!
!_COMPLEX       DO J = 1, NVAA+NLO
!_COMPLEX         DO I = 1, J
!_COMPLEX	   HS(I,J) = DCONJG(HS(I,J))
!_COMPLEX	 END DO
!_COMPLEX       END DO
#endif
!
!        backtransform eigenvectors to the original problem:
!           for H*x=(lambda)*S*x;
!           backtransform eigenvectors: x = inv(L)´*y or inv(U)*y
!
      CALL START_TIMER(time_backtransform)
#ifdef Parallel
      if(deschs(2).gt.0) then
!_REAL      CALL PDTRSM('L', 'U', 'N', 'N', NVAA+NLO, NE, 1.0D0,  &
!_REAL                 S, 1,1,DESCHS, Z, 1,1,DESCZ)
!_COMPLEX      CALL PZTRSM('L', 'U', 'N', 'N', NVAA+NLO, NE, (1.0D0,0.0D0),  &
!_COMPLEX                 S, 1,1,DESCHS, Z, 1,1,DESCZ)
      end if
#else
!_REAL      CALL DTRSM('L', 'U', 'N', 'N', NVAA+NLO, NE, 1.0D0,  &
!_REAL                 HS, HSROWS, Z, HSROWS)
!_COMPLEX      CALL ZTRSM('L', 'U', 'N', 'N', NVAA+NLO, NE, (1.0D0,0.0D0), &
!_COMPLEX                 HS, HSROWS, Z, HSROWS)
#endif

      CALL STOP_TIMER(time_backtransform)
      CALL BARRIER
#ifdef Parallel
      deallocate( ICLUSTR,GAP )
#endif
      deallocate( IFAIL, IWORK )
      deallocate( WORK )
!_COMPLEX      deallocate( RWORK )
      CALL BARRIER

!print
#ifndef Parallel
if(out) then
do i=1,ne
write(86,*) eigval(i)
write(89,1234) (Z(j,i),j=1,hsrows)
enddo
endif 
#endif

!
         Flop = DBLE(NVAA+NLO)*(NVAA+NLO)*(NVAA+NLO)/3.0D0
!_COMPLEX         Flop = Flop*4.0D0
         if(READ_CPU_TIME(time_cholesky).eq.0.d0) then
         WRITE (6,1001) 'Cholesky complete (CPU)' ,  &
                        READ_CPU_TIME(time_cholesky)
         else
         WRITE (6,1001) 'Cholesky complete (CPU)' ,  &
                        READ_CPU_TIME(time_cholesky) ,  &
                        1E-6* Flop/READ_CPU_TIME(time_cholesky)
         endif
#ifdef Parallel
         if(READ_WALL_TIME(time_cholesky).eq.0.d0) then
         WRITE (6,1001) 'Cholesky complete (WALL)' ,  &
                        READ_WALL_TIME(time_cholesky)
         else
         WRITE (6,1001) 'Cholesky complete (WALL)' ,  &
                        READ_WALL_TIME(time_cholesky) ,  &
                        1E-6* Flop/READ_WALL_TIME(time_cholesky)
         endif
#endif
         Flop = 3.0D0*DBLE(NVAA+NLO)*(NVAA+NLO)*(NVAA+NLO)/3.0D0
!_COMPLEX         Flop = Flop*4.0D0
         if(READ_CPU_TIME(time_transform).eq.0.d0) then
         WRITE (6,1001) 'Transform to eig.problem (CPU)' ,  &
                        READ_CPU_TIME(time_transform) 
         else
         WRITE (6,1001) 'Transform to eig.problem (CPU)' ,  &
                        READ_CPU_TIME(time_transform) ,  &
                        1E-6* Flop/READ_CPU_TIME(time_transform)
         endif
#ifdef Parallel
         if(READ_WALL_TIME(time_transform).eq.0.d0) then
         WRITE (6,1001) 'Transform to eig.problem (WALL)' ,  &
                        READ_WALL_TIME(time_transform) 
         else
         WRITE (6,1001) 'Transform to eig.problem (WALL)' ,  &
                        READ_WALL_TIME(time_transform) ,  &
                        1E-6* Flop/READ_WALL_TIME(time_transform)
         endif
#endif
         Flop = 4.0D0*DBLE(NVAA+NLO)*(NVAA+NLO)*(NVAA+NLO)/3.0D0
!_COMPLEX         Flop = Flop*4.0D0
         if(READ_CPU_TIME(time_sep).eq.0.d0) then
         WRITE (6,1001) 'Compute eigenvalues (CPU)' ,  &
                        READ_CPU_TIME(time_sep)
         else
         WRITE (6,1001) 'Compute eigenvalues (CPU)' ,  &
                        READ_CPU_TIME(time_sep),  &
                        1E-6* Flop/READ_CPU_TIME(time_sep)
         endif
#ifdef Parallel
         if(READ_WALL_TIME(time_sep).eq.0.d0) then
         WRITE (6,1001) 'Compute eigenvalues (WALL)' ,  &
                        READ_WALL_TIME(time_sep)
         else
         WRITE (6,1001) 'Compute eigenvalues (WALL)' ,  &
                        READ_WALL_TIME(time_sep) ,  &
                        1E-6* Flop/READ_WALL_TIME(time_sep)
         endif
#endif
         Flop = DBLE(NVAA+NLO)*(NVAA+NLO)*NE/4.0D0
!_COMPLEX         Flop = Flop*4.0D0
         if(READ_CPU_TIME(time_backtransform).eq.0.d0) then
         WRITE (6,1001) 'Backtransform (CPU)' ,  &
                        READ_CPU_TIME(time_backtransform) 
         else
         WRITE (6,1001) 'Backtransform (CPU)' ,  &
                        READ_CPU_TIME(time_backtransform) ,  &
                        1E-6* Flop/READ_CPU_TIME(time_backtransform)
         endif
#ifdef Parallel
         if(READ_WALL_TIME(time_backtransform).eq.0.d0) then
         WRITE (6,1001) 'Backtransform (WALL)' ,  &
                        READ_WALL_TIME(time_backtransform)
         else
         WRITE (6,1001) 'Backtransform (WALL)' ,  &
                        READ_WALL_TIME(time_backtransform) ,  &
                        1E-6* Flop/READ_WALL_TIME(time_backtransform)
         endif
#endif
 1001    FORMAT (1X,'Seclr4(',A,') :',t45,f10.3:t58,f10.3,' Mflops')

!
      RETURN
!
 6000 FORMAT (/,/,/,10X,'ERROR IN SCLRE2, IER1=',I5,'#E FOUND=',I5, &
              ' IER2=',I5)
!
!        End of 'SECLR4'
!
      END

