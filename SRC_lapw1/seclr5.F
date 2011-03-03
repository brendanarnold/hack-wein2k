      SUBROUTINE SECLR5
!
#ifdef Parallel
      use matrices, only : H, S, HSROWS, &
                           HPANEL, SPANEL, Z, EIGVAL
#else 
      use matrices, only : HS, HSDIAG, HSROWS, &
                           HPANEL, SPANEL, Z, EIGVAL
#endif
      use lapw_timer, only : START_TIMER, STOP_TIMER, time_iter,  &
                             READ_CPU_TIME, READ_WALL_TIME
      use parallel, only : MYROWHS, MYCOLHS, MYID, DESCLOCAL, DESCZ, &
                           ICTXTALL
      use lolog, only : nlo, loor, lapw, ilo
      use comi, only    : NV => NVAA, NE, NBLW => NBELW, NAT
      use coml, only    : SPRSWF
      use comr, only    : RKM
      IMPLICIT NONE
      INCLUDE 'param.inc'
#ifdef Parallel
      include 'mpif.h'
#endif
!
!     ..................................................................
!
!        interface between matrix set up in LAPW1 and iterative
!        diagonalization subroutines
!
!        Solve the secular equations H*Z = EIGVAL*S*Z
!        NVAA is the dimension of the secular matrices
!        NE is returned as the number of Eigenvalues between EL1 and EU1
!        except NE is limited to (NS/NVAA)
!        EL1 and EU1 are input as EIGVAL(1) and EIGVAL(2).
!        IF EIGVAL(2) .LT. EIGVAL(1) THEN no bounds are applied
!
!
!     ..................................................................
!
!        Parameters
!
      INTEGER            NUMIT
      PARAMETER          (NUMIT = 1)
!
!        Locals
!
      INTEGER            I, J
      INTEGER            IDUMMY, NVOLD, NEOLD
      DOUBLE PRECISION   EL1, EU1
      DOUBLE PRECISION   SXOLD, SYOLD, SZOLD, DDUMMY
      DOUBLE PRECISION   TFIRST, TEND
!_REAL      DOUBLE PRECISION   ZERO
!_COMPLEX      DOUBLE COMPLEX   ZERO
!_REAL      PARAMETER ( ZERO=0.0D0 )
!_COMPLEX      PARAMETER ( ZERO=(0.0D0,0.0D0) )
      INTEGER            J1, J2
      INTEGER            NBELW
      INTEGER            IERR
!_REAL      DOUBLE PRECISION, ALLOCATABLE :: ZI(:)
!_COMPLEX      COMPLEX*16, ALLOCATABLE :: ZI(:)
      DOUBLE PRECISION :: EIGVALI
!
!
      CHARACTER*10       C10DUM
      LOGICAL            FIRST,out,out1
      save first
      data out/.false./
      data out1/.false./
!
!        Intrinsics
!
       INTRINSIC         DCONJG
!
!        External Subroutines
!
!      EXTERNAL           SECLIT
!
!        Data statements
!
      DATA FIRST /.TRUE./
      
      read(5,'(L1)',err=111,end=111) out1
      if(out1)  out=out1
      print*, 'no out-error'
 111  continue
!      out=.true.
!      out=.false.
!
      CALL START_TIMER(time_iter)
      IF (EIGVAL(1) .GE. EIGVAL(2)) THEN
         EIGVAL(1) = -1.0D+70
         EIGVAL(2) =  1.0D+70
      ENDIF
!
      EL1 = EIGVAL(1)
      EU1 = EIGVAL(2)
!
      IF (FIRST .and. MYID.eq.0) THEN
!
!        skip vector-file header
!
         DO 10 I = 1, NAT
            READ (98,ERR=999,end=999) DDUMMY
            READ (98) DDUMMY
   10    CONTINUE
         FIRST = .FALSE.
      ENDIF
!
!        read eigenvectors of the previous LAPW-iteration
!
      IF(MYID.eq.0) THEN
        READ (98,ERR=999,end=999) SXOLD,SYOLD,SZOLD,C10DUM,NVOLD,NEOLD
        READ (98) IDUMMY
      END IF
#ifdef Parallel
      CALL MPI_BCAST(SXOLD, 1, MPI_DOUBLE_PRECISION, 0,  &
                     MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(SYOLD, 1, MPI_DOUBLE_PRECISION, 0,  &
                     MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(SZOLD, 1, MPI_DOUBLE_PRECISION, 0,  &
                     MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(C10DUM, 10, MPI_CHARACTER, 0,  &
                     MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(NVOLD, 1, MPI_INTEGER, 0,  &
                     MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(NEOLD, 1, MPI_INTEGER, 0,  &
                     MPI_COMM_WORLD, ierr)
#endif
      ALLOCATE(ZI(MAX(NVOLD, NV+NLO)))
      DO 20 I = 1, NEOLD
        IF(MYID.EQ.0) THEN
           READ(98) IDUMMY, EIGVALI
           READ(98) (ZI(J),J=1,NVOLD)
!_COMPLEX!....eigenvectors are stored complex conjugated
!_COMPLEX         DO 14 j = 1, NVold
!_COMPLEX           zI(j)=dconjg(zI(j))
!_COMPLEX  14     continue
!
!          promote old eigenvectors
!
           if(nvold.lt.nv+nlo) then
              do 16 j=nlo,1,-1
                j1=nv+j
                j2=nvold-nlo+j
                zI(j1)=zI(j2)
   16         continue    
              DO 15 J = NVOLD-nlo+1,NV
                ZI(J) = ZERO
   15         CONTINUE
!
!          or reduce old eigenvectors
!
           else if(nvold.gt.nv+nlo) then
              do 17 j=1,nlo
                j1=nv+j
                j2=nvold-nlo+j
                zI(j1)=zI(j2)
   17         continue    
           end if
        END IF
#ifdef Parallel
         CALL MPI_BCAST(EIGVALI, 1, MPI_DOUBLE_PRECISION, 0,  &
                        MPI_COMM_WORLD, ierr)
         CALL MPI_BCAST(EIGVALI, 1, MPI_DOUBLE_PRECISION, 0,  &
                        MPI_COMM_WORLD, ierr)
         EIGVAL(I) = EIGVALI
!_REAL         CALL PDGEMR2D(NV+NLO, 1, ZI, 1,1,DESCLOCAL, Z, 1, I, DESCZ, &
!_REAL                       ICTXTALL)
!_COMPLEX         CALL PZGEMR2D(NV+NLO, 1, ZI, 1,1,DESCLOCAL, Z, 1, I, DESCZ, &
!_COMPLEX                       ICTXTALL)
#else
         EIGVAL(I)=EIGVALI
!_REAL         CALL DCOPY(NV+NLO, ZI, 1, Z(1,I), 1)
!_COMPLEX         CALL ZCOPY(NV+NLO, ZI, 1, Z(1,I), 1)
#endif
   20 CONTINUE
      DEALLOCATE(ZI)
!
#ifndef Parallel
!_COMPLEX!       Conjugate Hamilton matrix      
!_COMPLEX        DO 23 J = 1, NV + NLO
!_COMPLEX          DO 22 I = J+1, NV + NLO
!_COMPLEX            HS( I, J ) = DCONJG( HS( I, J) )
!_COMPLEX   22     CONTINUE
!_COMPLEX   23   CONTINUE
#endif
!     
   24 continue
!print
#ifndef Parallel
if(out) then
write(87,*) HSrows,nv,nlo
write(88,1234) HS
write(87,1234) HSdiag
 1234  format(4e23.15)
endif
#endif
      
      DO 25 I = 1, NUMIT
         CALL SECLIT(NEOLD,NV+nlo,el1,eu1)
   25 CONTINUE
!print
#ifndef Parallel
if(out) then
do i=1,neold
write(86,*) eigval(i)
write(89,1234) (Z(j,i),j=1,hsrows)
enddo
endif 
#endif

      NE = NEOLD
      NBELW = 0
      DO 30 I = 1, NE
         IF (EIGVAL(I) .LT. EL1) NBELW = NBELW + 1
   30 CONTINUE
      CALL STOP_TIMER(time_iter)
      write(6,'(A,F10.3)') 'time in seclr5 (CPU) ', READ_CPU_TIME(time_iter)
#ifdef Parallel
      write(6,'(A,F10.3)') 'time in seclr5 (WALL)', READ_WALL_TIME(time_iter)
#endif
      RETURN
!
!.....full diagonalization
!
999   write(6,*) ' NO OLD EIGENVECTOR!!'
      call seclr4
!
!        End of 'SECLR5'
!
      END

