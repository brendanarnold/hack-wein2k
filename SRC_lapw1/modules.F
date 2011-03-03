      module parallel
!       ICTXT     - Context of processors holding matrices H and S
!       ICTXTALL  - Context including all processors
!       DESCHS(9) - BLACS array descriptor of H and S
!       DESCZ(9)  - BLACS array descriptor of Z
!       NPE       - Number of processing elements
!       MYID      - Number of current processor
!       NPROW     - Number of process rows
!       NPCOL     - Number of process columns 
!                   (NPROW x NPCOL chosen 'as square as possible')
!       MYROWHS   - Row number of current processor
!       MYCOLHS   - Column number of current processor
!       BLOCKSIZE - Blocksize for two dimensional distribution of H, S, and Z,
!                   used for both dimensions,
!                   also used as algorithmic block size in ScaLapack
!       LDHS      - Leading dimension of local portion of H, S, and Z
!       LCOLHS    - Number of local columns of H and S
!       LCOLZ     - Number of local columns of Z
        INTEGER :: ICTXT, ICTXTALL, ICTXTALLR
        INTEGER :: DESCHS(9), DESCZ(9), DESCPANEL(9)
        INTEGER :: DESCLOCAL(9), DESCR(9)
        INTEGER :: NPE, MYID, NPROW, NPCOL
        INTEGER :: MYROWHS, MYCOLHS
        INTEGER :: MYROWPANEL, MYCOLPANEL
        INTEGER :: BLOCKSIZE
        INTEGER :: LDHS, LCOLHS, LCOLZ
        INTEGER :: LLD, myrow, mycol, info

      CONTAINS
        SUBROUTINE INIT_PARALLEL
          IMPLICIT NONE
#ifdef Parallel
          include 'mpif.h'
          INTEGER :: IERR
          call MPI_INIT(IERR)
          call MPI_COMM_SIZE( MPI_COMM_WORLD, NPE, IERR)
          call MPI_COMM_RANK( MPI_COMM_WORLD, MYID, IERR)
          IF(MYID.EQ.0) write(6,*) 'Using ', NPE, ' processors, My ID = ', MYID
          CALL BARRIER
          CALL SL_INIT(ICTXTALL, 1, NPE)
	  call blacs_get(0,0,ictxtallr)
          call blacs_gridinit(ICTXTALLR,'R',npe,1)
#else
          NPE=1
          MYID=0
#endif
        END SUBROUTINE INIT_PARALLEL

        SUBROUTINE INIT_PARALLELMATRICES(N,NUME2)
          IMPLICIT NONE
          INTEGER :: N, NUME2
#ifdef Parallel
          include 'mpif.h'
          include 'param.inc'

          INTEGER :: IERR, NPROWHS, NPCOLHS, KONE1
          INTEGER :: NPROWPANEL, NPCOLPANEL
          INTEGER, EXTERNAL :: NUMROC

          BLOCKSIZE = 400
          BLOCKSIZE = 2
          BLOCKSIZE = 1
          BLOCKSIZE = 16
          BLOCKSIZE = 4
          BLOCKSIZE = 1000
          BLOCKSIZE = 80
          BLOCKSIZE = 64
! one-dimensional distribution for local matrix setup
          CALL BLACS_GRIDINFO(ICTXTALL, NPROWPANEL, NPCOLPANEL, &
                              MYROWPANEL, MYCOLPANEL)
          IF(MYROWPANEL.NE.-1) THEN
            CALL DESCINIT(DESCPANEL, N, NPCOLPANEL*BLOCKSIZE,  &
                          BLOCKSIZE, BLOCKSIZE, 0, 0, ICTXTALL, N, IERR)
            CALL DESCINIT(DESCLOCAL, N, N,  &
                          N, N,0,0,ICTXTALL, N, IERR)
          ELSE
            DESCPANEL(2) = -1
          END IF
! two-dimensional distribution for global matrix computation
!         NPCOL=SQRT(real(NPE))
!         NPROW=NPE/NPCOL
          if (npe.eq.1) then
            nprow=1
          elseif((npe.eq.2).or.(npe.eq.4)) then
            nprow=2
          else if((npe.eq.3).or.(npe.eq.12).or.(npe.eq.6).or.(npe.eq.9))then
            nprow=3
          elseif((npe.eq.8).or.(npe.eq.16)) then
            nprow=4
          elseif((npe.eq.25)) then
            nprow=5
          elseif((npe.eq.36)) then
            nprow=6
          elseif((npe.eq.49)) then
            nprow=7
          elseif((npe.eq.100)) then
            nprow=10
          elseif((npe.eq.32).or.(npe.eq.64)) then
            nprow=8
          else
            nprow=16
          endif
          if(nprow.gt.npe) stop 'number of processors not supported'
          npcol=npe/nprow

          write(6,*) 'Working with ',NPCOL,' x ',NPROW,' Process Grid'
          CALL SL_INIT(ICTXT, NPCOL, NPROW)
          call BLACS_GRIDINFO(ICTXT, NPROWHS, NPCOLHS, MYROWHS, MYCOLHS)
          LDHS = NUMROC(N, BLOCKSIZE, MYROWHS, 0, NPROWHS)
          LCOLHS = NUMROC(N, BLOCKSIZE, MYCOLHS, 0, NPCOLHS)
!          LCOLZ = NUMROC(NUME2, BLOCKSIZE, MYCOLHS, 0, NPCOLHS)
! for pdsyevx: Z has to be NxN
          LCOLZ = NUMROC(N, BLOCKSIZE, MYCOLHS, 0, NPCOLHS)
          IF(MYROWHS.NE.-1) THEN
            CALL DESCINIT(DESCHS, N, N, BLOCKSIZE, BLOCKSIZE, 0, 0,  &
                          ICTXT, LDHS, IERR)
!            CALL DESCINIT(DESCZ, N, MIN(N,NUME2), BLOCKSIZE, BLOCKSIZE, 
!     +                    0, 0, ICTXT, LDHS, IERR)
!       for pdsyevx: Z has to be NxN
            CALL DESCINIT(DESCZ, N, N, BLOCKSIZE, BLOCKSIZE,  &
                          0, 0, ICTXT, LDHS, IERR)
          ELSE
            DESCHS(2) = -1
            DESCZ(2) = -1
          END IF
          write(6,*) 'My ID = ', MYID
          write(6,*) 'DESCPANEL, DESCHS, DESCZ (2):',  &
                      DESCPANEL(2), DESCHS(2), DESCZ(2)


! for HNS: row-wise blockcyclic distribution over npe processors
!      call blacs_gridinfo(ICTXTALLR,npe, 1, myrow, mycol)
      KONE1=1
      call blacs_gridinfo(ICTXTALLR,npe, KONE1, myrow, mycol)
      lld = max(1,numroc(N, blocksize, myrow, 0, npe))
      call descinit(descr, N, NSLMAX*NSLMAX,         &
                    blocksize, NSLMAX*NSLMAX, 0, 0, ICTXTALLR, lld, INFO)
!      print *, myrow, '> nsl2 = ', nslmax*nslmax, ', lld = ', lld
!      print *, 'descr = ', descr

#else
        BLOCKSIZE = 64
        NPCOL = 1
        NPROW = 1
        ICTXT = -1
        LDHS = N
        LCOLHS = N
        LCOLZ = NUME2
        DESCHS(2) = -1
        DESCZ(2) = -1
#endif
        END SUBROUTINE INIT_PARALLELMATRICES

        SUBROUTINE BARRIER
#ifdef Parallel
          IMPLICIT NONE
          include 'mpif.h'
          INTEGER :: IERR
          CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
#endif
        END SUBROUTINE BARRIER


        SUBROUTINE CLOSE_PARALLEL
#ifdef Parallel
          IMPLICIT NONE
          include 'mpif.h'
          INTEGER :: IERR
          CALL MPI_FINALIZE(IERR)
#endif
        END SUBROUTINE CLOSE_PARALLEL

      end module parallel

!--------------------------------------------------------------

      module lapw_timer
!  64 Timers are available in ScaLapack
      double precision :: cpu_time(64)
      double precision :: old_cpu_time(64)
      double precision :: wall_time(64)
      double precision :: old_wall_time(64)
      integer, parameter :: time_total=1,  &
                   time_hamilt=10, time_hns=11, time_diag=12,  &
                 time_coor=15,  &
                 time_setwar=17, &
                 time_pkp=18, &
                 time_cholesky=20, time_transform=21, time_sep=22,  &
                 time_backtransform=23, &
                 time_pwpw=25, time_pwlo=26,  time_lolo=27,  time_hns30=28, &
                 time_iter=30, time_HxZ=31, time_SxZ=32, time_formF=33,  &
                 time_projectH=34, time_proj_gsep=35,  time_residual=36, &
                 time_backproject=37, time_SxF=38, time_HxF=39, &
                 time_albl=40, time_loop260=41, time_phase=42,  &
                 time_legendre=43, time_step_function=44, time_h=45,  &
                 time_us=46, time_loop210=47, time_loop230=48,  &
                 time_loop240=49, time_overlap=50, time_lo=51, &
                 time_g1=55, time_atpar=56
      logical    usecputime, usewalltime
      contains
        subroutine init_all_timer
#ifdef Parallel
          use parallel, only : ICTXTALL
          IMPLICIT NONE
          external slboot
          double precision  :: WTIME
          CALL SLBOOT
          CALL SLCOMBINE( ICTXTALL, 'All', '>', 'C', 1, 1, WTIME )
          if(WTIME.eq.-1.0D0) then
            usecputime=.false.
            CALL SLENABLE
          else
            usecputime=.true.
          end if
          CALL SLCOMBINE( ICTXTALL, 'All', '>', 'W', 1, 1, WTIME )
          if(WTIME.eq.-1.0D0) then
            usewalltime=.false.
            CALL SLENABLE
          else
            usewalltime=.true.
          end if
#endif
          cpu_time(:)=0.0D0
          wall_time(:)=0.0D0
        end subroutine init_all_timer

        subroutine init_timer(n)
          IMPLICIT NONE
          integer :: n
          cpu_time(n) = 0.0D0
          wall_time(n) = 0.0D0
        end subroutine init_timer

        subroutine start_timer(n)
          integer :: n
#ifdef Parallel
          external sltimer
          call SLTIMER(n)
          if(.not.usecputime) then
            call CPUTIM(old_cpu_time(n))
            cpu_time(n)=0.0D0
          end if
          if(.not.usewalltime) then
            call WALLTIM(old_wall_time(n))
            wall_time(n)=0.0D0
          end if
#else
          call CPUTIM(old_cpu_time(n))
          cpu_time(n)=0.0D0
          call WALLTIM(old_wall_time(n))
          wall_time(n)=0.0D0
#endif
        end subroutine start_timer

        subroutine stop_timer(n)
          integer :: n
          double precision :: current_time
#ifdef Parallel
          external sltimer
          call SLTIMER(n)
          if(.not.usecputime) then
            call CPUTIM(current_time)
            cpu_time(n) = cpu_time(n) + current_time - old_cpu_time(n)
          end if
          if(.not.usewalltime) then
            call WALLTIM(current_time)
            wall_time(n) = wall_time(n) + current_time - old_wall_time(n)
          end if
#else
          call CPUTIM(current_time)
          cpu_time(n) = cpu_time(n) + current_time - old_cpu_time(n)
          call WALLTIM(current_time)
          wall_time(n) = wall_time(n) + current_time - old_wall_time(n)
#endif
        end subroutine stop_timer
         
        double precision function read_cpu_time(n)
          use parallel, only : ICTXTALL
          integer :: n
          double precision :: WTIME
#ifdef Parallel
          external slcombine
          if(usecputime) then
            CALL SLCOMBINE( ICTXTALL, 'All', '>', 'C', 1, N, WTIME )
            read_cpu_time=WTIME
          else
            read_cpu_time=cpu_time(n)
          end if
#else
          read_cpu_time=cpu_time(n)
#endif
          return
        end function read_cpu_time

        double precision function read_wall_time(n)
          use parallel, only : ICTXTALL
          integer :: n
          double precision :: WTIME
#ifdef Parallel
          external slcombine
          if(usewalltime) then
            CALL SLCOMBINE( ICTXTALL, 'All', '>', 'W', 1, N, WTIME )
            read_wall_time=WTIME
          else
            read_wall_time=wall_time(n)
          end if
#else
          read_wall_time=wall_time(n)
#endif
        end function read_wall_time
      end module lapw_timer
          
!--------------------------------------------------------------

        module reallocate
          !     nur 1 (generischer) Name wird von aussen angesprochen
          interface doreallocate
            module procedure doreallocate_r8_d1
            module procedure doreallocate_r8_d2
            module procedure doreallocate_r8_d3
            module procedure doreallocate_i4_d1
            module procedure doreallocate_i4_d2
            module procedure doreallocate_warp_c16_d2x3
            module procedure doreallocate_c3_d1
            module procedure doreallocate_c10_d1
            module procedure hugo     !   ;)
          end interface
        contains

          !     leider sind mehrere subroutines notwendig fuer verschiedene Typen
          subroutine doreallocate_r8_d1(tf, newdimension)
            implicit none
            real*8, pointer :: hilfsfeld(:), tf(:)
            integer min1, newdimension
            allocate(hilfsfeld(newdimension))
            hilfsfeld=0.0D0
            !     nur 1 mal kopieren reicht
            !     auch fuer mehrdimensionale Felder schaut die Zuweisung gleich aus
            min1=min(newdimension,size(tf,1))
            hilfsfeld(1:min1)=tf(1:min1)
            deallocate(tf)
            !     der Zeiger wird nur auf das neue Feld umgebogen, nicht neu alloziert
            tf=>hilfsfeld
          end subroutine

          subroutine doreallocate_r8_d2(tf, newdimension1, newdimension2)
            implicit none
            real*8, pointer :: hilfsfeld(:,:), tf(:,:)
            integer min1, min2, newdimension1, newdimension2
            allocate(hilfsfeld(newdimension1,newdimension2))
            hilfsfeld=0.0D0
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
            hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine

          subroutine doreallocate_r8_d3(tf, newdimension1, newdimension2, newdimension3)
            implicit none
            real*8, pointer :: hilfsfeld(:,:,:), tf(:,:,:)
            integer min1, min2, min3, newdimension1, newdimension2, newdimension3
            allocate(hilfsfeld(newdimension1,newdimension2,newdimension3))
            hilfsfeld=0.0D0
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
            min3=min(newdimension3,size(tf,3))
            hilfsfeld(1:min1,1:min2,1:min3)=tf(1:min1,1:min2,1:min3)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine

          subroutine doreallocate_i4_d1(tf, newdimension)
            implicit none
            integer*4, pointer :: hilfsfeld(:), tf(:)
            integer min1, newdimension
            allocate(hilfsfeld(newdimension))
            hilfsfeld=0
            min1=min(newdimension,size(tf,1))
            hilfsfeld(1:min1)=tf(1:min1)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine

          subroutine doreallocate_i4_d2(tf, newdimension1, newdimension2)
            implicit none
            integer*4, pointer :: hilfsfeld(:,:), tf(:,:)
            integer min1, min2, newdimension1, newdimension2
            allocate(hilfsfeld(newdimension1,newdimension2))
            hilfsfeld=0
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
            hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine

          subroutine doreallocate_warp_c16_d2x3(tf, n1a, n1b, n2a, n2b, n3a, n3b)
            implicit none
            complex*16, pointer :: hilfsfeld(:,:,:), tf(:,:,:)
            integer min1, n1a, n1b, n2a, n2b, n3a, n3b
            allocate( hilfsfeld(-n1b:n1b, -n2b:n2b, -n3b:n3b) )
            hilfsfeld(-n1b:n1b, -n2b:n2b, -n3b:n3b)=(0.0D0,0.0D0)
            hilfsfeld(-n1a:n1a,-n2a:n2a,-n3a:n3a)=tf(-n1a:n1a,-n2a:n2a,-n3a:n3a)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine

          subroutine doreallocate_c3_d1(tf, newdimension1, dummy)
            implicit none
            character*3, pointer :: hilfsfeld(:), tf(:)
            integer :: min1, newdimension1
            double precision :: dummy
            allocate(hilfsfeld(newdimension1))
            min1=min(newdimension1,size(tf,1))
            hilfsfeld(1:min1)=tf(1:min1)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine

          subroutine doreallocate_c10_d1(tf, newdimension1)
            implicit none
            character*10, pointer :: hilfsfeld(:), tf(:)
            integer min1, newdimension1
            allocate(hilfsfeld(newdimension1))
            min1=min(newdimension1,size(tf,1))
            hilfsfeld(1:min1)=tf(1:min1)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine

        !     Es gibt auch Methoden, um das Programm unleserlich zu machen :-)
        !     das sollten wir besser vermeiden ;-)
          subroutine hugo
            write(*,*) " Hier ist Hugo"
          end subroutine
        end module

!--------------------------------------------------------------

      module matrices
        INCLUDE 'param.inc'
        integer hsrows
!
!      sequential code:
!        HS - Hamilton matrix (hermitian matrix, stored below diagonal 
!                                                in conjugated form)
!           - Overlap matrix  (hermitian matrix, stored above diagonal)
!        HSDIAG - Diagonal elements of Overlap matrix
!               - several times swapped with 
!                 diagonal elements of Hamilton matrix
!      parallel code:
!        H - Hamilton matrix (hermitian matrix, stored above diagonal)
!        S - Overlap matrix (hermitian matrix, stored above diagonal)
!      HS, HSDIAG, H, S are set up in HAMILT, HNS 
!                       and solved in SECLR4, SECLR5
!      HPANEL, SPANEL: hold panels (size HSROWS*BLOCKSIZE) of upper
!                      triangles of H, S, and HS, HSDIAG
!        
#ifdef Parallel
!_REAL      DOUBLE PRECISION, allocatable ::  S(:, :)
!_REAL      DOUBLE PRECISION, allocatable ::  H(:, :)
!_COMPLEX      DOUBLE COMPLEX, allocatable :: S(:, :)
!_COMPLEX      DOUBLE COMPLEX, allocatable :: H(:, :)
#else
!_REAL      DOUBLE PRECISION, allocatable ::  HS(:, :)
!_REAL      DOUBLE PRECISION, allocatable ::  HSDIAG(:)
!_COMPLEX      DOUBLE COMPLEX, allocatable :: HS(:, :)
!_COMPLEX      DOUBLE COMPLEX, allocatable :: HSDIAG(:)
#endif
!_REAL      DOUBLE PRECISION, allocatable ::  Z(:,:)
!_REAL      DOUBLE PRECISION, allocatable ::  HPANEL(:,:)
!_REAL      DOUBLE PRECISION, allocatable ::  SPANEL(:,:)
!_COMPLEX      DOUBLE COMPLEX, allocatable :: Z(:, :)
!_COMPLEX      DOUBLE COMPLEX, allocatable :: HPANEL(:, :)
!_COMPLEX      DOUBLE COMPLEX, allocatable :: SPANEL(:, :)
        DOUBLE PRECISION, allocatable ::  EIGVAL(:)
        DOUBLE PRECISION, pointer ::  XK(:), YK(:), ZK(:)
        DOUBLE PRECISION, pointer ::  RK(:)
        INTEGER, pointer :: KZZ(:,:)
        CONTAINS
        SUBROUTINE INIT_MATRICES( N, NUME2)
        use parallel, only : DESCHS, DESCZ, NPE, MYID, &
                             NPROW, NPCOL, MYROWHS, MYCOLHS, &
                             BLOCKSIZE, LDHS, LCOLHS, LCOLZ, &
                             INIT_PARALLELMATRICES
        use reallocate, only: doreallocate
        IMPLICIT NONE        
        INTEGER N, NUME2
        include 'param.inc'
        INTEGER INFO

        if( n.eq.0 .and. nume2.eq.0) then
          allocate( RK(NMATMAX+1) )
          allocate( XK(NMATMAX+1) )
          allocate( YK(NMATMAX+1) )
          allocate( ZK(NMATMAX+1) )
          allocate( KZZ(3,NMATMAX+1) )
          RK=0.0D0
          XK=0.0D0
          YK=0.0D0
          ZK=0.0D0
          KZZ=0.0D0
        else
          if(myid.eq.0) write(6,*) 'Matrix size ', N
          HSROWS = N
!          write(6,*) 'Initialize parallelization variables'
          CALL INIT_PARALLELMATRICES(N, NUME2)
          call doreallocate( RK,HSROWS+1 )
          call doreallocate( XK,HSROWS+1 )
          call doreallocate( YK,HSROWS+1 )
          call doreallocate( ZK,HSROWS+1 )
          call doreallocate( KZZ,3,HSROWS+1 )
#ifdef Parallel
          allocate( H(LDHS, LCOLHS))        
          allocate( S(LDHS, LCOLHS))        
          allocate( Z(LDHS, LCOLZ))
          allocate( HPANEL( HSROWS, BLOCKSIZE) )
          allocate( SPANEL( HSROWS, BLOCKSIZE) )
          H=0.0D0
          S=0.0D0
          Z=0.0D0
          HPANEL=0.0D0
          SPANEL=0.0D0
#else 
          allocate( HS(HSROWS, HSROWS) )
          allocate( HSDIAG( HSROWS ) )
          allocate( Z(HSROWS, NUME2) )
          allocate( HPANEL(HSROWS, BLOCKSIZE) )
          allocate( SPANEL(HSROWS, BLOCKSIZE) )
          HS=0.0D0
          HSDIAG=0.0D0
          Z=0.0D0
          HPANEL=0.0D0
          SPANEL=0.0D0
#endif          
          allocate( EIGVAL(HSROWS) )
          EIGVAL=0.0D0
        end if

        END SUBROUTINE INIT_MATRICES
        
        SUBROUTINE END_MATRICE1
          deallocate( RK, XK, YK, ZK )
#ifdef Parallel
          deallocate( H, S, HPANEL, SPANEL )
#else
          deallocate( HS, HSDIAG, HPANEL, SPANEL )
#endif  
        END SUBROUTINE END_MATRICE1

        SUBROUTINE END_MATRICES
          deallocate( KZZ )
#ifdef Parallel
          deallocate( Z )
#else
          deallocate( Z )
#endif  
          deallocate( EIGVAL )
        END SUBROUTINE END_MATRICES

#ifdef Parallel
        SUBROUTINE PAR_SYR2K(TYPE,NV,NLO,A1R,A1I,A2R,A2I,LINDEX,H)
          use parallel, only: LDHS, LCOLHS,DESCHS,lld
          include 'param.inc'
! TYPE=1: Update H with A1/B1 and A2/B2
!     =2: Update H(local orbital part only) with C1 and C2
          INTEGER TYPE,NV,NLO, LINDEX
          DOUBLE PRECISION A1R(:,:)
          DOUBLE PRECISION A2R(:,:)
          DOUBLE PRECISION A1I(:,:)
          DOUBLE PRECISION A2I(:,:)
          DOUBLE PRECISION HALF, ONE
          PARAMETER (ONE=1.0D0, HALF=0.5D0)
          DOUBLE COMPLEX CHALF, CONE
          PARAMETER (CONE=(1.0D0,0.0D0), CHALF=(0.5D0,0.0D0))
!_REAL          DOUBLE PRECISION H(LDHS, LCOLHS)
!_REAL          DOUBLE PRECISION, ALLOCATABLE :: ABCLM1(:,:)
!_REAL          DOUBLE PRECISION, ALLOCATABLE :: ABCLM2(:,:)
!_COMPLEX          DOUBLE COMPLEX H(LDHS, LCOLHS)
!_COMPLEX          DOUBLE COMPLEX, ALLOCATABLE :: ABCLM1(:,:)
!_COMPLEX          DOUBLE COMPLEX, ALLOCATABLE :: ABCLM2(:,:)
!_COMPLEX          DOUBLE COMPLEX, ALLOCATABLE :: A(:,:)
          INTEGER N
!  TYPE=1: Update H with A1/B1 and A2/B2
	  IF(TYPE.EQ.1) THEN
	    N=NV+NLO
	    INDEX1=1
	  ELSE
! TYPE=2: Update H(local orbital part only) with C1 and C2
	    N=NLO
	    INDEX1=NV+1
	  END IF
          allocate( ABCLM1(LDHS, size(A1R,2)) )
          allocate( ABCLM2(LDHS, size(A1R,2)) )
!_COMPLEX          allocate( A(lld, size(A1R,2)) )
!!_REAL          ABCLM1=0.0D0
!!_REAL          ABCLM2=0.0D0
!!_COMPLEX          ABCLM1=(0.0D0,0.0D0)
!!_COMPLEX          ABCLM2=(0.0D0,0.0D0)
!!_COMPLEX          A=(0.0D0,0.0D0)

!_REAL          CALL COPYTOABCLM(TYPE,LINDEX,NV,NLO,A1R,ABCLM1)
!_REAL          CALL COPYTOABCLM(TYPE,LINDEX,NV,NLO,A2R,ABCLM2)
!#define symmetric
#ifdef symmetric
!_REAL          IF(deschs(2).ge.0) &
!_REAL            CALL PDSYR2K('U','N',N,LINDEX, &
!_REAL                          HALF, ABCLM1,INDEX1,1,DESCHS, &
!_REAL                                ABCLM2,INDEX1,1,DESCHS, &
!_REAL                          ONE,  H,INDEX1,INDEX1, DESCHS)
#else
!_REAL          IF(deschs(2).ge.0) &
!_REAL            CALL PDSYR2M ('U','N','T',N,LINDEX, &
!_REAL                          ONE, ABCLM2,INDEX1,1,DESCHS, &
!_REAL                               ABCLM1,INDEX1,1,DESCHS, &
!_REAL                          ONE, H,INDEX1,INDEX1,DESCHS)
#endif
!_REAL          CALL COPYTOABCLM(TYPE,LINDEX,NV,NLO,A1I,ABCLM1)
!_REAL          CALL COPYTOABCLM(TYPE,LINDEX,NV,NLO,A2I,ABCLM2)
#ifdef symmetric
!_REAL          IF(deschs(2).ge.0) &
!_REAL            CALL PDSYR2K('U','N',N,LINDEX, &
!_REAL                          -HALF, ABCLM1,INDEX1,1,DESCHS, &
!_REAL                                 ABCLM2,INDEX1,1,DESCHS, &
!_REAL                          ONE,   H,INDEX1,INDEX1, DESCHS)
#else
!_REAL          IF(deschs(2).ge.0) &
!_REAL            CALL PDSYR2M('U','N','T',N,LINDEX, &
!_REAL                          -ONE, ABCLM2,INDEX1,1,DESCHS, &
!_REAL                                ABCLM1,INDEX1,1,DESCHS, &
!_REAL                          ONE,  H,INDEX1,INDEX1, DESCHS)
#endif
!_COMPLEX          A(1:lld,1:LINDEX)=DCMPLX(A1R(1:lld,1:LINDEX), &
!_COMPLEX                                     -A1I(1:lld,1:LINDEX))
!_COMPLEX          CALL COPYTOABCLM(TYPE,LINDEX,NV,NLO,A,ABCLM1)
!_COMPLEX          A(1:lld,1:LINDEX)=DCMPLX(A2R(1:lld,1:LINDEX), &
!_COMPLEX                                     -A2I(1:lld,1:LINDEX))
!_COMPLEX          CALL COPYTOABCLM(TYPE,LINDEX,NV,NLO,A,ABCLM2)
#ifdef symmetric
!            write(*,*) 'asdf: N, LINDEX, HSROWS:', N, LINDEX, HSROWS
!_COMPLEX          IF(deschs(2).ge.0) &
!_COMPLEX            CALL ZSYR2K('U','N',N,LINDEX, &
!_COMPLEX                          CHALF, ABCLM1,HSROWS, &
!_COMPLEX                                 ABCLM2,HSROWS, &
!_COMPLEX                          CONE,  H,HSROWS)
!!_COMPLEX          IF(deschs(2).ge.0) &
!!_COMPLEX            CALL PZSYR2K('U','N',N,LINDEX, &
!!_COMPLEX                          CHALF, ABCLM1,INDEX1,1,DESCHS, &
!!_COMPLEX                                 ABCLM2,INDEX1,1,DESCHS, &
!!_COMPLEX                          CONE,  H,INDEX1,INDEX1, DESCHS)
#else
!_COMPLEX          IF(deschs(2).ge.0) &
!_COMPLEX            CALL PZHER2M('D','U','N','T',N,LINDEX, &
!_COMPLEX                          CONE, ABCLM2,INDEX1,1,DESCHS, &
!_COMPLEX                                ABCLM1,INDEX1,1,DESCHS, &
!_COMPLEX                          CONE, H,INDEX1,INDEX1, DESCHS)
#endif
          deallocate( ABCLM1, ABCLM2 )
!_COMPLEX          deallocate( A )
        END SUBROUTINE PAR_SYR2K

        SUBROUTINE COPYTOABCLM(TYPE,LINDEX,NV,NLO,A,ABCLM)
        use parallel, only: LDHS, NPE, MYCOLPANEL, BLOCKSIZE, &
                            ICTXTALL, DESCHS, DESCR, &
			    ICTXTALLR
        IMPLICIT NONE
        INTEGER TYPE, LINDEX, NV, NLO
        include 'param.inc'
!_REAL        DOUBLE PRECISION A  (:,:)
!_REAL        DOUBLE PRECISION ABCLM(:,:)
!_COMPLEX        DOUBLE COMPLEX A  (:,:)
!_COMPLEX        DOUBLE COMPLEX ABCLM(:,:)

! TYPE=1: Update H with A1/B1 and A2/B2
        IF(TYPE.EQ.1) THEN
!_REAL        CALL PDGEMR2D(NV+NLO,LINDEX, A,1,1, DESCR,  &
!_REAL                                ABCLM, 1,1, DESCHS, ICTXTALL)
!_COMPLEX       CALL PZGEMR2D(NV+NLO,LINDEX, A,1,1, DESCR,  &
!_COMPLEX                                ABCLM, 1,1, DESCHS, ICTXTALL)
        ELSE IF (NLO.GT.0) THEN
! TYPE=2: Update H(local orbital part only) with C1 and C2
!_REAL        CALL PDGEMR2D(NLO,LINDEX, A,NV+1,1, DESCR,  &
!_REAL                                ABCLM, NV+1,1, DESCHS, ICTXTALL)
!_COMPLEX       CALL PZGEMR2D(NLO,LINDEX, A,NV+1,1, DESCR,  &
!_COMPLEX                                ABCLM, NV+1,1, DESCHS, ICTXTALL)
        END IF
        END SUBROUTINE COPYTOABCLM
#endif 

      end module matrices

!--------------------------------------------------------------

      module albl
!
!        AL(n,l,i) - coefficient al(kn) for atom i
!                    (including equivalent atoms)
!        BL(n,l,i) - coefficient bl(kn) for atom i
!                    (including equivalent atoms)
!
        DOUBLE PRECISION, allocatable ::   AL(:,:,:)
        DOUBLE PRECISION, allocatable ::   BL(:,:,:)
      contains
        subroutine init_albl(hsrows, lmax, ndif)
          IMPLICIT NONE
          integer hsrows, lmax, ndif
          allocate( AL(HSROWS, 0:LMAX-1, NDIF) )
          allocate( BL(HSROWS, 0:LMAX-1, NDIF) )
          AL=0.0D0
          BL=0.0D0
        end subroutine init_albl
        subroutine end_albl
          deallocate( AL, BL)
        end subroutine end_albl
      end module albl

!--------------------------------------------------------------

      module loabc
!
!        ALO(l,i)    - local orbital coefficient alm for atom i
!        BLO(l,i)    - local orbital coefficient blm for atom i
!        CLO(l,i)    - local orbital coefficient clm for atom i
!        DPELO(l,i)  - derivative of the energy derivative ul_dot(r,E)
!                      of the local orbital radial-function ul(r,E) with
!                      respect to r
!                      (evaluated at r = muffin tin radius of atom i)
!        DPLO(l,i)   - derivative of the local orbital radial-function
!                      ul(r,E) with respect to radius r
!                      (evaluated at r = muffin tin radius of atom i)
!        ELO(l,i)    - local orbital expansion energy El for atom i
!        PEILO(l,i)  - the norm of the energy derivative ul_dot(r,E)
!                      of the local orbital radial function ul(r,E1)
!                      integrated over the muffin tin sphere
!        PELO(l,i)   - ul_dot(r,E)
!                      (evaluated at r = muffin tin radius of atom i)
!        PLO(l,i)    - local orbital radial-function ul(r,E)
!                      (evaluated at r = muffin tin radius of atom i)
!        PI12LO(l,i) - the product of the local orbital radial functions
!                      ul(r,E1) and ul(r,E2) integrated over the muffin
!                      tin sphere (for atom i)
!        PE12LO(l,i) - the product of the derivative of ul(r,E1) with
!                      respect to r and ul(r,E2) integrated over the
!                      muffin tin sphere (for atom i)
!
        DOUBLE PRECISION, allocatable ::   ALO(:,:,:)
        DOUBLE PRECISION, allocatable ::   BLO(:,:,:)
        DOUBLE PRECISION, allocatable ::   CLO(:,:,:)
        DOUBLE PRECISION, allocatable ::   ELO(:,:,:)
        DOUBLE PRECISION, allocatable ::   PLO(:,:)
        DOUBLE PRECISION, allocatable ::   DPLO(:,:)
        DOUBLE PRECISION, allocatable ::   PELO(:,:)
        DOUBLE PRECISION, allocatable ::   DPELO(:,:)
        DOUBLE PRECISION, allocatable ::   PEILO(:,:)
        DOUBLE PRECISION, allocatable ::   PI12LO(:,:)
        DOUBLE PRECISION, allocatable ::   PE12LO(:,:)
      contains
        SUBROUTINE init_loabc(LOMAX, NLOAT, NATO)
          IMPLICIT NONE
          integer LOMAX, NLOAT, NATO
          allocate( ALO(0:LOMAX,nloat,NATO) )
          allocate( BLO(0:LOMAX,nloat,NATO) )
          allocate( CLO(0:LOMAX,nloat,NATO) )
          allocate( ELO(0:LOMAX,nloat,NATO) )
          allocate( PLO(0:LOMAX,NATO) )
          allocate( DPLO(0:LOMAX,NATO) )
          allocate( PELO(0:LOMAX,NATO) )
          allocate( DPELO(0:LOMAX,NATO) )
          allocate( PEILO(0:LOMAX,NATO) )
          allocate( PI12LO(0:LOMAX,NATO) )
          allocate( PE12LO(0:LOMAX,NATO) )
          ALO=0.0D0
          BLO=0.0D0
          CLO=0.0D0
          ELO=0.0D0
          PLO=0.0D0
          DPLO=0.0D0
          PELO=0.0D0
          DPELO=0.0D0
          PEILO=0.0D0
          PI12LO=0.0D0
          PE12LO=0.0D0
        END SUBROUTINE init_loabc
      end module loabc

!--------------------------------------------------------------

      module radial
        DOUBLE PRECISION, allocatable ::   A(:), A1(:,:)
        DOUBLE PRECISION, allocatable ::   A1LO(:,:)
        DOUBLE PRECISION, allocatable ::   AE(:), AE1(:,:)
        DOUBLE PRECISION, allocatable ::   AE1LO(:,:,:)
        DOUBLE PRECISION, allocatable ::   AP(:), B(:), B1(:,:)
        DOUBLE PRECISION, allocatable ::   B1LO(:,:)
        DOUBLE PRECISION, allocatable ::   BE(:), BE1(:,:)
        DOUBLE PRECISION, allocatable ::   BE1LO(:,:,:)
        DOUBLE PRECISION, allocatable ::   BP(:), VLM(:,:)
      contains
        SUBROUTINE init_radial(NRAD, LOMAX, NSLMAX, NLOAT, LMMX)
          IMPLICIT NONE
          integer NRAD, NSLMAX, LOMAX, NLOAT, LMMX
          allocate( A(NRAD) )
          allocate( A1(NRAD,NSLMAX) )
          allocate( A1LO(NRAD,0:LOMAX) )
          allocate( AE(NRAD) )
          allocate( AE1(NRAD,NSLMAX) )
          allocate( AE1LO(NRAD,nloat,0:LOMAX) )
          allocate( AP(NRAD) )
          allocate( B(NRAD) )
          allocate( B1(NRAD,NSLMAX) )
          allocate( B1LO(NRAD,0:LOMAX) )
          allocate( BE(NRAD) )
          allocate( BE1(NRAD,NSLMAX) )
          allocate( BE1LO(NRAD,nloat,0:LOMAX) )
          allocate( BP(NRAD) )
          allocate( VLM(NRAD,LMMX) )
          A=0.0D0
          A1=0.0D0
          A1LO=0.0D0
          AE=0.0D0
          AE1=0.0D0
          AP=0.0D0
          B=0.0D0
          B1=0.0D0
          B1LO=0.0D0
          BE=0.0D0
          BE1=0.0D0
          BE1LO=0.0D0
          BP=0.0D0
          VLM=0.0D0
        END SUBROUTINE init_radial
        SUBROUTINE end_radial
          deallocate(A, A1, A1LO, AE, AE1, AE1LO, AP, B, B1, B1LO, BE, BE1, BE1LO, BP, VLM)
        END SUBROUTINE end_radial
      end module radial

!--------------------------------------------------------------

      module orb
!   NMOD=1 LDA+U, NMOD=2 OP, NMOD=3 Bext calculation
!   NSP=1 spin up, NSP=-1 spin down
!   NATORB number of atoms for which orbital potential is included
!   Iatom are the indexes of atom types selected for orb.pot.
!   NLORB number of L's of given atom for which OP is considered
!   LORB L's of given atom
!   Bext is energy of ext. magn field muB*Bext in Hartree
!   VORB is orbital potential matrix created by ORB package
!       COMPLEX*16 VORB(NATO,0:3,-3:3,-3:3)
!       DOUBLE PRECISION BEXT
!       INTEGER NMOD, NSP, NATORB, IATOM(NATO), NLORB(NATO), LORB(3,NATO)
!       COMMON /ORB/  Bext,VORB, NMOD,NSP,Natorb,IATOM,NLORB,LORB
!       SAVE   /ORB/
        DOUBLE PRECISION            :: BEXT
        DOUBLE COMPLEX, allocatable :: VORB(:,:,:,:)
        INTEGER                     :: NMOD, NSP, NATORB
        INTEGER, allocatable        :: IATOM(:), NLORB(:), LORB(:,:)
      contains
        SUBROUTINE init_orb(NAT)
          IMPLICIT NONE
          integer NAT
          allocate( VORB(NAT,0:3,-3:3,-3:3) )
          allocate( IATOM(NAT) )
          allocate( NLORB(NAT) )
          allocate( LORB(3,NAT) )
          VORB=(0.0D0,0.0D0)
          IATOM=0
          NLORB=0
          LORB=0
        END SUBROUTINE init_orb
      end module orb

!--------------------------------------------------------------

      module lolog
!        LOOR(l,i) - selects which local orbitals should be used for
!                    atom i
!        NLO       - total number of used local orbitals
!
!      INTEGER            NLO
!      LOGICAL            LOOR(0:LOMAX,NATO),lapw(0:lmax-1,nato)
!      COMMON  /LOLOG/    NLO, LOOR, lapw, ilo(0:lomax,nato)
!      SAVE    /LOLOG/
        INTEGER              :: NLO
        LOGICAL, allocatable :: LOOR(:,:), LAPW(:,:)
        INTEGER, allocatable :: ILO(:,:)
      contains
        SUBROUTINE init_lolog(LOMAX, NATO, LMAX)
          IMPLICIT NONE
          integer LOMAX, NATO, LMAX
          allocate( LOOR(0:LOMAX, NATO) )
          allocate( LAPW(0:LMAX-1, NATO) )
          allocate( ILO(0:LOMAX, NATO) )
          LOOR=.true.
          LAPW=.true.
          ILO=0
        END SUBROUTINE init_lolog
      end module lolog

!--------------------------------------------------------------

      module assleg
!      DOUBLE PRECISION   YR(N,MAXDIM)
!      COMMON  /ASSLEG/   YR
!      SAVE    /ASSLEG/
        DOUBLE PRECISION, allocatable :: YR(:,:)
        INTEGER N, MAXDIM
        PARAMETER (N=6, MAXDIM=81)
      contains
        SUBROUTINE init_assleg
          IMPLICIT NONE
!          integer N,MAXDIM
          allocate( YR(N, MAXDIM) )
          YR=0.0D0
        END SUBROUTINE init_assleg
      end module assleg

!--------------------------------------------------------------

      module atspdt
!        DP(l,i)      - derivative of the radial-function ul(r,E) with
!                       respect to radius r
!                       (evaluated at r = muffin tin radius of atom i)
!        DPE(l,i)     - derivative of the energy derivative ul_dot(r,E)
!                       of the radial-function ul(r,E) with respect to
!                       radius r
!                       (evaluated at r = muffin tin radius of atom i)
!        GFAC(:,i)    - gaunt-factors for atom i
!        INS          - flag if nonspherical correction (warpin) should
!                       take place
!                       if (INS .EQ. 0) then do warpin otherwise don't
!        LM(:,j,i)    - L,M combination for lattice harmonics j for
!                       atom i
!                       LM(1,j,i) ... L for j-th gaunt factor
!                       LM(2,j,i) ... M for j-th gaunt factor
!        LMMAX(i)     - number of L,M combinations for the
!                       nonspherical contribution of atom i
!        LQIND(i)     - number of gaunt factors for atom i
!        LQNS(:,j,i)  - l,m,LM,l',m' combination corresponding to
!                       the gaunt-factor GFAC(j,i) for atom i
!                       LQNS(1,1:LQIND(i),i) ... l'+1 for atom i
!                       LQNS(2,1:LQIND(i),i) ... L
!                       LQNS(3,1:LQIND(i),i) ... l+1 for atom i
!                       LQNS(4,1:LQIND(i),i) ... index of LM for atom i
!                       LQNS(5,1:LQIND(i),i) ... m' for atom i
!                       LQNS(6,1:LQIND(i),i) ... m for atom i
!        E(l,i)       - expansion energy El for atom i
!        P(l,i)       - radial-function ul(r,E)
!                       (evaluated at r = muffin tin radius of atom i)
!        PE(l,i)      - derivative of the radial-function ul(r,E) with
!                       respect to energy E
!                       (evaluated at r = muffin tin radius of atom i)
!        PEI(l,i)     - norm of ul_dot(r,E) integrated over the muffin
!                       tin sphere
!        VNS1(l',L,l,i) - non-spherical matrix elements for atom i
!                         integral(u_{l'}*V_{L}*u_{l} dr)
!        VNS2(l',L,l,i) - non-spherical matrix elements for atom i
!                         integral(u_{l'}_dot*V_{L}*u_{l}_dot dr)
!        VNS3(l',L,l,i) - non-spherical matrix elements for atom i
!                         integral(u_{l'}*V_{L}*u_{l}_dot dr)
!                         VNS3(l,L,l') ...
!                            ... integral(u_{l'}_dot*V_{L}*u_{l} dr)
!
!      INTEGER            INS
!      INTEGER            LMMAX(NATO), LQIND(NATO)
!      INTEGER            LM(2,LMMX,NATO), LQNS(6,NGAU,NATO)
!      DOUBLE PRECISION   DP(LMAX,NATO), DPE(LMAX,NATO), E(LMAX,NATO)
!      DOUBLE PRECISION   P(LMAX,NATO), PE(LMAX,NATO), PEI(LMAX,NATO)
!      DOUBLE PRECISION   VNS1(NSLMAX,LMMX,NSLMAX,NATO)
!      DOUBLE PRECISION   VNS2(NSLMAX,LMMX,NSLMAX,NATO)
!      DOUBLE PRECISION   VNS3(NSLMAX,LMMX,NSLMAX,NATO)
!      COMPLEX*16         GFAC(NGAU,NATO)
!      COMMON  /ATSPDT/   GFAC, E, P, DP, PE, DPE, PEI, VNS1, VNS2, VNS3, &
!                         LMMAX, LM, INS, LQNS, LQIND
!      SAVE    /ATSPDT/
        INTEGER                       :: INS, NL
        INTEGER, allocatable          :: LMMAX(:), LQIND(:), LM(:,:,:), LQNS(:,:,:)
        DOUBLE PRECISION, allocatable :: DP(:,:), DPE(:,:), E(:,:), P(:,:), PE(:,:), PEI(:,:)
        DOUBLE PRECISION, allocatable :: VNS1(:,:,:,:), VNS2(:,:,:,:), VNS3(:,:,:,:)
        DOUBLE COMPLEX, allocatable   :: GFAC(:,:)
      contains
        SUBROUTINE init_atspdt(NATO, LMMX, NGAU, LMAX, NSLMAX)
          IMPLICIT NONE
          integer NATO, LMMX, NGAU, LMAX, NSLMAX
          allocate( LMMAX(NATO) )
          allocate( LQIND(NATO) )
          allocate( LM(2,LMMX,NATO) )
          allocate( LQNS(6,NGAU,NATO) )
          allocate( DP(LMAX,NATO) )
          allocate( DPE(LMAX,NATO) )
          allocate( E(LMAX,NATO) )
          allocate( P(LMAX,NATO) )
          allocate( PE(LMAX,NATO) )
          allocate( PEI(LMAX,NATO) )
          allocate( VNS1(NSLMAX,LMMX,NSLMAX,NATO) )
          allocate( VNS2(NSLMAX,LMMX,NSLMAX,NATO) )
          allocate( VNS3(NSLMAX,LMMX,NSLMAX,NATO) )
          allocate( GFAC(NGAU,NATO) )
          LMMAX=0
          LQIND=0
          LM=0
          LQNS=0
          DP=0.0D0
          DPE=0.0D0
          E=0.0D0
          P=0.0D0
          PE=0.0D0
          PEI=0.0D0
          VNS1=0.0D0
          VNS2=0.0D0
          VNS3=0.0D0
          GFAC=(0.0D0,0.0D0)
        END SUBROUTINE init_atspdt
      end module atspdt

!--------------------------------------------------------------

      module char
!        LATTIC  - lattice type
!                  'P' ...... primitive latice (cubic, tetragonal,
!                             orthorhombic, monoclinic, triclin)
!                  'FC' ..... face centered
!                  'BC' ..... body centered
!                  'HEX' .... hexagonal
!                  'CXY' .... c-base centered (only orthorombic)
!                  'CYZ' .... a-base centered (only orthorombic)
!                  'CXZ' .... b-base centered (only orthorombic)
!                  'R' ...... rhombohedral
!        NAME(i) - name of atom i in the unit-cell (e.g. 'Titanium')
!
!      CHARACTER*4        LATTIC
!      CHARACTER*10       NAME(NATO)
!      COMMON  /CHAR/     LATTIC, NAME
!      SAVE    /CHAR/
        CHARACTER*4 :: LATTIC
        CHARACTER*10, allocatable :: NAME(:)
      contains
        SUBROUTINE init_char(NAT)
          IMPLICIT NONE
          integer NAT
          allocate( NAME(NAT) )
        END SUBROUTINE init_char
      end module char

!--------------------------------------------------------------

      module comc
!        IPGR(j)  - point group of j-th k-point (only used for output
!                   (documentation) purposes)
!        KNAME(j) - name of j-th k-point (optional)
!
!      CHARACTER*3        IPGR(NKPT)
!      CHARACTER*10       KNAME(NKPT)
!      COMMON  /COMC/     KNAME, IPGR
!      SAVE    /COMC/
        CHARACTER*3, pointer :: IPGR(:)
        CHARACTER*10, pointer :: KNAME(:)
      contains
        SUBROUTINE init_comc(NKPT, I)
          use reallocate, only: doreallocate
          IMPLICIT NONE
          integer NKPT, I
          double precision :: dummy
          if(I.eq.0) THEN
            allocate( IPGR(NKPT) )
            allocate( KNAME(NKPT) )
          else
            call doreallocate(IPGR, NKPT, dummy)
            call doreallocate(KNAME, NKPT)
          end if
        END SUBROUTINE init_comc
      end module comc

!--------------------------------------------------------------

      module comi
!        LNSMAX - maximum considered azimuthal quantum number l of ul
!                 for non-spherical hamilton matrix contributions
!        NAT    - number of inequivalent atoms
!        NBELW  - number of Eigenvalues below given energy window
!        NE     - number of calculated Eigenvalues
!        NT     - maximum considered azimuthal quantum number + 1
!                 (l+1 of ul) for the spherical contribution
!        NVAA   - number of plane waves
!
!      INTEGER            NVAA, NE, NBELW, NAT, NT, LNSMAX
!      COMMON  /COMI/     NVAA, NE, NBELW, NAT, NT, LNSMAX
!      SAVE    /COMI/
        INTEGER :: NVAA, NE, NBELW, NAT, NT, LNSMAX
      end module comi

!--------------------------------------------------------------

      module coml
!        PRNTWF - if .TRUE. print wave functions to output-file
!        SPRSWF - if .TRUE. suppress wave function calculation
!        WFTAPE - if .TRUE. write wave functions to vector-file
!        REL    - if .TRUE. perform relativistic calculations
!        ITER   - if .TRUE. perform iterative diagonalization
!        NOHNS  - if .TRUE. supress nonspherical matrix elements
!
!      LOGICAL     SPRSWF,  PRNTWF,  WFTAPE,  REL, ITER,NOHNS,force
!      COMMON  /COML/     SPRSWF,  WFTAPE,  PRNTWF,  REL, ITER,NOHNS
!      SAVE    /COML/

        LOGICAL :: SPRSWF,  PRNTWF,  WFTAPE,  REL, ITER,NOHNS
      end module coml

!--------------------------------------------------------------

      module comr
!        EL        - lower boundary of energy window for Eigenvalues
!        EU        - upper boundary of energy window for Eigenvalues
!        EIGVAL(j) - calculated Eigenvalues
!        RKM       - Rmt*Kmax determines matrix size (convergence)
!                    Rmt .... muffin tin radius
!                    Kmax ... plane wave cut-off
!        WEIGHT(j) - weight of the j-th k-vector (order of group of k)
!
!      DOUBLE PRECISION   RKM,  EL,  EU
!      DOUBLE PRECISION   WEIGHT(NKPT)
!      COMMON  /COMR/     RKM,  WEIGHT,  EL,  EU
!      SAVE    /COMR/
        DOUBLE PRECISION :: RKM, EL, EU
        DOUBLE PRECISION, pointer :: WEIGHT(:)
      contains
        SUBROUTINE init_comr(NKPT, I)
          use reallocate, only: doreallocate
          IMPLICIT NONE
          integer NKPT, I
          if(I.eq.0) THEN
            allocate( WEIGHT(NKPT) )
            WEIGHT(1:NKPT)=0
          else
            call doreallocate(weight, NKPT)
          end if
        END SUBROUTINE init_comr
      end module comr

!--------------------------------------------------------------

      module gener
!
!        BR2(1:3,1:3) - reciprocal Bravais matrix
!
!      DOUBLE PRECISION   BR2(3,3)
!      COMMON  /GENER/    BR2
!      SAVE    /GENER/
        DOUBLE PRECISION ::  BR2(3,3)
      end module gener

!--------------------------------------------------------------

      module kpts
!        SX(j) - x-component of K-point j
!        SY(j) - y-component of K-point j
!        SZ(j) - z-component of K-point j
!
!      DOUBLE PRECISION   SX(NKPT),  SY(NKPT),  SZ(NKPT)
!      COMMON  /KPTS/     SX,  SY,  SZ
!      SAVE    /KPTS/
        DOUBLE PRECISION, pointer :: SX(:),  SY(:),  SZ(:)
      contains
        SUBROUTINE init_kpts(NKPT, I)
          use reallocate, only: doreallocate
          IMPLICIT NONE
          integer NKPT, I
          if(I.eq.0) THEN
            allocate( SX(NKPT) )
            allocate( SY(NKPT) )
            allocate( SZ(NKPT) )
            SX(1:NKPT)=0.0D0
            SY(1:NKPT)=0.0D0
            SZ(1:NKPT)=0.0D0
          else
            call doreallocate( SX, NKPT)
            call doreallocate( SY, NKPT)
            call doreallocate( SZ, NKPT)
          end if
        END SUBROUTINE init_kpts  
      end module kpts

!--------------------------------------------------------------

      module loint
!        VNS1LO(l',L,l) - local orbital non-spherical matrix elements
!                         integral(u_{l'}*V_{L}*u_{l} dr)
!        VNS2LO(l',L,l) - local orbital non-spherical matrix elements
!                         integral(u_{l'}_dot*V_{L}*u_{l}_dot dr)
!        VNS3LO(l',L,l) - local orbital non-spherical matrix elements
!                         integral(u_{l'}*V_{L}*u_{l}_dot dr)
!                         VNS3(l,L,l') ...
!                           ... integral(u_{l'}_dot*V_{L}*u_{l} dr)
!
!      DOUBLE PRECISION   VNS1LO(LOMAX+1,LMMX,NSLMAX,NATO)
!      DOUBLE PRECISION   VNS2LO(LOMAX+1,LMMX,NSLMAX,NATO)
!      DOUBLE PRECISION   VNS3LO(LOMAX+1,LMMX,NSLMAX,NATO)
!      COMMON  /LOINT/    VNS1LO, VNS2LO, VNS3LO
!      SAVE    /LOINT/
        DOUBLE PRECISION, allocatable :: VNS1LO(:,:,:,:), VNS2LO(:,:,:,:), VNS3LO(:,:,:,:)
      contains
        SUBROUTINE init_loint(LOMAX, LMMX, NSLMAX, NATO)
          IMPLICIT NONE
          integer LOMAX, LMMX, NSLMAX, NATO
          allocate( VNS1LO(LOMAX+1, LMMX, NSLMAX, NATO) )
          allocate( VNS2LO(LOMAX+1, LMMX, NSLMAX, NATO) )
          allocate( VNS3LO(LOMAX+1, LMMX, NSLMAX, NATO) )
          VNS1LO(:,:,:,:)=0.0D0
          VNS2LO(:,:,:,:)=0.0D0
          VNS3LO(:,:,:,:)=0.0D0
        END SUBROUTINE init_loint  
      end module loint

!--------------------------------------------------------------

      module lstapw
!        NV    - total number of considered reciprocal lattice vectors
!                (number of plane waves)
!        XK(j) - x-coordinate (carthesian) of j-th k-vector
!                (ordered by distance from origin)
!        YK(j) - y-coordinate (carthesian) of j-th k-vector
!                (ordered by distance from origin)
!        ZK(j) - z-coordinate (carthesian) of j-th k-vector
!                (ordered by distance from origin)
!        RKM       - Rmt*Kmax determines matrix size (input)
!
!      INTEGER            NV
!      DOUBLE PRECISION   RKMT
!      DOUBLE PRECISION   XK(HSROWS+1), YK(HSROWS+1), ZK(HSROWS+1)
!      COMMON  /LSTAPW/   RKMT, NV
!      SAVE    /LSTAPW/
        DOUBLE PRECISION :: RKMT
        INTEGER          :: NV
      end module lstapw

!--------------------------------------------------------------

      module onekpt
!        SX,SY,SZ - k-point in the irreducible Brillouin Zone
!
!      DOUBLE PRECISION   SX1, SY1, SZ1
!      COMMON  /ONEKPT/   SX1, SY1, SZ1
!      SAVE    /ONEKPT/
        DOUBLE PRECISION :: SX1, SY1, SZ1
      end module onekpt

!--------------------------------------------------------------

      module orth
!        ORTH - .TRUE. for orthogonal lattice
!               .FALSE. otherwise
!
!      LOGICAL            ORTHO
!      COMMON  /ORTH/     ORTHO
!      SAVE    /ORTH/
        LOGICAL          :: ORTHO
      end module orth

!--------------------------------------------------------------

      module out
!        IMAT(1:3,1:3,j) - matrix representation of (space group)
!                          symmetry operation j
!        IORD            - number of symmetry operations of space group
!        KKK(1:3,j)      - star vectors of reciprocal lattice
!                          vector KZZ(:)
!        KZZ1(1:3)       - reciprocal lattice vector
!        NKK             - number of vectors in KKK (size of the star)
!        TAU(1:3,j)      - non-primitive translation vector for symmetry
!                          operation j
!        TAUP(j)         - tauphase factor of vector KKK(:,j)
!        WARP(i,j,k)     - interstitial warpin for K vector (i,j,k)
!
!      INTEGER            IORD, NKK
!      INTEGER            KKK(3,NSYM), IMAT(3,3,NSYM)
!      DOUBLE PRECISION   TAU(3,NSYM)
!      COMPLEX*16         TAUP(NSYM)
!      COMPLEX*16         WARP(-KMAX1:KMAX1,-KMAX2:KMAX2,-KMAX3:KMAX3)
!      COMMON  /OUT/      WARP, TAU, TAUP, IORD, IMAT, KZZ1, KKK, NKK
!      SAVE    /OUT/
        INTEGER                       :: IORD, NKK
        INTEGER                       :: KZZ1(3)
        INTEGER, allocatable          :: KKK(:,:), IMAT(:,:,:)
        DOUBLE PRECISION, allocatable :: TAU(:,:)
        DOUBLE COMPLEX, allocatable   :: TAUP(:)
        DOUBLE COMPLEX, pointer       :: WARP(:,:,:)
        INTEGER                       :: WARP1,WARP2,WARP3
        LOGICAL                       :: kkk_allocated
        LOGICAL                       :: warp_allocated
      contains
        SUBROUTINE init_out(NSYM, KMAX1, KMAX2, KMAX3)
          IMPLICIT NONE
          integer NSYM, KMAX1, KMAX2, KMAX3
          DOUBLE COMPLEX :: CZERO
          PARAMETER ( CZERO=(0.0D0,0.0D0) )
          if(KMAX1.eq.0) THEN
            allocate( KKK(3,NSYM) )
            allocate( IMAT(3,3,NSYM) )
            allocate( TAU(3,NSYM) )
            allocate( TAUP(NSYM) )
            kkk_allocated=.true.
            warp_allocated=.false.
            KKK(1:3,1:nsym)=0
            IMAT(1:3,1:3,1:nsym)=0
            TAU(1:3,1:nsym)=0.0D0
            TAUP(1:nsym)=CZERO
          else
            allocate( WARP(-KMAX1:KMAX1,-KMAX2:KMAX2,-KMAX3:KMAX3) )
            WARP(-KMAX1:KMAX1,-KMAX2:KMAX2,-KMAX3:KMAX3)=CZERO
            WARP1=KMAX1
            WARP2=KMAX2
            WARP3=KMAX3
            warp_allocated=.true.
          end if
        END SUBROUTINE init_out
        SUBROUTINE end_out
          if(kkk_allocated) then
            deallocate( KKK, IMAT, TAU, TAUP )
            kkk_allocated=.false.
          end if
          if(warp_allocated) then
            deallocate( WARP )
            warp_allocated=.false.
          end if
        END SUBROUTINE end_out
      end module out

!--------------------------------------------------------------

      module potnlc
!        DH(i)   - step size of the logarithmic radial mesh for atom i
!        JRJ(i)  - number of radial (logarithmic) mesh points for atom i
!        RO(i)   - first radial mesh point for atom i
!        VR(j,i) - spherical part (l=0, m=0) of the total potential r*V
!                  at mesh point j for atom i
!
!      INTEGER            JRJ(NATO)
!      DOUBLE PRECISION   DH(NATO), RO(NATO), VR(NRAD,NATO)
!      COMMON  /POTNLC/   VR,RO,DH,JRJ
!      SAVE    /POTNLC/
        INTEGER, allocatable :: JRJ(:)
        DOUBLE PRECISION, allocatable :: DH(:), RO(:), VR(:,:)
      contains
        SUBROUTINE init_potnlc(NATO, NRAD)
          IMPLICIT NONE
          integer NATO, NRAD
          allocate( JRJ(NATO) )
          allocate( DH(NATO) )
          allocate( RO(NATO) )
          allocate( VR(NRAD,NATO) )
          JRJ(1:NATO)=0
          DH(1:NATO)=0.0D0
          RO(1:NATO)=0.0D0
          VR(1:NRAD,1:NATO)=0.0D0
        END SUBROUTINE init_potnlc
      end module potnlc

!--------------------------------------------------------------

      module rotmat
!
!        ROTIJ(1:3,1:3,i)  - rotate atom i (including equivalent ones)
!                            according to the symmetry. For atoms with
!                            only one equiv. atom this is always the
!                            identity matrix.
!        ROTLOC(1:3,1:3,i) - local rotation matrix for atom i (from
!                            global to local coordinate system)
!
!      DOUBLE PRECISION   ROTIJ(3,3,NDIF),ROTLOC(3,3,NATO)
!      COMMON  /ROTMAT/   ROTLOC, ROTIJ
!      SAVE    /ROTMAT/
        DOUBLE PRECISION, allocatable :: ROTLOC(:,:,:)
        DOUBLE PRECISION, pointer :: ROTIJ(:,:,:),tauij(:,:)
      contains
        SUBROUTINE init_rotmat(NDIF, NATO, I)
          use reallocate, only: doreallocate
          IMPLICIT NONE
          integer NDIF, NATO, I
          if(I.eq.0) THEN
            allocate( ROTLOC(3,3,NATO) )
            allocate( ROTIJ(3,3,NDIF) )
            allocate( tauIJ(3,NDIF) )
            ROTLOC(1:3,1:3,1:NATO)=0.0D0
            ROTIJ(1:3,1:3,1:NDIF)=0.0D0
            tauIJ(1:3,1:NDIF)=0.0D0
          else
            call doreallocate(ROTIJ, 3,3,NDIF)
            call doreallocate(tauIJ, 3,NDIF)
          end if
        END SUBROUTINE init_rotmat
      end module rotmat

!--------------------------------------------------------------

      module cut1
!      logical            cut
!      common /cut1/      cut
        LOGICAL :: CUT
      end module cut1

!--------------------------------------------------------------

      module totpot
!        JA(j)   - 1. component of the j-th representative star vector
!        JB(j)   - 2. component of the j-th representative star vector
!        JC(j)   - 3. component of the j-th representative star vector
!        POTK(j) - interstitial total potential for the j-th star
!
!      INTEGER            JA(NWAV),  JB(NWAV),  JC(NWAV)
!      COMPLEX*16         POTK(NWAV)
!      COMMON  /TOTPOT/   POTK,  JA,  JB,  JC
!      SAVE    /TOTPOT/
        INTEGER, allocatable        :: JA(:), JB(:), JC(:)
        DOUBLE COMPLEX, allocatable :: POTK(:)
      contains
        SUBROUTINE init_totpot(NWAV)
          IMPLICIT NONE
          integer NWAV
          allocate( JA(NWAV) )
          allocate( JB(NWAV) )
          allocate( JC(NWAV) )
          allocate( POTK(NWAV) )
          JA(1:NWAV)=0
          JB(1:NWAV)=0
          JC(1:NWAV)=0
          POTK(1:NWAV)=(0.0D0,0.0D0)
        END SUBROUTINE init_totpot
      end module totpot

!--------------------------------------------------------------

      module struk
!        NDF         - total number of atoms
!        multmax     - max multiplicity of an atom
!        ALAT(1:3)   - lattice constants (a,b,c)
!        ALPHA(1:3)  - angles of the unit-cell's unit-vectors
!        IATNR(i)    - atom index of (inequivalent) atom i
!                      also indicates cubic and non-cubic symmetries
!                      IATNR(i) .GT. 0 ... cubic symmetry
!                      IATNR(i) .LT. 0 ... non-cubic symmetry
!        MULT(i)     - number of equivalent atoms for inequiv. atom i
!        PIA(1:3)    - reciprocal lattice constants (2*pi/a, 2*pi/b,
!                      2*pi/c)
!        POS(1:3,i)  - position vector of atom i (including equivalent
!                      atoms) in the unit-cell
!        RMT(i)      - muffin tin radius of atom i
!        V(i)        - relative muffin tin spherevolume for atom i
!        VI          - inverse volume of the direct lattice unit-cell
!
!      INTEGER            IATNR(NATO), MULT(NATO)
!      DOUBLE PRECISION   VI
!      DOUBLE PRECISION   ALAT(3), ALPHA(3), PIA(3), POS(3,NDIF)
!      DOUBLE PRECISION   RMT(NATO), V(NATO)
!      COMMON  /STRUK/    POS, ALAT, ALPHA, RMT, V, PIA, VI, IATNR, MULT
!      SAVE    /STRUK/
        INTEGER              :: NDF, multmax
        INTEGER, allocatable :: IATNR(:), MULT(:)
        DOUBLE PRECISION     :: VI
        DOUBLE PRECISION     :: ALAT(3), ALPHA(3), PIA(3)
        DOUBLE PRECISION, allocatable :: RMT(:), V(:)
        DOUBLE PRECISION, pointer :: POS(:,:)
      contains
        SUBROUTINE init_struk(NATO, NDIF, I)
          use reallocate, only: doreallocate
          IMPLICIT NONE
          integer NATO, NDIF, I
          if(I.eq.0) THEN
            allocate( IATNR(NATO) )
            allocate( MULT(NATO) )
            allocate( RMT(NATO) )
            allocate( V(NATO) )
            allocate( POS(3, NDIF) )
            IATNR(1:NATO)=0
            MULT(1:NATO)=0
            RMT(1:NATO)=0.0D0
            V(1:NATO)=0.0D0
            POS(1:3,1:NDIF)=0.0D0
          else
            call doreallocate(POS, 3, NDIF)
          end if
          NDF = NDIF
        END SUBROUTINE init_struk
      end module struk

!--------------------------------------------------------------
