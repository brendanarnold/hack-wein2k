 MODULE defs
  REAL*8,PARAMETER             :: zero= 0.0d0
  COMPLEX*16,PARAMETER         :: zeroc= (0.0d0,0.0d0)
END MODULE defs

MODULE parallel
  INTEGER             :: ierr,npe,myid
  INTEGER,ALLOCATABLE :: nstart(:),nstop(:)
  INTEGER             :: comm_coul,comm_xc
  INTEGER             :: icoul
  INTEGER,ALLOCATABLE :: statusmpi(:)
  LOGICAL             :: coul,xc
  
  CONTAINS
  
  SUBROUTINE init_parallel
#ifdef Parallel
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    ALLOCATE(statusmpi(MPI_STATUS_SIZE))
    CALL MPI_INIT(ierr)
    IF(ierr.NE.MPI_SUCCESS) THEN
       PRINT*, 'MPI ini error'
       STOP
    ENDIF
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, npe, ierr)
    IF(ierr.NE.MPI_SUCCESS) THEN
       PRINT*, 'MPI_comm_size error'
       STOP
    ENDIF
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr)
#else
    myid=0
    npe=1
#endif
    ALLOCATE(nstart(0:npe-1),nstop(0:npe-1))
  END SUBROUTINE init_parallel

  SUBROUTINE begend(nat)
#ifdef Parallel
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER nat,modu,i,icolor1
    icolor1=MPI_UNDEFINED
    coul=.FALSE.
    xc=.FALSE.
    modu=(nat-1)/npe+1
    icoul=0
    icolor1=MPI_UNDEFINED
    DO i=0,npe-1
       nstart(i)=i*modu+1
       nstop(i)=MIN(nat,(i+1)*modu)
       IF(nstart(i).LE.nstop(i)) icoul=i
    ENDDO
    IF(nstart(myid).LE.nstop(myid)) THEN
       icolor1=1
       coul=.true.
       xc=.true.
    ENDIF
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, icolor1, myid, comm_coul,ierr)

!!$    icolor2=MPI_UNDEFINED
!!$    IF(myid.LT.nat) THEN
!!$       DO i=
!!$       nstart=myid*modu+1
!!$       nstop=MIN(nat,(myid+1)*modu)
!!$       icolor1=1
!!$       IF(nstop.EQ.nat) log_write=.TRUE.
!!$       IF(nstart.gt.nstop) THEN
!!$          icolor1=MPI_UNDEFINED
!!$          coul=.false.
!!$       ENDIF
!!$!       icolor2=1
!!$    endif
!!$    elseif(npe.ge.(2*nat)) then
!!$       icoul=nat-1
!!$       ixc=2*nat-1
!!$!       write(6,*) 'icoul,ixc',icoul,ixc
!!$       if(myid.lt.nat) then
!!$          modu=(nat-1)/(npe-nat)+1
!!$          nstart=myid*modu+1
!!$          nstop=MIN(nat,(myid+1)*modu)
!!$          icolor1=1
!!$       elseif(myid.ge.nat.and.myid.lt.(2*nat)) then
!!$          modu=(nat-1)/(npe-nat)+1
!!$          nstart=(myid-nat)*modu+1
!!$          nstop=MIN(nat,((myid-nat)+1)*modu)
!!$          icolor2=1
!!$      endif
!!$    endif
!    call MPI_COMM_SPLIT(MPI_COMM_WORLD, icolor2, myid, comm_xc,ierr)
#endif
  END SUBROUTINE begend

END MODULE parallel

MODULE ams
  REAL*8          :: atom_mass(103)

 CONTAINS
  SUBROUTINE init_ams
  REAL*8          :: a_m(103)
     DATA a_m /1.,4.,6.9,9.,10.8,12.,14.,16.,19.,20.2, &
          23.,24.3,27.,28.1,31.,32.,35.4,40.,39.1,40.,45., &
          47.9,50.9,52.,54.9,55.8,58.9,58.7,63.5,65.4,69.7, &
          72.6,74.9,79.,79.9,83.8,85.5,87.6,88.9,91.2,92.9, &
          95.9,98.,101.1,102.9,106.4,107.9,112.4,114.8, &
          118.7,121.8,127.6,126.9,131.3,132.9,137.3,138.9,140.1, &
          140.9,144.2,145.,150.4,152.,157.3,158.9,162.5, &
          164.9,167.3,168.9,173.,175.,178.5,180.9,183.8,186.2, &
          190.2,192.2,195.1,197.,200.6,204.4,207.2,209.,209., &
          210.,222.,223.,226.,227.,232.,231.,238.,237.,244.,243., &
          247.,247.,251.,252.,257.,258.,259.,262./     
     atom_mass(1:103) = a_m(1:103) 
   END SUBROUTINE init_ams
 END MODULE ams

module efeld
!!!  real*8  :: refeld
!!!  integer :: iefeld,jefeld(100)
  real*8  :: refeld, wefeld,fefeldpw(800)
  integer :: iefeld,jefeld(2000),kefeld,kebig
end module efeld

MODULE struct
  LOGICAL     :: ortho,inversion
  INTEGER     :: nat,nsym,ndif
  CHARACTER*4  :: lattic,irel,cform
  CHARACTER*80 :: title
  CHARACTER*10,ALLOCATABLE :: aname(:)
  REAL*8      :: br1(3,3),aa,bb,cc,br2(3,3),avec(3,3)
  REAL*8      :: a(3),alpha(3),pia(3),vol
  INTEGER,ALLOCATABLE :: iatnr(:),mult(:),jri(:),isplit(:)      !NATO
  REAL*8,ALLOCATABLE  :: r0(:),dx(:)                            !NATO   
  REAL*8,ALLOCATABLE  :: rotloc(:,:,:),rmt(:),zz(:),vsph(:)     !NATO   
  REAL*8,ALLOCATABLE  :: rotij(:,:,:),tauij(:,:)                !NDIF
  REAL*8,POINTER      :: pos(:,:)
END MODULE struct

MODULE symetr
  INTEGER   :: iord
  INTEGER,ALLOCATABLE :: iz(:,:,:),inum(:),kkk(:,:) !NSYM
  INTEGER,ALLOCATABLE :: kzz(:,:),inst(:)           !NVAW
  REAL*8,ALLOCATABLE  :: tau(:,:)                   !NSYM
END MODULE symetr

MODULE densit
  REAL*8,ALLOCATABLE :: clm(:,:,:),clmsp(:,:,:,:)
  COMPLEX*16,ALLOCATABLE :: rhok(:),rhoksp(:,:)
END MODULE densit

!!$module work
!!$  complex*16,allocatable :: YKA(:)
!!$  real*8,allocatable     :: V(:,:,:)
!!$  real*8,allocatable     :: r(:),vcoul2(:,:),value(:)
!!$end module work

module work
  complex*16,allocatable :: YKA(:)
  real*8,allocatable     :: V(:,:,:),exc(:,:,:,:),vxc(:,:,:,:)
  real*8,allocatable     :: r(:),vcoul2(:,:),value(:)
end module work

module vresp1 
  real*8,allocatable     :: vclmsp(:,:,:,:)
end module vresp1



