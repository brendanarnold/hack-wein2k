      module ldau
      real*8,allocatable ::   U(:,:),J(:,:)
      integer nldau
      end module ldau
      module struct
      real*8 AA,BB,CC,VOL
      real*8,allocatable :: RO(:),DXM(:),RMT(:),POS(:,:)
      character*10,allocatable :: ANA(:)
      integer,allocatable :: JRJ(:),NATOM(:),IATNR(:)
      integer NTYPE
      end module struct
      module opr
!.....allocatable with NATO
       complex*16,allocatable ::      dmat(:,:,:,:,:)

      real*8,allocatable :: VSP(:,:),EOP(:,:),zz(:), &
                    pop(:,:),al(:,:,:),alm(:,:), &
                    rotloc(:,:,:)

      integer,allocatable :: mult1(:),iatom(:),nlorb(:), &
                      lorb(:,:),NCALC(:)
      real*8  amix,Bexten
      integer nmod,natorb,nmodop
      end module opr
      module exex
       complex*16 vxcl(-3:3,-3:3)
       real*8 fxc, fupka(0:6),hybr,exc
       integer nexex,nsp
!!.....allocatable with NATO
       real*8,allocatable :: fxck(:,:),fupk(:,:),espxc(:)
       integer,allocatable ::  nlm1(:),lm1(:,:)
      end module exex
