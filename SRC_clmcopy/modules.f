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


