      module atomgrid
!.....allocatable with NATO
      real*8,allocatable :: POS(:,:)
      real*8,allocatable :: RMT(:),V(:),Rnot(:),DX(:)
      integer,allocatable :: IATNR(:),MULT(:),ISPLIT(:),jri(:)           
      real*8,allocatable :: ROTLOC(:,:,:),ROTIJ(:,:,:),TAUIJ(:,:) 
      end module atomgrid
