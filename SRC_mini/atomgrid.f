      module atomgrid
      integer ndif
!.....allocatable with NATO
      real*8,allocatable :: POS(:,:)
      real*8,allocatable :: RMT(:), V(:),zatom(:),rnot(:)
      integer,allocatable :: IATNR(:), MULT(:),isplit(:),JRI(:) 
      reaL*8,allocatable ::  ROTLOC(:,:,:),ROTIJ(:,:,:),TAUIJ(:,:)            
      end module atomgrid
