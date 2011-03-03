      module struct
      real*8, allocatable :: POS(:,:),ZNUC(:),RMT(:)
      integer,allocatable :: MULT(:),IATNR(:)
      end module struct
!
      module radgrd 
      real*8,allocatable :: RM(:,:),RNOT(:),DX(:)
      integer,allocatable :: JRI(:)
      end module radgrd
!
      module lolog 
      integer  NLO
      logical,allocatable :: LOOR(:,:),LORB(:)
      end module lolog
!
      module loabc
      real*8,allocatable :: ALO(:,:),BLO(:,:),CLO(:,:)  
      end module loabc
!
      module atspdt
      real*8,allocatable :: P(:,:),DP(:,:),PE(:,:),DPE(:,:) 
      end module atspdt
!
      module radfu
      real*8,allocatable ::  RRAD(:,:,:),RADE(:,:,:),RADL(:,:,:)    
      end module radfu
!
      module bessfu
      real*8,allocatable ::  FJ(:,:,:),DFJ(:,:,:),RAD(:)
      integer, allocatable :: IRAD(:)
      end module bessfu
!
      module work
      complex*16,allocatable :: aug(:,:,:)
      end module work
!
      module grid
      real*8,allocatable :: rgrid(:,:)
      integer,allocatable :: ireg(:),ilat(:,:),iri(:)
      integer npg
      end module grid
!
