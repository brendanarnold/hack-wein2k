      module felder
      INCLUDE 'param.inc'
      real*8, allocatable :: EBS(:,:), FC(:,:)     ! NKPT,NUME
      real*4, allocatable :: OPMAT(:,:,:) !NKPT,INUME,MG)
      real*8, allocatable :: MIMA(:,:)    !(NKPT,2)
      real*8, allocatable :: DENSTY(:,:,:) ! (INUMEden,MET,MG)
      integer ndif
      real*8  AA,BB,CC,ALPHA(3),PIA(3),VOL
      real*8, allocatable :: POS(:,:)   !   3,NDIF
      real*8, allocatable :: RMT(:),Vatom(:),R0(:) ! nato
      integer, allocatable :: IATNR(:),MULT(:),ISPLIT(:),JRI(:) ! nato
!      COMMON /STRUK/  POS(3,NDIF),AA,BB,CC,ALPHA(3),RMT(NATO),V(NATO), &
!                      PIA(3),VOL,IATNR(NATO),MULT(NATO),ISPLIT(NATO), &
!                      JRI(NATO),R0(NATO)
      CHARACTER*80 title
      CHARACTER*4  lattic
      CHARACTER*6  MODUS
      CHARACTER*10, allocatable:: ANAME(:)   !nato
!      COMMON /CHAR/   TITLE,LATTIC,MODUS,ANAME(NATO)
!      COMMON /ROTMAT/ ROTLOC(3,3,NATO),ROTIJ(3,3,NDIF),TAUIJ(3,NDIF)    
      real*8, allocatable :: ROTLOC(:,:,:)    ! (3,3,NATO)
      real*8, allocatable :: ROTIJ(:,:,:) ,TAUIJ(:,:) ! (3,3,NDIF),(3,NDIF) 
      end module felder
