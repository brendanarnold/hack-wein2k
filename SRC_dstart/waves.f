      module waves
!      include 'param.inc'
      INTEGER NWAVE
      INTEGER, pointer ::  KZZ(:,:)
      real*8,  pointer ::  ABSK(:)                            
      INTEGER, allocatable ::  INST(:)
      complex*16, allocatable ::  tauk(:)
!
      INTEGER, allocatable :: krcval(:,:,:)
      end module waves
