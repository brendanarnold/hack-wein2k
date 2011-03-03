      MODULE FELDER
      INCLUDE 'param.inc'
!
      COMPLEX*16,  ALLOCATABLE ::   A(:,:),B(:,:)
      REAL*8,      ALLOCATABLE ::   EE(:)
      INTEGER,     ALLOCATABLE ::   KV(:,:)
      COMPLEX*16,  ALLOCATABLE ::   PH(:,:)  ! NSYM,NMAT
      COMPLEX*16,  ALLOCATABLE ::   XM(:,:)  ! NSYM,NUME
      INTEGER,     ALLOCATABLE ::   L(:,:)   ! NSYM,NMAT
!
      END MODULE FELDER
!
