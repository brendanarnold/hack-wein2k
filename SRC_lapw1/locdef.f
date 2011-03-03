      SUBROUTINE LOCDEF(RBAS,GBAS,ROTLOC)
!
!    Redefines local rotation matix from struct-file with
!    unitary transformation  u(-1) * S * u  for non-orthogonal lattice
!
!      INTEGER j
      IMPLICIT NONE
      DOUBLE PRECISION rbas(3,3),gbas(3,3),ROTLOC(3,3)
      DOUBLE PRECISION B(3,3)
!      write(6,111) (rotloc(1,j),j=1,3)
!      write(6,111) (rotloc(2,j),j=1,3)
!      write(6,111) (rotloc(3,j),j=1,3)
             CALL MATMM(B,GBAS,ROTLOC)
!      write(6,111) (b(1,j),j=1,3)
!      write(6,111) (b(2,j),j=1,3)
!      write(6,111) (b(3,j),j=1,3)
             CALL MATMM(ROTLOC,B,RBAS)
!      write(6,111) (rotloc(1,j),j=1,3)
!      write(6,111) (rotloc(2,j),j=1,3)
!      write(6,111) (rotloc(3,j),j=1,3)
111   format(3f10.4)
      RETURN
      END
