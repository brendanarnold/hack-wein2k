!BOP
! !ROUTINE: DefineIndex
! !INTERFACE:
      SUBROUTINE DefineIndex
! !USES:
      use cross_dos, only : indexlm
! !DESCRIPTION:
!   Auxiliary routine for density of states.  The lm-DOS is stored in a
!   compact array.  The 'IndexLM' array, set up by DefineIndex, interfaces
!   to it and enables the user to use the lm-dos in a transparent way.
!   DOSLM(0) = s, DOSLM(1) = p, DOSLM(2) = d, DOSLM(3) = f
!   DOSLM(IndexLM(l, m)) = DOS\_lm

! !REVISION HISTORY:
!   Updated November 2004 (Kevin Jorissen)
!EOP
      implicit none

      IndexLM(0, 0) = 0
      IndexLM(1, -1) = 4
      IndexLM(1, 0) = 5
      IndexLM(1, 1) = 6
      IndexLM(2, -2) = 7
      IndexLM(2, -1) = 8
      IndexLM(2, 0) = 9
      IndexLM(2, 1) = 10
      IndexLM(2, 2) = 11
      IndexLM(3, -3) = 12
      IndexLM(3, -2) = 12
      IndexLM(3, -1) = 12
      IndexLM(3, 0) = 12
      IndexLM(3, 1) = 12
      IndexLM(3, 2) = 12
      IndexLM(3, 3) = 12

      RETURN
      END
      
