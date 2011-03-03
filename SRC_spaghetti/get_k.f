      subroutine get_k (aline,vkx,vky,vkz,kname)
!     ******************************************
!
      IMPLICIT REAL*8 (A-H,O-Z)
      character  aline*80,kname*12
      character  ram_file*80
!-----------------------------------------------------------------------
!
!.....COPY THE LINE TO AN INTERNAL FILE
      write(ram_file,'(a80)') aline
!
!.....READ THE K-POINT COORDINATES
      read(ram_file,'(7x,3f10.5,3x,a12)') vkx,vky,vkz,kname
      return
      end
