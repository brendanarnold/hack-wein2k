!BOP
! !ROUTINE: WriteDDlmLM
! !INTERFACE:
    subroutine WriteDDlmLM
! !USES:
    use energygrid
	use program_control,only : headers
! !DESCRIPTION:
!   Cross Density of States (i.e., products D(l,m) D(l',m')* , where D are coefficients in the
!   spherical harmonics - expansion of the wavefunction) are written to unit 56 for selected l,m
!   values and for each energy value considered.  Contrarily to the output of WriteXDOS, this output
!   is meant for human readers, not for computers ...
! !REVISION HISTORY:
!   Created November 2004 (Kevin Jorissen)
!EOP

      implicit none
           
      complex*16,external :: DDlmLM
      integer j
      
	  if (headers) then
      write(56,'(a)') '# Some cross DOS terms as a function of energy (eV).'
      write(56,'(a)') '# (direct DOS terms are in other files)'
      write(56,'(a)') '# energy, 00_1-1, 00_10, 00_11, 1-1_10, 1-1_11, 10_11, 00_2-2, 00_2-1, 00_20, 00_21, ',&
                      '00_22, 22_2-2, 22_21, 22_20, 21_2-1, 21_20'
	  endif

      do J = 1, IEMAX
        write(56,5555) ENE(J),DDlmLM(J,0,0,1,1,-1,1),DDlmLM(J,0,0,1,1,0,1),DDlmLM(J,0,0,1,1,1,1), &
        DDlmLM(J,1,-1,1,1,0,1),DDlmLM(J,1,-1,1,1,1,1),DDlmLM(J,1,0,1,1,1,1),DDlmLM(J,0,0,1,2,-2,1), &
        DDlmLM(J,0,0,1,2,-1,1),DDlmLM(J,0,0,1,2,0,1),DDlmLM(J,0,0,1,2,1,1),DDlmLM(J,0,0,1,2,2,1), &
        DDlmLM(J,2,2,1,2,-2,1),DDlmLM(J,2,2,1,2,1,1),DDlmLM(J,2,2,1,2,0,1),DDlmLM(J,2,1,1,2,-1,1),DDlmLM(J,2,1,1,2,0,1)
      enddo !J
      
5555  format(33F15.8)
      RETURN
      END



