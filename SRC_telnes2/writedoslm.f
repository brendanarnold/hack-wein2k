!BOP
! !ROUTINE: writedoslm
! !INTERFACE:
    subroutine WriteDOSLM
! !USES:
    use energygrid
    use cross_dos,only : doslm
	use program_control,only : headers
! !DESCRIPTION:
!   Partial density of states resolved in l and m quantum numbers is written
!   to file 57 as a function of energy.
! !REVISION HISTORY:
!   Created November 2004 (Kevin Jorissen)
!EOP

      implicit none
      
!     LOCAL STAFF :
      integer I,J

!     Write the file case.doslm

      if (headers) then
      write(57,'(a,a,a)') '##  Energy, DOS_0, DOS_1, DOS_2, DOS_3,', &
      ' DOS_(1-1), DOS_(1 0), DOS_(1 1), DOS_(2-2), DOS_(2-1), DOS_(2 0),', & 
      ' DOS_(2 1), DOS_(2 2).'
      write(57,4714) iemax
	  else
	  write(6,'(a)') ':INFO : Routine WriteDOSLM omits the header of case.dos.  You have to add it yourself ',&
                         'if you want to use the file for a new calculation.'
	  write(6,'(a)') 'This is the header (next *3* lines) :'
      write(6,'(a,a,a)') '##  Energy, DOS_0, DOS_1, DOS_2, DOS_3,', &
      ' DOS_(1-1), DOS_(1 0), DOS_(1 1), DOS_(2-2), DOS_(2-1), DOS_(2 0),', & 
      ' DOS_(2 1), DOS_(2 2).'
      write(6,4714) iemax
      endif

      do I=1, IEMAX
         write(57,1) ENE(I), (doslm(I,J), J=0, 12)
      enddo
      
1     format(f10.5,13(x,e12.5))
4714  format('#',37x,I5,/)  ! iemax can be read by our input routines, but will appear as a comment to gnuplot
 
      RETURN
      END



