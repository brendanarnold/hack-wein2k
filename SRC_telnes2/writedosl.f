!BOP
! !ROUTINE: writedosl
! !INTERFACE:
    subroutine WriteDOSL
! !USES:
    use energygrid
    use dimension_constants,only : lmax
    use densityofstates,only : dos
	use program_control,only : headers
! !DESCRIPTION:
!   Partial density of states resolved in l quantum number is written
!   to file 57 as a function of energy.  This file can be used as input
!   for a next run of TELNES2.
! !REVISION HISTORY:
!   Created November 2004 (Kevin Jorissen)
!EOP

      implicit none
      
!     LOCAL STAFF :
      integer I,J

!     Write the file case.dos:
      if (headers) then
      write(57,'(a)') '##  Energy, DOS_0, DOS_1, DOS_2, DOS_3'
	  write(57,4714) iemax
	  else
	  write(6,'(a)') ':INFO : Routine WriteDOSL omits the header of case.dos.  You have to add it yourself ',&
                         'if you want to use the file for a new calculation.'
	  write(6,'(a)') 'This is the header (next *3* lines) :'
      write(6,'(a)') '##  Energy, DOS_0, DOS_1, DOS_2, DOS_3'
	  write(6,4714) iemax
      endif
	  
      do I=1, iemax
         write(57,1) ENE(I), (DOS(I,J), J=0,lmax)
      enddo
      
1     format(f10.5,4(x,e12.5))
4714  format('#',37x,I5,/)  ! iemax can be read by our input routines, but will appear as a comment to gnuplot

      RETURN
      END


