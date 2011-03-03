!BOP
! !ROUTINE: writexdos
! !INTERFACE:
    subroutine WriteXDOS
! !USES:
    use energygrid
    use cross_dos,only : xdos
    use dimension_constants,only : nofcross
	use program_control,only : headers
! !DESCRIPTION:
!   Cross Density of States (i.e., products D(l,m) D(l',m')* , where D are coefficients in the
!   spherical harmonics - expansion of the wavefunction) are written to unit 58 for each l,m
!   and for each energy value considered.  This file can be used as input for the next run of
!   TELNES2.
! !REVISION HISTORY:
!   Created November 2004 (Kevin Jorissen)
!EOP

	  implicit none
      
!     LOCAL STAFF :
      integer I,J


!     Write the file case.dos:
      if (headers) then
      write(58,'(a,i5,a)') '##  Energy, all ',nofcross,' cross-dos components.'
      write(58,4714) iemax
	  else
	  write(6,'(a)') ':INFO : Routine WriteDOSL omits the header of case.dos.  You have to add it yourself ',&
                         'if you want to use the file for a new calculation.'
	  write(6,'(a)') 'This is the header (next *3* lines) :'
      write(6,'(a,i5,a)') '##  Energy, all ',nofcross,' cross-dos components.'
      write(6,4714) iemax
      endif

!     Write the file case.xdos:
      DO I=1, IEMAX
         WRITE(58,*) ene(i),(xdos(i,j),j=1,nofcross)
      ENDDO

 !1    format(f10.5,272(x,e12.5))      
 4714 FORMAT('#',37x,I5,/)  ! iemax can be read by our input routines, but will appear as a comment to gnuplot
      RETURN
      END

