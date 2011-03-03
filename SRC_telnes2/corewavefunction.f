!BOP
! !ROUTINE: Corewavefunction 
! !INTERFACE:
      SUBROUTINE Corewavefunction
! !USES:
      use input, only :              lc,corewffile,natom,nc,split
	  use dimension_constants,only : nrad
	  use program_control,only :     file2core,calcsplit,verbosity
	  use initialstate1
	  use initialstate2
      use crash
	  use constants,only :           ev2Ry
	  use struct,only :              jrj
! !DESCRIPTION:
!  This routine controls the acquirement of core wavefunctions.
!  Depending on input, it will either read them from a user specified file,
!  or calculate them by calling the wien lcore routine hfsd.
!
!  If the state does not have s-character, two edges are read/calculated :
!  for j=l+1/2 and for j=l-1/2.
!  If the states are calculated, also the energy splitting can be obtained
!  from the calculation.
!
!  If desired, the core state wavefunction is written to file 9.
!
!  Note that the program internally uses arrays containing wavefunction psi(r) * r.
!  However, in reading from/writing to file, the wavefunction psi(r) is used.
! !REVISION HISTORY:
!   Created November 2004 (Kevin Jorissen)
!EOP
	  implicit none
!  LOCAL VARIABLES
      real*8 energy1,energy2
	  integer i,kappa
	  
	  write(6,'(/,a)') 'Output from subroutine Corewavefunction :'
      call make_dp(nrad)      
      rewind(20)
	  rewind(7)
	  rewind(18)


      IF (file2core) THEN
!       Take the core wave function from file.  The file should contain the real wave function.
!       We multiply this by the radial coordinate, since this is what elnes traditionally uses.
	    write(6,'(a)') 'Reading core wave functions from file ...'
	    open(98,file=corewffile,form='formatted',status='old')
	    if(LC.gt.0) then
	      call make_dpdp(nrad)
          do i=1,NRAD
		    read(98,*) DR(i),DP(i),DQ(i),DPDP(i),DQDQ(i)
			dp(i)=dp(i)*dr(i)*dr(i)    !  The arrays as used by elnes already contain the metric factor r^2.
			dq(i)=dq(i)*dr(i)*dr(i)
			dpdp(i)=dpdp(i)*dr(i)*dr(i)
			dqdq(i)=dqdq(i)*dr(i)*dr(i)
		  enddo
        else   ! LC=0
          do i=1,NRAD
		    read(128,*) DR(i),DP(i),DQ(i)
			dp(i)=dp(i)*dr(i)*dr(i)    !  The arrays as used by elnes already contain the metric factor r^2.
			dq(i)=dq(i)*dr(i)*dr(i)
		  enddo
        endif ! LC
		close(98)
      ELSE

!     We calculate our core radials using the hfsd-Sub.
        if (lc.gt.0) then
	       call make_dpdp(nrad)
             KAPPA = LC    ! for the case j=l-1/2
             CALL ERRFLG(ERRFN,'Error in HFSD 1')
             CALL HFSD(NATOM,NC,LC,KAPPA,energy2)
!            save the results in DPdp and DQdq:
             DO I=1, NRAD+41
                DPdp(I) = DP(I)
                DQdq(I) = DQ(I)
             ENDDO
             REWIND(18)
             REWIND(7)
             REWIND(20)
        endif
!       hfsd(n,l,kappa) returns radial function: r*P(core,n',l') in DP(*) and DQ(*); 
      
        KAPPA = -LC-1    !   for the case j=l+1/2
        CALL ERRFLG(ERRFN,'Error in HFSD 2')
        CALL HFSD(NATOM,NC,LC,KAPPA,energy1)
      endif ! file2core


      if(.not.file2core.and.verbosity.ge.2) then
!     We will output the radial part of the core electron wave function :  u = DP / r^2 (DP = r^2 * u, ...)
!     We write it to to file case.corewavef.
!     To estimate the validity of the calculs, we can see what proportion
!     of the initial state electron lays out of the sphere.
        IF (LC.EQ.0) THEN
          DO I=1,jrj(natom)+40 !NRAD
            WRITE(9,4713) DR(I),DP(I) / DR(I) / DR(I),DQ(I) / DR(I) / DR(I)
          ENDDO   
        ELSE 
          DO I=1,jrj(natom)+40 !NRAD
            WRITE(9,4713) DR(I),DP(I) / DR(I) / DR(I),DQ(I) / DR(I) / DR(I),  &
                 DPdp(I)/DR(I)/DR(I),DQdq(I)/DR(I)/DR(I)
          ENDDO   
        ENDIF 
	  endif


      if(lc.gt.0.and.calcsplit) then
	    if(verbosity.ge.2) write(6,'(//)')
        write(6,'(a,f12.3,a)') 'Energy of j=l+1/2 core state : ',energy1/ev2Ry,' eV.'
		write(6,'(a,f12.3,a)') 'Energy of j=l-1/2 core state : ',energy2/ev2Ry,' eV.'
        split=dabs(energy1-energy2)/ev2Ry
		write(6,'(a,f12.3,a)') 'Split energy between the two edges : ',split,' eV.'
	  endif

      return
 4713 FORMAT(F10.5,7e14.5)
	  end   ! subroutine Corewavefunction

