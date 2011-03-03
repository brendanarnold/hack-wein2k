!BOP
! !ROUTINE: EnergyMesh
! !INTERFACE:
      SUBROUTINE EnergyMesh
! !USES:
      use energygrid
	  use input,only : emin,emax,split,lc
	  use dimension_constants
! !DESCRIPTION:
!     We initialise the energy scale, both in energy and in pixels.
!     The edge splitting (j=l+/-1/2) is translated into pixels.
! !REVISION HISTORY:
!   Created November 2004 (Kevin Jorissen)
!EOP
	  implicit none     
! LOCAL VARIABLES
      real*8 Difference
      integer i,j
	  logical test1, test2    


      write(6,'(/,a)') 'Output from subroutine InitEnergy:'
!     We take EMIN - SPLIT at the place of EMIN because the spectrum of the
!     second edge has to be calculated SPLIT eV before the other. It does not 
!     really increase the calculation-time because no calcul is done for E < 0.
      EMIN = EMIN - SPLIT

!     Test if EMIN and EMAX are within the given range of *.(x)dos; if not, set them appropriately.
      TEST1=.FALSE.
      TEST2=.FALSE.
      do J=1,IEMAX
         IF ((EMIN.LE.ENE(J)).AND.(.NOT.TEST1)) THEN
            JEMIN=J
            EMIN=ENE(J)
            TEST1=.TRUE.
         ENDIF
         IF ((EMAX.LE.ENE(J)).AND.(.NOT.TEST2)) THEN
            JEMAX=J
            EMAX=ENE(J)
            TEST2=.TRUE.
         ENDIF
      enddo
!     if everything else fails:
      IF (.NOT.TEST1) THEN 
         EMIN=ENE(1)
         JEMIN=1
      ENDIF
      IF (.NOT.TEST2) THEN 
         EMAX=ENE(IEMAX)
         JEMAX=IEMAX
      ENDIF

      write(6,'(a,a,i5,a,f8.3,a,i5,a,f8.3,a)') 'The calculation will use the energy interval from point ', &      
     'jemin = ', JEMIN, ' at energy emin = ',emin,' eV to point jemax = ',jemax,&
	 ' at energy emax = ',emax,' .'

!     initialise the energy scale for J < JEMIN
      Difference = ENE(JEMIN+1) - ENE(JEMIN)
      DO I = JEMIN, 0, -1
         ENE(I) =  ENE(I+1)-Difference 
      ENDDO      

      if (lc.gt.0) then
         JDecal = NINT(SPLIT/(ENE(JEMIN+1)-ENE(JEMIN)))
		 write(6,'(a,f7.3,a,i4,a)') 'Energy split is ',split,' eV, or ',jdecal,' pixels.'
	  endif
      
      RETURN
      END


