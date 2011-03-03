!BOP
! !ROUTINE: densityovstates
! !INTERFACE:
      SUBROUTINE densityovstates
! !USES:
      use input,only : natom,emin,de,emax,nspin
      use constants
	  use program_control,only : verbosity,mdos,averaged,wdos,xqtlalanovak
	  use dimension_constants
	  use struct,only : isplit
      use crash
	  use tetra_params
	  use energygrid
	  use cross_dos
	  use densityofstates
 ! !DESCRIPTION:
!   This important routine governs the handling of (cross) density of states.
!   It is able to either read from file or calculate from qtl-input.
!   It can write several types of DOS (l-resolved,lm-resolved, cross terms) to file.
!   If necessary, input for tetraforelnes is constructed.
!   A distinction is made between averaged and orientation sensitive calculations ;
!   in the first case, only l-resolved DOS is necessary;
!   in the second case, cross-DOS is necessary.
!   For the definition of tetraforelnes-input, the isplit parameter is taken into account.
!
!   Name of the routine is chosen to prevent confusion with the module
!   named densityofstates.
! !REVISION HISTORY:
!   Created November 2004 (Kevin Jorissen)
!EOP
      implicit none
!  local variables :
      integer j,i,el,em,i1
	  real*8 eminin,dein,emaxin

      write(6,'(/,a)') 'Output from subroutine DensityOfStates:'
 	  write(6,'(a,f8.3,a)') 'The Fermi level is at ',ef,' Rydberg.'
	  
! Forget the old ISPLIT values from case.struct
! Rather : use the ones from the qtl-file.  This is more consistent.
         REWIND 30 !KJ                                                      
         READ(30,1000)                                            !KJ
         READ(30,1010)                         !KJ
         READ(30,1020) j   !KJ	  
	           do I=1,j                                           
!              NATO=NSORT-1  BECAUSE OF THE INTERSTITIAL DOS            
               IF((I.EQ.j).and.(.not.xqtlalanovak)) cycle
               READ(30,1030) isplit(i)   ! KJ lun changed from 4 to 30                      
!!	       write(17,*) 'i,jatom,nequ,isplit,text :',i,jatom,nequ,isplit(i),text(i)
         enddo
	   1000 FORMAT(1X,/)                                                  
 1010 FORMAT(1X,15X,x,x,x,16X,x)                                    
 1020 FORMAT(1X,4x,10X,4x,8X,1x,7X,I3,8x,2x)                                  
 1030 FORMAT(6X,3x,7X,2x,9x,i2,2x,45x)          
!     Sanity check
      if((.not.averaged).and.mdos.and.isplit(natom).ne.99) then
! If we do not have to calculate the DOS, then we don't care about what's in the qtl-file.
! We only care about the DOS-files then, and if there's no case.xdos, we'll get a read error ...      
	    call errflg(errfn,'Orientation resolved calculation selected and DOS has to be calculated, but ISPLIT is not 99')
		stop 'Stupid User Alert - see error file'
	  endif


      IF (.not.mdos) THEN
!  read energy grid and (cross)dos from file

        if (.not.averaged) then
!!         if (isplit(natom).eq.99) then
!       read the cross DOS coefficients (and set in DOS the good values):
          call errflg(errfn,'Error in ReadCrossDOS')
          call ReadCrossDOS
          write(6,'(a)') 'Cross DOS read from wien2k case.xdos file.'
          if(verbosity.ge.1.or.wdos) then
            call errflg(errfn,'Error in WriteDOS')
            call WriteDOSLM
			if(verbosity.ge.2) then
			  call errflg(errfn,'Error in WriteDDlmLM')
			  call WriteDDlmLM
		    endif
          endif
        else
!       read the partial DOS :
          call errflg(errfn,'Error in ReadDOS')
          call ReadDOS
	    write(6,'(a)') 'DOS read from wien2k case.dos file.'
!       no need to write anything here ;-)
        endif 

      ELSE
!     calculate dos and crossdos ...

!       first initialize energy grid
        jemin=1
	    jemax=1+int(0.001+(emax-emin)/de)  ! 0.001 helps to avoid rounding errors
	    iemax=jemax
	    call make_ene(iemax)
	    do j=jemin,jemax
	      ene(j)=emin+dble(j-1)*de
	    enddo

        if (.not.averaged) then
!       make crossdos - we know that isplit=99

          call init_tetra_params(2*nofcross,isplit(natom))
		  do i=1,nofcross
			idos(2*i-1,2)=99+2*i
			idos(2*i,2)=100+2*i
		     if(.not.xqtlalanovak) then
!  the qtl-file contains information for all atoms and for the interstitial region		     
		       idos(2*i-1,1)=natom
		       idos(2*i,1)=natom
		     else
! the qtl-file contains information for only one atom		     
		       idos(2*i-1,1)=1
		       idos(2*i,1)=1
		     endif
		  enddo
		  eminin=emin*ev2Ry+ef; dein=de*ev2Ry; emaxin=emax*ev2Ry+ef
		  call make_crossdos(iemax,nofcross,lmax)
	      call errflg(errfn,'Error in Tetra')
          call tetraforelnes(eminin,dein,emaxin,dble(0),2*nofcross,'xqtl')  ! tetra works in Ry !!
		  eminin=(eminin-ef)/ev2Ry
		  emaxin=(emaxin-ef)/ev2Ry
		  if(dabs(eminin-emin).gt.0.00001) write(6,*) ':INFO : tetra changed emin = ',emin,'to ',eminin,' eV.'
		  if(dabs(emaxin-emax).gt.0.00001) write(6,*) ':INFO : tetra changed emax = ',emax,'to ',emaxin,' eV.'
          write(6,'(a)') 'cross DOS calculated'
          if(verbosity.ge.1.or.wdos) then
		    call errflg(errfn,'Error while calculating DOS from cross DOS')
!           Calculate the partial DOS:
            CALL DefineIndex
            DO I=1, IEMAX
               DO el = 0, lmax  ! used to be 3
                  DOSlm(I, el) = dble(0)         !!! doslm used to be called dos
                  DO em = -el, el
                     I1 = (el+1)**2-el+em
                     DOSlm(I, IndexLM(el, em)) = dble(xdos(i, ((I1+1)*I1)/2))  != DDR(I, ((I1+1)*I1)/2)
                     IF (el.NE.0) THEN
!                      if el = 0, IndexLM (0,0) = 0 = el ...
                       DOSlm(I, el) = DOSlm(I, el) + DOSlm(I, IndexLM(el,em))
                     ENDIF
                  ENDDO
                  IF (el.EQ.3) THEN   ! since apparently we still do not have charge splitting for f ??
				    do em=-el,el
                      DOSlm(I, IndexLM(3, em)) = DOSlm(I, 3) / dble(7)
                    enddo
                  ENDIF
               ENDDO
            ENDDO
            write(6,'(a)') 'l,m resolved DOS calculated from cross DOS.'
	        call errflg(errfn,'Error in WriteDosLM')
	        call WriteDosLM
	      endif
		  if(wdos.or.verbosity.ge.2) then
		    call errflg(errfn,'Error in WriteXDos')
			call WriteXDos
			call errflg(errfn,'Error in WriteDDlmLM')
			call WriteDDlmLM
		  endif

        elseif(averaged) then
!       only make s,p,d,f dos

          call init_tetra_params(4,isplit(natom))
		  do i=1,4
		     if(.not.xqtlalanovak) then
!  the qtl-file contains information for all atoms and for the interstitial region		     
		       idos(i,1)=natom
		     else
! the qtl-file contains information for only one atom		     
		       idos(i,1)=1
		     endif
			idos(i,2)=i+1
		  enddo
		  if(isplit(natom).eq.0) then
!         do nothing; this is the default
		  elseif(isplit(natom).eq.1) then
		    idos(3,2)=6
			idos(4,2)=7
		  elseif(isplit(natom).eq.2) then
		    idos(4,2)=7
		  elseif(isplit(natom).eq.3) then
		    idos(4,2)=8
		  elseif(isplit(natom).eq.4) then
		    idos(3,2)=6
			idos(4,2)=10
		  elseif(isplit(natom).eq.5) then
		    idos(4,2)=10
		  elseif(isplit(natom).eq.6) then
		    idos(3,2)=7
			idos(4,2)=8
		  elseif(isplit(natom).eq.8.or.isplit(natom).eq.15) then
		    idos(3,2)=7
			idos(4,2)=13
		  elseif(isplit(natom).eq.-2) then
		    idos(3,2)=6
			idos(4,2)=11
	          elseif(isplit(natom).eq.99) then
		    ! no action is required
		  else
		    write(6,'(a,i4,/,a)') ':WARNING - unrecognized value of isplit ',isplit(natom), &
			         'Continuing with settings for default isplit= 0.'
		  endif

		  call make_dos(iemax,lmax)
                  call make_crossdos(iemax,nofcross,lmax) !! why is this necessary?
	          call errflg(errfn,'Error in Tetra')
                  call tetraforelnes(emin*ev2Ry+ef,de*ev2Ry,emax*ev2Ry+ef,dble(0),4,'spdf')  ! tetra works in Ry !!
	          write(6,'(a)') 'DOS calculated'
	          if(verbosity.ge.1.or.wdos) then
	            call errflg(errfn,'Error in WriteDosL')
	            call WriteDosL
	          endif

	  endif  ! dos or cross dos ?

	ENDIF  ! from file or calculation?

	return
	end

