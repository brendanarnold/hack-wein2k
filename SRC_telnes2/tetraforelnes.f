!BOP
! !ROUTINE: tetraforelnes
! !INTERFACE:
      subroutine tetraforelnes(emin,de,emax,prev,nndos,modus)
! !USES:
      use reallocate
      use tetra_params    
      use dimension_constants,only : nofcross
	  use program_control, only : verbosity,xqtlalaNovak
	use cross_dos,only : xdos
	use densityofstates,only : dos

! !INPUT/OUTPUT PARAMETERS:
!   emin         :   lower limit of energy grid in Ry
!   de           :   step size of energy grid in Ry
!   emax         :   upper limit of energy grid in Ry
!   prev         :   broadening of spectrum as Gaussian half-width; not used, only here for compatibility
!   nndos        :   number of DOS cases to calculate
!   modus        :   distinguishes between regular dos ('spdf') and cross dos ('xqtl') calculation
!   xqtlalaNovak :   distinguishes between two kinds of input : qtl-file from x lapw2 -qtl ('false')
!                     and qtl-file from P. Novak's qtl - program ('true') [passed in module]
! !DESCRIPTION:
!   Calculates the density of states from the data in case.qtl by integration over
!   the Brillouin Zone.
!   This is a ripped version of the SRC\_tetra/tetra.f program.
!   Input is not taken from case.int, but passed through calling arguments and
!   the module tetra\_params.
!   Param.inc is removed and replaced by the module tetra\_params.
!   All variables are real*8 (some were real*4).
!   The routine writes close to zero output.
!
!   Lun 5 is changed to 3.  Lun 4 is changed to 30.
!
!   If the qtl-file does not contain data for the full range requested, this is
!   no problem : only the range 1:jend is copied to dos or xdos arrays for
!   elnes.  Since these arrays were initialized to zero in their make routines,
!   the spectrum will be zero at energies for which no qtl data is available.
!   The calling routine densityovstates will write a warning to the log file
!   case.outputelnes if this happens.
!   This behavior deviates from the old TELNES, where emax and jemax would be
!   changed if no qtl was available.
!   This scheme justifies the existence of jend, next to iemax.
! !REVISION HISTORY:
!   Taken from SRC\_tetra/tetra.f of wien2k\_04.10
!   Updated November 2004 (Kevin Jorissen)
!   Added xqtlalaNovak December 2004 (Kevin Jorissen)
!   Major cleanup, added nspin February 2005 (Kevin Jorissen)
!EOP


      implicit none
!   IN/OUTPUT
!     next variables added by KJ :
      character*4, intent(in) ::  modus
      real*8, intent(in) :: de,prev
	  real*8  emin,emax
	  integer, intent(in) :: nndos
!     end addition KJ

!   LOCALS
      CHARACTER *7  A7                                                  
      CHARACTER *9  A9,B9                                               
      CHARACTER *20 TITLE                                               
      CHARACTER*45,allocatable :: TEXT(:)                             
      CHARACTER *70 SYSTEM                                              
      CHARACTER*11  FORM,STATUS                                         
      CHARACTER*80  FNAME ,errfn                                              
                                                                       
      real*8,allocatable :: RNUMB(:,:), DENSTY(:,:)   ! KJ used to be real*4
      real*8  fermden(MG)      ! KJ used to be real*4
      COMMON /NCOUNT/ NNOC                                              
      real*8,pointer :: EBS(:,:), FC(:,:,:)     ! KJ used to be real*4
      real*8,allocatable :: QTL(:,:),xqtl(:,:)     ! KJ used to be real*4
      real*8,allocatable :: ehelp(:),denbr(:,:)      ! KJ used to be real*4

      real*8 alat1,alat2,alat3,eferm,spin,so,emax1,bmin,bmax,eval 
            integer nst,minwav,maxwav,kspin,nsort,kso,jcase,isort,jatom,nequ,isplit
	  integer nbtot,ny,nrband,i,isrt,j,i1,ndos,iat,jat,nymin,nymax,jend,nnoc

! KJ  : We do not want to open files.  We also do not want to read from case.int.
                                           
!     NST          number of k points                                   
!     NYMIN        LOWER BAND INDEX                                     
!     NYMAX        UPPER BAND INDEX                                     
!     ZEL          NUMBER OF ELECTRONS                                  
!     EMIN,DE,EMAX ENERGY GRID FOR DOS                                  
!                                                                       
!     DOS1  DOS: GES, S, P, D, EG, T2G, F, FOR ATOM NR.1                
!     DOS2  DOS: S, ... F, OUT FOR ATOM NR.2                            
!                                                                       
      read(3,*) nst
      rewind 3
      nymin = 1

!!           open(17,file='debug.tetra',form='formatted')
!KJ In next lines, lun 4 was changed to 30 (as in elnes.def)
         REWIND 30 !KJ                                                      
         READ(30,1000) SYSTEM                                           !KJ
         READ(30,1010) ALAT1,ALAT2,ALAT3,EFERM                        !KJ
         READ(30,1020) MINWAV,MAXWAV,KSPIN,NSORT,kso   !KJ

         allocate( TEXT(Nsort+1),QTL(Nsort+1,15),xqtl(nsort,lxdos2*(lxdos2+1)) )
         if(.not.xqtlalaNovak) &   !KJ - new qtl file doesn't contain interstitial dos - meaning also that last field of text and qtl arrays are not used.
!           INCREASE NSORT BY ONE FOR THE INTERSTITIAL DOS (OUT)        
            NSORT=NSORT+1                                               
         SPIN=dble(1)                                                        
         IF(KSPIN.EQ.2) SPIN=dble(0.5)                                        
         so=dble(1)
         if(kso.eq.2) then
               so=dble(2)
               spin=spin/dble(2)
         endif
!!	 write(17,*) 'minwav,maxmav,kspin,nsort,kso :',minwav,maxwav,kspin,nsort,kso
!!	 write(17,*) 'spin,so : ',spin,so
! in spin.pol SO case set total DOS to DOS of interstital (=normso)
         if(kso.eq.2.and.kspin.eq.2) then 
          DO  JCASE=1,NNDOS                                               
           if(IDOS(JCASE,1).eq.0) then
             IDOS(JCASE,1)=NSORT                                          
             IDOS(JCASE,2)=1                                           
           ENDIF                                                          
          enddo
         endif

!!          write(17,*) 'nsort :',nsort
!!          write(17,*) 'nndos :',nndos
!!	  do i=1,nndos
!!	    write(17,*) i,' : iat,jat = ',idos(i,1),idos(i,2)
!!	   enddo
                                                                  
         do ISORT=1,NSORT                                           
!              NATO=NSORT-1  BECAUSE OF THE INTERSTITIAL DOS            
               IF((ISORT.EQ.NSORT).and.(.not.xqtlalanovak)) cycle
               READ(30,1030) JATOM,NEQU,isplit,text(isort)   ! KJ lun changed from 4 to 30                      
!!	       write(17,*) 'isort,jatom,nequ,isplit,text(isort) :',isort,jatom,nequ,isplit,text(isort)
         enddo
         text(nsort)='                       '                               

         emax1=-999.d0
!         nbtot=1000  KJ too large for full k-meshes
         nbtot=50  ! KJ large enough for small systems (eg. graphite)
         allocate (EBS(Nst,NBTOT), FC(MG,Nst,NBTOT))
         ny=0

!KJ next block 'if verbosity' added by me
         if(verbosity.ge.2) then
           write(6,'(a)') 'mg,nst,nbtot,nndos,lxdos,idos ='
  		   write(6,'(100i5)') mg,nst,nbtot,nndos,lxdos,idos(:,:)
		 endif

 115     NY=NY+1
         if(ny.gt.nbtot) then !KJ I added the write(6 statements
		    write(6,*) ':INFO - tetra increases nbtot from ',nbtot,' to ',nbtot*2,' .'
			write(6,*) 'If you get crashes in tetra, you may not have enough memory : '
			write(6,*) 'Currently required is about ',nint(dble(8*mg*nst*2*nbtot/1024/1024)),' Mb.'
            nbtot=nbtot*2
            call doreallocate(ebs, nst, nbtot)                                
            call doreallocate(fc, mg, nst, nbtot)     
         endif                        
            READ(30,1040,end=141) NRBAND                        ! KJ lun changed from 4 to 30
	  !!  write(17,*) 'nrband : ',nrband
            BMIN=10000.                                                 
            BMAX=-10000.                                                
            do I=1,NST
	     !!  write(17,*) 'i = ',i,' / nst = ',nst
               do ISORT=1,NSORT                                     
                  READ(30,1050) EVAL,ISRT,( QTL(ISORT,J), J=1,13 )      ! KJ lun changed from 4 to 30
		!!  write(17,*) 'isort = ',isort,' / nsort = ',nsort
		 !! write(17,*) 'eval,isrt,qtl(isort,:)',eval,isrt,qtl(isort,:)
                  if((text(isort)(13:20).eq.'xdos(i,j').or.xqtlalaNovak) then   !KJ : added .or...., since xqtl will always be there in new format
                    if(lxdos.ne.3) stop 'LXDOS must be 3 for ISPLIT 99 or 88'
                    read(30,208)(xqtl(isort,i1),i1=1,(lxdos2)*(lxdos2+1))      ! KJ lun changed from 4 to 30
					if(xqtlalaNovak) read(30,*)                   ! KJ they repeat the eigenenergy, I don't know why, but they do it ;-)
					!KJ they also write an additional field for the total dos of the atom, but we just forget about that.  we don't need it.
		!!write(17,*) 'xqtl(isort,:) ',xqtl(isort,:)
                  else if(text(isort)(13:20).eq.'xdos(i,i') then
                    read(30,208)(xqtl(isort,i1),i1=1,lxdos2)      ! KJ lun changed from 4 to 30
                  endif
               enddo            
               EBS(I,NY)=EVAL                                           
               BMIN=DMIN1(EVAL,BMIN)  ! careful : amin1 was for real*4                              
               BMAX=DMAX1(EVAL,BMAX)
	     !!  write(17,*) 'ebs(i,ny),bmin,bmax : ',ebs(i,ny),bmin,bmax
               do NDOS=1,NNDOS                                        
                 IAT=IDOS(NDOS,1)                                               
                 JAT=IDOS(NDOS,2)                                               
                 if(iat.gt.NSORT) then
		  write(6,*) 'nndos,ndos,iat,jat,nsort :',nndos,ndos,iat,jat,nsort
		  stop 'iat in input too large'
		 endif
                 IF(IAT.EQ.0) FC(NDOS,I,NY)=1.d0*SPIN      
                 IF(IAT.GT.0.and.jat.le.14)then
                   FC(NDOS,I,NY)=QTL(IAT,JAT)*SPIN*so    
                 else if(IAT.GT.0) then
                   FC(NDOS,I,NY)=xQTL(IAT,JAT-100)*SPIN
                 endif
		!write(17,*) 'ndos = ',ndos,' / nndos = ',nndos
		!write(17,*) 'fc(ndos,i,ny) = ',fc(ndos,i,ny)
               enddo
	       !!if(nrband.eq.1) write(125,'(13g13.5)') qtl(iat,1:13)
            enddo                                                    

            emax1=dmax1(bmin,emax1)   
!	    stop !band 1 finished
 140        goto 115
 141        nymax=ny-1
            if(emax.gt.emax1) then
              write(6,'(a)') ' EMAX reduced due to lower HIGHEST BAND-minimum'
              emax=emax1
            endif

      JEND=1 + int(0.001+(EMAX-EMIN)/DE    )  !KJ this to avoid rounding errors.
      allocate (RNUMB(jend,MG), DENSTY(jend,MG))
      allocate (ehelp(jend),denbr(jend,mg))
      DENSTY(:,:)=dble(0)
      RNUMB(:,:)=dble(0)
      write(6,55) EMIN,DE,EMAX 
     !! write(122,'(10g13.5)') fc(1:4,:,1) !:200)
     !! write(123,'(10g13.5)') ebs(:,1) !:200)
      CALL ARBDOS (0,NNDOS,1,JEND,EMIN,EMAX,NYMIN,NYMAX,DE,ebs,fc,nst,rnumb,densty)

     !! write(124,'(10g13.5)') densty(:,1:4)
!KJ : I do not wish to write files for the moment, or to apply broadening.
      if(modus.eq.'spdf') then
	    do j=1,4
	      dos(1:jend,j-1)=densty(1:jend,j)/dble(13.6058)*dble(2)  ! from /Ry to /eV ; I don't know what the *2 is for.
	    enddo
	  elseif(modus.eq.'xqtl') then
	    do j=1,nofcross
	      xdos(1:jend,j)=dcmplx(densty(1:jend,2*j-1),densty(1:jend,2*j))/dble(13.6058)*dble(2)
	    enddo
	  else
	    stop 'strange modus is confusing tetra.'
	  endif
	  rewind(30)
      deallocate(rnumb,densty,ehelp,denbr,ebs,fc,text,qtl,xqtl)  !KJ added KJ
      return    ! end addition by KJ
!                                                                       
      STOP ' LEGAL END TETRA'                                           
!                                                                       
  55   format(1X,'EMIN, DE, EMAX',14X,':',3F10.5,/)                    
  208 FORMAT((10F10.5))                             
 1000 FORMAT(1X,A70,/)                                                  
 1010 FORMAT(1X,15X,3F8.5,16X,F10.5)                                    
 1020 FORMAT(1X,I4,10X,I4,8X,I1,7X,I3,8x,i2)                                  
 1030 FORMAT(6X,I3,7X,I2,9x,i2,2x,a45)                                       
 1040 FORMAT(7X,I3)                                                     
 1050 FORMAT(F10.5,I3,F8.5,3X,12F8.5)                                
      END                                                               
