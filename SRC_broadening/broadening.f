      PROGRAM broadening                              
	
!     The BROADENING program is a series of broadening routines collected in 
!     such a way as to be useful for WIEN2k TELNES and XSPEC calculations.
!     It performs the following :
!      - initial state (core hole) broadening (convolution with Lorentz)
!      - final state (valence) broadening  (convolution with Lorentz of energy dependent width)
!      - experimental (spectrometer) broadening (convolution with a Gaussian)
!
!     BROADENING applies these three steps in slightly different ways according to the calculation modus :
!     - XAS
!     - XES
!     - ELNES
!                            
!     HISTORY :
!     $Id: broadening.f,v 1.5 2004/11/18 07:32:51 jluitz Exp $
!     Broadening of theoretical spectra
!     (C) 2002 by Joachim Luitz

!     new version based on the original "lorentz" program
!     (C)1996 by Joachim Luitz
!      with modifications by H. Enkisch, 18.01.1999
!     and many other additions

!     Update November 2004 Kevin Jorissen :
!         implicit none
!         separate broadening subroutines created
!         dynamical array allocation
!         elnes switch
!         partial spectra output also

!     Final touch November 2005 Kevin Jorissen

      implicit none

!  For file handling :
      CHARACTER*80 fname
      CHARACTER *11 STATUS,FORM
      INTEGER iunit,IRECL
      CHARACTER*80 DEFFN,errfn
!  Input from case.inb
      LOGICAL absorb,emis,elnes
	real*8 xint1,xint2,s,wshift,ga,gb,split
	integer w,ncol,col1,col2
	real*8 e0,e1,e2
!  Local variables
	real*8 ef,emin,delta,eadd
	integer isplit,istep,iadd,header
	real*8 gamma0,gamma
      character*67 errmsg
      character*4 modus
      real*8,allocatable :: X(:),Y(:),ya(:),yb(:),yenda(:),yendb(:),yend(:),ytemp(:,:)
      real*8 inta,intb
      real*8,parameter ::  pi=3.14159265359	
	integer i1,i2,i,nimax


! **************** OPENING FILES ***************************************

!     read filename for .def file
      call GTFNAM(DEFFN,errfn)
! windows :
!      deffn(1:14)='broadening.def'
!	errfn(1:16)='broadening.error'
!  end windows

      call ERRFLG(errfn,'Error in BROADENING')

      open (1,FILE=DEFFN,STATUS='OLD',err=910)
 8003 continue
      read (1,*,end=8001,err=911) iunit,fname,STATUS,FORM,IRECL
      open (iunit,FILE=fname,STATUS=STATUS,FORM=FORM,err=912)
      goto 8003
 8001 continue

! **************** readING INPUT *************************************

!     Set defaults
      emis=.false.
      absorb=.false.
	  elnes=.false.
      wshift=dble(0)
      ef=dble(0)

!     Read main input file case.inb, to get parameters for broadening the spectrum.

!     Read title into fname and forget it.
      read(5,'(a)',end=914,err=914) fname
	  write(6,'(a,/,a)') 'Name of our calculation : ',fname
!     Are we doing XES, ELNES, or XAS?  Set logical switch.
      read(5,'(4a)',end=914,err=914) modus
      if(modus.eq.'emis'.or.modus.eq.'EMIS')then
         emis=.true.
         write(6,*) 'We are going to do emission spectra.'
      elseif(modus.eq.'elne'.or.modus.eq.'ELNE') then
         elnes=.true.
         write(6,*) 'We are going to do ELNES spectra!'
      else  ! default is XAS
         absorb=.true.
         write(6,*) 'We are going to do absorption spectra.'
      endif


      if(elnes) then
!     how many columns in the file, and which two to use for broadening
!     (this allows us to do different types of files.)
        read(5,*,end=914,err=914) ncol,col1,col2
	  if(col1.gt.ncol.or.col2.gt.ncol) stop 'Stupid User Alert : Check column indices!'
	  write(6,*) 'Input file contains ',ncol,' spectra.  We use those in columns ', &
	                col1,' and ',col2,' .'
	endif
!     For XAS and XES, broadening.def specifies the correct file, I think.

!     Two edges are used.  The second one will be shifted by energy split.  The two edges
!     will also be rescaled by the respective weights xint1 and xint2.	     
      read(5,*,end=914,err=914) split,xint1,xint2
	if (xint1.lt.dble(0).or.xint2.lt.dble(0)) stop 's.U.E. : negative weights given!'
      write(6,*) 'Weights for rescaling edges : ',xint1,' and ',xint2,' .'
	if((xint1+xint2).eq.dble(0)) write(6,*) 'Ehh, this means your result will be ... Hmm ...  Zero!'

!     Widths for core hole lorentzian broadening :
      read(5,*,end=915,err=915) ga,gb
	if(ga.lt.dble(0).or.gb.lt.dble(0)) stop 's.U.E. : and what is negative broadening supposed to mean??'
	write(6,*) 'Core hole widths are ',ga,' and ',gb,' eV.'

!     Now obtain information for valence broadening :which type is used?
      read(5,*,err=9016,end=915) W,wshift
 9016 continue
      if(w.lt.0.or.w.gt.3) stop 'Invalid w in input file for broadening.'
      write(6,*) 'Valence broadening modus ',w,' and edge offset is ',wshift,' eV.'

!     Get width of Gaussian for experimental broadening.
      read(5,*,end=915,err=915) s
	  if(s.lt.dble(0)) stop 's.U.E. : and what is negative broadening supposed to mean??'
      write(6,*) 'Spectrometer FWHM resolution is ',s,' eV.'

!     Following parameters are only required for quadratic valence broadening.
      if(emis.OR.W.EQ.2) THEN
!KJ         read(5,'(4a)',end=914,err=914)modus
         read(5,*,end=915,err=915) E0
         read(5,*,end=915,err=915) E1
         read(5,*,end=915,err=915) E2
	   write(6,*) 'Parameters for quadratic energy dependent broadening :'
	   write(6,*) 'E0 = ',E0,', E1 = ',E1,', E2 = ',E2,' eV.'
      ELSE
         E0=dble(0)
         E1=dble(0)
         E2=dble(0)
      endif

!     Close the input file - we've read everything in it.
	close(5) 
      

! ******************* READING THE SPECTRUM ***************************************

!   first scan the file to count the header lines :
      nimax=0
	header=0
  1   read(47,'(a4)',end=2) modus
      if (modus(1:1).eq.'#'.or.modus(2:2).eq.'#'.or.modus.eq.'    ') then
	  header=header+1  ! comment or blank line
	else
	  nimax=nimax+1    ! 'real' line containing a spectrum
	endif
	goto 1
  2   continue
!   now position on first line with spectrum data (read past header) :
      rewind(47)
	if(header.gt.0) then
	  do i=1,header
	    read(47,*)
	  enddo
	endif

      isplit=0
      istep=0
      iadd=1  !KJ for elnes
!   allocate and initialize first set of arrays :
      allocate(x(nimax),y(nimax))
	x(:)=dble(0)
	y(:)=dble(0)


      
!     read the theoretical spectrum
      if (elnes) then
!       allocate remaining arrays :
	  allocate(ytemp(max(2,ncol),nimax))
        allocate(ya(nimax),yb(nimax),yend(nimax),yenda(nimax),yendb(nimax))
	  ytemp(:,:)=dble(0)
	  ya(:)=dble(0)
	  yb(:)=dble(0)
	  yend(:)=dble(0)
	  yenda(:)=dble(0)
	  yendb(:)=dble(0)
        do i=1,nimax
           read(47,*) x(i),(ytemp(i1,i),i1=1,ncol)
        enddo

	  if(col1.gt.0) ya(1:nimax)=ytemp(col1,1:nimax) ! if not, array remains zero.
	  if(col2.gt.0) yb(1:nimax)=ytemp(col2,1:nimax)
	  deallocate(ytemp)

	else  ! for xas and xes
        do i=1,nimax
           read(47,3000) x(i),y(i)
        enddo
 3000   format(f10.5,e14.5,e14.5)
      endif	  

      close(47)


!     get energy step width
      delta=(X(2)-X(1))/2
      isplit=ABS(INT((nimax-1)/(X(nimax)-X(1))*split)) ! L2-L3 split in pixels
      istep=ABS(INT((nimax-1)/(X(nimax)-X(1))*(wshift-x(1)))) ! fermi energy in pixels

      write(6,*) 'Spectra are split by',isplit,' channels (',split,' eV).'

      if(.not.elnes) then  ! elnes has already split the spectra itself.
!       how far should we go to the left + right?
        emin=x(1)
        eadd=dble(10)
        if (split.gt.5) then
           eadd=split+dble(5)
        endif
        iadd=ABS(INT((nimax-1)/(X(nimax)-X(1))*eadd))
        write(6,*) 'channels to add: ',iadd, ' (',eadd,' eV)'

!       enlarge x and y arrays first :
        allocate(ytemp(2,nimax))
	  ytemp(1,:)=x
	  ytemp(2,:)=y
	  deallocate(x,y)
	  allocate(x(nimax+2*iadd),y(nimax+2*iadd))
	  x(1:nimax)=ytemp(1,1:nimax)
	  y(1:nimax)=ytemp(2,1:nimax)
	  deallocate(ytemp)


!       moving data by added channels
        do i=nimax,1,-1
           x(i+iadd)=x(i)
           y(i+iadd)=y(i)
        enddo

!       adding 0's at left
        do i=1,iadd
           x(i)=emin+delta*2*(i-1-iadd)
           y(i)=dble(0)
        enddo

!       adding 0's at right
        do i=nimax+1+iadd, nimax+2*iadd
           x(i)=emin+delta*2*(i-1-iadd)
           y(i)=dble(0)
        enddo

!       now we have a new range...
        nimax=nimax+2*iadd

!       allocate remaining arrays
        allocate(ya(nimax),yb(nimax),yend(nimax),yenda(nimax),yendb(nimax))
	  ya(:)=dble(0)
	  yb(:)=dble(0)
	  yend(:)=dble(0)
	  yenda(:)=dble(0)
	  yendb(:)=dble(0)
      endif  ! not elnes




! ************************ VALENCE BROADENING ***********************

	if(.not.elnes) then
!     first we do valence broadening: it's the same for L3 and L2, and L3 and L2 are the same
        call ValenceBroadening(X,Y,yend,w,absorb,istep,wshift,E0,E1,E2,ef,delta,nimax)
        ya=yend
!       this doesn't hurt, even if we only have a K spectrum
        yb=yend
      else   ! elnes
!     L2 and L3 are different, and L2 is already shifted by split w.r.t. L3.
        call ValenceBroadening(x,ya,yend,w,.true.,istep,wshift,E0,E1,E2,ef,delta,nimax)
        ya=yend
        call ValenceBroadening(x,yb,yend,w,.true.,istep+isplit,wshift+split,E0,E1,E2,ef,delta,nimax)
        yb=yend
      endif  ! elnes ?

      
! ******************** CORE LIFE TIME BROADENING ************************************

!     first for L3:  ya -> yenda
      ga=ga/dble(2)
      call CoreBroadening(x,ya,yenda,nimax,delta,ga)
!     now we do L2:  yb -> yendb
      gb=gb/dble(2)
	call CoreBroadening(x,yb,yendb,nimax,delta,gb)

! ******************* SPECTROMETER BROADENING *****************************************

! KJ : We could construct total spectrum first and call SpecBro only for total spectrum.
! KJ : But I'd rather have all spectra.
      call SpectrometerBroadening(x,yenda,ya,nimax,delta,s)  ! yenda -> ya
      call SpectrometerBroadening(x,yendb,yb,nimax,delta,s)  ! yendb -> yb           

! ******************* COMPUTE TOTAL SPECTRUM *****************************************
!           ya,yb -> y

      if(elnes) then
	   y= xint1*ya+xint2*yb
	elseif(split.eq.0) then
!        -> forget the second part, hehehe
         write(6,*)'no split -> restore'
	   y=ya*(xint1+xint2)
	else
!        we did split :-(
!        oh, oh, now we need some brains.....

!        fortunately we have that routine already in the original txspec
!        thanks to some clever students. Pooh, and I was afraid that
!        I might have to think today :-)=

!        a few days later; I hardly slept...
!        ...
!        lots of stuff deleted after hours or painful editing
!        ...

!        Obviously the "clever students" were not so clever

         do i=1,nimax
!        don't blame me if your result is zero when you set xint1 and xint2 to zero.
            y(i)=xint1*ya(i)+xint2*yb(i-isplit)
         enddo
      endif

! ****************************** WRITE OUTPUT *****************************	
       
!     write header
      if(elnes) then
	  write(46,'(a)') '# ELNES absorption spectra'
	elseif(absorb) then
	  write(46,'(a)') '# XANES absorption spectra'
	elseif(emis) then
	  write(46,'(a)') '# XES emission spectra'
	endif
      write(46,*) '# split= ',split
      write(46,*) '# xint1,xint2= ',xint1,xint2
      write(46,*) '# ga,gb= ',ga,gb
      write(46,*) '# W,wshift= ',W,wshift
      write(46,*) '# s= ',s
      
      if (emis.or.w.eq.2) then
         write(46,*) '#W=',W
         write(46,*) '#E0=',E0
         write(46,*) '#E1=',E1
         write(46,*) '#E2=',E2
      endif

!     another header ...
      if(elnes) then
	   write(46,'(/,a,/)') '# Energy        Total spectrum      First edge        Second edge'
	else
	   write(46,'(/,a,/)') '# Energy        Total spectrum'
	endif

!     write out broadened spectrum
      do i=iadd,nimax-iadd
         if(elnes) then
	      write(46,'(4G18.10)') x(i),y(i),ya(i),yb(i)
	   else
            write(46,3000) x(i),y(i)
	   endif
      enddo


! ***************************** END OF PROGRAM *****************************

      deallocate(x,y,ya,yb,yend,yenda,yendb)

      call errclr(errfn)
      stop 'BROADENING done'

!     Errors

!     broadening.def couldn't be opened
 910  write(errmsg,9000) fname
      call outerr('BROADENING',errmsg)
      goto 999

!     broadening.def is corrupt
 911  write(errmsg,9001)fname
      call outerr('BROADENING',errmsg)
      goto 999

!     Error opening a unit
 912  write (errmsg,9010) iunit
      call outerr('BROADENING',errmsg)
      write (errmsg,9020) fname
      call outerr('BROADENING',errmsg)
      write (errmsg,9030) STATUS, FORM
      call outerr('BROADENING',errmsg)
      goto 999

 914  write(errmsg,9040)
      call outerr('BROADENING',errmsg)
      goto 999

 915  write(errmsg,9050)
      call outerr('BROADENING',errmsg)
      goto 999

 999  stop 'BROADENING - Error'


 9000 format('can''t open definition file ',A40)
 9001 format('can''t read definition file ',A40)
 9010 format('can''t open unit: ',I2)
 9020 format('       filename: ',A50)
 9030 format('         status: ',A,'  form: ',A)
 9040 format('input file  has wrong format or is corrupt!')
 9050 format('Not enough broadening parameters specified!')
      end

