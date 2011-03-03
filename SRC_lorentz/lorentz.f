      PROGRAM lorentz                                                     
!
!     Apply lorentz broadening to theoretical spectra
!
!     (C)1996 by Joachim Luitz
!
!     =====  Modifyed by H. Enkisch, 18.01.1999  =======================v
!     
!     Modifications are:
!
!     1. Consideration of the finite life time of the intermediate state
!        in the absorption spectra ==> Gamma also in the ABS case. 
!
!     2. Consideration of the spectrometer resolution by means of a
!        Gaussian instead of a Lorentzian.
!
!     ===================================================================^ 
!
      implicit real*8 (a-h,o-z)
      INCLUDE 'param.inc'

      EXTERNAL OUTERR

      CHARACTER*80 FNAME
      CHARACTER *11 STATUS,FORM

      CHARACTER*80 DEFFN,ERRFN
      CHARACTER*67 ERRMSG
      CHARACTER*4 DUMMY
      LOGICAL EMIS,ABS
      INTEGER IUNIT,IRECL
      DIMENSION X(IEMAX0),Y(IEMAX0)
      dimension yend(IEMAX0)
      data  pi/3.14159/

!     read filename for .def file

      call GTFNAM(DEFFN,ERRFN)
      CALL ERRFLG(ERRFN,'Error in LORENTZ')

      OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
 8003 CONTINUE
       READ (1,*,END=8001,ERR=911) IUNIT,FNAME,STATUS,FORM,IRECL
       IF(IUNIT.NE.6) then
       OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=912)
       END IF
      GOTO 8003
 8001 CONTINUE

!     read the theoretical spectrum
      DO i=1,IEMAX0+1
         READ(47,3000,ERR=2000,END=8004)X(i),Y(i)
      ENDDO
 8004 CONTINUE
      nimax=i
 3000 FORMAT(f10.5,e14.5)
 3001 FORMAT(f10.3,f10.3)
!     get parameters for broadening the spectrum
!     From our *.inxs we get the crucial data,
!     READ title into FNAME  and forget it
      READ(5,'(a)',END=914,ERR=914)FNAME
      READ(5,*,END=914,ERR=914)NATOM
      READ(5,*,END=914,ERR=914)NC
      READ(5,*,END=914,ERR=914)LC
      READ(5,*,END=914,ERR=914)SPLIT,XINT1,XINT2

!      READ(5,*,END=914,ERR=914)Kappa
      READ(5,*,END=914,ERR=914)EMIN,DE,EMAX
      READ(5,'(4a)',END=914,ERR=914)DUMMY
      IF(DUMMY.EQ.'EMIS')then
        EMIS=.TRUE.
      ELSE
        EMIS=.FALSE.
      ENDIF
!
      IF(DUMMY(1:3).EQ.'ABS')then
        ABS=.TRUE.
      ELSE
        ABS=.FALSE.
      ENDIF
!

        READ(5,*,END=915,ERR=915)S
        READ(5,*,END=915,ERR=915)G0
        EF=0.0d0
        IF(EMIS) THEN
           READ(5,*,END=915,ERR=915)W
           READ(5,'(4a)',END=914,ERR=914)DUMMY
           READ(5,*,END=915,ERR=915)E0
           READ(5,*,END=915,ERR=915)E1
           READ(5,*,END=915,ERR=915)E2
        ELSE
           E0=EMAX
           E1=EMAX
           E2=EMAX
        ENDIF
   
!     clear data, to be on the safe side
 202  do 300 i1=1,nimax
        yend(i1)=0
 300  continue

!     get energy step width
      delta=(X(2)-X(1))/2

!     Factor of 2 for core lifetime!!!
      g0 = g0/2
      
!     EMIS first:      
!     go into lifetime broadening, if g0.gt.0 and W.gt.0

      if(EMIS) then
         if ((g0.gt.0).or.(W.gt.0)) then

!      write(*,*)'Lifetime broadening'

!     now calculate Lifetime broadening using
!     Lorentz equation

            if(E0.NE.E2) then 
               write (6,*)"Doing 3-Parameter Range lifetime broadening"
            else
               write (6,*)"Doing 1-Parameter Range lifetime broadening"
            endif
            
            do 100 i1=1,nimax
               do 100 i2=1,nimax
                  if(E0.NE.E2) then 
                     if (X(i1).gt.E0) then
                        gamma=g0+W*(1-((X(i1)-E0)/(EF-E0)))**2
                     elseif (X(i1).gt.E1) then
                        gamma=g0+W
                     else
                        gamma=g0+W+W*(1-((X(i1)-E2)/(E1-E2)))**2
                     end if
                  else
                     gamma=g0+W*(1-((X(i1)-E0)/(EF-E0)))**2
                  endif
!     if (i2.eq.1) then
!     write(99,*)X(i1), gamma
!     endif

!               yend(i1)=yend(i1)+y(i2)/pi*
! Hartmut  Enkisch reported that in this line i1 and i2 are
! in wrong order and thus leading to wrong integral intensities.... 
! Thanks!
               yend(i2)=yend(i2)+y(i1)/pi* &
                    (atan((X(i1)-X(i2)+delta)/gamma) &
                    -(atan((X(i1)-X(i2)-delta)/gamma)))
          
 100        continue
      
!     write(*,*)'done'

!     broadened xspex goes to y for 2nd broadening
            do i=1,nimax
               y(i)=yend(i)
            enddo
            
         end if
      endif

!     now select for ABS
!     (this is the old spectrometer broadening!)
     
      if (ABS) then
        do 102 i1=1,nimax
          do 102 i2=1,nimax
           yend(i1)=yend(i1)+y(i2)/pi* &
                (atan((X(i1)-X(i2)+delta)/g0) &
               -(atan((X(i1)-X(i2)-delta)/g0)))
 102  continue

!     broadened xspex goes to y for 2nd broadening
            do i=1,nimax
               y(i)=yend(i)
            enddo

      endif
!
!     write(*,*)'Spectro. broadening'
!     now we do spectrometer broadening, if selected
      if (S.gt.0) then
!
!     S in case.inxs is FWHM, transform S into sigma
!
        S = S / 2.35

        do 101 i1=1,nimax
          yend(i1)=0
          do 101 i2=1,nimax
           yend(i1)=yend(i1)+y(i2)* &
              2*delta/( SQRT(2*pi) * S )* &
              EXP( ((X(i1)-X(i2))**2)/(-2.0*S**2) ) 
   
!


 101  continue
!      write(*,*)'done'
      end if

!     write header
      write(46,*) '#S=',S
      write(46,*) '#GAMMA=',G0
      if (EMIS) then
         write(46,*) '#W=',W
         write(46,*) '#E0=',E0
         write(46,*) '#E1=',E1
         write(46,*) '#E2=',E2
      endif
!     write out broadened spectrum
      DO i=1,nimax-1
         write(46,3000)X(i),YEND(I)
      ENDDO

      CALL ERRCLR(ERRFN)
      STOP 'Lorentz done'

!     Errors

!     xspec.def couldn't be opened
 910  WRITE(ERRMSG,9000) FNAME
      CALL OUTERR('LORENTZ',ERRMSG)
      GOTO 999

!     xspec.def is corrupt
 911  WRITE(ERRMSG,9001)FNAME
      CALL OUTERR('LORENTZ',ERRMSG)
      GOTO 999

!     Error opening a unit
 912  WRITE (ERRMSG,9010) IUNIT
      CALL OUTERR('LORENTZ',ERRMSG)
      WRITE (ERRMSG,9020) FNAME
      CALL OUTERR('LORENTZ',ERRMSG)
      WRITE (ERRMSG,9030) STATUS, FORM
      CALL OUTERR('LORENTZ',ERRMSG)
      GOTO 999

 914  WRITE(ERRMSG,9040)
      CALL OUTERR('LORENTZ',ERRMSG)
      GOTO 999

 915  WRITE(ERRMSG,9050)
      CALL OUTERR('LORENTZ',ERRMSG)
      GOTO 999

!     Too many data lines in Input file! Redimension IEMAX0!)
 2000 WRITE(ERRMSG,9090)
      CALL OUTERR('LORENTZ',ERRMSG)
      GOTO 999

 999  STOP 'LORENTZ - Error'


 9000 FORMAT('can''t open definition file ',A40)
 9001 FORMAT('can''t read definition file ',A40)
 9010 FORMAT('can''t open unit: ',I2)
 9020 FORMAT('       filename: ',A50)
 9030 FORMAT('         status: ',A,'  form: ',A)
 9040 FORMAT('CASE.INXS has wrong format or is corrupt!')
 9050 FORMAT('Not enough broadening parameters specified!')
 9090 FORMAT('Too many data lines in Input file! Redimension IEMAX0!')
      end


