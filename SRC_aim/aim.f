      PROGRAM AIM                                                     
!
!     AIM (1999-2000) by Javier D. Fuhr and Jorge O. Sofo
!
!     Instituto Balseiro and Centro Atomico Bariloche
!     S. C. de Bariloche - Rio Negro
!     Argentina
!     e-mail: fuhr@cab.cnea.gov.ar
!             sofo@cab.cnea.gov.ar
!
!     Implementation of Bader's Atom In Molecules analysis
!       of electron density.
!
!     AIM can:
!       - Find critical points.
!       - Follow a gradient line starting from any point.
!       - Obtain the atomic surfaces.
!       - Calculate the radious of the atomic surface for any given
!         angles theta and phi.
!       - Integrate the electron density inside the atomic surface.
!       - Calculate the dipole inside the atomic 
!         surface (not well tested).
!
!     TAPE 5  = case.inb       (input) 
!     TAPE 6  = case.outputaim (output)
!     TAPE 8  = case.struct    (input)
!     TAPE 9  = case.clmsum    (input)
!     TAPE 21 = case.surf      (output/input)
!     TAPE 22 = case.crit      (output)
!     TAPE 23 = case.surfin    (input)
!

      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*80      FNAME
      CHARACTER*11      STATUS,FORM
      iarg=iargc()
      call getarg(iarg,fname)
!      if(fname.eq.'      ') call getarg(1,fname)
      OPEN(1,FILE=fname,STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM, &
         ERR=8002)
      GOTO 8003
 8000 WRITE(*,*) ' ERROR IN OPENING AIM.DEF !!!!'
      STOP 'AIM.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, &
      '  FORM:',FORM
      STOP 'OPEN FAILED'
 8001 CONTINUE
        call cputim(t1)   
      CALL MAIN1()
        call cputim(t2)
      write(6,*) 'total cputime:', t2-t1   
      STOP                                                              
      END                                                               
