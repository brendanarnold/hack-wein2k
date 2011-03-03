     PROGRAM  KRAM

!
!     original program by Peter Vogl
!     adapted for LAPW optic package by Robert Abt
!
!     LAST UPDATES: August   1999  Claudia Ambrosch-Draxl
!                   December 1999  CAD & Zhenya Sherman
!                   December 2002  CAD: implementation of relectivity and complex refractive index

      implicit real*8 (a-h,o-z)
!ad
      include 'param.inc'
!ad
      parameter (pi    = 3.141592654d0)
      parameter (h     = 6.6262d-34)
      parameter (hbar  = h/2/pi)
      parameter (eps0  = 8.854d-12)
      parameter (c     = 2.997925d+8)
      parameter (eV2J  = 1.602177d-19)
      parameter (sigf  = eV2J/h/2 *1.d-15)
      parameter (sigf1 = eV2J/hbar*eps0*1.d-2)

!ad

      REAL*8  xe(maxde),ffx(maxde,mpol),absk(maxde,mpol)
      REAL*8  wpl(mpol),gampl(mpol)
      REAL*8  sumr1(mpol),sumr2(mpol),sumr3(mpol)
      REAL*8  eloss(maxde,mpol),reflect(maxde,mpol)
      REAL*8  cn(maxde,mpol),ck(maxde,mpol)
      REAL*8  fhelp1(maxde),fhelp2(maxde),diag(mpol)
!ad
      COMPLEX*16 eps(maxde,mpol),fx(maxde,mpol),sig(maxde,mpol)
      COMPLEX*16 iimag,ione,chelp
!ad
      CHARACTER*11   STATUS, FORM
      CHARACTER*80   FNAME,DEFFN,ERRFN,SYSTEM
      CHARACTER*5    ecv
      CHARACTER*6    hloss
      CHARACTER*9    header(mpol),hsim,hsre
      CHARACTER*8    hcn,hck,hrf,habs
      CHARACTER*2    hpol(mpol)
      CHARACTER*7    hdim,hdre
!ad
      absf  = 1.0d-6*dsqrt(2.d0)*eV2J/hbar/c
      iimag=(0D0,1D0)
      ione=(1D0,0D0)
!*                                                               *  
!*     absk=dsqrt(2)*omega/(hbar c) * dsqrt(-e1+dsqrt(e1^2+e2^2))   *
!*     if omega is given in eV !                                 *
!*     -->  abkfac = eV2J / hbar / c [SI] * 10^-6                *
!*     -->  absk = [ 10^4 cm^-1]                                 *
!*     sigma in (Ohm cm)^-1  (case.absorp)                       *
!*     sigma in 10^15 sec^-1 (case.sigmak)                       *
!*                                                               *
!ad
!ad
!ad ____________________________ OPEN FILES ____________________________
!ad
!ad
      CALL GTFNAM(DEFFN,ERRFN)
      OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
   10 CONTINUE
         READ (1,*,END=20,ERR=920) IUNIT,FNAME,STATUS,FORM,IRECL
         OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
      GOTO 10
   20 CONTINUE
      CLOSE (1)
!ad
!ad ____________________________ READ INPUT ____________________________
!ad
!ad
      read(5,*) gamma
      read(5,*) esh
      read(5,*) iadd
!ad
!ad _________________________ READ case.JOINT __________________________
!ad
!ad ___________________________ read headers ___________________________
!ad
!ad
      read(10,330) NCOL,VOL
      read(10,331) ecv,(header(j),j=1,ncol)

      do j=1,ncol

      hpol(j)=header(j)(8:9)
      write(*,*) hpol(j)
      diag(j)=0.0d0
      if(hpol(j).eq.'xx') diag(j)=1.0d0
      if(hpol(j).eq.'yy') diag(j)=1.0d0
      if(hpol(j).eq.'zz') diag(j)=1.0d0
!ad
      if(header(j)(1:2).eq.'Re'.and.gamma.lt.1.d-15) then
      gamma=0.01
      write(6,332) gamma
      endif
!ad
      if(header(j)(1:2).eq.'Re'.and.abs(esh).gt.1.d-15) then
      write(6,333) 
      endif

      enddo

      hdim='Im_eps_'
      hdre='Re_eps_'
      hsim='Im_sigma_'
      hsre='Re_sigma_'
      habs=' absorp_ '
      hloss='eloss_'
      hcn='ref_ind_'
      hck='extinct_'
      hrf='reflect_'
!ad
!ad ________________________ define energy units _______________________
!ad
!ad
      write(*,*) ' Energy units: ',ecv
!ad
      if(ecv.eq.' [eV]') then
      efactor=1.0d0
      ec=13.6057d0
      elseif(ecv.eq.' [Ry]') then
      efactor=13.6057d0
      ec=1.0d0
      else
      write(*,*) ' No valid energy units!'
      stop
      endif
      sigfac=efactor*sigf
      absfac=efactor*absf
      sumfac=8*pi*pi/VOL*ec*ec
!ad
      write(*,*) ' Lorentzian broadening with gamma: ',gamma,ecv
!ad
!ad
!ad ___________________________ read columns ___________________________
!ad
!ad
      read(10,*)
      i=1
 101  read(10,15,END=100)  xe(i),(ffx(i,j),j=1,ncol)
      i=i+1
      if(i.gt.maxde)  then
      write(6,312) i-1
      stop
      endif
      goto 101
 100  continue
      write(*,*) i-1,' data points'
      close(10)
!ad
      negrid = i-1
!ad
!ad
!ad ______________________ apply scissors operator _____________________
!ad
!ad
      if (abs(esh).gt.0.d0) call eshift(xe,ffx,negrid,esh,ncol,ecv)
      egrid  = xe(2)-xe(1)
      write(*,*) ' ENERGY INCREMENT: ',egrid  
!ad
!ad
!ad ______________________ KRAMERS KRONIG ANALYSIS _____________________
!ad
!ad  prepare intraband contributions for calculating the loss function
!ad
      if (iadd.eq.1) then
      read(5,*) (wpl(icol),icol=1,ncol)
      read(5,*) (gampl(icol),icol=1,ncol)
      write(*,*) 'intraband contributions prepared'
      endif
!ad 
!ad _________________________ without broadening _______________________
!ad
      if (gamma.lt.1D-15) then
!ad   
!ad   initialize sum rule check
!ad
      do j=1,ncol
      sumr1(j)=0.d0
      sumr2(j)=0.d0
      sumr3(j)=0.d0
      enddo
!ad

       do i = 2,negrid-2,2
       do j = 1,ncol
         xx = xe(i)
         call green(fer1,ffx(1,j),xe,egrid,negrid,i)
         call green_reg(fer2,ffx(1,j),xe,egrid,negrid,-xx)
         epskk = diag(j) - 1./pi * (fer1 + fer2)
!ad
       el1help=0.d0
       el2help=0.d0
!ad
!ad    add intraband contributions 
!ad
       if(iadd.eq.1) then
       wpl2=wpl(j)**2
       el1help=eps1(wpl2,gampl(j),xe(i))-1
       el2help=eps2(wpl2,gampl(j),xe(i))
       endif
!
!cad   eps(i,j)=dcmplx(epskk,ffx(i,j))
!
       el2=ffx(i,j)+el2help
       el1=epskk   +el1help
       eps(i,j)=dcmplx(el1,el2)
!
       sig(i,j)=-iimag*sigfac*eps(i,j)*xe(i)
       absk(i,j)=absfac*xe(i)*dsqrt(-dble(eps(i,j))+abs(eps(i,j)))

       abse2=el2**2+el1**2
       abse=sqrt(abse2)
       cnh=sqrt(0.5*(abse+el1))
       ckh=sqrt(0.5*(abse-el1))
       rh=((cnh-1)**2+ckh**2)/((cnh+1)**2+ckh**2)
!
       cn(i,j)=cnh
       ck(i,j)=ckh
       reflect(i,j)=rh
!
       eloss(i,j)=el2/abse2
       sumr1(j)=sumr1(j)+el2*xe(i)*2*egrid/sumfac
       sumr2(j)=sumr2(j)+eloss(i,j)*xe(i)*2*egrid/sumfac
       sumr3(j)=sumr3(j)+eloss(i,j)/xe(i)*2*egrid
!cad
       enddo
       write(77,77) xe(i),(sumr1(j),sumr2(j),sumr3(j),j=1,ncol)
       enddo
       write(*,*) 'sum rule 1: Int(sigma)dw   ',(sumr1(j),j=1,ncol)
       write(*,*) 'sum rule 2: Int(eloss.w)dw ',(sumr2(j),j=1,ncol)
       write(*,*) 'sum rule 3: Int(eloss/w)dw ',(sumr3(j),j=1,ncol)
       write(*,*) 'KK without broadening done'
!ad
!ad
!ad __________________________ with broadening _________________________
!ad
!ad   initialize sum rule check
!ad
      do j=1,ncol
      sumr1(j)=0.d0
      sumr2(j)=0.d0
      sumr3(j)=0.d0
      enddo
!ad
!ad __________________ distinguish real and imaginary parts ____________
!ad
!ad
      else
!ad
      do j=1,ncol
      IF (header(j)(1:2).eq.'Re') THEN
!ad
      do i=1,negrid
      fhelp1(i)=ffx(i,j)
      fhelp2(i)=0.0d0
      enddo

      write(*,*) 'losmo2 called to perform KK for Re to Im'
      call losmo2(gamma,xe,fhelp1,fhelp2,negrid)
      call losmo(gamma,xe,fhelp1,ffx(1,j),negrid)
!ad
      ELSE IF (header(j)(1:2).eq.'Im') THEN
!ad
      do i=1,negrid
      fhelp1(i)=0.0d0
      fhelp2(i)=ffx(i,j)
      enddo

      write(*,*) 'losmo1 called to perform KK for Im to Re'
      call losmo1(gamma,xe,fhelp1,fhelp2,negrid,diag(j))
      call losmo(gamma,xe,fhelp2,ffx(1,j),negrid)
!ad
      ELSE
      write(*,109) j
 
      END IF

      do i=2,negrid-2,2
      eps(i,j)=dcmplx(fhelp1(i),fhelp2(i))
      sig(i,j)=-iimag*sigfac*eps(i,j)*xe(i)
      absk(i,j)=absfac*xe(i)*dsqrt(-dble(eps(i,j))+abs(eps(i,j)))
!cad
      abse2=fhelp1(i)**2+fhelp2(i)**2
      abse=sqrt(abse2)
      cnh=sqrt(0.5*(abse+fhelp1(i)))
      ckh=sqrt(0.5*(abse-fhelp1(i)))
      rh=((cnh-1)**2+ckh**2)/((cnh+1)**2+ckh**2)
!
      cn(i,j)=cnh
      ck(i,j)=ckh
      reflect(i,j)=rh
!

      
      eloss(i,j)=fhelp2(i)/(fhelp2(i)**2+fhelp1(i)**2)
      enddo

      enddo
!ad
!ad    add intraband contributions for loss function
!ad
       if(iadd.eq.1) then
       do i = 2,negrid-2,2
       do j = 1,ncol
!ad
!ad
       el1help=0.d0
       el2help=0.d0
       wpl2=wpl(j)**2
       el1help=eps1(wpl2,gampl(j),xe(i))-1
       el2help=eps2(wpl2,gampl(j),xe(i))
!ad
       el2=aimag(eps(i,j))+el2help
       el1=dble(eps(i,j))+el1help

       abse2=el2**2+el1**2
       abse=sqrt(abse2)
       cnh=sqrt(0.5*(abse+el1))
       ckh=sqrt(0.5*(abse-el1))
       rh=((cnh-1)**2+ckh**2)/((cnh+1)**2+ckh**2)
!
       cn(i,j)=cnh
       ck(i,j)=ckh
       reflect(i,j)=rh
!
       eloss(i,j)=el2/abse2
       sumr1(j)=sumr1(j)+el2*xe(i)*2*egrid/sumfac
       sumr2(j)=sumr2(j)+eloss(i,j)*xe(i)*2*egrid/sumfac
       sumr3(j)=sumr3(j)+eloss(i,j)/xe(i)*2*egrid
!cad
!cad   add intraband contributions also to epsilon tensor
!cad
       eps(i,j)=dcmplx(el1,el2)
       sig(i,j)=-iimag*sigfac*eps(i,j)*xe(i)
!ad
!ad
       enddo
       write(77,77) xe(i),(sumr1(j),sumr2(j),sumr3(j),j=1,ncol)
       enddo
       write(*,*) 'sum rule 1: Int(sigma)dw   ',(sumr1(j),j=1,ncol)
       write(*,*) 'sum rule 2: Int(eloss.w)dw ',(sumr2(j),j=1,ncol)
       write(*,*) 'sum rule 3: Int(eloss/w)dw ',(sumr3(j),j=1,ncol)
       write(*,*) 'KK with broadening done'
       endif
!ad
      endif
!ad
!ad ___________________________ WRITE OUTPUT ___________________________
!ad
!ad ___________________ header for dielectric tensor __________________
!ad
!ad

         system=' '
         write(12,110) system,gamma,ecv,esh,ecv
         if(iadd.ne.1) write(12,114)
         if(iadd.eq.1) then 
           write(12,115) (wpl(icol),icol=1,ncol)
           write(12,116) (gampl(icol),icol=1,ncol)
         endif
         write(12,111) ecv,(hdre,hpol(j),hdim,hpol(j),j=1,ncol)
         write(12,13)
!ad
!ad ___________________ header for optical conductivity ________________
!ad
!ad
         write(13,110) system,gamma,ecv,esh,ecv
         if(iadd.ne.1) write(13,114)
         if(iadd.eq.1) then
           write(13,115) (wpl(icol),icol=1,ncol)
           write(13,116) (gampl(icol),icol=1,ncol)
         endif
         write(13,210) 
         write(13,112) ecv,(hsre,hpol(j),hsim,hpol(j),j=1,ncol)
         write(13,13) 
!ad
!ad __________________ header for absorption coefficient _______________
!ad
!ad
         write(14,110) system,gamma,ecv,esh,ecv
         if(iadd.ne.1) write(14,114)
         if(iadd.eq.1) then
           write(14,115) (wpl(icol),icol=1,ncol)
           write(14,116) (gampl(icol),icol=1,ncol)
         endif
         write(14,211)
         write(14,112) ecv,(hsre,hpol(j),j=1,ncol),(habs,hpol(j),j=1,ncol)
         write(14,13)
!ad
!ad ____________ header for loss function and reflectivity _____________
!ad
!ad
         write(15,110) system,gamma,ecv,esh,ecv
         if(iadd.ne.1) write(15,114)
         if(iadd.eq.1) then
           write(15,115) (wpl(icol),icol=1,ncol)
           write(15,116) (gampl(icol),icol=1,ncol)
         endif
         write(15,112) ecv,(hloss,hpol(j),j=1,ncol)
         write(15,13)
!ad
!ad ______________ header for complex refractive index  ________________
!ad
!ad
         write(16,110) system,gamma,ecv,esh,ecv
         if(iadd.ne.1) write(16,114)
         if(iadd.eq.1) then
           write(16,115) (wpl(icol),icol=1,ncol)
           write(16,116) (gampl(icol),icol=1,ncol)
         endif
         write(16,112) ecv,(hcn,hpol(j),j=1,ncol),(hck,hpol(j),j=1,ncol)
         write(16,13)
!ad
!ad ____________ header for loss function and reflectivity _____________
!ad
!ad
         write(17,110) system,gamma,ecv,esh,ecv
         if(iadd.ne.1) write(17,114)
         if(iadd.eq.1) then
           write(17,115) (wpl(icol),icol=1,ncol)
           write(17,116) (gampl(icol),icol=1,ncol)
         endif
         write(17,112) ecv,(hrf,hpol(j),j=1,ncol)
         write(17,13)
!ad
!ad _________________________ write out RESULTS ________________________
!ad
!ad
      do 700 i=2,negrid-2,2
      write(12,113) xe(i),(DBLE(eps(i,j)),aimag(eps(i,j)),j=1,ncol)
      write(13,113) xe(i),(DBLE(sig(i,j)),aimag(sig(i,j)),j=1,ncol)
      write(14,113) xe(i),(DBLE(sig(i,j))*sigf1/sigf,j=1,ncol), &
                              (absk(i,j),            j=1,ncol)
      write(15,113) xe(i),(eloss(i,j),j=1,ncol)
      write(16,113) xe(i),(cn(i,j),j=1,ncol),(ck(i,j),j=1,ncol)
      write(17,113) xe(i),(reflect(i,j),j=1,ncol)
  700 continue
!ad
      close(12)
      close(13)
      close(14)
      close(15)
!ad
      goto 999
  910 write(*,*) ' kram.def could not be opened'
      goto 999
  920 write(*,*) fname, ' could not be opened - check def-file'
      goto 999 
!     stop ' Kramers-Kronig analysis done'
  999 continue
!ad
!ad ______________________________ FORMATS _____________________________
!ad
!ad
   11 format(9x,a6)
   13 format('#')
   15 format(f13.5,9(1x,e18.8))
  109 format(1x,'Error reading header of column ',I2)
  110 format('#',A80,/,'# Lorentzian broadening with gamma= ',&
      F8.6,1x,A5,/,'# Im(epsilon) shifted by ',F8.4,1x,a6)
  111 format('# Energy',A5,1x,9(A7,A2,5x,A7,A2,5x))
  112 format('# Energy',A5,1x,9(A9,A2,3x,A9,A2,3x))
  113 format(f11.6,18e14.6)
  114 format('# No intraband contributions added',/,'#')
  115 format('# Intraband contributions added: w_p=',6f8.3) 
  116 format('#                              Gamma=',6f8.3)   
  210 format('# optical conductivity sigma in [10^15 / sec]',/,'#')
  211 format('# optical conductivity sigma in [1 / (Ohm cm)]', &
             6x,' absorption in [10^4 / cm)]',/,'#')
  212 format('# loss function')
  312 format(' Parameter MAXDE should be greater than ',i8)
  330 format(1x,I1,36x,F18.10)
  331 format(10x,A5,5x,8(A9,10x),A9)
  332 format(/,' WARNING: GAMMA set to ',f6.4,/)
  333 format(/,' WARNING:',/, & 
          ' Scissors operator should not be applied to real parts!',/)
 5000 format(A80)
   77 format(f13.6,2x,4(3f9.4,2x))

      end



