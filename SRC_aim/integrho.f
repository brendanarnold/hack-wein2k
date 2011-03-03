      SUBROUTINE INTEGRHO()
!
!     AIM (1999-2000) by Javier D. Fuhr and Jorge O. Sofo
!
!     Instituto Balseiro and Centro Atomico Bariloche
!     S. C. de Bariloche - Rio Negro
!     Argentina
!     e-mail: fuhr@cab.cnea.gov.ar
!             sofo@cab.cnea.gov.ar
!
!     This routine integrates the electron density inside the
!       atomic surface already calculated.
!
!     Major changes, L. D. Marks, January 2006
!       1) Replaced radial integral by a quadrature
!       2) Added two alternative integrals for checking
!       3) Added error analysis code
!
      use sphe
      use zsphe
      implicit none
      include 'param.inc'

      real*8 cth,th,ph,wcth,wph,rs,volume,sumtheta
      real*8 shift,pi,tpi,rmin,rmft,rmedian,sumpw
      real*8 themin,themax,phimin,phimax,theta,phi,rmax,weight,frmin
      real*8 rsmin,rsmax,ct1,ct2,sum,sum2,f1,ErrorSum,ErrorAbs

      integer nth,nph,npr_use
      integer npr,index,i,j,ii

      character*4 switch

      logical srho,sgrho,shrho,deb
      logical weit,gaus
      
      COMMON /DEBUG/ deb
      COMMON /INTEG/ cth(NGAUSS),th(NGAUSS),ph(NGAUSS),wcth(NGAUSS), &
         wph(NGAUSS),rs(NGAUSS,NGAUSS),nth,nph
!      COMMON /SPHE/ POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO), &
!         RNOT(NATO),ndat,IATNR(NDIF),IOP(NDIF),MULT(NATO)
!      COMMON /ZSPHE/ ZZ(NDIF)
      dimension shift(3)
!     For sorting radial distances
      integer,allocatable :: tag(:)
      real*8,allocatable  :: radii(:)
!
!     For alternative integrations
      real*8 t1, t2,tt,sumint,sumpar,tol,t,wsum
      integer ierr,jj,nsize
      real*8 ans, result
      real*8,allocatable :: work(:),ftab(:),xtab(:)
      real*8,allocatable :: ftheta(:), xtheta(:),y(:,:),e(:),fthet2(:)

      pi=acos(-1.d0)
      tpi=2.0*pi

      read(5,112) switch
 112  format(a4)
      if(switch.eq.'WEIT') then
        weit=.true.
      else
        weit=.false.
      endif
      gaus=.false.
      if(.not.weit) then
        read(5,*) switch
        if(switch.eq.'GAUS') then
          gaus=.true.
        else
          gaus=.false.
        endif
      endif

      read(5,*) npr
      write(6,*) 'npr = ',npr
      srho=.true.
      sgrho=.false.
      shrho=.false.
      
      rewind(21)
      read(21,*) index,shift
      read(21,*) nth,themin,themax,frmin
      read(21,*) nph,phimin,phimax
      allocate (radii(nth*nph+20),tag(nth*nph+20))
!     Fix numerical accuracy errors in angles
      tol=1d-5
      do i=-8,8
        if(i.eq.0)then
                if(abs(themin) .lt. tol)themin=0.D0
                if(abs(themax) .lt. tol)themax=0.D0
                if(abs(phimin) .lt. tol)phimin=0.D0
                if(abs(phimax) .lt. tol)phimax=0.D0
        else
                t=tpi/dble(i)
                if(abs(themin-t) .lt. tol)themin=t
                if(abs(themax-t) .lt. tol)themax=t
                if(abs(phimin-t) .lt. tol)phimin=t
                if(abs(phimax-t) .lt. tol)phimax=t
        endif
      enddo
      write(6,*) 'NTH NPH ',nth,nph
!     Nsize is more than we really need, but it's still small
      nsize=max(nth,nph)+20
      allocate (y(3,nsize), e(nsize), work(nsize),ftab(nsize),xtab(nsize))
      allocate (ftheta(nsize), xtheta(nsize),fthet2(nsize) )

      ii=0
      do i=1,nth
        do j=1,nph
          if(weit) then
            read(21,*) theta,phi,rmax,weight            
          else
            read(21,*) theta,phi,rmax
          end if
          ii=ii+1
          radii(ii)=rmax
          tag(ii)=ii
        end do
      end do
!     Find the median rmax
      call SORTAG (radii, nth*nph, TAG)
      rmedian=radii(nth*nph/2)
      deallocate(radii)
      read(21,*) rsmin,rsmax

      rmft=rmt(iabs(iatnr(index)))

      rmin=rmft

      if(rmin.gt.rsmin) rmin=rsmin

      rewind(21)
      read(21,*) index,shift
      read(21,*) nth,themin,themax
      read(21,*) nph,phimin,phimax
      tol=1d-5
      do i=-8,8
        if(i.eq.0)then
                if(abs(themin) .lt. tol)themin=0.D0
                if(abs(themax) .lt. tol)themax=0.D0
                if(abs(phimin) .lt. tol)phimin=0.D0
                if(abs(phimax) .lt. tol)phimax=0.D0
        else
                t=tpi/dble(i)
                if(abs(themin-t) .lt. tol)themin=t
                if(abs(themax-t) .lt. tol)themax=t
                if(abs(phimin-t) .lt. tol)phimin=t
                if(abs(phimax-t) .lt. tol)phimax=t
        endif
      enddo
      if(gaus) then
        ct1=cos(themin)
        ct2=cos(themax)
        call gauleg(ct1,ct2,cth,wcth,nth)
!        call gauleg(phimin,phimax,ph,wph,nph)
        call gaulep(phimin,phimax,ph,wph,nph)
      endif
      sum2=0.0D0
!      sumpw=0.0D0
      call intrhospher(f1,index,rmin)
      sum2=sum2+f1
      write(6,*) 'RSMIN: ',rsmin,'  RMT:',rmft
      if (rsmin.lt.rmft*frmin) then
      write(6,*) ':ERROR    The Bader surface is within the atomic sphere, '
      write(6,*) ':ERROR    but your input parameter frmin',frmin,' is too big'
      write(6,*) ':ERROR    lower frmin and run again'
      stop 'Error: frmin too big'
      endif
      write(6,*) ':RHOSPHE ',rmin,sum2
      write(6,901)
901   format(/'Results for phi integrations:'/ &
        '    Theta       Quadrature        Spline    ',&
        '     Parabolic          Weight')
      sum=0.0
      volume=0.D0

      do i=1,nth
        sumtheta=0.D0
        wsum=0.D0
        do j=1,nph
          if(weit) then
            read(21,*) theta,phi,rmax,weight
          else
            read(21,*) theta,phi,rmax
            if(gaus) then
              weight=wcth(i)*wph(j)
            else
              weight=1.d0
            end if
          end if
!
!         Different numbers here
          npr_use=nint(npr*rmax/rmedian)
          npr_use=max(npr/2,npr_use)
!
!       Check, integrate PW's in RMT
!          call intPWonly(f1,theta,phi,rmft,index,npr)
!          sumpw=sumpw+weight*f1
          call intrhoint(f1,theta,phi,rmin,rmax,index,npr_use)
!          write(6,*) ':R ',theta,phi,rmin,rmax,f1,weight
          sumtheta=sumtheta+weight*f1
          volume=volume+rmax*rmax*rmax*weight/3.D0
          wsum=wsum+weight
!         For numerical integration
          ftab(j+3)=f1
          xtab(j+3)=phi
        end do
!       Reflect about origin
        if (abs (xtab(4)-phimin) .lt. 1D-8)then
            do jj=1,3
                ftab(4-jj)=ftab(4+jj)
                xtab(4-jj)=2.D0*phimin-xtab(4+jj)
            enddo
        else
            do jj=1,3
                ftab(4-jj)=ftab(3+jj)
                xtab(4-jj)=2.D0*phimin-xtab(3+jj)
            enddo
        endif
        if(abs (xtab(nph+3)-phimax) .lt. 1D-8)then
            do jj=1,3
                ftab(nph+3+jj)=ftab(nph+3-jj)
                xtab(nph+3+jj)=2.D0*phimax-xtab(nph+3-jj)
            enddo
        else
            do jj=1,3
                ftab(nph+3+jj)=ftab(nph+4-jj)
                xtab(nph+3+jj)=2.D0*phimax-xtab(nph+4-jj)
            enddo
        endif
!       Numerical integration via cubic spline
        call cspint ( ftab, xtab, nph+6, phimin, phimax, y, e, work, result )
!       Parabolic fit integration
        call davint ( xtab, ftab, nph+6, phimin, phimax, ans, ierr)
!       For second integration
        ftheta(i+3)=result
        xtheta(i+3)=cos(theta)
        fthet2(i+3)=ans
!       Rescaling adjustment
        write(6,410)theta,sumtheta*(phimax-phimin)/wsum,result,ans, &
                abs(wsum/(phimax-phimin))
410     format(F12.6,3D16.8,'  ',D16.8)
        sum=sum+sumtheta
      end do

!
        write(6,*)' '
!       Theta numerical integration
!       Reflect about origin
        t1=cos(themax)
        t2=cos(themin)
        if (abs (xtheta(4)-t2) .lt. 1D-8)then
            do jj=1,3
                ftheta(4-jj)=ftheta(4+jj)
                fthet2(4-jj)=fthet2(4+jj)
                xtheta(4-jj)=2.D0*t2-xtheta(4+jj)
            enddo
        else
            do jj=1,3
                ftheta(4-jj)=ftheta(3+jj)
                fthet2(4-jj)=fthet2(3+jj)
                xtheta(4-jj)=2.D0*t2-xtheta(3+jj)
            enddo
        endif
        if(abs (xtheta(nth+3)-t1) .lt. 1D-8)then
            do jj=1,3
                ftheta(nth+3+jj)=ftheta(nth+3-jj)
                fthet2(nth+3+jj)=fthet2(nth+3-jj)
                xtheta(nth+3+jj)=2.D0*t1-xtheta(nth+3-jj)
            enddo
        else
            do jj=1,3
                ftheta(nth+3+jj)=ftheta(nth+4-jj)
                fthet2(nth+3+jj)=fthet2(nth+4-jj)
                xtheta(nth+3+jj)=2.D0*t1-xtheta(nth+4-jj)
            enddo
        endif
!       Make sure that the values are in order, probably reversed
        do j=1,nth+6
                tag(j)=j
                work(j)=ftheta(j)
                ftab(j)=fthet2(j)
!                write(6,*)j,xtheta(j),acos(xtheta(j))
        enddo
        call sortag(xtheta,nth+6,tag)
        do j=1,nth+6
                jj=tag(j)
                ftheta(j)=work(jj)
                fthet2(j)=ftab(jj)
!                write(6,999)'Debug ',j,jj,xtheta(j),ftheta(j),fthet2(j)
!999             format(a,2i3,3F16.10)
        enddo
!       Make sure we are giving the limits right
        if(t2 .gt. t1 )then
                tt=t1
                t1=t2
                t2=tt
        endif
!       Spline integral
        call cspint ( ftheta, xtheta, nth+6, t1, t2, y, e, work, result )
!       Parabolic integral
        call davint ( xtheta, fthet2, nth+6, t2, t1, ans, ierr)
        result=result*(pi/(themin-themax))*(tpi/(phimax-phimin))
        ans   =-ans  *(pi/(themin-themax))*(tpi/(phimax-phimin))

      if (gaus.or.weit) then
        sum=sum*(pi/(themin-themax))*(tpi/(phimax-phimin))
        volume=volume*(pi/(themin-themax))*(tpi/(phimax-phimin))
!        sumpw=sumpw*(pi/(themin-themax))*(tpi/(phimax-phimin))
      else
        sum=sum/(nth*nph)*2.0*tpi
        volume=volume/(nth*nph)*2.0*tpi
!        sumpw=sumpw/(nth*nph)*2.0*tpi
      endif
!
!     Gradient checking
!     This one looks wrong, fails symmetry tests
!      call CheckBspline(ErrorSum,ErrorAbs)
      call CheckCspline(ErrorSum,ErrorAbs,weit)
      write(6,*)' '
!
600   continue
!      write(6,303) 'Plane component inside sphere',sumpw
!      write(6,*)' '
      write(6,303) 'Integrated interstial charges:'
      write(6,303) ':RHO from Cubic Spline ',result
      write(6,303) ':RHO from Parabolic    ',ans
      write(6,303) ':RHOINTE (Quadrature)  ',sum,' :VOLUME ',volume
      write(6,*)' '
      sumint=sum2+result
      sumpar=sum2+ans
      sum2=sum2+sum
      write (6,404) index,zz(index),sumpar,zz(index)-sumpar
      write (6,403) index,zz(index),sumint,zz(index)-sumint

      write (6,302) index,zz(index),sum2,zz(index)-sum2

      return

 200  format(3F10.8,2F15.10)
 201  format(':RSUR ',3F10.8,2F15.10)
 202  format(2F15.10)
 300  FORMAT(i3,3F15.10)
 301  FORMAT(i3,2F15.10)
 302  FORMAT(':RHOTOT     for IND-ATOM',i4,'  Z=',f5.1,'  CHARGE:',2F13.8)
 303  format(a,F13.8,a,f12.6)
 403  FORMAT(':CUBESPLINE for IND-ATOM',i4,'  Z=',f5.1,'  CHARGE:',2F13.8)
 404  FORMAT(':PARABOLIC  for IND-ATOM',i4,'  Z=',f5.1,'  CHARGE:',2F13.8)
      end
!
      subroutine intrhospher(f1,index,rmin)
      use rad
      use sphe
      use blank
      implicit none
      include 'param.inc'

      real*8 f1,rmin
      real*8 pi,sum

      integer index
      integer jatom,ilm,j

      logical deb

      COMMON /DEBUG/ deb
!      COMMON /RAD/ RM(NRAD,NATO),JRI(NATO)                       
!      COMMON /SPHE/ POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO), &
!         RNOT(NATO),ndat,IATNR(NDIF),IOP(NDIF),MULT(NATO)
!      COMMON // CLM(NRAD,NCOM,NATO),CLM2(NRAD,NCOM,NATO), &
!         LM(2,NCOM,NATO),LMMAX(NATO)      

      pi=acos(-1.0D0)                     
      jatom=iabs(iatnr(index))
      sum=0.0D0

      ilm=0
      j=1

!$$$      do i=1,lmmax(jatom)
!$$$        l=lm(1,i,jatom)
!$$$        if(l.eq.0) then
!$$$          ilm=i
!$$$          goto 25
!$$$        end if
!$$$      end do

      ilm=1

 25   if (ilm.ne.0) then
        
        do while((.not.(rm(j+2,jatom).gt.rmin)).and. &
           ((j+1).lt.jri(jatom)))
          sum=sum+clm(j,ilm,jatom)*rm(j,jatom)+ &
             4.0D0*clm(j+1,ilm,jatom)*rm(j+1,jatom)+ &
             clm(j+2,ilm,jatom)*rm(j+2,jatom)
          j=j+2
        end do

        if(j.gt.jri(jatom)) j=jri(jatom)
        rmin=rm(j,jatom)
        sum=sum*dx(jatom)/3.0D0
        sum=sum*2.0D0*sqrt(pi)
      end if
        
!$$$      sum = 0.0                                                        
!$$$      npoint = jri(jatom)
!$$$      js = 2-mod(npoint,2)
!$$$      jf = jri(jatom)-2
!$$$
!$$$      do  j=js,jf,2                                                    
!$$$        sum = sum+clm(j,ilm,jatom)*rm(j,jatom)+4.*clm(j+1,ilm,jatom)*
!$$$     $     rm(j+1,jatom)+clm(j+2,ilm,jatom)*rm(j+2,jatom)
!$$$      enddo
!$$$      sum = sum*dx(jatom)/3.0                                              
!$$$      sum = sum + (1.-mod(npoint,2))*0.5*
!$$$     $   (clm(1,ilm,jatom)+clm(js,ilm,jatom))  
!$$$     $   *(rm(js,jatom)-rm(1,jatom))       
      
!       write(*,*) 'jatom,index,jri,rnot,rmt,rmin',jatom,index,jri(jatom),rm(1,jatom),rm(jri(jatom),jatom),rmin
      f1=sum

      return 
      end
!
      subroutine intrhoint(f1,theta,phi,r0,rmax,index,npr)
      use rad
      use sphe
      implicit none
      include 'param.inc'

      real*8 f1,theta,phi,r0,rmax
      real*8 v,rho,grho,hrho,cutoff,ct,st,cf,sf,dr,r,sum,rr
      
      integer index,npr,ibest,ipbest
      integer i,ipos,iat,ist

      logical srho,sgrho,shrho,deb

      COMMON /DEBUG/ deb
!      COMMON /RAD/ RM(NRAD,NATO),JRI(NATO)                       
!      COMMON /SPHE/ POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO), &
!         RNOT(NATO),ndat,IATNR(NDIF),IOP(NDIF),MULT(NATO)
      DIMENSION v(3),grho(3),hrho(3,3)
      real*8,allocatable :: rwts(:),rvals(:)
      allocate (rwts(npr),rvals(npr))

      srho=.true.
      sgrho=.false.
      shrho=.false.
      cutoff=0.0D0

      ct=cos(theta)
      st=sin(theta)
      cf=cos(phi)
      sf=sin(phi)

      if(rmax.lt.r0) then
        f1=0.0D0
      else
!       New version, radial quadrature
        call gauleg(r0,rmax,rvals,rwts,npr)
        do i=1,npr
               r=rvals(i)
               v(1)=pos(1,index)+r*st*cf
               v(2)=pos(2,index)+r*st*sf
               v(3)=pos(3,index)+r*ct
               call vgh_rho(v,rho,grho,hrho,srho,sgrho,shrho,rr,iat,ipos,ist,ibest,ipbest)
               sum=sum+rho*r*r*rwts(i)
        enddo
        f1=sum
        deallocate (rwts,rvals)
!
!        Old version
!        dr=(rmax-r0)/npr
!        r=r0+dr*0.5D0  !/2.0
!        sum=0.0D0
!        do i=1,npr
!          v(1)=pos(1,index)+r*st*cf
!          v(2)=pos(2,index)+r*st*sf
!          v(3)=pos(3,index)+r*ct
!          call vgh_rho(v,rho,grho,hrho,srho,sgrho,shrho,rr,iat,ipos,ist,ibest,ipbest)
!          sum=sum+rho*r*r
!          r=r+dr
!        end do
!        f1=sum/npr*(rmax-r0)
      end if

      return
      end
      subroutine intPWonly(f1,theta,phi,myrmt,index,npr)
!     For integrating PW's only
      use rad
      use sphe
      implicit none
      include 'param.inc'

      real*8 f1,theta,phi,r0,myrmt
      real*8 v,rho,grho,hrho,cutoff,ct,st,cf,sf,dr,r,sum,rr
      
      integer index,npr,ibest,ipbest
      integer i,ipos,iat,ist

      logical srho,sgrho,shrho,deb

      COMMON /DEBUG/ deb
!      COMMON /RAD/ RM(NRAD,NATO),JRI(NATO)                       
!      COMMON /SPHE/ POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO), &
!         RNOT(NATO),ndat,IATNR(NDIF),IOP(NDIF),MULT(NATO)
      DIMENSION v(3),grho(3),hrho(3,3)
      real*8,allocatable :: rwts(:),rvals(:)
      allocate (rwts(npr),rvals(npr))

      srho=.true.
      sgrho=.false.
      shrho=.false.
      cutoff=0.0D0

      ct=cos(theta)
      st=sin(theta)
      cf=cos(phi)
      sf=sin(phi)
        call gauleg(0.0D0,myrmt,rvals,rwts,npr)
        do i=1,npr
               r=rvals(i)
               v(1)=pos(1,index)+r*st*cf
               v(2)=pos(2,index)+r*st*sf
               v(3)=pos(3,index)+r*ct
               call INTERST(v,rho,grho,hrho,srho,sgrho,shrho)
               sum=sum+rho*r*r*rwts(i)
        enddo
        f1=sum
        deallocate (rwts,rvals)
      return
      end

