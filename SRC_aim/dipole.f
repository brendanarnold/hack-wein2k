      SUBROUTINE DIPOLE()

!
!     AIM (1999-2000) by Javier D. Fuhr and Jorge O. Sofo
!
!     Instituto Balseiro and Centro Atomico Bariloche
!     S. C. de Bariloche - Rio Negro
!     Argentina
!     e-mail: fuhr@cab.cnea.gov.ar
!             sofo@cab.cnea.gov.ar
!
!     Calculates dipole from the atomic surface and
!       the electron density.
!
!     It needs the atomic surface already calculated or read
!       from case.surf. 
!     Small changes by L. D. Marks, January 2006
!
      use sphe
      implicit none
      include 'param.inc'

      real*8 cth,th,ph,wcth,wph,rs
      real*8 shift,pi,tpi,rmin,rmft
      real*8 themin,themax,phimin,phimax,theta,phi,rmax,weight
      real*8 rsmin,rsmax,ct1,ct2,dipx,dipy,dipz
      real*8 sumx,sumy,sumz,spth(PARTMAX),spfi(PARTMAX)
      real*8 sum2x,sum2y,sum2z,the,fi,dipr,dth,dfi

      integer nth,nph,ipth,ipfi,ibest,ipbest
      integer npr,index,i,j,npth,npfi

      character*4 switch

      logical srho,sgrho,shrho,deb
      logical weit,gaus
      
      COMMON /DEBUG/ deb
      COMMON /INTEG/ cth(NGAUSS),th(NGAUSS),ph(NGAUSS),wcth(NGAUSS), &
         wph(NGAUSS),rs(NGAUSS,NGAUSS),nth,nph
!      COMMON /SPHE/ POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO), &
!         RNOT(NATO),ndat,IATNR(NDIF),IOP(NDIF),MULT(NATO)
      dimension shift(3)

      pi=acos(-1.d0)
      tpi=2.0*pi

      read(5,*) switch
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
      read(21,*) nth,themin,themax
      read(21,*) nph,phimin,phimax

      write(6,*) 'NTH NPH ',nth,nph

      rewind(23)
      read(23,*) npth
      write(6,*) 'NPTH=',npth
       if(npth.gt.partmax) stop 'npth'
      do i=1,npth
        read(23,*) spth(i)
        write(6,*) 'SPTH(',i,')=',spth(i)
      enddo
      read(23,*) npfi
      write(6,*) 'NPFI=',npfi
       if(npfi.gt.partmax) stop 'npfi'
      do i=1,npfi
        read(23,*) spfi(i)
        write(6,*) 'SPFI(',i,')=',spfi(i)
      enddo


      do i=1,nth
        do j=1,nph
          if(weit) then
            read(21,*) theta,phi,rmax,weight
          else
            read(21,*) theta,phi,rmax
          end if
        end do
      end do
      read(21,*) rsmin,rsmax

      rmft=rmt(iabs(iatnr(index)))

      rmin=rmft

      if(rmin.gt.rsmin) rmin=rsmin

      rewind(21)
      read(21,*) index,shift
      read(21,*) nth,themin,themax
      read(21,*) nph,phimin,phimax

      if(gaus) then
        ct1=cos(themin)
        ct2=cos(themax)
        call gauleg(ct1,ct2,cth,wcth,nth)
!        call gauleg(phimin,phimax,ph,wph,nph)
        call gaulep(phimin,phimax,ph,wph,nph)
      endif
      
      sumx=0.0
      sumy=0.0
      sumz=0.0
      call dipolespher(dipx,dipy,dipz,index,rmin)
      sumx=sumx+dipx
      sumy=sumy+dipy
      sumz=sumz+dipz
      write(6,*) ':DIPOLESPHE ',rmin,sumx,sumy,sumz

      sum2x=0.0
      sum2y=0.0
      sum2z=0.0
      do i=1,nth
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
          call dipoleint(dipr,theta,phi,rmin,rmax,index,npr)
          sum2x=sum2x+dipr*sin(theta)*cos(phi)*weight
          sum2y=sum2y+dipr*sin(theta)*sin(phi)*weight
          sum2z=sum2z+dipr*cos(theta)*weight
          dfi=0.0
          do ipfi=1,npfi
            dfi=dfi+(phimax-phimin)
            fi=spfi(ipfi)*phi+dfi+0.5*(1.-spfi(ipfi))*(phimax-phimin)
            sum2x=sum2x+dipr*sin(theta)*cos(fi)*weight
            sum2y=sum2y+dipr*sin(theta)*sin(fi)*weight
            sum2z=sum2z+dipr*cos(theta)*weight
          enddo
          dth=0.0D0
          do ipth=1,npth
            dth=dth+(themax-themin)
            the=spth(ipth)*theta+dth+0.5*(1.-spth(ipth))*(themax-themin)
            sum2x=sum2x+dipr*sin(the)*cos(phi)*weight
            sum2y=sum2y+dipr*sin(the)*sin(phi)*weight
            sum2z=sum2z+dipr*cos(the)*weight
            dfi=0.0
            do ipfi=1,npfi
              dfi=dfi+(phimax-phimin)
              fi=spfi(ipfi)*phi+dfi+0.5*(1.-spfi(ipfi))*(phimax-phimin)
              sum2x=sum2x+dipr*sin(the)*cos(fi)*weight
              sum2y=sum2y+dipr*sin(the)*sin(fi)*weight
              sum2z=sum2z+dipr*cos(the)*weight
            enddo
          enddo
        enddo
      enddo

      sumx=sumx+sum2x
      sumy=sumy+sum2y
      sumz=sumz+sum2z
      write(6,*) ':DIPOLE INTE ',sum2x,sum2y,sum2z
      write(6,*) ':DIPOTOT ',index,sumx,sumy,sumz

      return

 200  format(3F10.8,2F15.10)
 201  format(':RSUR ',3F10.8,2F15.10)
 202  format(2F15.10)
 300  FORMAT(i3,3F15.10)
 301  FORMAT(i3,2F15.10)
      end


      subroutine dipolespher(dipx,dipy,dipz,index,rmin)
      use sphe
      use rad
      use blank
      implicit none
      include 'param.inc'

      complex*16 y,imag,imag1

      real*8 rmin
      real*8 pi,cth(NGAUSS),fi(NGAUSS),wcth(NGAUSS),wfi(NGAUSS)
      real*8 sumx,sumy,sumz,cfi,sfi,sth,sum
      real*8 intang(NCOM,3),minu,sqin2,dipx,dipy,dipz,ang

      integer index
      integer jatom,ilm,i,j,lmmx,nilmm1
      integer l,m,idm,idp,npg

      logical deb

      COMMON /DEBUG/ deb
!      COMMON /RAD/ RM(NRAD,NATO),JRI(NATO)                       
!      COMMON /SPHE/ POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO), &
!         RNOT(NATO),ndat,IATNR(NDIF),IOP(NDIF),MULT(NATO)
!      COMMON // CLM(NRAD,NCOM,NATO),CLM2(NRAD,NCOM,NATO), &
!         LM(2,NCOM,NATO),LMMAX(NATO)   
      DIMENSION y((lmax2+2)*(lmax2+2))
      PARAMETER (NPG = 60)
!      DATA IMAG/(0.0,1.0)/,SQIN2/0.7071067812/      
      parameter (imag=(0.d0,1.d0) , sqin2=0.707106781186547D0)
      pi=acos(-1.0D0)                     

      jatom=iabs(iatnr(index))
      
      nilmm1=0
      dipx=0.0
      dipy=0.0
      dipz=0.0
      sumx=0.0
      sumy=0.0
      sumz=0.0
      lmmx = lmmax(jatom)

      call gauleg(-1.d0,1.d0,cth,wcth,NPG)
!      call gauleg(0.d0,2.d0*pi,fi,wfi,NPG)
      call gaulep(0.d0,2.d0*pi,fi,wfi,NPG)
      do i=1,nilmm1
        intang(i,1)=0.0
        intang(i,2)=0.0
        intang(i,3)=0.0
      enddo

      do i=1,NPG
        do j=1,NPG
          sth=sqrt(1.0-cth(i)*cth(i))
          cfi=cos(fi(j))
          sfi=sin(fi(j))
          call ylm(cth(i),sth,cfi,sfi,lmax2+1,y)
          do ilm=1,lmmx
            l=iabs(lm(1,ilm,jatom))
            m=lm(2,ilm,jatom)
            minu=1.
            imag1=(1.D0,0.D0)
            if(lm(1,ilm,jatom).lt.0) then
              imag1=-imag
              minu=-1.D0
            endif
            if(mod(m,2).eq.1) then
              imag1=-imag1
              minu=-minu
            endif
            if(m.eq.0) then
              idm=l*(l+1)+1
              ang=y(idm)
            else
              idm=l*(l+1)+1
              idp=idm+m
              idm=idm-m
              ang=(y(idp)+minu*y(idm))*sqin2*imag1
            endif
            intang(ilm,1)=intang(ilm,1)+sth*cfi*ang*wcth(i)*wfi(j)
            intang(ilm,2)=intang(ilm,2)+sth*sfi*ang*wcth(i)*wfi(j)
            intang(ilm,3)=intang(ilm,3)+cth(i)*ang*wcth(i)*wfi(j)
          enddo
        enddo
      enddo

!      write(6,*) '**** INTANG ****'
!      do k=1,lmmx
!        write(6,*) k,lm(1,k,jatom),lm(2,k,jatom),
!     $     intang(k,1),intang(k,2),intang(k,3)
!      enddo

!      write(6,*) '**** DIPOLE ****'
      do ilm=1,lmmx
        sum=0.d0
        j=1
        do while((.not.(rm(j+2,jatom).gt.rmin)).and. &
           ((j+1).lt.jri(jatom)))
          sum=sum+clm(j,ilm,jatom)*rm(j,jatom)*rm(j,jatom)+ &
             4.0*clm(j+1,ilm,jatom)*rm(j+1,jatom)*rm(j+1,jatom)+ &
             clm(j+2,ilm,jatom)*rm(j+2,jatom)*rm(j+2,jatom)
          j=j+2
        end do
        if(j.gt.jri(jatom)) j=jri(jatom)
        rmin=rm(j,jatom)
        sum=sum*dx(jatom)/3.0
        if(lm(1,ilm,jatom).eq.0) then
          sum=sum/sqrt(4.0*pi)
        endif
        sumx=sumx+sum*intang(ilm,1)
        sumy=sumy+sum*intang(ilm,2)
        sumz=sumz+sum*intang(ilm,3)
!        write(6,*) ilm,sum,sumx,sumy,sumz
      enddo

      write(6,*) '**** DIPOLE SPHERE ****'
      write(6,*) sumx,sumy,sumz

      dipx=sumx
      dipy=sumy
      dipz=sumz
        
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
      

      return 
      end
!
!
!
      subroutine dipoleint(f1,theta,phi,r0,rmax,index,npr)
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

      srho=.true.
      sgrho=.false.
      shrho=.false.
      cutoff=0.0

      ct=cos(theta)
      st=sin(theta)
      cf=cos(phi)
      sf=sin(phi)

      if(rmax.lt.r0) then
        f1=0.0
      else
        dr=(rmax-r0)/(1.d0*npr)
        r=r0+0.5D0*dr
        sum=0.0D0
        
        do i=1,npr
          v(1)=pos(1,index)+r*st*cf
          v(2)=pos(2,index)+r*st*sf
          v(3)=pos(3,index)+r*ct
          call vgh_rho(v,rho,grho,hrho,srho,sgrho,shrho,rr,iat,ipos,ist,ibest,ipbest)
          sum=sum+rho*r*r*r
          r=r+dr
        end do
        f1=sum/npr*(rmax-r0)
      end if

      return
      end

