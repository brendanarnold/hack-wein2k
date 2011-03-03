      subroutine grhosphe(grho,ir,r,jatom,xy,cth,sth,cfi,sfi,iatnr)

!
!     AIM (1999-2000) by Javier D. Fuhr and Jorge O. Sofo
!
!     Instituto Balseiro and Centro Atomico Bariloche
!     S. C. de Bariloche - Rio Negro
!     Argentina
!     e-mail: fuhr@cab.cnea.gov.ar
!             sofo@cab.cnea.gov.ar
!
!     This routine calculates the gradient of the electron density
!       for a point inside a MT.
!
!     GRHO: returns the gradient
!     IR:   index of radial point
!     R:    radius
!     JATOM: index of atom
!     XY:    sqrt(x^2+y^2) (in local coordinates)
!     CTH,STH,CFI,SFI: cos(theta),sin(theta),cos(phi),sin(phi)
!     IATNR
!

      use blank
      use rad
      implicit none
      include 'param.inc'

      complex*16 yl,dtyl,dtdtyl,dfdtyl
      complex*16 imag,imag1

      real*8 rho,ang,grho,dscrho,drrho,dtang,dfang,change
      real*8 minu,sqin2,xy,r,yl1,tstpi,snull
      real*8 rhopp,rhor2,drhor2
      real*8 cth,sth,cfi,sfi,pi,tpi,fpi,sqpi,sqtpi

      integer lmmx,ilm,ir,jatom,l,m,idx,idm,idp,i,j
      integer iatnr,jrj

      logical deb
      logical sy,syp,sypp
      
!     lmax2+2 needed for the derivative of ylm
      COMMON /DEBUG/ deb
!      COMMON // CLM(NRAD,NCOM,NATO),CLM2(NRAD,NCOM,NATO), &
!         LM(2,NCOM,NATO),LMMX(NATO)      
!      COMMON /RAD/ RM(NRAD,NATO),JRI(NATO)
      COMMON /YLMS/ yl((lmax2+4)*(lmax2+4)),dtyl((lmax2+3)*(lmax2+3)), &
         dtdtyl((lmax2+2)*(lmax2+2)),dfdtyl((lmax2+2)*(lmax2+2))
      COMMON /CTES/ pi,tpi,fpi,sqpi,sqtpi

      dimension rho(ncom),ang(ncom),grho(3),dscrho(3)
      dimension drrho(ncom),dtang(ncom),dfang(ncom),change(3,3)
!      data imag/(0.0,1.0)/,sqin2/0.7071067812/,snull/1.d-8/
      parameter (imag=(0.d0,1.d0), sqin2=0.707106781186547D0)
      parameter (snull=1.d-8)     
      jrj=jri(jatom)
      sy=.true.
      syp=.true.
      sypp=.false.

      if(r.lt.snull) then
        do i=1,3
          grho(i)=0.d0
        enddo
        goto 25
      else
!$$$        call ylm(v,lmax2+1,yl)
!$$$        call dtylm(v,lmax2,yl,dtyl)
        call gen_change(change,cth,sth,cfi,sfi)
        lmmx=lmmax(jatom)        
        
        if (xy.eq.0.d0) then
          tstpi=2.0*sqtpi
          do i=1,3
            grho(i)=0.d0
          enddo
          do ilm=1,lmmx
            l=iabs(lm(1,ilm,jatom))
            m=lm(2,ilm,jatom)
            call splint(rm(1,jatom),clm(1,ilm,jatom), &
               clm2(1,ilm,jatom), &
               jrj,r,ir,rhor2,drhor2,rhopp,sy,syp,sypp)
            rho(ilm)=rhor2/(r*r)
            drrho(ilm)=drhor2/(r*r)-2.D0*rho(ilm)/r
            if(m.eq.0) then
              idx=l*(l+1)+1
              if(cth.lt.0.) then
                grho(3)=grho(3)-drrho(ilm)*yl(idx)
              else
                grho(3)=grho(3)+drrho(ilm)*yl(idx)
              endif
            else
              if(m.eq.1) then
                yl1=sqrt(l*(l+1.D0)*(2.D0*l+1.D0))/tstpi
                if(lm(1,ilm,jatom).lt.0) then
                  grho(2)=grho(2)+yl1*rho(ilm)/r
                else
                  grho(1)=grho(1)+yl1*rho(ilm)/r
                endif
              endif
            endif
          enddo
        else
          do ilm=1,lmmx                                                  
            call splint(rm(1,jatom),clm(1,ilm,jatom), &
               clm2(1,ilm,jatom), &
               jrj,r,ir,rhor2,drhor2,rhopp,sy,syp,sypp)
            rho(ilm)=rhor2/(r*r)
            drrho(ilm)=drhor2/(r*r)-2.D0*rho(ilm)/r
            l=iabs(lm(1,ilm,jatom))
            m=lm(2,ilm,jatom)
            minu=1.
            imag1=(1.,0.)
            if(lm(1,ilm,jatom).lt.0) then
              imag1=-imag
              minu=-1.
            end if
            if(mod(m,2).eq.1) then
              imag1=-imag1
              minu=-minu
            end if
            if(m.eq.0) then                                                
              idx=l*(l+1)+1
              ang(ilm)=yl(idx)
              dtang(ilm)=dtyl(idx)
              dfang(ilm)=0.D0
            else                                                           
              idm=l*(l+1)+1
              idp=idm+m
              idm=idm-m
              ang(ilm)=(yl(idp)+minu*yl(idm))*sqin2*imag1
              dtang(ilm)=(dtyl(idp)+minu*dtyl(idm))*sqin2*imag1
              dfang(ilm)=(yl(idp)-minu*yl(idm))*sqin2*imag1*imag*m 
            endif
          enddo
          
!$$$          if(deb) write(6,*) 'ang: ',(ang(ilm),ilm=1,lmmx)
!$$$          if(deb) write(6,*) 'dtang: ',(dtang(ilm),ilm=1,lmmx)
!$$$          if(deb) write(6,*) 'dfang: ',(dfang(ilm),ilm=1,lmmx)

!          write(6,*) 'ang: ',(ang(ilm),ilm=1,lmmx)
!          write(6,*) 'dtang: ',(dtang(ilm),ilm=1,lmmx)
!          write(6,*) 'dfang: ',(dfang(ilm),ilm=1,lmmx)

          if(iatnr.gt.0) then
            call sumd(rho,drrho,ang,dtang,dfang,dscrho,lmmx,lm,jatom)
            if(deb) then
              do ilm=1,lmmx 
                write(6,*)'grhosphe: ',ilm,rho(ilm),drrho(ilm),ang(ilm)
              enddo
            endif
          else
            dscrho(1)=0.D0
            dscrho(2)=0.D0
            dscrho(3)=0.D0
            do ilm=1,lmmx 
              if(deb) then
                write(6,*)'grhosphe: ',ilm,ang(ilm),dtang(ilm), &
                 dfang(ilm),rho(ilm),drrho(ilm)
              endif
              dscrho(1) = dscrho(1)+drrho(ilm)*ang(ilm)
              dscrho(2) = dscrho(2)+rho(ilm)*dtang(ilm)
              dscrho(3) = dscrho(3)+rho(ilm)*dfang(ilm)
            enddo
          endif
!          write(6,*) 'grhosphere ',dscrho(1),dscrho(2),dscrho(3)
          dscrho(2) = dscrho(2)/r
          dscrho(3) = dscrho(3)/(r*sth)

!     Now we go from spherical coordinates(dscrho) 
!     to Cartesian ones(grho)
          do i=1,3
            grho(i)=0.D0
            do j=1,3
              if(deb) then
                write(6,*)'grhosphe: ',i,j,dscrho(j),change(i,j)
              endif
              grho(i) = grho(i)+dscrho(j)*change(j,i)
            enddo
          enddo
        endif
      endif

 25   return                                                            
      end                                                               

