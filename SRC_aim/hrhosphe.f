      SUBROUTINE HRHOSPHE(hrho,ir,r,jatom,xy,cth,sth,cfi,sfi,iatnr,ist)

!
!     AIM (1999-2000) by Javier D. Fuhr and Jorge O. Sofo
!
!     Instituto Balseiro and Centro Atomico Bariloche
!     S. C. de Bariloche - Rio Negro
!     Argentina
!     e-mail: fuhr@cab.cnea.gov.ar
!             sofo@cab.cnea.gov.ar
!
!     This routine calculates the Hessian of the electron density for
!       point inside a MT.
!
!     NOTE: It is not yet implement the limit for sqrt(x^2+y^2) (in
!           local coordinates) equal 0.
!
      use blank
      use rad
      implicit none
      include 'param.inc'

      complex*16 yl,dtyl,dtdtyl,dfdtyl
      complex*16 imag,imag1

      real*8 rho,ang,hrho,drrho,drdrrho
      real*8 dtang,dfang,dtdtang,dfdtang,dfdfang,change
      real*8 minu,sqin2,xy,r,snull
      real*8 rhor2,drhor2,ddrhor2
      real*8 cth,sth,cfi,sfi,pi,tpi,fpi,sqpi,sqtpi,invr,invr2
      real*8 mat

      integer lmmx,ilm,ir,jatom,l,m,idx,idm,idp,i,j
      integer iatnr,jrj,ist

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

      dimension rho(ncom),ang(ncom),hrho(3,3),mat(3,3)
      dimension drrho(ncom),dtang(ncom),dfang(ncom),change(3,3)
      dimension drdrrho(ncom),dtdtang(ncom),dfdtang(ncom),dfdfang(ncom)
!      data imag/(0.0,1.0)/,sqin2/0.7071067812/,snull/1.e-8/
      parameter (imag=(0.d0,1.d0), sqin2=0.707106781186547D0 )
      parameter (snull=1.d-8)     
      jrj=jri(jatom)
      ist=0

      if(r.lt.snull) then
        do i=1,3
          do j=1,3
            hrho(i,j)=0.D0
          enddo
        enddo
        ist=2
        goto 25
      else
        call gen_change(change,cth,sth,cfi,sfi)
!        call gen_invchange(change,invchange)
        lmmx=lmmax(jatom)        
        
        if (xy.eq.0.d0) then
          do i=1,3
            do j=1,3
              hrho(i,j)=0.D0
            enddo
          enddo
          ist=3
          goto 25
        else
          sy=.true.
          syp=.true.
          sypp=.true.

          do ilm=1,lmmx                                                  
            call splint(rm(1,jatom),clm(1,ilm,jatom), &
               clm2(1,ilm,jatom), &
               jrj,r,ir,rhor2,drhor2,ddrhor2,sy,syp,sypp)
            rho(ilm)=rhor2/(r*r)
            drrho(ilm)=drhor2/(r*r)-2.*rho(ilm)/r
            drdrrho(ilm)=ddrhor2/(r*r)-4.D0*drhor2/(r*r*r)+ &
               6.D0*rho(ilm)/(r*r)
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
              dtdtang(ilm)=dtdtyl(idx)
              dfdtang(ilm)=0.D0
              dfang(ilm)=0.D0
              dfdfang(ilm)=0.D0
            else                                                           
              idm=l*(l+1)+1
              idp=idm+m
              idm=idm-m
              ang(ilm)=(yl(idp)+minu*yl(idm))*sqin2*imag1
              dtang(ilm)=(dtyl(idp)+minu*dtyl(idm))*sqin2*imag1
              dtdtang(ilm)=(dtdtyl(idp)+minu*dtdtyl(idm))*sqin2*imag1
!              dfdtang(ilm)=dtang(ilm)*imag*m
              dfdtang(ilm)=(dfdtyl(idp)+minu*dfdtyl(idm))*sqin2*imag1
              dfang(ilm)=(yl(idp)-minu*yl(idm))*sqin2*imag1*imag*m 
              dfdfang(ilm)=-(yl(idp)+minu*yl(idm))*sqin2*imag1*m*m 
            endif
          enddo
          
          if (deb) write(6,*) 'ang: ',(ang(ilm),ilm=1,lmmx)
          if (deb) write(6,*) 'dtang: ',(dtang(ilm),ilm=1,lmmx)
          if (deb) write(6,*) 'dtdtang: ',(dtdtang(ilm),ilm=1,lmmx)
          if (deb) write(6,*) 'dfdtang: ',(dfdtang(ilm),ilm=1,lmmx)
          if (deb) write(6,*) 'dfang: ',(dfang(ilm),ilm=1,lmmx)
          if (deb) write(6,*) 'dfdfang: ',(dfdfang(ilm),ilm=1,lmmx)
          if(iatnr.gt.0) then
            call sumdd(r,cth,sth,rho,drrho,drdrrho,ang, &
               dtang,dtdtang,dfdtang,dfang,dfdfang,hrho,lmmx,lm,jatom)
          else
            do i=1,3
              do j=1,3
                hrho(i,j)=0.
              enddo
            enddo
            invr=1./r
            invr2=invr*invr
            do ilm=1,lmmx 
              hrho(1,1) = hrho(1,1)+drdrrho(ilm)*ang(ilm)
              hrho(1,2) = hrho(1,2)-rho(ilm)*dtang(ilm)*invr2+ &
                 drrho(ilm)*dtang(ilm)*invr
              hrho(1,3) = hrho(1,3)+(-rho(ilm)*dfang(ilm)*invr2+ &
                 drrho(ilm)*dfang(ilm)*invr)/sth
              hrho(2,2) = hrho(2,2)+rho(ilm)*dtdtang(ilm)*invr2+ &
                 drrho(ilm)*ang(ilm)*invr
              hrho(2,3) = hrho(2,3)+rho(ilm)* &
                 (-cth*dfang(ilm)+sth*dfdtang(ilm))/(sth*sth)*invr2
              hrho(3,3) = hrho(3,3)+rho(ilm)*(dfdfang(ilm)+ &
                 cth*sth*dtang(ilm))/(sth*sth)*invr2+ &
                 drrho(ilm)*ang(ilm)*invr
            enddo
            hrho(2,1)=hrho(1,2)
            hrho(3,1)=hrho(1,3)
            hrho(3,2)=hrho(2,3)
          endif
          if (deb) then
            write (6,*) 'HRHOSPHE HRHO = ',(hrho(1,i),i=1,3)
            write (6,*) '                 ',(hrho(2,i),i=1,3)
            write (6,*) '                 ',(hrho(3,i),i=1,3)
          endif

          call mat3prod(mat,change,hrho,.true.,.false.)
          call mat3prod(hrho,mat,change,.false.,.false.)
!          call mat3prod(mat,change,hrho,.false.,.false.)
!          call mat3prod(hrho,mat,change,.false.,.true.)

        endif
      endif
 25   return                                                            
      end
                                         
