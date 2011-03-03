      SUBROUTINE ONESTEP(v,chg,grho,h,np,npmax,deltar,dir)

!
!     AIM (1999-2000) by Javier D. Fuhr and Jorge O. Sofo
!
!     Instituto Balseiro and Centro Atomico Bariloche
!     S. C. de Bariloche - Rio Negro
!     Argentina
!     e-mail: fuhr@cab.cnea.gov.ar
!             sofo@cab.cnea.gov.ar
!
!     Advance one step following the gradient from V(3). 
!       It returns a
!       new point in V(3) and the value and gradient of the
!       electron density at this point in CHG and GRHO(3)
!
!     H*GRHO is the step
!     NPMAX: maximum number of divisions
!     NP: returns the number of divisions needed
!     DELTAR: return the length of the step 
!
!     Small changes by L. D. Marks, January 2006

      implicit none
      include 'param.inc'

      real*8 v,chg,grho,h,deltar,vinter,vk,vkold,hrho
      real*8 dt,tdt,r,drmin,dir
      real*8 stpos,stgrho,stngrho,stlen
      
      integer np,npmax,i,j,iat,ipos,ist,ibest,ipbest

      logical srho,sgrho,shrho
      logical deb,ldeb

      COMMON /DEBUG/ deb
      COMMON /FOLL/ stpos(3,MAXFOLL),stgrho(3,MAXFOLL), &
           stngrho(MAXFOLL),stlen(MAXFOLL),ldeb

      dimension vinter(3,200)
      DIMENSION V(3),vk(3),vkold(3)
      dimension grho(3),hrho(3,3)
!      data drmin/1.e-4/
      parameter (drmin = 1D-4)

      sgrho=.true.
      srho=.false.
      shrho=.false.

      dt=h*dir
      np=1
      deltar=1.0D0
      do i=1,3
        vk(i)=v(i)
      end do

      do while((np.lt.3).or. &
         ((np.le.npmax).and.(deltar.gt.drmin)))
         np=np*2
         tdt=dt
         dt=dt*0.5d0
          vkold(1)=vk(1)
          vkold(2)=vk(2)
          vkold(3)=vk(3)
        call vgh_rho(vk,chg,grho,hrho,srho,sgrho,shrho,r, &
           iat,ipos,ist,ibest,ipbest)
        if (ist.eq.1) then
          np=npmax+5
          return
        endif

          vinter(1,1)=v(1)+dt*grho(1)
          vinter(2,1)=v(2)+dt*grho(2)
          vinter(3,1)=v(3)+dt*grho(3)

        do j=2,np
          call vgh_rho(vinter(1,j-1),chg,grho,hrho,srho,sgrho,shrho,r, &
             iat,ipos,ist,ibest,ipbest)
          if (ist.eq.1) then
            np=npmax+5
            return
          endif
          if(j.eq.2) then
              vinter(1,2)=v(1)+tdt*grho(1)
              vinter(2,2)=v(2)+tdt*grho(2)
              vinter(3,2)=v(3)+tdt*grho(3)
          else
              vinter(1,j)=vinter(1,j-2)+tdt*grho(1)
              vinter(2,j)=vinter(2,j-2)+tdt*grho(2)
              vinter(3,j)=vinter(3,j-2)+tdt*grho(3)
          endif
        end do

        call vgh_rho(vinter(1,np),chg,grho,hrho,srho,sgrho,shrho,r, &
           iat,ipos,ist,ibest,ipbest)
        if (ist.eq.1) then
          np=npmax+5
          return
        endif

          vinter(1,np+1)=vinter(1,np-1)+dt*grho(1)
          vinter(2,np+1)=vinter(2,np-1)+dt*grho(2)
          vinter(3,np+1)=vinter(3,np-1)+dt*grho(3)
        
        deltar=0.d0
        do i=1,3
          vk(i)=(vinter(i,np)+vinter(i,np+1))*0.5d0
          deltar=deltar+(vkold(i)-vk(i))*(vkold(i)-vk(i))
        end do
      end do
      
      deltar=0.d0
      do i=1,3
        deltar=deltar+(vk(i)-v(i))*(vk(i)-v(i))
      enddo
      deltar=sqrt(deltar)

        v(1)=vk(1)
        v(2)=vk(2)
        v(3)=vk(3)

      srho=.true.
      call vgh_rho(v,chg,grho,hrho,srho,sgrho,shrho,r, &
         iat,ipos,ist,ibest,ipbest)
      if (ist.eq.1) then
        np=npmax+5
        return
      endif
      if(deb.or.ldeb) write(6,*) ':VKf ',np,vk

 25   return 
      end
      
