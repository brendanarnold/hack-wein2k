      subroutine vgh_rho(v,rho,grho,hrho,srho,sgrho,shrho,r, &
         iat,ipos,ist,ibest,ipbest)
!
!     AIM (1999-2000) by Javier D. Fuhr and Jorge O. Sofo
!
!     Instituto Balseiro and Centro Atomico Bariloche
!     S. C. de Bariloche - Rio Negro
!     Argentina
!     e-mail: fuhr@cab.cnea.gov.ar
!             sofo@cab.cnea.gov.ar
!
!     Calculate value, gradient, and Hessian of rho depending
!     on the logical switchs srho, sgrho, and shrho.
!
!     The values are return in rho, grho(3), hrho(3,3)
!     
!     iat:  index of atom           (==0 if v is in the interstitial)
!     ipos: index of the unit cell  (==0 if v is in the interstitial)
!     r:    distance to the position of the atom      
!
!

      implicit none
      INCLUDE 'param.inc'
      
      real*8 v(3),rho,grho(3),hrho(3,3),r
      integer iat,ipos,ist,ibest,ipbest
      logical srho,sgrho,shrho,inter
      
      inter=.true.
      ist=0
      call search_at(v,r,inter,iat,ipos,ibest,ipbest)
      if (inter) then
        call interst(v,rho,grho,hrho,srho,sgrho,shrho)
      else
        call sphere(v,r,rho,grho,hrho,iat,srho,sgrho,shrho,ist)
      endif

      return
      end
      
      subroutine search_at(v,r,inter,iat,ipos,ibest,ipbest)
!     Small changes by L. D. Marks, January 2006
!     Added return of closest R, and closest atom (may be buggy)
      use sphe
      use atpos
      IMPLICIT NONE
      include 'param.inc'
      logical inter,deb,deb2
      real*8 br1,br2,br3,br4,r,rjmin,vtmp
      real*8 vt,v(3),vt1(3),drmin,rbest
      integer ipos,iat,ii,jj,jatom,ibest,ipbest
      COMMON /DEBUG/ deb
!      COMMON /BRAV/ br1(3,3),br2(3,3),br3(3,3),br4(3,3)
      parameter (drmin=1d-5)      
      save deb2
      data deb2/.true./
      inter=.true.

        ibest=0
        rbest=1e5

!.... LOOP OVER DIFFERENT ATOMS
      do iat=1,ndat 
        jatom=iabs(iatnr(iat))
        rjmin=(rmt(jatom)*rmt(jatom)-drmin)
!.... LOOP OVER CELL ORIGINS
        do ipos=1,npos
          r=0.d0
          do ii=1,3
            vt=v(ii)-(ATP(II,IPOS)+POS(II,IAT))
            r=r+vt*vt
          enddo
          if(r.lt.rjmin)then
            inter=.false.
            r=sqrt(r)
            return
          else if(r.lt.rbest)then
            rbest=r
            ibest=iat
            ipbest=ipos
          endif
        enddo
      enddo
      iat=0
      ipos=0
      r=sqrt(rbest)
 25   return
      end

