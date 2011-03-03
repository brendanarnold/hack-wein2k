!***********************************************************************
      subroutine pairdis(rr,np,gg,nsh,pp,rmta,dlambda)
      implicit double precision (a-h,o-z)
      integer np,nsh,i,ip1,ip2,ish1,ish2,ish3,i1,i2
      real*8 rr(3,np),gg(3,3)
      real*8 pp(3,np),rmta(np),rij(3),pij(3)
      real*8 rmtaij2,rij2,rpij, dlam, dlambda,eps,pij2,rad
      real*8 r12(3),p12(3),rmt12, rhelp
      parameter (eps = 1.d-6 )
!
       write(*,*)' in pairdis  ---------- ' 
       write(*,*)'  rr = ' 
      write(*,100)((rr(i,j),i=1,3),j=1,np)
      write(*,*)'  gg = ' 
      write(*,100)gg
      write(*,*)'  pp = ' 
      write(*,100)((pp(i,j),i=1,3),j=1,np)
      write(*,*)'  rmta = '
      write(*,100)(rmta(j),j=1,np)
      write(*,*)'  nsh = ' ,nsh
 110  format(30i3)
 100  format(3f16.8)
      do 10 ip1=1,np
         do 20 ip2 = ip1+1,np
            rmtaij2 = ( rmta(ip1)+rmta(ip2) )**2
            rhelp =  rmta(ip1)+rmta(ip2) 
            do 30 i=1,3
               pij(i) = pp(i,ip1) - pp(i,ip2)
 30         CONTINUE
            if( dabs(pij(1)) .gt. eps .or. &
                dabs(pij(2)) .gt. eps .or. &
                dabs(pij(3)) .gt. eps )then
               pij2 = 1.d0/( pij(1)**2 + pij(2)**2 + pij(3)**2 )
               rmtaij2 = rmtaij2 * pij2
               do 40 ish1=-nsh,nsh
               do 40 ish2=-nsh,nsh
               do 40 ish3=-nsh,nsh
                  rpij = 0.d0
                  rij2 = 0.d0
                  do 50 i =1,3
!
                     rij(i) = rr(i,ip1) - rr(i,ip2) +  &
                     dble(ish1) * gg(i,1) +  &
                     dble(ish2) * gg(i,2) +  &
                     dble(ish3) * gg(i,3) 
!
                     rij2 = rij2 + rij(i)**2
!
                     rpij = rpij + rij(i) * pij(i)
!
 50               CONTINUE
                  if( dabs( rij2 ) .gt. eps .and.  &
                      rpij .lt. -eps  )then
                     rpij = rpij * pij2
                     rij2 = rij2 * pij2
                     rad = rmtaij2 + rpij**2 - rij2
!
                     if( rad .gt. 0.d0 )then
                        dlam = dmin1( - rpij,-rpij - dsqrt(rad) )
                        if( dlam .lt. dlambda )then
                           dlambda = dlam
                           do 60 i=1,3
                              r12(i) = rij(i)
                              p12(i) = pij(i)
 60                        CONTINUE
                           rmt12 = rhelp 
                           i1 = ip1
                           i2 = ip2
                          if ( ip1.eq.1 .and. ip2.eq.2 ) &
             write(*,200)i1,i2,r12,p12,rmt12,dlambda
                        endif
                     endif
                  endif
 40            CONTINUE
            endif
 20      CONTINUE
 10   CONTINUE
      write(*,*)'   output  test '
      write(*,200)i1,i2,r12,p12,rmt12,dlambda
200   format('  teilchen    ',2i5 / &
             '    r12       ',3f10.5/ &
             '    p12       ',3f10.5/ &
             '    rmt12     ',f10.5/ &
             '    dlambda   ',f10.5)
      dlambda =  dlambda * 0.99d0
      return
      end
