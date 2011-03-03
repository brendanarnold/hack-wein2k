      subroutine mixin(xwert,iter,n,nat)
      implicit double precision (a-h,o-z)
      include 'param.inc'
      DIMENSION  xwert(mxpm,maxit),alph(nato+1,4)
      COMMON /STRUK/  POS(3,NDIF),AA,BB,CC,ALPHA(3),RMT(NATO),V(NATO), &
                      PIA(3),VOL,IATNR(NATO),MULT(NATO),ISPLIT(NATO)
      tiny=1.d-6
      iterho=0
      itecon=1
      read(13,*,end=2000) iterho
      if (iterho.lt.1) goto 2000
      do 10 i=1,n,3
         itec=0
         help=dble(i+2)/3.0d0
         help=dint(help)
         iatom=int(help)
         xa=xwert(i,iter)
         xb=xwert(i,iter-1)
         xc=xwert(i,iter-2)
         xd=xwert(i,iter-3)
         xe=xwert(i,iter+1)
         ya=xwert(i+1,iter)
         yb=xwert(i+1,iter-1)
         yc=xwert(i+1,iter-2)
         yd=xwert(i+1,iter-3)
         ye=xwert(i+1,iter+1)
         za=xwert(i+2,iter)
         zb=xwert(i+2,iter-1)
         zc=xwert(i+2,iter-2)
         zd=xwert(i+2,iter-3)
         ze=xwert(i+2,iter+1)
         ichoic=0
         call choice(xa,xb,xe,ya,yb,ye,za,zb,ze,ichoic)
         write(6,*) 'iatom,ichoic',iatom,ichoic
         if (ichoic.eq.0) then
            alph(iatom,1)=1.0d0
            alph(iatom,2)=0.0d0
            alph(iatom,3)=0.0d0
            alph(iatom,4)=0.0d0
         endif
         if (ichoic.eq.1) then
            itec=3
         alph(iatom,1)=1 + (ya - ye)/(-ya + yb) +  &
      (((xa - xe)*(-ya + yb) - (-xa + xb)*(ya - ye))* &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc) -  &
      ((-xa + xc)*(-ya + yb) - (-xa + xb)*(-ya + yc))* &
      (-((ya - ye)*(-za + zb)) + (-ya + yb)*(za - ze)))/ &
      (((-xa + xd)*(-ya + yb) - (-xa + xb)*(-ya + yd))* &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc) -  &
      ((-xa + xc)*(-ya + yb) - (-xa + xb)*(-ya + yc))* &
      (-((-ya + yd)*(-za + zb)) + (-ya + yb)*(-za + zd))) -  &
      (-ya + yd)*(((xa - xe)*(-ya + yb) - (-xa + xb)*(ya - ye))* &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc) -  &
      ((-xa + xc)*(-ya + yb) - (-xa + xb)*(-ya + yc))* &
      (-((ya - ye)*(-za + zb)) + (-ya + yb)*(za - ze)))/ &
      ((-ya + yb)*(((-xa + xd)*(-ya + yb) - (-xa + xb)*(-ya + yd))* &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc) -  &
      ((-xa + xc)*(-ya + yb) - (-xa + xb)*(-ya + yc))* &
      (-((-ya + yd)*(-za + zb)) + (-ya + yb)*(-za + zd)))) -  &
      (-((-ya + yd)*(-za + zb)) + (-ya + yb)*(-za + zd))* &
      (((xa - xe)*(-ya + yb) - (-xa + xb)*(ya - ye))* &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc) -  &
      ((-xa + xc)*(-ya + yb) - (-xa + xb)*(-ya + yc))* &
      (-((ya - ye)*(-za + zb)) + (-ya + yb)*(za - ze)))/ &
      ((-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc)* &
      (((-xa + xd)*(-ya + yb) - (-xa + xb)*(-ya + yd))* &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc) -  &
      ((-xa + xc)*(-ya + yb) - (-xa + xb)*(-ya + yc))* &
      (-((-ya + yd)*(-za + zb)) + (-ya + yb)*(-za + zd)))) -  &
      (-ya + yc)*(-((-((-ya + yd)*(-za + zb)) + (-ya + yb)*(-za + zd))* &
      (((xa - xe)*(-ya + yb) - (-xa + xb)*(ya - ye))* &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc) -  &
      ((-xa + xc)*(-ya + yb) - (-xa + xb)*(-ya + yc))* &
      (-((ya - ye)*(-za + zb)) + (-ya + yb)*(za - ze)))/ &
      ((-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc)* &
      (((-xa + xd)*(-ya + yb) - (-xa + xb)*(-ya + yd))* &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc) -  &
      ((-xa + xc)*(-ya + yb) - (-xa + xb)*(-ya + yc))* &
      (-((-ya + yd)*(-za + zb)) + (-ya + yb)*(-za + zd))))) +  &
      (-((ya - ye)*(-za + zb)) + (-ya + yb)*(za - ze))/ &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc))/(-ya + yb) +  &
      (-((ya - ye)*(-za + zb)) + (-ya + yb)*(za - ze))/ &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc)
          alph(iatom,2)=-((ya - ye)/(-ya + yb)) +  &
      (-ya + yd)*(((xa - xe)*(-ya + yb) - (-xa + xb)*(ya - ye))* &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc) -  &
      ((-xa + xc)*(-ya + yb) - (-xa + xb)*(-ya + yc))* &
      (-((ya - ye)*(-za + zb)) + (-ya + yb)*(za - ze)))/ &
      ((-ya + yb)*(((-xa + xd)*(-ya + yb) - (-xa + xb)*(-ya + yd))* &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc) -  &
      ((-xa + xc)*(-ya + yb) - (-xa + xb)*(-ya + yc))* &
      (-((-ya + yd)*(-za + zb)) + (-ya + yb)*(-za + zd)))) +  &
      (-ya + yc)*(-((-((-ya + yd)*(-za + zb)) + (-ya + yb)*(-za + zd))* &
      (((xa - xe)*(-ya + yb) - (-xa + xb)*(ya - ye))* &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc) -  &
      ((-xa + xc)*(-ya + yb) - (-xa + xb)*(-ya + yc))* &
      (-((ya - ye)*(-za + zb)) + (-ya + yb)*(za - ze)))/ &
      ((-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc)* &
      (((-xa + xd)*(-ya + yb) - (-xa + xb)*(-ya + yd))* &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc) -  &
      ((-xa + xc)*(-ya + yb) - (-xa + xb)*(-ya + yc))* &
      (-((-ya + yd)*(-za + zb)) + (-ya + yb)*(-za + zd))))) +  &
      (-((ya - ye)*(-za + zb)) + (-ya + yb)*(za - ze))/ &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc))/(-ya + yb)
          alph(iatom,3)=(-((-ya + yd)*(-za + zb)) +  &
      (-ya + yb)*(-za + zd))* &
      (((xa - xe)*(-ya + yb) - (-xa + xb)*(ya - ye))* &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc) -  &
      ((-xa + xc)*(-ya + yb) - (-xa + xb)*(-ya + yc))* &
      (-((ya - ye)*(-za + zb)) + (-ya + yb)*(za - ze)))/ &
      ((-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc)* &
      (((-xa + xd)*(-ya + yb) - (-xa + xb)*(-ya + yd))* &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc) -  &
      ((-xa + xc)*(-ya + yb) - (-xa + xb)*(-ya + yc))* &
      (-((-ya + yd)*(-za + zb)) + (-ya + yb)*(-za + zd)))) -  &
      (-((ya - ye)*(-za + zb)) + (-ya + yb)*(za - ze))/ &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc)
          alph(iatom,4)=-((((xa - xe)*(-ya + yb)  &
      - (-xa + xb)*(ya - ye))* &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc) -  &
      ((-xa + xc)*(-ya + yb) - (-xa + xb)*(-ya + yc))* &
      (-((ya - ye)*(-za + zb)) + (-ya + yb)*(za - ze)))/ &
      (((-xa + xd)*(-ya + yb) - (-xa + xb)*(-ya + yd))* &
      (-(yb*za) + yc*za + ya*zb - yc*zb - ya*zc + yb*zc) -  &
      ((-xa + xc)*(-ya + yb) - (-xa + xb)*(-ya + yc))* &
      (-((-ya + yd)*(-za + zb)) + (-ya + yb)*(-za + zd))))
         endif
         if (ichoic.eq.2) then
            itec=1
            alph(iatom,1)=1+(xa-xe)/(-xa+xb)
            alph(iatom,2)=1.0d0-alph(iatom,1)
            alph(iatom,3)=0.0d0
            alph(iatom,4)=0.0d0
         endif
         if (ichoic.eq.3) then
            itec=1
            alph(iatom,1)=1+(ya-ye)/(-ya+yb)
            alph(iatom,2)=1.0d0-alph(iatom,1)
            alph(iatom,3)=0.0d0
            alph(iatom,4)=0.0d0
         endif
         if (ichoic.eq.4) then
            itec=1
            alph(iatom,1)=1+(za-ze)/(-za+zb)
            alph(iatom,2)=1.0d0-alph(iatom,1)
            alph(iatom,3)=0.0d0
            alph(iatom,4)=0.0d0
         endif
         if (ichoic.eq.5) then
            itec=2
            alph(iatom,1)=1 + ((xa - xe)* &
      (-ya + yb) - (-xa + xb)*(ya - ye))/ &
      (xb*ya - xc*ya - xa*yb + xc*yb + xa*yc - xb*yc) - &
      (-ya + yc)*((xa - xe)*(-ya + yb) - (-xa + xb)*(ya - ye))/ &
      ((-ya + yb)*(xb*ya - xc*ya - xa*yb + xc*yb + xa*yc - xb*yc)) + &
      (ya - ye)/(-ya + yb)
            alph(iatom,2)=(-ya + yc)*((xa - xe)*(-ya + yb) -  &
      (-xa + xb)*(ya - ye))/ &
      ((-ya + yb)*(xb*ya - xc*ya - xa*yb + xc*yb + xa*yc - xb*yc)) - &
      (ya - ye)/(-ya + yb)
            alph(iatom,3)=-(((xa - xe)*(-ya + yb) -  &
      (-xa + xb)*(ya - ye))/ &
      (xb*ya - xc*ya - xa*yb + xc*yb + xa*yc - xb*yc))
            alph(iatom,4)=0.0d0
         endif
         if (ichoic.eq.6) then
            itec=2
            alph(iatom,1)=1 + ((xa - xe)* &
      (-za + zb) - (-xa + xb)*(za - ze))/ &
      (xb*za - xc*za - xa*zb + xc*zb + xa*zc - xb*zc) - &
      (-za + zc)*((xa - xe)*(-za + zb) - (-xa + xb)*(za - ze))/ &
      ((-za + zb)*(xb*za - xc*za - xa*zb + xc*zb + xa*zc - xb*zc)) + &
      (za - ze)/(-za + zb)
            alph(iatom,2)=(-za + zc)*((xa - xe)*(-za + zb) -  &
      (-xa + xb)*(za - ze))/ &
      ((-za + zb)*(xb*za - xc*za - xa*zb + xc*zb + xa*zc - xb*zc)) - &
      (za - ze)/(-za + zb)
            alph(iatom,3)=-(((xa - xe)*(-za + zb) -  &
      (-xa + xb)*(za - ze))/ &
      (xb*za - xc*za - xa*zb + xc*zb + xa*zc - xb*zc))
            alph(iatom,4)=0.0d0
         endif
         if (ichoic.eq.7) then
            itec=2
            alph(iatom,1)=1 + ((ya - ye)* &
      (-za + zb) - (-ya + yb)*(za - ze))/ &
      (yb*za - yc*za - ya*zb + yc*zb + ya*zc - yb*zc) - &
      (-za + zc)*((ya - ye)*(-za + zb) - (-ya + yb)*(za - ze))/ &
      ((-za + zb)*(yb*za - yc*za - ya*zb + yc*zb + ya*zc - yb*zc)) + &
      (za - ze)/(-za + zb)
            alph(iatom,2)=(-za + zc)*((ya - ye)*(-za + zb) -  &
      (-ya + yb)*(za - ze))/ &
      ((-za + zb)*(yb*za - yc*za - ya*zb + yc*zb + ya*zc - yb*zc)) - &
      (za - ze)/(-za + zb)
            alph(iatom,3)=-(((ya - ye)*(-za + zb) -  &
      (-ya + yb)*(za - ze))/ &
      (yb*za - yc*za - ya*zb + yc*zb + ya*zc - yb*zc))
            alph(iatom,4)=0.0d0
         endif
!         write(6,*) 'xpl',xwert(i,iter+1),xwert(i,iter),xwert(i,iter-1)
!         frag=xwert(i,iter)*xwert(i+2,iter-1)-
!     *        xwert(i,iter-2)*xwert(i+2,iter-1)+
!     *        xwert(i,iter-2)*xwert(i+2,iter)-
!     *        xwert(i,iter)*xwert(i+2,iter-2)+
!     *        xwert(i,iter-1)*xwert(i+2,iter-2)-
!     *        xwert(i,iter-1)*xwert(i+2,iter)
!         al=(-xwert(i,iter+1)*xwert(i+2,iter-2)+
!     *       xwert(i,iter-2)*xwert(i+2,iter+1)+
!     *       xwert(i,iter+1)*xwert(i+2,iter-1)-
!     *       xwert(i,iter-1)*xwert(i+2,iter+1)-
!     *       xwert(i,iter-2)*xwert(i+2,iter-1)+
!     *       xwert(i,iter-1)*xwert(i+2,iter-2))/
!     *       frag
!         bl=-(xwert(i,iter+1)*xwert(i+2,iter)-
!     *       xwert(i,iter-2)*xwert(i+2,iter)+
!     *       xwert(i,iter)*xwert(i+2,iter-2)-
!     *       xwert(i,iter+1)*xwert(i+2,iter-2)+
!     *       xwert(i,iter-2)*xwert(i+2,iter+1)-
!     *       xwert(i,iter)*xwert(i+2,iter+1))/
!     *       frag
!         cl=(xwert(i,iter)*xwert(i+2,iter-1)-
!     *       xwert(i,iter+1)*xwert(i+2,iter-1)+
!     *       xwert(i,iter+1)*xwert(i+2,iter)-
!     *       xwert(i,iter)*xwert(i+2,iter+1)-
!     *       xwert(i,iter-1)*xwert(i+2,iter)+
!     *       xwert(i,iter-1)*xwert(i+2,iter+1))/
!     *       frag
!         write(6,*) 'al,bl,cl',al,bl,cl
!         xtest1=xwert(i,iter)*al+
!     *          bl*xwert(i,iter-1)+cl*xwert(i,iter-2)
!     *          al*(xwert(i,iter)-xwert(i,iter-1))-
!     *          bl*(xwert(i,iter-1)-xwert(i,iter-2))
!         xtest3=al*xwert(i+2,iter)+
!     *          bl*xwert(i+2,iter-1)+cl*xwert(i+2,iter-2)
!         write(6,*) 'xtest',xtest1,xwert(i,iter+1)
!         write(6,*) 'ztest',xtest3,xwert(i+2,iter+1)
         write(6,*) 'alph',alph(iatom,1),alph(iatom,2),alph(iatom,3), &
                     alph(iatom,4)
         xtest1=xwert(i,iter)*alph(iatom,1)+ &
          alph(iatom,2)*xwert(i,iter-1)+alph(iatom,3)*xwert(i,iter-2) &
               +alph(iatom,4)*xwert(i,iter-3)
         xtest3=xwert(i+2,iter)*alph(iatom,1)+ &
                alph(iatom,2)*xwert(i+2,iter-1)+ &
                alph(iatom,3)*xwert(i+2,iter-2) &
               +alph(iatom,4)*xwert(i+2,iter-3)
         write(6,*) 'xtest',xtest1,xwert(i,iter+1)
         write(6,*) 'ztest',xtest3,xwert(i+2,iter+1)
!         alph(iatom,1)=al
!         alph(iatom,2)=bl
!         alph(iatom,3)=cl
!         alph(iatom,4)=0.0d0
          if (itecon.lt.itec) itecon=itec
  10  CONTINUE
      do 30 j=1,4
         alph(iatom+1,j)=0.d0
         do 20 i=1,iatom
            alph(iatom+1,j)=alph(iatom+1,j)+alph(i,j)
  20     CONTINUE
         alph(iatom+1,j)=alph(iatom+1,j)/iatom
  30  CONTINUE
      write(6,1000)'INTER'
      do 40 i=1,iatom+1
         write(6,1010)(alph(i,j),j=1,4)
  40  CONTINUE
! 2000 write(6,*)'mixin',alph,iterho,iatom
 2000 call inter(alph,iterho,nat,itecon)
      return
 1000 format(a5)
 1010 format(4f10.4)
      end
