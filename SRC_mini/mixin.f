      subroutine mixin(q,xwert,iter,n,nat,jspin)
      implicit double precision (a-h,o-z)
      include 'param.inc'
      DIMENSION  q(nato,2,maxit),xwert(mxpm,maxit),alph(nato+1,4,2)
      COMMON /STRUK/  POS(3,NDIF),AA,BB,CC,ALPHA(3),RMT(NATO),V(NATO), &
                      PIA(3),VOL,IATNR(NATO),MULT(NATO),ISPLIT(NATO)
      tiny=1.d-6
      iterho=0
      itecon=1
      read(13,*,end=2000) iterho
      if (iterho.lt.1) goto 2000
      if (iter.gt.30) then 
	 call lsqfit(xwert,q,iter,jspin,nat)
         do 100 ispin=1,jspin
         do 100 iatom=1,nat
            itec=1
            alph(iatom,1,ispin)=1+ &
                      (q(iatom,ispin,iter)-q(iatom,ispin,iter+1))/ &
                      (-q(iatom,ispin,iter)+q(iatom,ispin,iter-1))
            alph(iatom,2,ispin)=1.0d0-alph(iatom,1,ispin)
            alph(iatom,3,ispin)=0.0d0
            alph(iatom,4,ispin)=0.0d0
 100     continue
      else 
      jspin=2
      do 10 ispin=1,jspin
      do 10 i=1,n,3
         itec=0
         help=dble(i+2)/3.0d0
         iatom=nint(help)
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
         if(ichoic.eq.1) ichoic=6
         if (ichoic.eq.0) then
            alph(iatom,1,ispin)=1.0d0
            alph(iatom,2,ispin)=0.0d0
            alph(iatom,3,ispin)=0.0d0
            alph(iatom,4,ispin)=0.0d0
         endif
         if (ichoic.eq.1) then
            itec=3
         alph(iatom,1,ispin)=1 + (ya - ye)/(-ya + yb) +  &
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
          alph(iatom,2,ispin)=-((ya - ye)/(-ya + yb)) +  &
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
          alph(iatom,3,ispin)=(-((-ya + yd)*(-za + zb)) +  &
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
          alph(iatom,4,ispin)=-((((xa - xe)*(-ya + yb)  &
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
            alph(iatom,1,ispin)=1+(xa-xe)/(-xa+xb)
            alph(iatom,2,ispin)=1.0d0-alph(iatom,1,ispin)
            alph(iatom,3,ispin)=0.0d0
            alph(iatom,4,ispin)=0.0d0
         endif
         if (ichoic.eq.3) then
            itec=1
            alph(iatom,1,ispin)=1+(ya-ye)/(-ya+yb)
            alph(iatom,2,ispin)=1.0d0-alph(iatom,1,ispin)
            alph(iatom,3,ispin)=0.0d0
            alph(iatom,4,ispin)=0.0d0
         endif
         if (ichoic.eq.4) then
            itec=1
            alph(iatom,1,ispin)=1+(za-ze)/(-za+zb)
            alph(iatom,2,ispin)=1.0d0-alph(iatom,1,ispin)
            alph(iatom,3,ispin)=0.0d0
            alph(iatom,4,ispin)=0.0d0
         endif
         if (ichoic.eq.5) then
            itec=2
            alph(iatom,1,ispin)=1 + ((xa - xe)* &
      (-ya + yb) - (-xa + xb)*(ya - ye))/ &
      (xb*ya - xc*ya - xa*yb + xc*yb + xa*yc - xb*yc) - &
      (-ya + yc)*((xa - xe)*(-ya + yb) - (-xa + xb)*(ya - ye))/ &
      ((-ya + yb)*(xb*ya - xc*ya - xa*yb + xc*yb + xa*yc - xb*yc)) + &
      (ya - ye)/(-ya + yb)
            alph(iatom,2,ispin)=(-ya + yc)*((xa - xe)*(-ya + yb) -  &
      (-xa + xb)*(ya - ye))/ &
      ((-ya + yb)*(xb*ya - xc*ya - xa*yb + xc*yb + xa*yc - xb*yc)) - &
      (ya - ye)/(-ya + yb)
            alph(iatom,3,ispin)=-(((xa - xe)*(-ya + yb) -  &
      (-xa + xb)*(ya - ye))/ &
      (xb*ya - xc*ya - xa*yb + xc*yb + xa*yc - xb*yc))
            alph(iatom,4,ispin)=0.0d0
         endif
         if (ichoic.eq.6) then
            itec=2
            alph(iatom,1,ispin)=1 + ((xa - xe)* &
      (-za + zb) - (-xa + xb)*(za - ze))/ &
      (xb*za - xc*za - xa*zb + xc*zb + xa*zc - xb*zc) - &
      (-za + zc)*((xa - xe)*(-za + zb) - (-xa + xb)*(za - ze))/ &
      ((-za + zb)*(xb*za - xc*za - xa*zb + xc*zb + xa*zc - xb*zc)) + &
      (za - ze)/(-za + zb)
            alph(iatom,2,ispin)=(-za + zc)*((xa - xe)*(-za + zb) -  &
      (-xa + xb)*(za - ze))/ &
      ((-za + zb)*(xb*za - xc*za - xa*zb + xc*zb + xa*zc - xb*zc)) - &
      (za - ze)/(-za + zb)
            alph(iatom,3,ispin)=-(((xa - xe)*(-za + zb) -  &
      (-xa + xb)*(za - ze))/ &
      (xb*za - xc*za - xa*zb + xc*zb + xa*zc - xb*zc))
            alph(iatom,4,ispin)=0.0d0
         endif
         if (ichoic.eq.7) then
            itec=2
            alph(iatom,1,ispin)=1 + ((ya - ye)* &
      (-za + zb) - (-ya + yb)*(za - ze))/ &
      (yb*za - yc*za - ya*zb + yc*zb + ya*zc - yb*zc) - &
      (-za + zc)*((ya - ye)*(-za + zb) - (-ya + yb)*(za - ze))/ &
      ((-za + zb)*(yb*za - yc*za - ya*zb + yc*zb + ya*zc - yb*zc)) + &
      (za - ze)/(-za + zb)
            alph(iatom,2,ispin)=(-za + zc)*((ya - ye)*(-za + zb) -  &
      (-ya + yb)*(za - ze))/ &
      ((-za + zb)*(yb*za - yc*za - ya*zb + yc*zb + ya*zc - yb*zc)) - &
      (za - ze)/(-za + zb)
            alph(iatom,3,ispin)=-(((ya - ye)*(-za + zb) -  &
      (-ya + yb)*(za - ze))/ &
      (yb*za - yc*za - ya*zb + yc*zb + ya*zc - yb*zc))
            alph(iatom,4,ispin)=0.0d0
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
         write(6,*) 'alph',alph(iatom,1,ispin),alph(iatom,2,ispin), &
         alph(iatom,3,ispin),alph(iatom,4,ispin)
         xtest1=xwert(i,iter)*alph(iatom,1,ispin)+ &
                alph(iatom,2,ispin)*xwert(i,iter-1)+ &
                alph(iatom,3,ispin)*xwert(i,iter-2) &
               +alph(iatom,4,ispin)*xwert(i,iter-3)
         xtest3=xwert(i+2,iter)*alph(iatom,1,ispin)+ &
                alph(iatom,2,ispin)*xwert(i+2,iter-1)+ &
                alph(iatom,3,ispin)*xwert(i+2,iter-2) &
               +alph(iatom,4,ispin)*xwert(i+2,iter-3)
         write(6,*) 'xtest',xtest1,xwert(i,iter+1)
         write(6,*) 'ztest',xtest3,xwert(i+2,iter+1)
!         alph(iatom,1)=al
!         alph(iatom,2)=bl
!         alph(iatom,3)=cl
!         alph(iatom,4)=0.0d0
          if (itecon.lt.itec) itecon=itec
  10  CONTINUE
      endif
      do 30 ispin=1,jspin
      do 30 j=1,4
         alph(nat+1,j,ispin)=0.d0
         do 20 i=1,nat
            alph(nat+1,j,ispin)=alph(nat+1,j,ispin)+alph(i,j,ispin)
  20     CONTINUE
         alph(nat+1,j,ispin)=alph(nat+1,j,ispin)/iatom
  30  CONTINUE
      write(6,1000)'INTER'
      do 40 ispin=1,jspin
      do 40 i=1,nat+1
         write(6,1010)ispin,i,(alph(i,j,ispin),j=1,4)
  40  CONTINUE
! 2000 write(6,*)'mixin',alph,iterho,iatom
 2000 call inter(alph,iterho,nat,itecon)
      return
 1000 format(a5)
 1010 format(2i2,4f10.4)
      end
