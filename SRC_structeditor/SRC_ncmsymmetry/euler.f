subroutine euler(mod,rot_new,a,b,c)
  implicit real*8 (a-h,o-z)
  
  integer mod
  real tmp

  dimension rot_old(3,3),rot_new(3,3),z(3),zz(3)
  dimension y(3),yy(3),yyy(3),pom(3),x(3),xx(3)
  
  zero=1.0d-20
  pi=acos(-1d0)
  
  do i=1,3
     do j=1,3
        rot_old(i,j)=0
        if (i.eq.j) rot_old(i,i)=1
     end do
  end do
3 format(20x,3(f10.8,2x))
  
  do j=1,3
     y(j)=rot_old(j,2)
     yyy(j)=rot_new(j,2)
     z(j)=rot_old(j,3)
     zz(j)=rot_new(j,3)
  end do
  
  call vecprod(z,zz,yy)
  y_norm=dsqrt(dot(yy,yy))
  
  if (y_norm.lt.1d-10) then
     
     if (abs(dot(y,yyy)).gt.1d0) then
        aa=dot(y,yyy)/abs(dot(y,yyy))
        a=acos(aa)
     else
        a=acos(dot(y,yyy))
     end if
     
     if (dot(z,zz).gt.zero) then
        c=zero
        b=zero
        if (yyy(1).gt.zero) a=2*pi-a
     else
        c=a
        a=zero
        b=pi
        if (yyy(1).lt.zero) c=2*pi-c
     end if
     
  else
     
     do j=1,3
        yy(j)=yy(j)/y_norm
     end do
     
     aa=dot(y,yy)
     bb=dot(z,zz)
     cc=dot(yy,yyy)
     if (abs(aa).gt.1d0) aa=aa/abs(aa)
     if (abs(bb).gt.1d0) bb=bb/abs(bb)
     if (abs(cc).gt.1d0) cc=cc/abs(cc)
     b=acos(bb)
     a=acos(aa)
     c=acos(cc)
     if (yy(1).gt.zero) a=2*pi-a
     call vecprod(yy,yyy,pom)
     if (dot(pom,zz).lt.-zero) c=2*pi-c
  end if

 if (mod.eq.1) then
!   space fixed axes (exchange a,c)
    tmp=a
    a=c
    c=tmp
 endif 


end subroutine euler


