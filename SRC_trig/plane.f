!Definition of a plane for lapw5 by 3 atoms.
!Specify lattice parameters (works so far only for orthogonal systems or H lat)
!        3 atoms               First atom=center of plot,
!                              second atom below
!                              third atom to the left
!        x-,y-length of plot (bohr).
implicit real*8 (a-h,o-z)
dimension a(3),x0(3),x1(3),x2(3),xn(3),xd1(3),xd2(3),xdhelp(3)
character*1 lat
lat=' '
read(*,*) a
read(*,*) x0
read(*,*) x1
read(*,*) x2
read(*,*) alaenge, blaenge
read(*,*,end=1) lat
 1     continue
!           convert from fractions to bohr
x0=x0*a
x1=x1*a
x2=x2*a
if(lat.eq.'H') then
x0(1)=x0(1)*sqrt(3.d0)/2.d0
x1(1)=x1(1)*sqrt(3.d0)/2.d0
x2(1)=x2(1)*sqrt(3.d0)/2.d0
x0(2)=-x0(1)/sqrt(3.d0)+x0(2)
x1(2)=-x1(1)/sqrt(3.d0)+x1(2)
x2(2)=-x2(1)/sqrt(3.d0)+x2(2)
endif
!
xd1=x0-x1
xd2=x0-x2


!
call vecprod(xd1,xd2,xn)
xdhelp=xd2
call vecprod(xd1,xn,xd2)

xlaenge=flaenge(xd1)
xd1=xd1*blaenge/xlaenge

xlaenge=flaenge(xd2)
xd2=xd2*alaenge/xlaenge   

cosalpha=dot_product(xd2,xdhelp)/flaenge(xd2)/flaenge(xdhelp)

print*, 'cos alpha',cosalpha

x0=x0 - xd1/2. - xd2/2.*sign(1.d0,cosalpha)
x1=x0 + xd1
x2=x0 + xd2*sign(1.d0,cosalpha)

!         convert from bohr to fractions
x0=x0/a
x1=x1/a
x2=x2/a
if(lat.eq.'H') then
x0(1)=x0(1)*2.d0/sqrt(3.d0)
x1(1)=x1(1)*2.d0/sqrt(3.d0)
x2(1)=x2(1)*2.d0/sqrt(3.d0)
!x0(2)=x0(1)/sqrt(3.d0)+x0(2)
!x1(2)=x1(1)/sqrt(3.d0)+x1(2)
!x2(2)=x2(1)/sqrt(3.d0)+x2(2)
x0(2)=x0(1)/2.d0+x0(2)
x1(2)=x1(1)/2.d0+x1(2)
x2(2)=x2(1)/2.d0+x2(2)
endif

write(*,*) nint(x0(1)*100000),nint(x0(2)*100000),nint(x0(3)*100000),100000
write(*,*) nint(x2(1)*100000),nint(x2(2)*100000),nint(x2(3)*100000),100000
write(*,*) nint(x1(1)*100000),nint(x1(2)*100000),nint(x1(3)*100000),100000

end

subroutine vecprod(x,y,z)
implicit real*8 (a-h,o-z)
dimension x(3),y(3),z(3)
z(1)=x(2)*y(3)-x(3)*y(2)
z(2)=-(x(1)*y(3)-x(3)*y(1))
z(3)=x(1)*y(2)-x(2)*y(1)
end

function flaenge(x)
implicit real*8 (a-h,o-z)
dimension x(3)
flaenge=sqrt(x(1)**2+x(2)**2+x(3)**2)
end
