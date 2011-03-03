      subroutine drho(jri,clmsp,dcsp,ddcsp,r)
!
      implicit real*8 (a-h,o-z)
      INCLUDE 'param.inc'
!
!
!
      dimension clmsp(npt,2),dcsp(npt,2),ddcsp(npt,2),r(npt), &
       fdx(npt),x(npt),fx(npt),b(npt),c(npt),d(npt) 
!
      n=jri
      do 100 j=1,n
 100  x(j)=dlog(r(j))
!
       do 1 ispin=1,2
!
        do 2 j=1,jri
!        write (6,*) 'clmsp',clmsp(j,ispin)
 2      fx(j)=clmsp(j,ispin)
        call spline(n,x,fx,b,c,d)
        do 3 j=1,jri
         u=x(j)
         fdx(j)=sevald(n,u,x,b,c,d)
         dcsp(j,ispin)=(fdx(j)-2.d0* &
                        clmsp(j,ispin))/r(j)/r(j)/r(j)
 3      continue
        call spline(n,x,fdx,b,c,d)
        do 4 j=1,jri
         u=x(j)
         fx(j)=sevald(n,u,x,b,c,d)
         ddcsp(j,ispin)=((fx(j)-4.d0*clmsp(j,ispin))/ &
                            r(j)/r(j)/r(j)-5.d0*dcsp(j,ispin))/ &
                            r(j)
 4      continue
 1    continue
!      write(6,*) 'r',(r(ir),ir=371,jri)
!      write(6,*) 'clmsp',(clmsp(ir,1),ir=371,jri)
!      write(6,*) 'dcsp',(dcsp(ir,1),ir=371,jri)
!      write(6,*) 'ddcsp',(ddcsp(ir,1),ir=371,jri)

      return
      end      

