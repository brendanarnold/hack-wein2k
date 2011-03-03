      subroutine drho(jatom,jrj,lmmax, &
                      clmsp,dcsp,ddcsp)
!
      use struct
      use work
      implicit real*8 (a-h,o-z)
!
      INCLUDE 'param.inc'
!
      real*8 :: clmsp(nrad,ncom,nat,2), &
                dcsp(nrad,ncom,2),ddcsp(nrad,ncom,2), &
       fdx(nrad),x(nrad),fx(nrad),b(nrad),c(nrad),d(nrad) 
!
      n=jrj
      do 100 j=1,n
 100  x(j)=dlog(r(j))
!
      do 1 lm1=1,lmmax
       do 1 ispin=1,2
!
        do 2 j=1,jrj
 2      fx(j)=clmsp(j,lm1,jatom,ispin)
        call spline(n,x,fx,b,c,d)
        do 3 j=1,jrj
         u=x(j)
         fdx(j)=sevald(n,u,x,b,c,d)
         dcsp(j,lm1,ispin)=(fdx(j)-2.d0* &
                            clmsp(j,lm1,jatom,ispin))/r(j)/r(j)/r(j)
 3      continue
        call spline(n,x,fdx,b,c,d)
        do 4 j=1,jrj
         u=x(j)
         fx(j)=sevald(n,u,x,b,c,d)
         ddcsp(j,lm1,ispin)=((fx(j)-4.d0*clmsp(j,lm1,jatom,ispin))/ &
                            r(j)/r(j)/r(j)-5.d0*dcsp(j,lm1,ispin))/ &
                            r(j)
 4      continue
 1    continue
      return
      end      
