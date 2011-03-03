      subroutine dfrad(r,f,df,n)

      use param
      implicit real*8 (a-h,o-z)
      
      dimension r(nrad)
      dimension f(nrad),df(nrad)
      dimension fdx(nrad),x(nrad),fx(nrad)
      dimension b(nrad),c(nrad),d(nrad) 
!---
      IF (N .GT. NRAD) WRITE (*,*) 'In dfrad: N > NRAD'
!---
      
      DO 5000 j=1,n
         x(j)=dlog(r(j))
         fx(j)=f(j)
 5000 CONTINUE
      
      call spline(n,x,fx,b,c,d)

      DO 5001 j=1,n
         u=x(j)
         fdx(j)=sevald(n,u,x,b,c,d)
         df(j)=fdx(j)/r(j)
 5001 CONTINUE
      
      return
      end      
