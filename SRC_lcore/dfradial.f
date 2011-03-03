      subroutine dfradial(r,f,df,jri)

      implicit real*8 (a-h,o-z)
      INCLUDE 'param.inc'

      dimension r(nrad)
      dimension f(nrad),df(nrad)
      dimension fdx(nrad),x(nrad),fx(nrad)
      dimension b(nrad),c(nrad),d(nrad) 
!
      n=jri
      DO 5005 j=1,n
        x(j)=dlog(r(j))
        fx(j)=f(j)
 5005 CONTINUE
     
      call spline(n,x,fx,b,c,d)
      DO 5006 j=1,jri
        u=x(j)
        fdx(j)=sevald(n,u,x,b,c,d)
        df(j)=fdx(j)/r(j)
 5006 CONTINUE

      return
      end      
