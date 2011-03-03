!MASS-Lib http://www.rs6000.ibm.com/resource/technology/MASS
      subroutine vrec(y,x,n)
      implicit none
      real*8 x(*),y(*)
      integer n,j
      do 10 j=1,n
      y(j)=1.d0/x(j)
   10 continue
      return
      end
      subroutine vsincos(x,y,z,n)
      implicit none
      real*8 x(*),y(*),z(*)
      integer n, j
      do 10 j=1,n
      x(j)=sin(z(j))
      y(j)=cos(z(j))
   10 continue
      return
      end
      subroutine vcos(y,x,n)
      implicit none
      real*8 x(*),y(*)
      integer n,j
      do 10 j=1,n
      y(j)=cos(x(j))
   10 continue
      return
      end
      subroutine vcosisin(y,x,n)
      implicit none
      complex*16 y(*)
      real*8 x(*)
      integer n,j
      do 10 j=1,n
      y(j)=dcmplx(cos(x(j)),sin(x(j)))
   10 continue
      return
      end
