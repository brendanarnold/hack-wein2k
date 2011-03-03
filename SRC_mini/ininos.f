      subroutine ininos(xmin,x0,delta,ekin,n,weight,mult)
      implicit double precision (a-h,o-z)
      include 'param.inc'
      DIMENSION  xmin(n),x0(n)
      COMMON /STRUK/  AA,BB,CC,ALPHA(3),PIA(3),VOL
      DIMENSION  weight(*),mult(*)
      ekin=0.0d0
      do 10 i=1,n,3
         help=dble(i+2)/3.0d0
         iatom=nint(help)
      call distxy(x0(i),x0(i+1),x0(i+2), &
      xmin(i),xmin(i+1),xmin(i+2),dists)
         write (6,*) 'distance', dists
         do 20 j=1,mult(iatom)
            ekin=ekin+weight(iatom)*((dists/delta)**2)/2.0d0
  20     CONTINUE
  10  CONTINUE
      return
      end
