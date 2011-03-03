      subroutine mold(xmin1,x0,xplus1,n,fwert,delta,weight,temp,mult)
      implicit double precision (a-h,o-z)
      include 'param.inc'
      COMMON /STRUK/  AA,BB,CC,ALPHA(3),PIA(3),VOL
      DIMENSION  xmin1(n),x0(n),xplus1(n),fwert(n)
      DIMENSION  weight(*),mult(*)
!
!     force and gradient in Ry/bohr 
!
      iatom=1
      ihelp=0
      do 100 i=1,n
         ihelp=ihelp+1
         write(6,*) 'iatom',iatom
         xplus1(i)=2.0d0*x0(i)-xmin1(i)+ &
                   fwert(i)/2.0d0*(delta**2)/weight(iatom)
         if (ihelp.eq.3) then
            ihelp=0
            iatom=iatom+1
         endif
 100  CONTINUE
!
! calculate temperature
!
      boltzk=3.166678911d-6
      free=-6.0d0
      do 10 i=1,n,3
         help=dble(i+2)/3.0d0
         iatom=nint(help)
         do 20 j=1,mult(iatom)
            free=free+3.0d0
  20     CONTINUE
  10  CONTINUE
      if (free.lt.0.1d0) free=1.0d0
      call ininos(x0,xplus1,delta,ekin,n,weight,mult)
      temp=2.0d0*ekin/boltzk/free
      write(6,*)'  new position'
      write(6,*)(xplus1(i),i=1,n)
      write(6,*)'fwert',(fwert(i),i=1,n)
      write(6,*)'  temperature',temp
550   format( 20x,3f12.4 )
      return
      end
