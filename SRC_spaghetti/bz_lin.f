      subroutine bz_lin(v,nv,lines,nlines,x,nbreak,break)
!     *****************************************
!
      IMPLICIT REAL*8 (A-H,O-Z)
!      INCLUDE 'param.inc'
!
      dimension  v(3,*)
      dimension  lines(*)
      dimension  x(*)
!
      dimension  v0(3),vdir(3)
      logical  break(*)
!
      data toler  /1.d-05/
!-----------------------------------------------------------------------
!
!.....INITIALIZE LINE-CHECK;  the 1. and the 2. k-point always build
!     a Brillouin-Zone line
      nlines=1
      lines(1)=1
      v0(1)  =v(1,1)
      v0(2)  =v(2,1)
      v0(3)  =v(3,1)
      vdir(1)=v(1,2) - v(1,1)
      vdir(2)=v(2,2) - v(2,1)
      vdir(3)=v(3,2) - v(3,1)
      xsum=sqrt( vdir(1)**2 + vdir(2)**2 + vdir(3)**2 )
      x(1)=0.
      x(2)=xsum
      nbreak=0
      break(1)=.false.
      break(2)=.false.
      dmax=10.*xsum
      dbreak=2.d0*xsum          
!     dmax indicates gap between 2 lines
      jkp=2
!
!.....START SEARCH LOOP FOR ALL OTHER K-POINTS
 10   continue
      jkp=jkp+1
         if(jkp.gt.nv) goto 100
         d=sqrt( (v(1,jkp)-v(1,jkp-1))**2 + (v(2,jkp)-v(2,jkp-1))**2 &
                  + (v(3,jkp)-v(3,jkp-1))**2 )
         if(d.gt.dbreak) then
	    d=dmax
	    xsum=xsum + d
            x(jkp)=xsum
	    nbreak=nbreak+1
	    break(jkp)=.true.
            nlines=nlines+1
            lines(nlines)=jkp-1
            nlines=nlines+1
            lines(nlines)=jkp
            v0(1)=v(1,jkp)
            v0(2)=v(2,jkp)
            v0(3)=v(3,jkp)
            vdir(1)=v(1,jkp+1)-v(1,jkp)
            vdir(2)=v(2,jkp+1)-v(2,jkp)
            vdir(3)=v(3,jkp+1)-v(3,jkp)
            dbreak=2.d0*sqrt( vdir(1)**2 + vdir(2)**2 + vdir(3)**2 )
            goto 10
         else
	    break(jkp)=.false.
         endif
	 xsum=xsum + d
         x(jkp)=xsum
         eps1=(v(1,jkp)-v0(1))*vdir(2) - (v(2,jkp)-v0(2))*vdir(1)
         eps2=(v(2,jkp)-v0(2))*vdir(3) - (v(3,jkp)-v0(3))*vdir(2)
         eps3=(v(1,jkp)-v0(1))*vdir(3) - (v(3,jkp)-v0(3))*vdir(1)
         if (abs(eps1).gt.toler .or. abs(eps2).gt.toler &
             .or. abs(eps3).gt.toler)  then
            nlines=nlines+1
            lines(nlines)=jkp-1
            v0(1)=v(1,jkp-1)
            v0(2)=v(2,jkp-1)
            v0(3)=v(3,jkp-1)
            vdir(1)=v(1,jkp)-v(1,jkp-1)
            vdir(2)=v(2,jkp)-v(2,jkp-1)
            vdir(3)=v(3,jkp)-v(3,jkp-1)
            dbreak=2.d0*sqrt( vdir(1)**2 + vdir(2)**2 + vdir(3)**2 )
         endif
      goto 10
!
!.....READY
 100  continue
      lines(nlines+1)=nv
      return
      end
