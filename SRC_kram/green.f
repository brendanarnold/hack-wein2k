      subroutine green(fer,de,xe,grid,ngrid,ising)
!ad
!ad  routine written by Peter Vogl
!ad
! --------------------------------------------------------
!    this routine computes the hilbert transform with
!    simpson's rule. fer= p sum(de(e')/e-e')
!    xe..energy grid
!    grid...grid-size for de(e')
!    ngrid...nr of gridpoints for de
!    ising nr of gridpoint e, where singularity occurs 

      implicit double precision (a-h,o-z)
      dimension de(1),xe(1)

      gridb3=grid/3.
      es=xe(ising)
      fer=0.
      il=ising-1
      ir=ngrid-ising
      imax=max0(il,ir)
      do 300 inc=2,imax,2
      add1=0.
      add2=0.
      add3=0.
      incl=ising-inc
      incr=ising+inc
      if(inc.gt.2) goto 100
      dl1=de(incl)
      dl2=de(incl+1)
      dr1=de(incr)
      dr2=de(incr-1)
!ad
!
!    grid: .     ..*..      .
! prgr.prt:/ 100/ here/ 200 /
!    sum(e-2h,e+2h)(de(e')/e-e')=-4h*d'(e)
!    where h is increment (quadratic expansion)
!    d'(e)=(d(e-2h)-d(e+2h)+8d(e+h)-8d(e-h))/12h + o(d''''')
!
!ad
      add1= - ( dl1-dr1 + 8.*(dr2-dl2) )/3.
      goto 250
100   continue

      if(incl.lt.1) goto 200

!ad
!ad    simpson left side of singularity
!ad

      dl1=de(incl)/(es-xe(incl))
      dl2=de(incl+1)/(es-xe(incl+1))
      dl3=de(incl+2)/(es-xe(incl+2))
      add2=gridb3*(dl1+4.*dl2+dl3)
200   if(incr.gt.ngrid) goto 250

!ad
!ad    simpson right side of singularity
!ad
 
      dr1=de(incr)/(es-xe(incr))
      dr2=de(incr-1)/(es-xe(incr-1))
      dr3=de(incr-2)/(es-xe(incr-2))
      add3=gridb3*(dr1+4.*dr2+dr3)
250   fer=fer+add1+add2+add3
300   continue

!ad
      return
      end
