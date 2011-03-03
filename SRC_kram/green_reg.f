      subroutine green_reg(fer,de,xe,grid,ngrid,xsing)
!ad
!ad  routine written by Peter Vogl
!ad
! --------------------------------------------------------
!    this routine computes the hilbert transform with
!    simpson's rule. fer= p sum(de(e')/e-e')
!    xe..energy grid
!    grid...grid-size for de(e')
!    ngrid...nr of gridpoints for de
!    xsing.. energy in same units as xe, where singularity occurs.
!            must be outside of integration

      implicit double precision (a-h,o-z)
      dimension de(1),xe(1)

      gridb3=grid/3.
      es=xsing
      if( (es .lt. xe(ngrid)) .and. (es .gt. xe(1)) ) &
        stop 'error in green_reg: singiularity not outside'
!
      fer=0.
!
      do 300 incl=1,ngrid-2,2
         dl1=de(incl)/(es-xe(incl))
         dl2=de(incl+1)/(es-xe(incl+1))
         dl3=de(incl+2)/(es-xe(incl+2))

         fer = fer + gridb3*(dl1+4.*dl2+dl3)
!
300   continue
      return
      end
