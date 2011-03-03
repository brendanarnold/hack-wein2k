      SUBROUTINE GEAEX(D,S,U,V,EX,VX)
!
!      GEA
!
!  INPUT:  D = DENSITY
!          S = ABS(GRAD D)/(2*KF*D)
!          U = (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KF)**3)
!          V = (LAPLACIAN D)/(D*(2*KF)**2)
!  OUTPUT:  EX = EXCHANGE ENERGY PER ELECTRON (in atomic units)
!           VX = EXCHANGE POTENTIAL (in atomic units)
!engel
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER( cvxlda = - 0.9847450218426966d0 )
      PARAMETER( cexlda = - 0.738558766d0 )
      PARAMETER( a =  0.12345679 )

      if( d.gt.0.d0 ) then

      droot = d**0.3333333333333333d0
      xi = s**2

      f   = 1.d0 + a * xi 
      fp  =        a 
      vx = cvxlda*droot * ( f - 1.5d0*v*fp )
      ex = cexlda*droot * f

      else

        vx = 0.d0
        ex = 0.d0

      end if

      return
      end

