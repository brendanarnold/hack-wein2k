      SUBROUTINE vxi35a(D,S,U,V,EX,VX)
!engel
!  This subroutine evaluates the exchange-only energy per electron 
!  and exchange-only potential for the GGA presented in the paper
!  E. Engel, S. H. Vosko, Phys. Rev. B47, 13164 (1993), Eq.(21).
!  Its structure and use are completely equivalent to the subroutine
!  EXCH of the PW91 subroutine package available from J. P. Perdew.
!  Thus this subroutine may replace EXCH in order to be used together
!  with the PW91 form for the correlation energy functional in the 
!  local density approximation, i.e. eliminating Perdew's subroutine
!  CORGGA.
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
      PARAMETER( a1 =  1.64712732976569d0 )
      PARAMETER( a2 =  .9801182732515644d0 )
      PARAMETER( a3 =  1.739938440729993d-02 )
      PARAMETER( b1 =  1.52367053964223d0 )
      PARAMETER( b2 =  .3672289494300232d0 )
      PARAMETER( b3 =  1.128174995137635d-02 )

      if( d.gt.0.d0 ) then

      droot = d**0.3333333333333333d0
      xi = s**2

      a   = 1.d0 + a1 * xi +        a2 * xi**2 +        a3 * xi**3
      ap  =        a1      + 2.d0 * a2 * xi    + 3.d0 * a3 * xi**2
      app =                  2.d0 * a2         + 6.d0 * a3 * xi
      b   = 1.d0 + b1 * xi +        b2 * xi**2 +        b3 * xi**3
      bp  =        b1      + 2.d0 * b2 * xi    + 3.d0 * b3 * xi**2
      bpp =                  2.d0 * b2         + 6.d0 * b3 * xi

      ra = ap / a
      rap = app / a - ra * ra
      rb = bp /b
      rbp = bpp / b - rb * rb
      f = a / b
      fp = f * ( ra - rb )
      fpp = fp * ( ra - rb ) + f * ( rap - rbp )

      vx = cvxlda*droot * ( f - 1.5d0*v*fp - (3.d0*u-4.d0*s**3)*s*fpp )
      ex = cexlda*droot * f

      vx=v
      ex=xi

      else

        vx = 0.d0
        ex = 0.d0

      end if

      return
      end


