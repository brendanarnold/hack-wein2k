      subroutine FX_tpss(rho,gr,tau,EX)
      implicit double precision (a-h,o-z)

      parameter ( pi = 3.1415926535897932384626433832795d0 )
      parameter ( THRD = 1.d0/3.d0 )
      parameter ( ee = 1.537d0 )
      parameter ( cc = 0.39774d0 )
      parameter ( dk = 0.804d0 )
      parameter ( bb = 0.40d0 )

!c gr = the magnitude of the gradient of the electron density
!c tau = the kinetic energy density
!c FXMGGA = the enhancement factor of the spin-unpolarized density
      PARAMETER(AX=-0.738558766382022405884230032680836D0)
      EXUNIF = AX*RHO**(THRD)
  
      tau0 = 0.3d0*(3.d0*pi*pi)**(2.d0*THRD)*rho**(5.d0*THRD)
         if (rho.gt.1.0d-18) then
            tauw = 0.125d0*gr*gr/rho
         else 
            tauw = 0.0d0
         endif
      if(tauw.gt.tau.and.rho.lt.1.d-1) write(6,*)'tauw.gt.tau(fx_tpss)',tauw,tau
             if(tauw.gt.tau) tau=tauw
!      tauw = 0.125d0*gr*gr/rho
      p = gr*gr/( 4.d0*(3.d0*pi*pi)**(2.d0/3.d0)*rho**(8.d0/3.d0))
      alp = (tau - tauw)/tau0
      alp = dabs(alp)
      alpterm = (9.d0/20.d0)*(alp - 1.d0)/dsqrt(1.d0 + bb*alp*(alp - 1.d0))
      qtil = alpterm + 2.d0*p/3.d0

      xx = tauw/tau
      xx2 = xx*xx

      dmu = 10.d0/81.d0
      pbemu = 0.21951d0
      d1 = 146.d0/2025.d0
      d2 = -73.d0/405.d0
      effec = cc*4.d0*xx2/(1.d0 + xx2)**2.d0
      hh = 2.d0*dsqrt(ee)*dmu*(3.d0*xx/5.d0)**2.d0
      ptil = dsqrt(0.5*(0.6d0*xx)**2.d0 + 0.5d0*p**2.d0)
      
      xtil = ((dmu + effec)*p + d1*qtil*qtil +              &
             d2*qtil*ptil + (dmu**2.d0/dk)*p**2.d0 + hh +   &
             ee*pbemu*p**3.d0)/(1.d0 + dsqrt(ee)*p)**2.d0
      FXMGGA = 1.d0 + dk - dk/(1.d0 + xtil/dk)
      EX = EXUNIF*FXMGGA
      return
      end
