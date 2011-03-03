!-------------------------------------------------------------------
!     exchange-correlation functional of Van Voorhis and Scuseria
!     (VS98)
!
! Reference: T. Van Voorhis, G.E. Scuseria
!            ``A novel form of the exchange-correlation 
!              energy functional''
!            J. Chem. Phys. 109, p.400 (1998)
!
! Input:  up: spin-up density
!         dn: spin-down density
!         grup: | grad up |
!         grdn: | grad dn |
!         tup: 0.5*sum_{i,occ} | grad psi_{i,up} |^2
!              kinetic energy density of up-electrons
!         tdn: 0.5*sum_{i,occ} | grad psi_{i,down} |^2
!              kinetic energy density of down-electrons
! Output: exvs:  VS98 exchange energy per particle
!         ecvs:  VS98 correlation energy per particle
!
! Requires subroutine corlsd to calculate the correlation energy
! per particle of the uniform electron gas using the parametrization
! of Perdew and Wang
!
!-------------------------------------------------------------------
      subroutine excvs98(up,dn,grup,grdn,tup,tdn,exvs,ecvs)
      implicit real*8 (a-h,o-z)
!
         if (up.gt.1.0d-18) then
            tauwup = 0.125d0*grup*grup/up
         else 
            tauwup = 0.0d0
         endif
         if (dn.gt.1.0d-18) then
            tauwdn = 0.125d0*grdn*grdn/dn
         else 
            tauwdn = 0.0d0
         endif
         tauw = tauwup + tauwdn
         tau = tup + tdn
       if(tauw.gt.tau.and.up.lt.1.d-1) write(6,*) 'tauw.gt.tau',tauw,tau
             if(tauwup.gt.tup) tup=tauwup
             if(tauwdn.gt.tdn) tdn=tauwdn

      exvs = 0.0d0
      ecvs = 0.0d0

      exup = 0.0d0
      exdn =0.0d0

!     spin-up exchange energy

      if (up.gt.1.0d-18) then
         call exvs98(up,grup,tup,exup)
         exup = exup*up/(up+dn)
      endif

!     spin-down exchange energy

      if (dn.gt.1.0d-18) then
         call exvs98(dn,grdn,tdn,exdn)
         exdn = exdn*dn/(up+dn)
      endif

      exvs = exup + exdn
      
      if ((up+dn).gt.1.0d-18) then
         call ecvs98(up,dn,grup,grdn,tup,tdn,ecvs)
      endif

      return
      end
!-------------------------------------------------------------------
      subroutine exvs98(rho,gr,tau,exvs)

      implicit real*8 (a-h,o-z)

      parameter ( PI = 3.1415926535897932384626433832795d0 )
      parameter ( THRD = 0.33333333333333333333333333333d0 )
      parameter ( FOURTHRD = 4.0d0*thrd)
      parameter ( FAC = 3.09366772628013593097d0 )
!     FAC = (3*PI**2)**(1/3)
      parameter ( ZERO = 0.0d0 )

      parameter ( ALP = 0.001867d0 )
      parameter ( AA = -0.9800d0 )
      parameter ( BB = -0.003557d0 )
      parameter ( CC = 0.006250d0 )
      parameter ( DD = -0.00002354d0 )
      parameter ( EE = -0.0001283d0 )
      parameter ( FF = 0.0003575d0 )

      cf = 0.6d0*FAC*FAC*(2.0d0**(2.0d0*THRD))
      ttau = 2.0d0*tau

      xx = gr/(rho**FOURTHRD)
      zz = ttau/(rho**(5.0d0*THRD)) - cf
      ga = 1.0d0 + ALP*(xx*xx + zz)

      ter = (DD*xx**4.0d0 + EE*xx*xx*zz + FF*zz*zz)/(ga**3.0d0)
      ter = ter + (BB*xx*xx + CC*zz)/(ga*ga) + AA/ga
      exvs = ter*rho**THRD

      return
      end

!-------------------------------------------------------------------

      subroutine ecvs98(up,dn,grup,grdn,tup,tdn,ecvs)
      implicit real*8 (a-h,o-z)

      parameter ( PI = 3.1415926535897932384626433832795d0 )
      parameter ( THRD = 0.33333333333333333333333333333d0 )
      parameter ( FOURTHRD = 4.0d0*thrd)
      parameter ( FAC = 3.09366772628013593097d0 )
!     FAC = (3*PI**2)**(1/3)

!     parameters for opposite spin correlation 

      parameter ( ALP1 = 0.003050d0 )
      parameter ( AA1 = 0.7035d0 )
      parameter ( BB1 = 0.007695d0 )
      parameter ( CC1 = 0.051530d0 )
      parameter ( DD1 = 0.00003394d0 )
      parameter ( EE1 = -0.001269d0 )
      parameter ( FF1 = 0.001296d0 )

!     parameters for same spin correlation 

      parameter ( ALP2 = 0.005151d0 )
      parameter ( AA2 = 0.3271d0 )
      parameter ( BB2 = -0.03229d0 )
      parameter ( CC2 = -0.02942d0 )
      parameter ( DD2 = 0.002134d0 )
      parameter ( EE2 = -0.005452d0 )
      parameter ( FF2 = 0.01578d0 )

      ttup = 2.0d0*tup
      ttdn = 2.0d0*tdn

      rho = up + dn
      zeta = (up - dn)/rho
      rs = (3.0d0/(4.0d0*PI*rho))**THRD
      call corlsd(rs,zeta,ec1,vc1,vc2,dum1,dum2,dum3)

      if (up.gt.1.0d-18) then
         rs = (3.0d0/(4.0d0*PI*up))**THRD
         zeta = 1.0d0
         call corlsd(rs,zeta,ec0up,vc1,vc2,dum1,dum2,dum3)
      else
         ec0up = 0.0d0
      endif
 
      if (dn.gt.1.0d-18) then 
         rs = (3.0d0/(4.0d0*PI*dn))**THRD
         zeta = 1.0d0
         call corlsd(rs,zeta,ec0dn,vc1,vc2,dum1,dum2,dum3)
      else
         ec0dn = 0.0d0
      endif

!     opposite spin correlation in LDA

      ecopp = ec1 - up*ec0up/rho - dn*ec0dn/rho

!     same spin correlation in LDA
!     spin-up
 
      ecupl = up*ec0up/rho

!     spin-down
   
      ecdnl = dn*ec0dn/rho

!     prepare for calculation of VS correlation

      cf = 0.6d0*FAC*FAC*(2.0d0**(2.0d0*THRD))
      if (up.gt.1.0d-18) then 
         xxup = grup/(up**FOURTHRD)
         zzup = ttup/(up**(5.0d0*THRD)) - cf
      else 
         xxup = 0.0d0
         zzup = 0.0d0
      endif
      if (dn.gt.1.0d-18) then 
         xxdn = grdn/(dn**FOURTHRD)
         zzdn = ttdn/(dn**(5.0d0*THRD)) - cf
      else
         xxdn = 0.0d0
         zzdn = 0.0d0
      endif
      xx = dsqrt(xxup*xxup + xxdn*xxdn)
      zz = zzup + zzdn

!     opposite spin correlation 

      ga = 1.0d0 + ALP1*(xx*xx + zz)
      ter = (DD1*xx**4.0d0 + EE1*xx*xx*zz + FF1*zz*zz)/(ga**3.0d0)
      ter = ter + (BB1*xx*xx + CC1*zz)/(ga*ga) + AA1/ga
      ecopvs = ter*ecopp 

!     same spin correlation
!     spin-up
      
      if (up.gt.1.0d-18) then
         ga = 1.0d0 + ALP2*(xxup*xxup + zzup)
         ter = (DD2*xxup**4.0d0 + EE2*xxup*xxup*zzup +  &
                 FF2*zzup*zzup)/(ga**3.0d0) 
         ter = ter + (BB2*xxup*xxup + CC2*zzup)/(ga*ga)
         ter = ter + AA2/ga
         ddup = 1.0d0 - xxup*xxup/(4.0d0*(zzup + cf))
         ecs1vs = ter*ddup*ecupl
      else
         ecs1vs = 0.0d0
      endif

!     spin-down

      if (dn.gt.1.0d-18) then
         ga = 1.0d0 + ALP2*(xxdn*xxdn + zzdn)
         ter = (DD2*xxdn**4.0d0 + EE2*xxdn*xxdn*zzdn + &
                 FF2*zzdn*zzdn)/(ga**3.0d0)
         ter = ter + (BB2*xxdn*xxdn + CC2*zzdn)/(ga*ga)
         ter = ter + AA2/ga
         dddn = 1.0d0 - xxdn*xxdn/(4.0d0*(zzdn + cf))
         ecs2vs = ter*dddn*ecdnl
      else
         ecs2vs = 0.0d0
      endif

      ecvs = ecopvs + ecs1vs + ecs2vs 
      return 
      end

!-------------------------------------------------------------------
!     end of subroutines for VS98 exchange-correlation functional 
!-------------------------------------------------------------------
