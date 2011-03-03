      SUBROUTINE CORPBE_tpss(RS,ZETA,TT,ECC,HH,ECPBE)

!C  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!C       : ZETA=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!C       : tt=ABS(GRAD rho)/(rho*2.*KS*GG)  -- only needed for PBE
!C       : HH,ecgga,vcgga
!c----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
!c numbers for use in LSD energy spin-interpolation formula.
!c      GAM= 2^(4/3)-2
!c      FZZ=f''(0)= 8/(9*GAM)
!c numbers for construction of PBE
!c      gamma=(1-log(2))/pi^2
!c      betamb=coefficient in gradient expansion for correlation.
      parameter( pi = 3.1415926535897932384626433832795d0 )
      parameter(thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd)
      parameter(sixthm=thrdm/2.d0)
      parameter(thrd4=4.d0*thrd)
      parameter(GAM=0.5198420997897463295344212145565d0)
      parameter(fzz=8.d0/(9.d0*GAM))
      parameter(gamma=0.03109069086965489503494086371273d0)
      parameter(betamb=0.06672455060314922d0,DELT=betamb/gamma)
      parameter ( FAC = 3.09366772628013593097d0 )
!c     FAC = (3*PI**2)**(1/3)
      parameter ( ALPHA = 1.91915829267751300662482032624669d0 )
      cc = (3.0d0*pi*pi/16.0d0)**thrd

!c EC0=unpolarized LSD correlation energy
!c EC0RS=dEC0/drs
!c EC1=fully polarized LSD correlation energy
!c EC1RS=dEC1/drs
!c ALFM=-spin stiffness
!c ALFMRS=-dalphac/drs
!c F=spin-scaling factor
!c construct LSD ec
      F = ((1.D0+ZETA)**THRD4+(1.D0-ZETA)**THRD4-2.D0)/GAM
      CALL GCOR(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0, &
         0.49294D0,1.00D0,RS,EC0,EC0RS)
      CALL GCOR(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0, &
         0.62517D0,1.00D0,RS,EC1,EC1RS)
      CALL GCOR(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0, &
         0.49671D0,1.00D0,RS,ALFM,ALFMRS)
!C  ALFM IS MINUS THE SPIN STIFFNESS ALFC
      ALFC = -ALFM
      Z4 = ZETA**4.0d0
      ECC = EC0*(1.D0-F*Z4)+EC1*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
!c----------------------------------------------------------------------
!c PBE correlation energy 
!c GG=phi(zeta)
!c DELT=betamb/gamma
!c B=A
      gg = 0.5d0*((1.0d0+zeta)**thrd2 + (1.0d0-zeta)**thrd2)
      GG3 = GG**3.0d0
      PON=-ECC/(GG3*gamma)
      B = DELT/(DEXP(PON)-1.D0)
      B2 = B*B
      TT2 = TT*TT
      TT3 = TT2*TT
      TT4 = TT2*TT2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*TT2
      Q5 = 1.D0+B*TT2+B2*TT4
      HH = GG3*gamma*DLOG(1.D0+DELT*TT2*Q4/Q5)
      ECPBE = ECC + HH
      return
      end

