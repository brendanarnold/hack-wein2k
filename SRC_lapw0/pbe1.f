!######################################################################
!----------------------------------------------------------------------
!pb      subroutine easypbe(up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap,
      subroutine pbe1(up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap, &
                 agr,delgr,lcor,lpot, &
          exuplsd,exdnlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd, &
          exuppw91,exdnpw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91, &
          exuppbe,exdnpbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)
!pb     1           exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
!pb     1           expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91,
!pb     1           expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! EASYPBE is a driver for the PBE subroutines, using simple inputs
! K. Burke, May 14, 1996.
! inputs: up=up density
!	: agrup=|grad up|
!	: delgrup=(grad up).(grad |grad up|) 
!	: uplap=grad^2 up=Laplacian of up
!	: dn,agrdn,delgrdn,dnlap=corresponding down quantities
!	: agr=|grad rho|
!	: delgr=(grad rho).(grad |grad rho|) 
!	: lcor=flag to do correlation(=0=>don't)
!	: lpot=flag to do potential(=0=>don't)
! outputs: exlsd=LSD exchange energy density, so that
!		ExLSD=int d^3r rho(r) exlsd(r)
!	 : vxuplsd=up LSD exchange potential
!	 : vxdnlsd=down LSD exchange potential
!        : exclsd=LSD exchange-correlation energy density
!	 : vxcuplsd=up LSD exchange-correlation potential
!	 : vxcdnlsd=down LSD exchange-correlation potential
!        : expw91,vxuppw91,vxdnpw91,ecpw91,etc.=PW91 quantities
!        : expbe,vxuppbe,vxdnpbe,ecpbe,etc.=PBE quantities
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! needed constants:
! pi32=3 pi**2
! alpha=(9pi/4)**thrd
      implicit real*8(a-h,o-z)
      parameter(thrd=1.d0/3.d0,thrd2=2.d0*thrd)
      parameter(pi32=29.608813203268075856503472999628d0)
      parameter(pi=3.1415926535897932384626433832795d0)
      parameter(alpha=1.91915829267751300662482032624669d0)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! PBE exchange
! use  Ex[up,dn]=0.5*(Ex[2*up]+Ex[2*dn]) (i.e., exact spin-scaling)
! do up exchange
! fk=local Fermi wavevector for 2*up=(3 pi^2 (2up))^(1/3) 
! s=dimensionless density gradient=|grad rho|/ (2*fk*rho)_(rho=2*up)
! u=delgrad/(rho^2*(2*fk)**3)_(rho=2*up)
! v=Laplacian/(rho*(2*fk)**2)_(rho=2*up)
      rho2=2.d0*up
      if(rho2.gt.1d-18)then
        fk=(pi32*rho2)**thrd
        s=2.d0*agrup/(2.d0*fk*rho2)
        u=4.d0*delgrup/(rho2*rho2*(2.d0*fk)**3)
        v=2.d0*uplap/(rho2*(2.d0*fk)**2)
        call exchpbe(rho2,s,u,v,0,lpot,exuplsd,vxuplsd)
        call exchpw91(rho2,s,u,v,exuppw91,vxuppw91)
        call exchpbe(rho2,s,u,v,1,lpot,exuppbe,vxuppbe)
      else
	exuplsd=0.d0
	vxuplsd=0.d0
	exuppw91=0.d0
	vxuppw91=0.d0
	exuppbe=0.d0
	vxuppbe=0.d0
      endif
! repeat for down
      rho2=2.d0*dn
      if(rho2.gt.1d-18)then
        fk=(pi32*rho2)**thrd
        s=2.d0*agrdn/(2.d0*fk*rho2)
        u=4.d0*delgrdn/(rho2*rho2*(2.d0*fk)**3)
        v=2.d0*dnlap/(rho2*(2.d0*fk)**2)
        call exchpbe(rho2,s,u,v,0,lpot,exdnlsd,vxdnlsd)
        call exchpw91(rho2,s,u,v,exdnpw91,vxdnpw91)
        call exchpbe(rho2,s,u,v,1,lpot,exdnpbe,vxdnpbe)
      else
	exdnlsd=0.d0
	vxdnlsd=0.d0
	exdnpw91=0.d0
	vxdnpw91=0.d0
	exdnpbe=0.d0
	vxdnpbe=0.d0
      endif
10    continue 
! construct total density and contribution to ex
      rho=up+dn
      exlsd=(exuplsd*up+exdnlsd*dn)/rho
      expw91=(exuppw91*up+exdnpw91*dn)/rho
      expbe=(exuppbe*up+exdnpbe*dn)/rho
      if(lcor.eq.0)return
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Now do correlation
! zet=(up-dn)/rho
! g=phi(zeta)
! rs=(3/(4pi*rho))^(1/3)=local Seitz radius=alpha/fk
! sk=Ks=Thomas-Fermi screening wavevector=sqrt(4fk/pi)
! twoksg=2*Ks*phi
! t=correlation dimensionless gradient=|grad rho|/(2*Ks*phi*rho)
! uu=delgrad/(rho^2*twoksg^3)
! rholap=Laplacian
! vv=Laplacian/(rho*twoksg^2)
! ww=(|grad up|^2-|grad dn|^2-zet*|grad rho|^2)/(rho*twoksg)^2
! ec=lsd correlation energy
! vcup=lsd up correlation potential
! vcdn=lsd down correlation potential
! h=gradient correction to correlation energy
! dvcup=gradient correction to up correlation potential
! dvcdn=gradient correction to down correlation potential
      if(rho.lt.1.d-18)return
      zet=(up-dn)/rho
      g=((1.d0+zet)**thrd2+(1.d0-zet)**thrd2)/2.d0
      fk=(pi32*rho)**thrd
      rs=alpha/fk
      sk=sqrt(4.d0*fk/pi)
      twoksg=2.d0*sk*g
      t=agr/(twoksg*rho)
      uu=delgr/(rho*rho*twoksg**3)
      rholap=uplap+dnlap
      vv=rholap/(rho*twoksg**2)
      ww=(agrup**2-agrdn**2-zet*agr**2)/(rho*rho*twoksg**2)
      call CORPBE(RS,ZET,T,UU,VV,WW,1,lpot,ec,vcup,vcdn, &
                        H,DVCUP,DVCDN)
      eclsd=ec
      ecpbe=ec+h
      vcuplsd=vcup
      vcdnlsd=vcdn
      vcuppbe=vcup+dvcup
      vcdnpbe=vcdn+dvcdn
      call CORLSD1(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
      call CORPW91(RS,ZET,G,EC,ECRS,ECZET,T,UU,VV,WW,H,DVCUP,DVCDN)
      ecpw91=ec+h
      vcuppw91=vcup+dvcup
      vcdnpw91=vcdn+dvcdn
      return
      end
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE EXCHPBE(rho,S,U,V,lgga,lpot,EX,VX)
!----------------------------------------------------------------------
!  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
!  K Burke's modification of PW91 codes, May 14, 1996
!  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  INPUT rho : DENSITY
!  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
!  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
!  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
!   (for U,V, see PW86(24))
!  input lgga:  (=0=>don't put in gradient corrections, just LDA)
!  input lpot:  (=0=>don't get potential and don't need U and V)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
!     {\bf 40},  3399  (1989) (E).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Formulas:
!   	e_x[unif]=ax*rho^(4/3)  [LDA]
! ax = -0.75*(3/pi)^(1/3)
!	e_x[PBE]=e_x[unif]*FxPBE(s)
!	FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
! uk, ul defined after [a](13) 
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(thrd=1.d0/3.d0,thrd4=4.d0/3.d0)
      parameter(pi=3.14159265358979323846264338327950d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
      parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! construct LDA exchange energy density
      exunif = AX*rho**THRD
      if(lgga.eq.0)then
	ex=exunif
        vx=ex*thrd4
	return
      endif
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! construct PBE enhancement factor
      S2 = S*S
      P0=1.d0+ul*S2
      FxPBE = 1d0+uk-uk/P0
      EX = exunif*FxPBE
      if(lpot.eq.0)return
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first and second derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
!  Fss=d Fs/ds
      Fs=2.d0*uk*ul/(P0*P0)
      Fss=-4.d0*ul*S*Fs/P0
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! calculate potential from [b](24) 
      VX = exunif*(THRD4*FxPBE-(U-THRD4*S2*s)*FSS-V*FS)
      RETURN
      END
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE CORPBE(RS,ZET,T,UU,VV,WW,lgga,lpot,ec,vcup,vcdn, &
                        H,DVCUP,DVCDN)
!----------------------------------------------------------------------
!  Official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
!       : UU=(GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
!       : VV=(LAPLACIAN rho)/(rho * (2*KS*G)**2)
!       : WW=(GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2
!       :  UU,VV,WW, only needed for PBE potential
!       : lgga=flag to do gga (0=>LSD only)
!       : lpot=flag to do potential (0=>energy only)
!  output: ec=lsd correlation energy from [a]
!        : vcup=lsd up correlation potential
!        : vcdn=lsd dn correlation potential
!        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
!        : dvcup=nonlocal correction to vcup
!        : dvcdn=nonlocal correction to vcdn
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, 
!     {\sl Generalized gradient approximation made simple}, sub.
!     to Phys. Rev.Lett. May 1996.
! [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
!     construction of a generalized gradient approximation:  The PW91
!     density functional}, submitted to Phys. Rev. B, Feb. 1996.
! [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      parameter(thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd)
      parameter(sixthm=thrdm/2.d0)
      parameter(thrd4=4.d0*thrd)
      parameter(GAM=0.5198420997897463295344212145565d0)
      parameter(fzz=8.d0/(9.d0*GAM))
      parameter(gamma=0.03109069086965489503494086371273d0)
      parameter(bet=0.06672455060314922d0,delt=bet/gamma)
      parameter(eta=1.d-12)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor2(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0, &
          0.49294D0,rtrs,EU,EURS)
      CALL gcor2(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0, &
          0.62517D0,rtRS,EP,EPRS)
      CALL gcor2(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0, &
          0.49671D0,rtRS,ALFM,ALFRSM)
      ALFC = -ALFM
      Z4 = ZET**4
      F=((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZET=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
              -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      if(lgga.eq.0)return
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! PBE correlation energy
! G=phi(zeta), given after [a](3)
! DELT=bet/gamma
! B=A of [a](8)
      G=((1.d0+ZET)**thrd2+(1.d0-ZET)**thrd2)/2.d0
      G3 = G**3
      PON=-EC/(G3*gamma)
      B = DELT/(DEXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      H = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
      if(lpot.eq.0)return
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
      G4 = G3*G
      T6 = T4*T2
      RSTHRD = RS/3.D0
      GZ=(((1.d0+zet)**2+eta)**sixthm- &
      ((1.d0-zet)**2+eta)**sixthm)/3.d0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      hB = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      hBT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      hRST = RSTHRD*T2*hBT*BEC*ECRS
      hZ = 3.D0*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
      hT = 2.d0*BET*G3*Q9/Q8
      hZT = 3.D0*GZ*hT/G+hBT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      hTT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
      RETURN
      END
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE GCOR2(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
! slimmed down version of GCOR used in PW91 routines, to interpolate
! LSD correlation energy, as given by (10) of
! J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
! K. Burke, May 11, 1996.
      IMPLICIT REAL*8 (A-H,O-Z)
      Q0 = -2.D0*A*(1.D0+A1*rtrs*rtrs)
      Q1 = 2.D0*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2.D0*B2+rtrs*(3.D0*B3+4.D0*B4*rtrs))
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1*(1.d0+Q1))
      RETURN
      END
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE EXCHPW91(D,S,U,V,EX,VX)
!  GGA91 EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
!  INPUT D : DENSITY
!  INPUT S:  ABS(GRAD D)/(2*KF*D)
!  INPUT U:  (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KF)**3)
!  INPUT V: (LAPLACIAN D)/(D*(2*KF)**2)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(a1=0.19645D0,a2=0.27430D0,a3=0.15084D0,a4=100.d0)
      parameter(ax=-0.7385588D0,a=7.7956D0,b1=0.004d0)
      parameter(thrd=0.333333333333D0,thrd4=1.33333333333D0)
! for Becke exchange, set a3=b1=0
      FAC = AX*D**THRD
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
      P0 = 1.D0/DSQRT(1.D0+A*A*S2)
      P1 = DLOG(A*S+1.D0/P0)
      P2 = DEXP(-A4*S2)
      P3 = 1.D0/(1.D0+A1*S*P1+B1*S4)
      P4 = 1.D0+A1*S*P1+(A2-A3*P2)*S2
      F = P3*P4
      EX = FAC*F
!  LOCAL EXCHANGE OPTION
!     EX = FAC
!  ENERGY DONE. NOW THE POTENTIAL:
      P5 = B1*S2-(A2-A3*P2)
      P6 = A1*S*(P1+A*S*P0)
      P7 = 2.D0*(A2-A3*P2)+2.D0*A3*A4*S2*P2-4.D0*B1*S2*F
      FS = P3*(P3*P5*P6+P7)
      P8 = 2.D0*S*(B1-A3*A4*P2)
      P9 = A1*P1+A*A1*S*P0*(3.D0-A*A*S2*P0*P0)
      P10 = 4.D0*A3*A4*S*P2*(2.D0-A4*S2)-8.D0*B1*S*F-4.D0*B1*S3*FS
      P11 = -P3*P3*(A1*P1+A*A1*S*P0+4.D0*B1*S3)
      FSS = P3*P3*(P5*P9+P6*P8)+2.D0*P3*P5*P6*P11+P3*P10+P7*P11
      VX = FAC*(THRD4*F-(U-THRD4*S3)*FSS-V*FS)
!  LOCAL EXCHANGE OPTION:
!     VX = FAC*THRD4
      RETURN
      END
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE CORLSD1(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
!  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
!  INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
!  OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
!     POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET (ECZET)
!  OUTPUT: CORRELATION CONTRIBUTION (ALFC) TO THE SPIN STIFFNESS
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(gam=0.5198421D0,fzz=1.709921D0)
      parameter(thrd=0.333333333333D0,thrd4=1.333333333333D0)
      F = ((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
      CALL GCOR1(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0, &
          0.49294D0,1.00D0,RS,EU,EURS)
      CALL GCOR1(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0, &
          0.62517D0,1.00D0,RS,EP,EPRS)
      CALL GCOR1(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0, &
          0.49671D0,1.00D0,RS,ALFM,ALFRSM)
!  ALFM IS MINUS THE SPIN STIFFNESS ALFC
      ALFC = -ALFM
      Z4 = ZET**4
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
!  ENERGY DONE. NOW THE POTENTIAL:
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
              -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      RETURN
      END
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE GCOR1(A,A1,B1,B2,B3,B4,P,RS,GG,GGRS)
!  CALLED BY SUBROUTINE CORLSD
      IMPLICIT REAL*8 (A-H,O-Z)
      P1 = P + 1.D0
      Q0 = -2.D0*A*(1.D0+A1*RS)
      RS12 = DSQRT(RS)
      RS32 = RS12**3
      RSP = RS**P
      Q1 = 2.D0*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RS12+2.D0*B2+3.D0*B3*RS12+2.D0*B4*P1*RSP)
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
      RETURN
      END
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE CORpw91(RS,ZET,G,EC,ECRS,ECZET,T,UU,VV,WW,H, &
                         DVCUP,DVCDN)
!  pw91 CORRELATION, modified by K. Burke to put all arguments 
!  as variables in calling statement, rather than in common block
!  May, 1996.
!  INPUT RS: SEITZ RADIUS
!  INPUT ZET: RELATIVE SPIN POLARIZATION
!  INPUT T: ABS(GRAD D)/(D*2.*KS*G)
!  INPUT UU: (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KS*G)**3)
!  INPUT VV: (LAPLACIAN D)/(D * (2*KS*G)**2)
!  INPUT WW:  (GRAD D)*(GRAD ZET)/(D * (2*KS*G)**2
!  OUTPUT H: NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
!  OUTPUT DVCUP,DVCDN:  NONLOCAL PARTS OF CORRELATION POTENTIALS
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(xnu=15.75592D0,cc0=0.004235D0,cx=-0.001667212D0)
      parameter(alf=0.09D0)
      parameter(c1=0.002568D0,c2=0.023266D0,c3=7.389D-6,c4=8.723D0)
      parameter(c5=0.472D0,c6=7.389D-2,a4=100.D0)
      parameter(thrdm=-0.333333333333D0,thrd2=0.666666666667D0)
      BET = XNU*CC0
      DELT = 2.D0*ALF/BET
      G3 = G**3
      G4 = G3*G
      PON = -DELT*EC/(G3*BET)
      B = DELT/(DEXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      T6 = T4*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      Q6 = C1+C2*RS+C3*RS2
      Q7 = 1.D0+C4*RS+C5*RS2+C6*RS3
      CC = -CX + Q6/Q7
      R0 = 0.663436444d0*rs
      R1 = A4*R0*G4
      COEFF = CC-CC0-3.D0*CX/7.D0
      R2 = XNU*COEFF*G3
      R3 = DEXP(-R1*T2)
      H0 = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
      H1 = R3*R2*T2
      H = H0+H1
!  LOCAL CORRELATION OPTION:
!     H = 0.0D0
!  ENERGY DONE. NOW THE POTENTIAL:
      CCRS = (C2+2.*C3*RS)/Q7 - Q6*(C4+2.*C5*RS+3.*C6*RS2)/Q7**2
      RSTHRD = RS/3.D0
      R4 = RSTHRD*CCRS/COEFF
      GZ = ((1.D0+ZET)**THRDM - (1.D0-ZET)**THRDM)/3.D0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      H0B = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      H0RS = -RSTHRD*H0B*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      H0BT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      H0RST = RSTHRD*T2*H0BT*BEC*ECRS
      H0Z = 3.D0*GZ*H0/G + H0B*(BG*GZ+BEC*ECZET)
      H0T = 2.*BET*G3*Q9/Q8
      H0ZT = 3.D0*GZ*H0T/G+H0BT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      H0TT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      H1RS = R3*R2*T2*(-R4+R1*T2/3.D0)
      FACT4 = 2.D0-R1*T2
      H1RST = R3*R2*T2*(2.D0*R4*(1.D0-R1*T2)-THRD2*R1*T2*FACT4)
      H1Z = GZ*R3*R2*T2*(3.D0-4.D0*R1*T2)/G
      H1T = 2.D0*R3*R2*(1.D0-R1*T2)
      H1ZT = 2.D0*GZ*R3*R2*(3.D0-11.D0*R1*T2+4.D0*R1*R1*T4)/G
      H1TT = 4.D0*R3*R2*R1*T*(-2.D0+R1*T2)
      HRS = H0RS+H1RS
      HRST = H0RST+H1RST
      HT = H0T+H1T
      HTT = H0TT+H1TT
      HZ = H0Z+H1Z
      HZT = H0ZT+H1ZT
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
!  LOCAL CORRELATION OPTION:
!     DVCUP = 0.0D0
!     DVCDN = 0.0D0
      RETURN
      END
!----------------------------------------------------------------------
