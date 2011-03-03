      SUBROUTINE VXCLM2(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
                FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
                exu,exd,GDGGT,tauup,taudn,igrad,ir) 
      IMPLICIT REAL*8(A-H,O-Z)
      integer index
      data index/0/
      save index
!      
!      
! selects the GGA functional by the number given in input0,
! adding 50 means radial gradient only
!
!    
!
!.....for low density use LDA to avoid oszillations
      if((fu.le.1.d-6).or.(fd.le.1.d-6))then
        if(fu.lt.1.d-8) fu=1.d-8
        if(fd.lt.1.d-8) fd=1.d-8
      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
                 t,gdggt,1,1, &
             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
!      exu=2.d0*(exuppb+ecpbe)
!      exd=2.d0*(exdnpb+ecpbe)
!      vxu=2.d0*(vxuppb+vcuppb)
!      vxd=2.d0*(vxdnpb+vcdnpb)
      exu=2.d0*(exupls  +  eclsd)
      exd=2.d0*(exdnls  +  eclsd)
      vxu=2.d0*(vxupls+vcupls)
      vxd=2.d0*(vxdnls+vcdnls)
!      vxu=0.d0
!      vxd=0.d0
!      exu=0.d0
!      exd=0.d0
      return
      endif
!
!      igd=mod(igrad,10)
!
! potential option 11: GGA-WC for exchange + GGA-PBE for correlation
!
      IF (IGRAD .EQ. 11) THEN

      GDGGU=GXU*GGXU+GYU*GGYU+GZU*GGZU
      GDGGD=GXD*GGXD+GYD*GGYD+GZD*GGZD
      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)

      CALL WC05(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
                exupwc05,exdnwc05,vxupwc05,vxdnwc05)
      exu = 2D0*EXUPWC05
      exd = 2D0*EXDNWC05
      vxu = 2D0*VXUPWC05
      vxd = 2D0*VXDNWC05

      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
                 t,gdggt,1,1, &
             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
      exu = exu + 2D0*ecpbe
      exd = exd + 2D0*ecpbe
      vxu = vxu + 2D0*vcuppb
      vxd = vxd + 2D0*vcdnpb

      RETURN
      ENDIF
!
! potential option 12: HPBE GGA-XC-Energy (PKZB), PBE potential
!
      IF(igrad.EQ.12) THEN
      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
                 t,gdggt,1,1, &
             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
      vxu=2.d0*(vxuppb+vcuppb)
      vxd=2.d0*(vxdnpb+vcdnpb)
      tauup=(tauup+2.d0*g2u*0.25d0-vxu*fu)*0.5d0
      taudn=(taudn+2.d0*g2d*0.25d0-vxd*fd)*0.5d0
      if(ir.ne.index) then
      if(tauup.lt.0.0d0)  write(6,'(5e13.5)')  &
      tauup,2.d0*tauup-2.d0*g2u*.25d0+vxu*fu,2.d0*g2u*.25d0,-vxu*fu,fu
      endif
      index=ir
      if(tauup.lt.0.0d0) tauup=1.d-8
      if(taudn.lt.0.0d0) taudn=1.d-8
      call dfxhpbe(fu,gmagu,tauup,fd,gmagd,taudn,exu,exd)
!      if(ir.eq.410)  write(6,'(a11,5e13.5)') 'g2,-vx,tau',g2u,-vxu,tauup
!      if(ir.eq.781)  write(6,'(a11,5e13.5)') 'g2,-vx,tau',g2u,-vxu,tauup
!
!
!      if(ir.eq.781) write(*,'(a2,i4,5e13.5)') 'ir',ir,fu,fd,exu,exd
      exu=2.0d0*exu
      exd=2.0d0*exd
!
!  HPBE-C GGA
!
      call sicpbe(fu,fd,gmagu,gmagd,t,tauup,taudn,ecsicpbe)
!      if(ir.ne.index)write(*,'(2HIR,i4,f14.5,2e15.5)')ir,fu,exu,ecsicpbe
!      index=ir
!      if(ir.eq.781) write(*,'(a2,i4,5e13.5)') 'ir',ir,fu,fd,ecsicpbe
      exu=exu+2.0d0*ecsicpbe
      exd=exd+2.0d0*ecsicpbe
!
      RETURN
      ENDIF
!
! potential option 13: PBE (Perdew-Burke-Ernzerhof) 96 GGA
!
      IF(igrad.EQ.13) THEN
      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
                 t,gdggt,1,1, &
             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
      exu=2.d0*(exuppb+ecpbe)
      exd=2.d0*(exdnpb+ecpbe)
      vxu=2.d0*(vxuppb+vcuppb)
      vxd=2.d0*(vxdnpb+vcdnpb)
      RETURN
      ENDIF
!
! potential option 14: Perdew-Wang 91 GGA 
!
      IF(igrad.EQ.14)THEN
      CALL VPW91(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
                FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
                exu,exd,GDGGT)
      RETURN
      ENDIF
!
! potential option 15: Engel-Vosko 93 GGA exchange with 
!                      Perdew-Wang 91 GGA correlation term 
!
      IF(igrad.EQ.15)THEN
      CALL VXC15(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
                FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
                exu,exd,GDGGT)
      RETURN
      ENDIF
!
! potential option 16: GGA-RPBE for exchange + GGA-PBE for correlation
!
      IF (IGRAD .EQ. 16) THEN

      GDGGU=GXU*GGXU+GYU*GGYU+GZU*GGZU
      GDGGD=GXD*GGXD+GYD*GGYD+GZD*GGZD
      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)

      CALL RPBE(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
                exuprpbe,exdnrpbe,vxuprpbe,vxdnrpbe)
      exu = 2D0*EXUPRPBE
      exd = 2D0*EXDNRPBE
      vxu = 2D0*VXUPRPBE
      vxd = 2D0*VXDNRPBE

      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
                 t,gdggt,1,1, &
             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
      exu = exu + 2D0*ECPBE
      exd = exd + 2D0*ECPBE
      vxu = vxu + 2D0*VCUPPB
      vxd = vxd + 2D0*VCDNPB

      RETURN
      ENDIF
!
! potential option 17: BPW91 (GGA-Becke88 for exchange + GGA-PW91 for correlation)
!
      IF (IGRAD .EQ. 17) THEN

      GDGGU=GXU*GGXU+GYU*GGYU+GZU*GGZU
      GDGGD=GXD*GGXD+GYD*GGYD+GZD*GGZD
      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)

      CALL B88(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
               exupb88,exdnb88,vxupb88,vxdnb88)
      exu = 2D0*EXUPB88
      exd = 2D0*EXDNB88
      vxu = 2D0*VXUPB88
      vxd = 2D0*VXDNB88

      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
                 t,gdggt,1,1, &
             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
      exu = exu + 2D0*ECPW
      exd = exd + 2D0*ECPW
      vxu = vxu + 2D0*VCUPPW
      vxd = vxd + 2D0*VCDNPW

      RETURN
      ENDIF
!
! potential option 18: LDA and GGA parts of 3-parameters hybrid functional B3PW91
!
      IF (IGRAD .EQ. 18) THEN

      PI=4D0*ATAN(1D0)
      FT=FU+FD
      RS=(0.75D0/(PI*FT))**(1D0/3D0)
      ZET=(FU-FD)/FT
      GDGGU=GXU*GGXU+GYU*GGYU+GZU*GGZU
      GDGGD=GXD*GGXD+GYD*GGYD+GZD*GGZD
      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)

      AXLDA=0.28D0
      exu=AXLDA*XCENER(FT,FU,PI,5)
      exd=AXLDA*XCENER(FT,FD,PI,5)
      vxu=AXLDA*XCPOT (FT,FU,PI,5)
      vxd=AXLDA*XCPOT (FT,FD,PI,5)

      CALL CORLSD(RS,ZET,ECLSD,VLSDU,VLSDD,ECRS,ECZET,ALFC)
      ACLDA=0.19D0
      exu=exu+ACLDA*2D0*ECLSD
      exd=exd+ACLDA*2D0*ECLSD
      vxu=vxu+ACLDA*2D0*VLSDU
      vxd=vxd+ACLDA*2D0*VLSDD

      CALL B88(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
               exupb88,exdnb88,vxupb88,vxdnb88)
      AXGGA=0.72D0
      exu=exu+AXGGA*2D0*EXUPB88
      exd=exd+AXGGA*2D0*EXDNB88
      vxu=vxu+AXGGA*2D0*VXUPB88
      vxd=vxd+AXGGA*2D0*VXDNB88

      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
                 t,gdggt,1,1, &
             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
      ACGGA=0.81D0
      exu=exu+ACGGA*2D0*ECPW
      exd=exd+ACGGA*2D0*ECPW
      vxu=vxu+ACGGA*2D0*VCUPPW
      vxd=vxd+ACGGA*2D0*VCDNPW

      RETURN
      ENDIF
!
! potential option 25: Engel-Vosko 93 GGA exchange with
!                      LSDA correlation term 
!
      IF(IGRAD.EQ.25)THEN
      CALL EV92(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
                FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
                exu,exd,GDGGT)
      RETURN
      ENDIF
!
! potential option 26: Engel-Vosko 93 GGA exchange potential with
!                      LSDA potential correlation term combined with 
!                      Perdew-Wang 91 GGA exchange energy with
!                      Perdew-Wang 91 GGA correlation term 
!
      IF(IGRAD.EQ.26)THEN
      CALL VXC26(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
                FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
                exu,exd,GDGGT)
      RETURN
      ENDIF
!
! potential option 27: TPSS meta GGA-XC-Energy, PBE potential
!
      IF(igrad.EQ.27) THEN
      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
                 t,gdggt,1,1, &
             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
      vxu=2.d0*(vxuppb+vcuppb)
      vxd=2.d0*(vxdnpb+vcdnpb)
      tauup=(tauup+2.d0*g2u*0.25d0-vxu*fu)*0.5d0
      taudn=(taudn+2.d0*g2d*0.25d0-vxd*fd)*0.5d0
      if(ir.ne.index) then
      if(tauup.lt.0.0d0)  write(6,'(5e13.5)')  &
      tauup,2.d0*tauup-2.d0*g2u*.25d0+vxu*fu,2.d0*g2u*.25d0,-vxu*fu,fu
      endif
      index=ir
      if(tauup.lt.0.0d0) tauup=1.d-8
      if(taudn.lt.0.0d0) taudn=1.d-8
      call dfxtpss(fu,gmagu,tauup,fd,gmagd,taudn,t,exu,exd)
!      if(ir.eq.410)  write(6,'(a11,5e13.5)') 'g2,-vx,tau',g2u,-vxu,tauup
!      if(ir.eq.781)  write(6,'(a11,5e13.5)') 'g2,-vx,tau',g2u,-vxu,tauup
!
!
!      if(ir.eq.781) write(*,'(a2,i3,5e13.5)') 'ir',ir,fu,fd,exu,exd
      exu=2.0d0*exu
      exd=2.0d0*exd
!
!  TPSS-C GGA
!
      grupdn=GXU*GXD+GYU*GYD+GZU*GZD
      call sicpbe_tpss(fu,fd,gmagu,gmagd,t,tauup,taudn,grupdn,ecsicpbe)
!      if(ir.ne.index)write(*,'(2HIR,i4,f14.5,2e15.5)')ir,fu,exu,ecsicpbe
!      index=ir
!      if(ir.eq.781) write(*,'(a2,i3,5e13.5)') 'ir',ir,fu,fd,ecsicpbe
      exu=exu+2.0d0*ecsicpbe
      exd=exd+2.0d0*ecsicpbe
!
      RETURN
      ENDIF
!
! potential option 30:  EX = 1.d0, E-TOT should give number of electrons
!
      IF(igrad.EQ.30) THEN
      exu=1.d0
      exd=1.d0
      vxu=1.d0
      vxd=1.d0
      RETURN
      ENDIF
!
! potential option 33:  EX-PBE only
!
      IF(igrad.EQ.33) THEN
      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
                 t,gdggt,1,1, &
             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
      exu=2.d0*(exuppb)
      exd=2.d0*(exdnpb)
      vxu=2.d0*(vxuppb)
      vxd=2.d0*(vxdnpb)
!  modified for hybrid functionals of EECE
!     vxu=2.d0*(vxuppb+vcuppb)
!     vxd=2.d0*(vxdnpb+vcdnpb)
      RETURN
      ENDIF
!
! potential option 35: for analysing xi and eta
!                      case.vns and case.vsp contains xi and eta in terms of YLM
!
      IF(IGRAD.EQ.35)THEN
      CALL VXi35(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
                FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
                exu,exd,GDGGT)
      RETURN
      ENDIF
!
! potential option 36: for analysing (laplace rho)
!                      case.vns and case.vsp contains laplace rho in terms of YLM
!
      IF(IGRAD.EQ.36) THEN
      vxu=g2u
      vxd=g2d
      exu=0.d0
      exd=0.d0
      RETURN
      ENDIF
!
! potential option 40: HPBE GGA test
!
      IF(igrad.EQ.40) THEN
      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
                 t,gdggt,1,1, &
             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
      exu=2.0d0*g2u
      exd=2.0d0*g2d
      vxu=2.d0*(vxuppb+vcuppb)
      vxd=2.d0*(vxdnpb+vcdnpb)
      return
      endif
!
! potential option 41: HPBE GGA (X only)
!
      IF(igrad.EQ.41) THEN
      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
                 t,gdggt,1,1, &
             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
      vxu=2.d0*(vxuppb+vcuppb)
      vxd=2.d0*(vxdnpb+vcdnpb)
      tsaveup=tauup
      tsavedn=taudn
      tauup=(tauup+2.d0*g2u*0.25d0-vxu*fu)*0.5d0
      taudn=(taudn+2.d0*g2d*0.25d0-vxd*fd)*0.5d0
      tauwup=gmagu**2/8.d0/fu 
      tauwdn=gmagd**2/8.d0/fd 
!cc      if(tauwup.gt.tauup) tauup=tauwup
!cc      if(tauwdn.gt.taudn) taudn=tauwdn
      if((ir.ne.index)) then
      if(tauwup.gt.tauup) write(6,'(a14,5e13.5)') 'tauwup.gt.tauw',tauwup,tauup,fu
!      if(tauup.lt.0.0d0)  write(6,'(5e13.5)')  &
!      tauwup,tauup,fu
      endif
      index=ir
      if(ir.eq.1) then
         if(tauwup.gt.tauup) write(6,'(a14,5e13.5)')'tauwup.gt.tauw',tauwup,tauup,fu
         if(tauwup.lt.tauup) write(6,'(a17,5e13.5)')'tauwup.lt.tauw,ok',tauwup,tauup,fu
      endif
      call dfxhpbe(fu,gmagu,tauup,fd,gmagd,taudn,exu,exd)
!      if(ir.eq.410)  write(6,'(a11,5e13.5)') 'g2,-vx,tau',g2u,-vxu,tauup
!      if(ir.eq.781)  write(6,'(a11,5e13.5)') 'g2,-vx,tau',g2u,-vxu,tauup
!     *tauup,2.d0*tauup-1.d0*g2u*.25d0+vxu*fu,1.d0*g2u*.25d0,-vxu*fu,fu
!cc      call dfxhpbe(fu,gmagu,tauup,fd,gmagd,taudn,exu,exd)
!
      if(tauwup.gt.tauup) tauup=tauwup
      if(tauwdn.gt.taudn) taudn=tauwdn
      vxu=tauup
      vxd=tauup
!
      exu=2.0d0*exu
      exd=2.0d0*exd
!
      return
      endif
!
! potential option 42: HPBE-C GGA
!
      IF(igrad.EQ.42) THEN
      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
                 t,gdggt,1,1, &
             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
      vxu=2.d0*(vxuppb+vcuppb)
      vxd=2.d0*(vxdnpb+vcdnpb)
      ecu1=ecpbe
      ecd1=ecpbe
      T=SQRT((GXU)**2+(GYU)**2+(GZU)**2)
      gdggt1=gdggu
      zero=0.d0
      zero1=1.d-14
      zero2=1.d-14
      zero3=1.d-14
      zero4=1.d-14
      call pbe1(fu,gmagu,gdggu,g2u,zero1,zero2,zero3,zero4, &
                t,gdggt1,1,1, &
             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
      ecu2=ecpbe
      T=SQRT((GXd)**2+(GYd)**2+(GZd)**2)
      gdggt1=gdggd
      call pbe1(0.d0,0.d0,0.d0,0.d0,fd,gmagd,gdggd,g2d, &
                 t,gdggt1,1,1, &
             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
      ecd2=ecpbe
!
      tsaveup=tauup
      tsavedn=taudn
      tauup=(tauup+2.d0*g2u*0.25d0-vxu*fu)*0.5d0
      taudn=(taudn+2.d0*g2d*0.25d0-vxd*fd)*0.5d0
      tauwup=gmagu**2/8.d0/fu 
      tauwdn=gmagd**2/8.d0/fd 
!
!cc      if(tauwup.gt.tauup) tauup=tauwup
!cc      if(tauwdn.gt.taudn) taudn=tauwdn
!
      f1=1.d0 + 0.53d0 *((tauwup+tauwdn)/(tauup+taudn))**2
      f2up=(1.d0+ 0.53d0) *(tauwup/tauup)**2 * ecu2
      f2dn=(1.d0+ 0.53d0) *(tauwdn/taudn)**2 * ecd2
!      f2up=ecu2 + 0.53d0 *(tauwup/tauup)**2 * ecu2
!      f2dn=ecd2 + 0.53d0 *(tauwdn/taudn)**2 * ecd2
!
      if(tauup.lt.0.0d0)  write(6,'(5e13.5)')  &
      tauup,2.d0*tauup-2.d0*g2u*.25d0+vxu*fu,2.d0*g2u*.25d0,-vxu*fu,fu
      if(ir.eq.410)  write(6,'(a11,5e13.5)') 'tau,tauw',tauup,tauwup
      if(ir.eq.781)  write(6,'(a11,5e13.5)') 'tau,tauw',tauup,tauwup
!     *tauup,2.d0*tauup-1.d0*g2u*.25d0+vxu*fu,1.d0*g2u*.25d0,-vxu*fu,fu
!cc      call dfxhpbe(fu,gmagu,tauup,fd,gmagd,taudn,exu,exd)
!
      vxu=tauwup/max(tauup,tauwup)
      vxd=tauwdn/max(taudn,tauwdn)
!
      exu=2.0d0*(f1*ecu1-f2up)
      exd=2.0d0*(f1*ecd1-f2dn)
!
      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
      call sicpbe(fu,fd,gmagu,gmagd,t,tauup,taudn,ecsicpbe)
      exu=2.0d0*ecsicpbe
      exd=2.0d0*ecsicpbe
      if(ir.eq.371) write(6,'(a10,9e12.4,/,7e12.4)') 'esic,ec',exu, &
          2.0d0*(f1*ecu1-f2up),f1,ecu1,f2up,tauwup,tauup,ecu2
      if(ir.eq.371) write(6,'(6e12.4)')fu,gmagu,gdggu,g2u,zero1, &
                zero2,zero3,zero4, &
                t,gdggt1,1,1, &
             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb
!      vxu=2.0d0*(f1*ecu1-f2u)
!      vxd=2.0d0*(f1*ecd1-f2d)
!c       vxu=tauwup/tauup
!c       vxd=tauwdn/taudn
!
      RETURN
      ENDIF
!
! potential option 43: Van Voorhis+G.Scuseria M-GGA (VSXC) (PBE potential)
!
      IF(igrad.EQ.43) THEN
!
      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
                 t,gdggt,1,1, &
             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
      vxu=2.d0*(vxuppb+vcuppb)
      vxd=2.d0*(vxdnpb+vcdnpb)
      tauup=(tauup+2.d0*g2u*0.25d0-vxu*fu)*0.5d0
      taudn=(taudn+2.d0*g2d*0.25d0-vxd*fd)*0.5d0
      if(tauup.lt.0.0d0)  write(6,'(5e13.5)')  &
      tauup,2.d0*tauup-2.d0*g2u*.25d0+vxu*fu,2.d0*g2u*.25d0,-vxu*fu,fu
      if(tauup.lt.0.0d0) tauup=1.d-8
      if(taudn.lt.0.0d0) taudn=1.d-8
      call excvs98(fu,fd,gmagu,gmagd,tauup,taudn,exvs,ecvs)
      exu=2.D0*(exvs+ecvs)
      exd=2.D0*(exvs+ecvs)
      if(ir.ne.index)write(*,'(2HIR,i4,f14.5,2e15.5)') &
      ir,fu,2.d0*exvs,2.d0*ecvs
      index=ir
!      if(ir.eq.410)  write(6,'(a11,5e13.5)') 'g2,-vx,tau',g2u,-vxu,tauup
!      if(ir.eq.781)  write(6,'(a11,5e13.5)') 'g2,-vx,tau',g2u,-vxu,tauup
!
      RETURN
      ENDIF
!
!
!
!
!
!      IF(igrad.EQ.11)THEN
!      CALL VXCPW2(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU,
!     *          FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD,
!     *          exu,exd,GDGGT)
!      RETURN
!      ENDIF
!
!     IF(igrad.EQ.12)THEN
!     CALL VBECKE(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU,
!    *          FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD,GDGGT)
!     RETURN
!     ENDIF
!
! potential option 11: HPBE GGA-X-Energy, PBE-C-Energy, PBE potential
!
!      IF(igrad.EQ.11) THEN
!      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
!      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
!      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
!      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
!                 t,gdggt,1,1, &
!             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
!             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
!             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
!      vxu=2.d0*(vxuppb+vcuppb)
!      vxd=2.d0*(vxdnpb+vcdnpb)
!      tauup=(tauup+2.d0*g2u*0.25d0-vxu*fu)*0.5d0
!      taudn=(taudn+2.d0*g2d*0.25d0-vxd*fd)*0.5d0
!      if(tauup.lt.0.0d0)  write(6,'(a10,5e13.5)') 'tau-up neg',  &
!      tauup,2.d0*tauup-2.d0*g2u*.25d0+vxu*fu,2.d0*g2u*.25d0,-vxu*fu,fu
!      if(tauup.lt.0.0d0) tauup=1.d-8
!      if(taudn.lt.0.0d0) taudn=1.d-8
!      call dfxhpbe(fu,gmagu,tauup,fd,gmagd,taudn,exu,exd)
!      if(ir.eq.410)  write(6,'(a11,5e13.5)') 'g2,-vx,tau',g2u,-vxu,tauup
!      if(ir.eq.781)  write(6,'(a11,5e13.5)') 'g2,-vx,tau',g2u,-vxu,tauup
!
!
!      exu=2.0d0*(exu+ecpbe)
!      exd=2.0d0*(exd+ecpbe)
!
!      RETURN
!      ENDIF
!
! potential option 17: Lembarki exchange potential with 
!                      Perdew-Wang 91 GGA correlation term 
!
!      IF(igrad.EQ.17)THEN
!      CALL VXC17(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
!                FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
!                exu,exd,GDGGT)
!      RETURN
!      ENDIF
!
! potential option 18: TH1 (Tozer-Handy 98) GGA
!
!      IF(igrad.EQ.18) THEN
!      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
!      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
!      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
!      call th1(fu,gxu,gyu,gzu,gmagu,gdggu,g2u,fd,gxd,gyd,gzd, &
!               gmagd,gdggd,g2d,t,gdggt,exc,vxcu,vxcd)
!      if(fu.lt.1.0d0) write(*,*) fu,gxu,gxd,exc,vxcu,vxcd
!      exu=2.D0*exc
!      exd=2.D0*exc
!      vxu=2.D0*vxcu
!      vxd=2.D0*vxcd
!      RETURN
!      ENDIF
!
! potential option 19: TH2 (Tozer-Handy 98) GGA
!
!      IF(igrad.EQ.19) THEN
!      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
!      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
!      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
!      call th2(fu,gxu,gyu,gzu,gmagu,gdggu,g2u,fd,gxd,gyd,gzd, &
!               gmagd,gdggd,g2d,t,gdggt,exc,vxcu,vxcd)
!      if(fu.lt.1.0d0) write(*,*) fu,gxu,gxd,exc,vxcu,vxcd
!      exu=2.D0*exc
!      exd=2.D0*exc
!      vxu=2.D0*vxcu
!      vxd=2.D0*vxcd
!      if(fu.lt.1.0d0.and.fu.gt.0.999d0) write(*,*) 
!     * fu,gxu,gyu,gzu,gmagu,gdggu,g2u,fd,gxd,gyd,gzd,
!     *         gmagd,gdggd,g2d,t,gdggt,exc,vxcu,vxcd
!      RETURN
!      ENDIF
!
! potential option 20: HCTH (Tozer-Handy 98) GGA
!
!      IF(igrad.EQ.20) THEN
!      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
!      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
!      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
!      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
!                 t,gdggt,1,1, &
!             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
!             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
!             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
!      vxu=2.d0*(vxuppb+vcuppb)
!      vxd=2.d0*(vxdnpb+vcdnpb)
!
!      zab=gxu*gxd+gyu*gyd+gzu*gzd
!      call hcth(dfdra, dfdza, dfdrb, dfdzb, dfdzab,  &
!                        fu, fd, gmagu, gmagd, zab, .true., exc)
!      exu=2.D0*exc
!      exd=2.D0*exc
!      exu=2.D0*exc/(fu+fd)
!      exd=2.D0*exc/(fu+fd)
! ????????
!      RETURN
!      ENDIF
!
! potential option 21: TH2 (Tozer-Handy 98)-V, HCTH-Energy
!
!      IF(igrad.EQ.21) THEN
!      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
!      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
!      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
!      call th2(fu,gxu,gyu,gzu,gmagu,gdggu,g2u,fd,gxd,gyd,gzd, &
!               gmagd,gdggd,g2d,t,gdggt,exc,vxcu,vxcd)
!      if(fu.lt.1.0d0) write(*,*) fu,gxu,gxd,exc,vxcu,vxcd
!      exu=2.D0*exc
!      exd=2.D0*exc
!      vxu=2.D0*vxcu
!      vxd=2.D0*vxcd
!      zab=gxu*gxd+gyu*gyd+gzu*gzd
!      call hcth(dfdra, dfdza, dfdrb, dfdzb, dfdzab,  &
!                        fu, fd, gmagu, gmagd, zab, .true., exc)
!      exu=2.D0*exc/(fu+fd)
!      exd=2.D0*exc/(fu+fd)
!      RETURN
!      ENDIF
!
! potential option 22: Filatov-Thiel GGA (PBE potential)
!
!      IF(igrad.EQ.22) THEN
!      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
!      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
!      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
!      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
!                 t,gdggt,1,1, &
!             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
!             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
!             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
!      vxu=2.d0*(vxuppb+vcuppb)
!      vxd=2.d0*(vxdnpb+vcdnpb)
!
!      call excft98(fu,fd,gmagu, gmagd, g2u,g2d, exft,ecft)
!      exu=2.D0*(exft+ecft)
!      exd=2.D0*(exft+ecft)
!      RETURN
!      ENDIF
!
! potential option 23: BLYP GGA (PBE potential)
!
!      IF(igrad.EQ.23) THEN
!      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
!      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
!      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
!      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
!                 t,gdggt,1,1, &
!             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
!             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
!             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
!      vxu=2.d0*(vxuppb+vcuppb)
!      vxd=2.d0*(vxdnpb+vcdnpb)
!
!      call excblyp(fu,fd,gmagu, gmagd, t, exb88,eclyp)
!      exu=2.D0*(exb88+eclyp)
!      exd=2.D0*(exb88+eclyp)
!      RETURN
!      ENDIF
!
! potential option 24: Revised PBE-X + PBE-C (pbe-potential)
!
!      IF(IGRAD.EQ.24)THEN
!      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
!      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
!      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
!      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
!                 t,gdggt,1,1, &
!             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
!             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
!             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
!      vxu=2.d0*(vxuppb+vcuppb)
!      vxd=2.d0*(vxdnpb+vcdnpb)
!      call exrevpbe(fu,fd,gmagu,gmagd,exrpbe)
!      exu=2.d0*(exrpbe+ecpbe)
!      exd=2.d0*(exrpbe+ecpbe)
!
!      RETURN
!      ENDIF
!
! potential option 29: GEA (with coefficient of Engel-Vosko 93 GGA) exchange with
!                      LSDA correlation term 
!
!      IF(IGRAD.EQ.29)THEN
!      CALL GEA(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
!               FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
!               exu,exd,GDGGT)
!      RETURN
!      ENDIF
!
! potential option 31:  EX-LDA only
!
!      IF(igrad.EQ.31) THEN
!      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
!      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
!      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
!      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
!                 t,gdggt,1,1, &
!             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
!             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
!             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
!      exu=2.d0*(exupls)
!      exd=2.d0*(exdnls)
!  modified for hybrid functionls of EECE
!      vxu=2.d0*(vxupls)
!      vxd=2.d0*(vxdnls)
!      vxu=2.d0*(vxuppb+vcuppb)
!      vxd=2.d0*(vxdnpb+vcdnpb)
!      RETURN
!      ENDIF
!
! potential option 32:  EC-LDA only
!
!      IF(igrad.EQ.32) THEN
!      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
!      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
!      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
!      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
!                 t,gdggt,1,1, &
!             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
!             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
!             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
!      exu=2.d0*(eclsd)
!      exd=2.d0*(eclsd)
!      vxu=2.d0*(vxuppb+vcuppb)
!      vxd=2.d0*(vxdnpb+vcdnpb)
!      RETURN
!      ENDIF
!
! potential option 34:  EC-PBE only
!
!      IF(igrad.EQ.34) THEN
!      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
!      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
!      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
!      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
!                 t,gdggt,1,1, &
!             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
!             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
!             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
!      exu=2.d0*(ecpbe)
!      exd=2.d0*(ecpbe)
!      vxu=2.d0*(vxuppb+vcuppb)
!      vxd=2.d0*(vxdnpb+vcdnpb)
!      RETURN
!      ENDIF
!
! potential option 44: HPBE GGA-X-Energy + KCIS-C energy, PBE potential
!
!      IF(igrad.EQ.44) THEN
!      GDGGu=GXU*GGXU+GYU*GGYU+GZU*GGZU
!      GDGGd=GXD*GGXD+GYD*GGYD+GZD*GGZD
!      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)
!      call pbe1(fu,gmagu,gdggu,g2u,fd,gmagd,gdggd,g2d, &
!                 t,gdggt,1,1, &
!             exupls,exdnls,vxupls,vxdnls,eclsd,vcupls,vcdnls, &
!             exuppw,exdnpw,vxuppw,vxdnpw,ecpw,vcuppw,vcdnpw, &
!             exuppb,exdnpb,vxuppb,vxdnpb,ecpbe,vcuppb,vcdnpb)
!      vxu=2.d0*(vxuppb+vcuppb)
!      vxd=2.d0*(vxdnpb+vcdnpb)
!      tauup=(tauup+2.d0*g2u*0.25d0-vxu*fu)*0.5d0
!      taudn=(taudn+2.d0*g2d*0.25d0-vxd*fd)*0.5d0
!      if(tauup.lt.0.0d0)  write(6,'(5e13.5)')  &
!      tauup,2.d0*tauup-2.d0*g2u*.25d0+vxu*fu,2.d0*g2u*.25d0,-vxu*fu,fu
!      if(tauup.lt.0.0d0) tauup=1.d-8
!      if(taudn.lt.0.0d0) taudn=1.d-8
!      call dfxhpbe(fu,gmagu,tauup,fd,gmagd,taudn,exu,exd)
!      if(ir.eq.410)  write(6,'(a11,5e13.5)') 'g2,-vx,tau',g2u,-vxu,tauup
!      if(ir.eq.781)  write(6,'(a11,5e13.5)') 'g2,-vx,tau',g2u,-vxu,tauup
!
!
!      exu=2.0d0*exu
!      exd=2.0d0*exd
!
!  HPBE-C GGA
!
!      call kciscor(fu,fd,gmagu,gmagd,t,tauup,taudn,ecsic)
!      if(ir.ne.index)write(*,'(2HIR,i4,f14.5,2e15.5)') ir,fu,exu,ecsic
!      index=ir
!      exu=exu+2.0d0*ecsic
!      exd=exd+2.0d0*ecsic
!
!      RETURN
!      ENDIF
!
! potential option 45: Leeuwen exchange potential with 
!                      Perdew-Wang 91 GGA correlation term 
!
!      IF(igrad.EQ.45)THEN
!      CALL VXC16(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
!                FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
!                exu,exd,GDGGT)
!      RETURN
!      ENDIF
!
      RETURN
      END
