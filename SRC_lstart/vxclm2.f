      SUBROUTINE VXCLM2(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
                FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
                exu,exd,GDGGT,igrad) 
      IMPLICIT REAL*8(A-H,O-Z)
!
      if((fu.le.1.d-18).or.(fd.le.1.d-18))then
      vxu=0.d0
      vxd=0.d0
      return
      endif
!
!      igd=mod(igrad,10)
!      IF(igrad.EQ.11)THEN
!      CALL VXCPW2(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
!                FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
!                exu,exd,GDGGT)
!      RETURN
!      ENDIF
!     IF(igrad.EQ.12)THEN
!     CALL VBECKE(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU,
!    *          FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD,GDGGT)
!     RETURN
!     ENDIF
! potential option 13: PBE (Perdew-Burke-Enzenhofer) 96 GGA
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
! potential option 11: GGA-WC05 for exchange + GGA-PBE for correlation
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
      IF(igrad.EQ.14)THEN
      CALL VPW91(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
                FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
                exu,exd,GDGGT)
      RETURN
      ENDIF
!
      IF(igrad.EQ.15)THEN
      CALL VXC15(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
                FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
                exu,exd,GDGGT)
      RETURN
      ENDIF
!
      IF(igrad.EQ.17)THEN
      CALL VXC17(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
                FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
                exu,exd,GDGGT)
      RETURN
      ENDIF
!
      IF(IGRAD.EQ.24)THEN
      CALL VXC24(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
                FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
                exu,exd,GDGGT)
      RETURN
      ENDIF
!
      IF(IGRAD.EQ.25)THEN
      CALL EV92(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
                FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
                exu,exd,GDGGT)
      RETURN
      ENDIF

      IF(IGRAD.EQ.29)THEN
      CALL GEA(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
               FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
               exu,exd,GDGGT)
      RETURN
      ENDIF
      RETURN
      END
