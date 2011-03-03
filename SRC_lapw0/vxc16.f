      SUBROUTINE VXC16(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
                FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
                exu,exd,GDGGT)
      IMPLICIT REAL*8(A-H,O-Z)
      save
      COMMON /GAS/ FK,SK,G,EC,ECRS,ECZET
      DATA THRD,FTDS/0.3333333333333333D0,1.3333333333333333D0/
      DATA TTDS/0.6666666666666667/
      read(5,*,end=10) rhokup,rhokdn
      write(21,9000)rhokup,rhokdn
 9000 format(':krh  :',2f10.7)
 10   PI=4.D0*ATAN(1.D0)
      PI2=PI*PI
! EXCHANGE - UP
      D=2.D0*FU
      GDGG=GXU*GGXU+GYU*GGYU+GZU*GGZU
      TKF=2.D0*(3.D0*PI2*D)**THRD
      S=GMAGU/(TKF*FU)
      U=GDGG/(TKF*(FU*TKF)**2)
      V=G2U/(FU*TKF**2)
      VXU=rhokup*vxu/fu
      CALL EXCH(D,S,U,V,EXUP,VXUP)
      VXUP=vxu+2.0d0*exup
! EXCHANGE - DN
      D=2.D0*FD
      GDGG=GXD*GGXD+GYD*GGYD+GZD*GGZD
      TKF=2.D0*(3.D0*PI2*D)**THRD
      S=GMAGD/(TKF*FD)
      U=GDGG/(TKF*(FD*TKF)**2)
      V=G2D/(FD*TKF**2)
      VXD=rhokdn*vxd/fd
      CALL EXCH(D,S,U,V,EXDN,VXDN)
      VXDN=VXD+2.0d0*exdn
! LOCAL CORRELATION
      D=FU+FD
      RS=(0.75D0/(PI*D))**THRD
      ZET=(FU-FD)/D
      CALL CORLSD(RS,ZET,ECLSD,VLSDU,VLSDD,ECRS,ECZET,ALFC)
      EC=ECLSD
! GRADIENT CORRECTION TO CORRELATION
      G=0.5D0*((1.D0+ZET)**TTDS+(1.D0-ZET)**TTDS)
      GP=GXU*GXD+GYU*GYD+GZU*GZD
      GDGZET=((1.D0-ZET)*GMAGU**2-(1.D0+ZET)*GMAGD**2-2.D0*ZET*GP)/D
      FK=(3.D0*PI2*D)**THRD
      SK=SQRT(4.D0*FK/PI)
      TKSG=2.D0*G*SK
      T=SQRT((GXU+GXD)**2+(GYU+GYD)**2+(GZU+GZD)**2)/(D*TKSG)
      U=GDGGT/(TKSG*(TKSG*D)**2)
      V=(G2U+G2D)/(D*TKSG**2)
      W=GDGZET/(D*TKSG**2)
      CALL CORGGA(RS,ZET,T,U,V,W,ECGGA,VGGAU,VGGAD)
!
!      VXU=2.D0*(VXUP)
!      VXD=2.D0*(VXDN)
      VXU=2.D0*(VXUP+VLSDU+VGGAU)
      VXD=2.D0*(VXDN+VLSDD+VGGAD)
      exu=2.d0*(exup+eclsd+ecgga)
      exd=2.d0*(exdn+eclsd+ecgga)
!      exu=2.d0*(exup+eclsd)
!      exd=2.d0*(exdn+eclsd)
!
      RETURN
      END
