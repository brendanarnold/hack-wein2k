      SUBROUTINE PWXAD5(iex,I1,I2,I3,tvec,am, &
       fft,tauup,taudn,ifft1,ifft2,ifft3,vxu,vxd,exu,exd)
!
!
      IMPLICIT REAL*8(A-H,O-Z)
      complex*16 fft(ifft1,ifft2,ifft3,2)
!
      DIMENSION TVEC(10,27), &
                     BVT(10),BV1(10),BV2(10),AM(125)
!
!   EXCHANGE-CORRELATION POTENTIAL SECTION
!
      DO 72 LX1=1,10
      BV1(LX1)=0.D0
   72 BV2(LX1)=0.D0
      IND=0
      DO 73 I=I1-1,I1+1
       DO 73 J=I2-1,I2+1
        DO 73 K=I3-1,I3+1
         IND=IND+1
         ia=i
         ja=j
         ka=k
         if(i.lt.1) ia=i+ifft1
         if(j.lt.1) ja=j+ifft2
         if(k.lt.1) ka=k+ifft3
         if(i.gt.ifft1) ia=i-ifft1
         if(j.gt.ifft2) ja=j-ifft2
         if(k.gt.ifft3) ka=k-ifft3
         DUP=fft(ia,ja,ka,1)
         DDN=fft(ia,ja,ka,2)
          DO 373 LX2=1,10
          BV1(LX2)=BV1(LX2)+DUP*TVEC(LX2,IND)
  373     BV2(LX2)=BV2(LX2)+DDN*TVEC(LX2,IND)
      if(i1.eq.1.and.i2.eq.1.and.i3.eq.1)  &
       write(6,1819) i,j,k,dup,bv1
 1819 format(3i3,11f10.7)
   73 CONTINUE
!
      CALL CHSLV(AM,10,BV1)
!      if(i1.eq.1.and.i2.eq.1.and.i3.eq.1) 
!     * write(6,*) bv1
      CALL CHSLV(AM,10,BV2)
      DO 274 LX2=1,10
  274 BVT(LX2)=BV1(LX2)+BV2(LX2)
      GXU=BV1(2)
      GYU=BV1(3)
      GZU=BV1(4)
      GMAGU=dSQRT(GXU*GXU+GYU*GYU+GZU*GZU)
      G2U=2.D0*(BV1(5)+BV1(6)+BV1(7))
      if(abs(gmagu).le.1.d-10) then
          ggxu =2.d0*bv1(5)
          ggyu =2.d0*bv1(6)
          ggzu =2.d0*bv1(7)
          gdgu=g2u
      else
          FACT=1.D0/GMAGU
          GGXU=FACT*(2.D0*BV1(5)*GXU+BV1(10)*GYU+BV1(9)*GZU)
          GGYU=FACT*(BV1(10)*GXU+2.D0*BV1(6)*GYU+BV1(8)*GZU)
          GGZU=FACT*(BV1(9)*GXU+BV1(8)*GYU+2.D0*BV1(7)*GZU)
      end if
      GX=BVT(2)
      GY=BVT(3)
      GZ=BVT(4)
      gmagt=dSQRT(GX*GX+GY*GY+GZ*GZ)
      if(abs(gmagt).le.1.d-10) then
          ggx =2.d0*bv1(5)
          ggy =2.d0*bv1(6)
          ggz =2.d0*bv1(7)
          gdggt=2.D0*(BVt(5)+BVt(6)+BVt(7))
      else
          FACT=1.D0/dSQRT(GX*GX+GY*GY+GZ*GZ)
          GGX=FACT*(2.D0*BVT(5)*GX+BVT(10)*GY+BVT(9)*GZ)
          GGY=FACT*(BVT(10)*GX+2.D0*BVT(6)*GY+BVT(8)*GZ)
          GGZ=FACT*(BVT(9)*GX+BVT(8)*GY+2.D0*BVT(7)*GZ)
          GDGGT=GX*GGX+GY*GGY+GZ*GGZ
      end if
      GXD=BV2(2)
      GYD=BV2(3)
      GZD=BV2(4)
      GMAGD=dSQRT(GXD*GXD+GYD*GYD+GZD*GZD)
      G2D=2.D0*(BV2(5)+BV2(6)+BV2(7))
      if(abs(gmagd).le.1.d-10) then
          ggxd =2.d0*bv2(5)
          ggyd =2.d0*bv2(6)
          ggzd =2.d0*bv2(7)
          gdgg=g2d
      else
          FACT=1.D0/GMAGD
          GGXD=FACT*(2.D0*BV2(5)*GXD+BV2(10)*GYD+BV2(9)*GZD)
          GGYD=FACT*(BV2(10)*GXD+2.D0*BV2(6)*GYD+BV2(8)*GZD)
          GGZD=FACT*(BV2(9)*GXD+BV2(8)*GYD+2.D0*BV2(7)*GZD)
      endif
!      gxu=0.
!      gyu=0.
!      gzu=0.
!      gmagu=0.
!      g2u=0.
!      ggxu=0.
!      ggyu=0.
!      ggzu=0.
!      gxd=0.
!      gyd=0.
!      gzd=0.
!      gmagd=0.
!      g2d=0.
!      ggxd=0.
!      ggyd=0.
!      ggzd=0.
      DUP=fft(i1,i2,i3,1)
      DDN=fft(i1,i2,i3,2)
          if(i1.eq.1.and.i2.eq.1.and.i3.eq.1) then
      write(6,8886)i3,dup,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU
      write(6,8886)i3,ddn,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD
          end if
 8886     format (1(1x,i3),4(3x,f10.6,1x,f10.6))
      CALL VXCLM2(dup,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
                  ddn,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD, &
                  vxu,vxd,exu,exd,GDGGT,tauup, &
                  taudn,mod(iex,50),i3+i2*1000+i1*1000000)
      if(i1.eq.1.and.i2.eq.1.and.i3.eq.1) write(6,*)vxu,vxd,gdggt
!
       RETURN
       END

