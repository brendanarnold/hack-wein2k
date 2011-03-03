!
      SUBROUTINE VXCPW2(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
                FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD,VXU,VXD, &
                exu,exd,GDGGT)
!
      IMPLICIT REAL*8(A-H,K,O-Z)
      data mist/0/
!
      DATA PI,TPI2/3.1415926535897931D0,29.60881320326807D0/
      DATA B,C,XM/14.D0,0.2D0,0.06666666666666667D0/
      DATA FTDS,THRD/1.3333333333333333D0,0.33333333333333333D0/
      DATA CAGU,CAB1U,CAB2U,CAAU,CABU,CACU,CADU &
       /-.1423D0,1.0529D0,.3334D0,.0311D0,-.048D0,.0020D0,-0.0116D0/
      DATA CAGP,CAB1P,CAB2P,CAAP,CABP,CACP,CADP &
       /-.0843D0,1.3981D0,.2611D0,.01555D0,-.0269D0,.0007D0,-0.0048D0/
      FGGA(S)=(1.D0+S*S*(0.0864D0/XM+S*S*(B+C*S*S)))**XM
      SIDF(S)=XM*(0.1728D0/XM+S*S*(4.D0*B+6.D0*C*S*S)) &
                   /(1.D0+S*S*(0.0864D0/XM+S*S*(B+C*S*S)))**(1.D0-XM)
      DSIDF(S)=XM*(XM-1.D0)*S*(0.1728D0/XM+S*S*(4.D0*B+6.D0*C*S*S))**2 &
                   /(1.D0+S*S*(0.0864D0/XM+S*S*(B+C*S*S)))**(2.D0-XM) &
           +  XM*(24.D0*C*S*S+8.D0*B)*S &
                   /(1.D0+S*S*(0.0864D0/XM+S*S*(B+C*S*S)))**(1.D0-XM)
      KF(RHO)=(TPI2*RHO)**THRD
      RS(RHO)=(0.75D0/(PI*RHO))**THRD
      PHI(GM,RHO,CFR)=8.1290825D-4*GM/(CFR*RHO**1.1666666666667D0)
      AX=-0.75D0*(3.D0/PI)**THRD
!
      IF((FU.LE.0.D0).OR.(FD.LE.0.D0))THEN
      VXU=0.D0
      VXD=0.D0
      WRITE( 6,'('' ---> NMTNEW: NEGATIVE DENSITY:'',2D11.3)')FU,FD
      WRITE(16,'('' ---> NMTNEW: NEGATIVE DENSITY:'',2D11.3)')FU,FD
      RETURN
      ENDIF
!
      itot=0
      FT=FU+FD
      GXT=GXU+GXD
      GYT=GYU+GYD
      GZT=GZU+GZD
      GMAGT=SQRT(GXT*GXT+GYT*GYT+GZT*GZT)
      if (gmagt.le.2.d-10) itot=1
      G2T=G2U+G2D
!     GGXT=GGXU+GGXD
!     GGYT=GGYU+GGYD
!     GGZT=GGZU+GGZD
!     GDGGT=GXT*GGXT+GYT*GGYT+GZT*GGZT
      GDGGU=GXU*GGXU+GYU*GGYU+GZU*GGZU
      GDGGD=GXD*GGXD+GYD*GGYD+GZD*GGZD
!
!     EXCHANGE
!
      KFU=KF(2.D0*FU)
      KFD=KF(2.D0*FD)
      SU=GMAGU/(2.D0*KFU*FU)
      SD=GMAGD/(2.D0*KFD*FD)
      EXCUX=AX*FGGA(SU)*(2.D0*FU)**THRD
      EXCDX=AX*FGGA(SD)*(2.D0*FD)**THRD
      TU=G2U/(FU*(2.D0*KFU)**2)
      TD=G2D/(FD*(2.D0*KFD)**2)
      UU=GDGGU/(FU**2*(2.D0*KFU)**3)
      UD=GDGGD/(FD**2*(2.D0*KFD)**3)
      VXCUX=AX*(2.D0*FU)**THRD*(FTDS*FGGA(SU)-TU*SIDF(SU)- &
                                       (UU-FTDS*SU**3)*DSIDF(SU))
      VXCDX=AX*(2.D0*FD)**THRD*(FTDS*FGGA(SD)-TD*SIDF(SD)- &
                                       (UD-FTDS*SD**3)*DSIDF(SD))
!
!     CEPERLY ALDER CORRELATION
!
      CHI=(FU-FD)/(FT)
      FCHI=((1.D0+CHI)**FTDS+(1.D0-CHI)**FTDS-2.D0)/(2.D0**FTDS-2.D0)
      DFCHI=FTDS*((1.D0+CHI)**THRD-(1.D0-CHI)**THRD)/(2.D0**FTDS-2.D0)
      RSTOT=RS(FT)
      IF(RSTOT.GT.1.D0)THEN
      DEN=1.D0/(1.D0+CAB1U*SQRT(RSTOT)+CAB2U*RSTOT)
      ECU=CAGU*DEN
      VCU=ECU*DEN*(1.D0+1.1666666666666666D0*CAB1U*SQRT(RSTOT)+ &
                                                  FTDS*CAB2U*RSTOT)
      DEN=1.D0/(1.D0+CAB1P*SQRT(RSTOT)+CAB2P*RSTOT)
      ECP=CAGP*DEN
      VCP=ECP*DEN*(1.D0+1.1666666666666666D0*CAB1P*SQRT(RSTOT)+ &
                                                  FTDS*CAB2P*RSTOT)
      ELSE
      DLRS=LOG(RSTOT)
      ECU=(CAAU+CACU*RSTOT)*DLRS+CABU+CADU*RSTOT
      VCU=(CAAU+0.66666666666666666D0*CACU*RSTOT)*DLRS+CABU+ &
           THRD*(RSTOT*(2.D0*CADU-CACU)-CAAU)
      ECP=(CAAP+CACP*RSTOT)*DLRS+CABP+CADP*RSTOT
      VCP=(CAAP+0.66666666666666666D0*CACP*RSTOT)*DLRS+CABP+ &
           THRD*(RSTOT*(2.D0*CADP-CACP)-CAAP)
      ENDIF
      EXCUC=(ECU+FCHI*(ECP-ECU))
      VXCUC=(VCU+FCHI*(VCP-VCU)+(ECP-ECU)*( 1.D0-CHI)*DFCHI)
      VXCDC=(VCU+FCHI*(VCP-VCU)+(ECP-ECU)*(-1.D0-CHI)*DFCHI)
!
!     GRADIENT CORRECTION TO CORRELATION
!
      FACT1=((1.D0+CHI)/2.D0)**(5.D0/3.D0)
      FACT2=((1.D0-CHI)/2.D0)**(5.D0/3.D0)
      DPARM=2.D0**THRD*SQRT(FACT1+FACT2)
      CFP=0.001667D0+(0.002568D0+RSTOT*(0.023266D0+7.389D-6*RSTOT))/ &
          (1.D0+RSTOT*(8.723D0+RSTOT*(0.472D0+7.389D-2*RSTOT)))
      DCFP=(-RSTOT/(3.D0*FT))*( (0.023266D0+1.4778D-5*RSTOT)* &
          (1.D0+RSTOT*(8.723D0+RSTOT*(0.472D0+7.389D-2*RSTOT))) &
        -(0.002568D0+RSTOT*(0.023266D0+7.389D-6*RSTOT))* &
         (8.723D0+RSTOT*(0.944D0+0.22167D0*RSTOT)) &
        )/(1.D0+RSTOT*(8.723D0+RSTOT*(0.472D0+7.389D-2*RSTOT)))**2
      PHIP=PHI(GMAGT,FT,CFP)
      EXCUC=EXCUC+EXP(-PHIP)*CFP*GMAGT*GMAGT/ &
                                       (DPARM*(FT)**(7.D0/3.D0))
      EXCDC=EXCUC
      PREF=-EXP(-PHIP)*CFP/(DPARM*FT**THRD)
!
      if(itot.eq.0) then
      TERM1=PREF*( (2.D0-PHIP)*G2T/FT &
       -(FTDS+PHIP*(-3.666666666666667D0+PHIP*1.166666666666667D0))* &
                                                         (GMAGT/FT)**2 &
        +PHIP*(PHIP-3.D0)*GDGGT/(FT*GMAGT) )
      else
      term1=pref*( phip*phip-4.d0*phip+2.d0)*g2t/ft
      endif
!
      PREF=-PREF*5.D0*FT**THRD/(6.D0*(FT**2*DPARM)**2)
      TF2=2.D0**(2.D0/3.D0)
      TF1=TF2*(1.D0-PHIP)
      TF2=TF2*(2.D0-PHIP)
      TERM1=TERM1+EXP(-PHIP)*(GMAGT**2)*(PHIP*(PHIP-1.D0)-1.D0)*DCFP/ &
                                                       (DPARM*FT**FTDS)
      TERMU=PREF* &
           (FU**(2.D0/3.D0)-FD**(2.D0/3.D0))*(TF1*FD*GMAGT**2-TF2*FT* &
           (GXD*GXT+GYD*GYT+GZD*GZT))
      TERMD=PREF* &
           (FD**(2.D0/3.D0)-FU**(2.D0/3.D0))*(TF1*FU*GMAGT**2-TF2*FT* &
           (GXU*GXT+GYU*GYT+GZU*GZT))
      VXCUC=VXCUC+(TERM1+TERMU)
      VXCDC=VXCDC+(TERM1+TERMD)
!
!     HARTREES TO RY
!
      VXU=2.D0*(VXCUX+VXCUC)
      VXD=2.D0*(VXCDX+VXCDC)
      EXU=2.D0*(EXCUX+EXCUC)
      EXD=2.D0*(EXCDX+EXCDC)
      if(mist.eq.0) then
      mist=mist+1
      write(6,*) vxu,vxcux,vxcuc,term1,termu
      end if
!
      RETURN
      END
