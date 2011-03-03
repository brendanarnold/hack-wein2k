      SUBROUTINE XCPOT1 (LM,JATOM,LMMAX,LR2V,iex,JSPIN,sigmaj,imfield,bextene)
!!!      use efeld
      use struct
      use densit
      use work
      use vresp1
!                                                                       
!     XCPOT1 CALCULATES THE EXCHANGE-CORRELATION POTENTIAL INSIDE       
!     THE SPHERES.                                                      
!     AFTER SUMATION OF  CHARGE DENSITY  BY CLM(R) FOR TEN  THETA       
!     VALUES PER MESH-POINT, THE ( CUBIC ) HARMONIC-EXPANSION  OF       
!     EXCHANGE-POTENTIAL IS GENERATED AND STORED IN ARRAY VEX.
!                                                                       
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
!

      REAL*8         :: sigmaj(2)
      COMPLEX*16           YL,dyt,dyt2,dyp,dyp2,dytp         
      LOGICAL           LR2V, robfit
!
      DIMENSION         YY(nspd,NCOM)                                    
      DIMENSION         bbu(NCOM),bbbu(NCOM),bbd(NCOM),          &
                        bbbd(NCOM) ! , AA2(NCOM,NCOM)
!      DIMENSION         spt_vxcu(nspd),spt_excu(nspd),spt_vxcd(nspd), &
!                        spt_excd(nspd)
!                                                                       
      DIMENSION  YL((lmax2+1)*(lmax2+1))
      DIMENSION  LM(2,NCOM+3,NAT),LMMAX(NAT)         
      DIMENSION  A2(3),ROTA1(3)                                          
!
!
      dimension spt(3,nspd),weight(nspd), &
                vxcsu(nspd),vxcsd(nspd), &
                vxu(nspd),vxd(nspd),exu(nspd),exd(nspd), &
                dcsp(nrad,ncom,2),ddcsp(nrad,ncom,2), &
                dyt((lmax2+1)*(lmax2+1)),dyt2((lmax2+1)*(lmax2+1)), &
                dyp((lmax2+1)*(lmax2+1)),dyp2((lmax2+1)*(lmax2+1)), &
                dytp((lmax2+1)*(lmax2+1)), &
                dydt(nspd,ncom),dydt2(nspd,ncom), &
                dydp(nspd,ncom),dydp2(nspd,ncom),dydtp(nspd,ncom), &
                qrhou(ncom),qrhod(ncom),qyy(ncom), &
                qsin(nspd),qcos(nspd)
      COMMON /phil/ dcsp,ddcsp

!      DATA THRD/0.3333333333333333D0/
       parameter (thrd=1.D0/3.D0 )
!                                                                       
!---------------------------------------------------------------------  
!.....DEFINE CONSTANTS, INITIALIZE VARIABLES 
!c       iex1=iex
!c       iex=1
      IF(.NOT.ALLOCATED(vclmsp)) allocate(VCLMSP(1:NRAD,1:NCOM,1:NAT,1:2))
         pi=acos(-1.d0)
         rtpi=sqrt(4.D0*pi)
         pisq=pi*pi
         wnorm=4.D0*pi
       1980 FORMAT(3X)                                                        
        iatest=-1
     INDEX=1
     DO J=1,JATOM-1
       INDEX=INDEX+MULT(J)
     ENDDO

 2000 FORMAT(16X,I2//)                                                  
 2010 FORMAT(16X,I2,5X,I2/)                                             
 2020 FORMAT(3X,5E13.7)                                                 
 2021 FORMAT(3X,4E19.12)                                                 
 2030 FORMAT(///)                                                       
 2031 FORMAT(/)                                                         
!                                                                       
!                                                                       

      NCOMX=NCOM                                                        
      PI=ACOS(-1.D0)                                                    
      TPI=2.D0*PI                                                       
      SQFP=SQRT(4.D0*PI)                                                
      PIST=PI/36.D0                                                     
      PIN=PI/9.D0-1.D-10
      PI4=PI*0.25D0                                                  
!                                                                       
      LLMM=LMMAX(JATOM)                                                 
!************************************************
!
!     generates Gauss-Legendre points on sphere
!
      LUSE=0
      do L=1,LLMM
        LUSE=MAX(LUSE,abs(lm(1,L,JATOM)))
      enddo                                                                       
      call gpoint(spt,weight,npt1,LUSE)
!-----------
      do 1 ind=1,npt1
         a2(1)=spt(1,ind)
         a2(2)=spt(2,ind)
         a2(3)=spt(3,ind)
         CALL ROTATE (A2,ROTLOC(1,1,JATOM),ROTA1)                        
         CALL YLM (ROTA1,lmax2,YL)                                          
!
         if(iex.gt.10) then
         call dylm(rota1,lmax2,yl,dyt,dyt2,dyp,dyp2,dytp,qsink,qcosk)
         qsin(ind)=qsink
         qcos(ind)=qcosk
         endif
!
         LM3=LMMAX(JATOM)                                                  
         DO 3 LM1=1,LM3                                                    
          CALL SUML(LM1,YK,YL,LM(1,1,JATOM))                             
          YY(IND,LM1)=YK                                                 
 3       CONTINUE                                                          
!
         if(iex.gt.10) then
         do 1001 lm1=1,lm3
          call suml(lm1,yk,dyt,lm(1,1,jatom))
          dydt(ind,lm1)=yk
          call suml(lm1,yk,dyt2,lm(1,1,jatom))
          dydt2(ind,lm1)=yk
          call suml(lm1,yk,dyp,lm(1,1,jatom))
          dydp(ind,lm1)=yk
          call suml(lm1,yk,dyp2,lm(1,1,jatom))
          dydp2(ind,lm1)=yk
          call suml(lm1,yk,dytp,lm(1,1,jatom))
          dydtp(ind,lm1)=yk
 1001    continue
         endif       
!        
         IF(IATNR(JATOM).GT.0.AND.LMMAX(JATOM).GT.1) THEN                  
!         LLMM=LMMAX(JATOM)-2
         CALL SUMFAC(llmm,YY,IND,LM(1,1,jatom),LMMAX(JATOM))
!
         if(iex.gt.10) then
         call sumfac(llmm,dydt,ind,LM(1,1,jatom),LMMAX(JATOM))
         call sumfac(llmm,dydt2,ind,LM(1,1,jatom),LMMAX(JATOM))
         call sumfac(llmm,dydp,ind,LM(1,1,jatom),LMMAX(JATOM))
         call sumfac(llmm,dydp2,ind,LM(1,1,jatom),LMMAX(JATOM))
         call sumfac(llmm,dydtp,ind,LM(1,1,jatom),LMMAX(JATOM))
         endif
!        
         ENDIF                                                             
!                                                                       
  1   CONTINUE
!                                                                       
!                                                                       
      WRITE(6,1530) JATOM,ANAME(JATOM),1.D0 !,RCOND
!                            
!========================================================
!         derivatives dCLM/dR and ddCLM/ddR
!========================================================
      if(iex.gt.10) then
        call drho(jatom,jri(jatom),lmmax(jatom), &
                              clmsp,dcsp,ddcsp)
      end if
!========================================================
!
!     LOOP OVER ALL RADIAL MESH POINTS                                  
!
      init=1
      rp=r(1)
!
      DO 34 ir=1,JRI(JATOM)                                               
!
       do 10 lm1=1,llmm
        bbu(lm1)=0.0D0
        bbd(lm1)=0.0D0
        bbbu(lm1)=0.0D0
        bbbd(lm1)=0.0D0
 10    continue
!
       DO 77 k=1,npt1
!       spt_vxcu(k)=0
!       spt_vxcd(k)=0
!       spt_excu(k)=0
!       spt_excd(k)=0
!
        do 1002 lm1=1,lmmax(jatom)
         qrhou(lm1)=clmsp(ir,lm1,jatom,1)
         qrhod(lm1)=clmsp(ir,lm1,jatom,2)
         qyy(lm1)=yy(k,lm1)
 1002   continue
        call srolyl(qrhou,qrhod,qyy,sqfp,iatnr(jatom),llmm, &
                    fu,fd,ft,LM(1,1,jatom),LMMAX(JATOM))
        RIR=1.D0/(r(ir)*r(ir))
        ft=ft*RIR
        fu=fu*RIR
        fd=fd*RIR
        if(fu.lt.1.d-8) then
          write(6,*)' xcpot1: ir,k',ir,k,' rhoup negative:',fu
          fu=1.d-8
          ft=fu+fd
        end if
        if(fd.lt.1.d-8) then
          write(6,*)' xcpot1: ir,k',ir,k,' rhodn negative:',fd
          fd=1.d-8
          ft=fu+fd
        end if
!
!.......iex 41: calc part of tau: rhoe and Vcoul_rho
        if(iex.eq.41.or.iex.eq.42.or.iex.eq.43.or.iex.eq.44 &
                    .or.iex.eq.12.or.iex.eq.27) then
        do  lm1=1,lmmax(jatom)
         qrhou(lm1)=vclmsp(ir,lm1,jatom,1)
         qrhod(lm1)=vclmsp(ir,lm1,jatom,2)
         qyy(lm1)=yy(k,lm1)
        enddo
        call srolyl(qrhou,qrhod,qyy,sqfp,iatnr(jatom),llmm, &
                    tauup,taudn,tautot,LM(1,1,jatom),LMMAX(JATOM))
        tautot=-tautot*RIR
        tauup=-tauup*RIR
        taudn=-taudn*RIR
        do  lm1=1,lmmax(jatom)
         qrhou(lm1)=v(ir,lm1,jatom)
         qrhod(lm1)=v(ir,lm1,jatom)
         qyy(lm1)=yy(k,lm1)
        enddo
         qrhou(1)=qrhou(1)*RTPI
         qrhod(1)=qrhod(1)*RTPI
        call srolyl(qrhou,qrhod,qyy,sqfp,iatnr(jatom),llmm, &
                    vhup,vhdn,vhtot,LM(1,1,jatom),LMMAX(JATOM))
!.......v was not multiplied with r**2
        tautot=tautot-vhtot*ft
        tauup=tauup -vhup*fu
        taudn=taudn -vhdn*fd
        endif
!
!--------------------------------
        if(iex.gt.10) then
!
        do 1003 lm1=1,lmmax(jatom)
         qrhou(lm1)=dcsp(ir,lm1,1)
         qrhod(lm1)=dcsp(ir,lm1,2)
 1003   continue
        call srolyl(qrhou,qrhod,qyy,sqfp,iatnr(jatom),llmm, &
                    dndru,dndrd,dndrt,LM(1,1,jatom),LMMAX(JATOM))
!
        do 1004 lm1=1,lmmax(jatom)
         qrhou(lm1)=ddcsp(ir,lm1,1)
         qrhod(lm1)=ddcsp(ir,lm1,2)
 1004   continue
        call srolyl(qrhou,qrhod,qyy,sqfp,iatnr(jatom),llmm, &
                    d2ndru,d2ndrd,d2ndrt,LM(1,1,jatom),LMMAX(JATOM))
!
        rpp=1.D0/(rp*rp)
        do 1005 lm1=1,lmmax(jatom)
         qrhou(lm1)=clmsp(ir,lm1,jatom,1)*rpp
         qrhod(lm1)=clmsp(ir,lm1,jatom,2)*rpp
         qyy(lm1)=dydt(k,lm1)
         if(iex.gt.60) qyy(lm1)=0.0d0
 1005   continue
        call srolyl(qrhou,qrhod,qyy,sqfp,iatnr(jatom),llmm, &
                    dndtu,dndtd,dndtt,LM(1,1,jatom),LMMAX(JATOM))
!
        do 1006 lm1=1,lmmax(jatom)
         qyy(lm1)=dydt2(k,lm1)
         if(iex.gt.60) qyy(lm1)=0.0d0
 1006   continue
        call srolyl(qrhou,qrhod,qyy,sqfp,iatnr(jatom),llmm, &
                    d2ndtu,d2ndtd,d2ndtt,LM(1,1,jatom),LMMAX(JATOM))
!
        do 1007 lm1=1,lmmax(jatom)
         qyy(lm1)=dydp(k,lm1)
         if(iex.gt.60) qyy(lm1)=0.0d0
 1007   continue
        call srolyl(qrhou,qrhod,qyy,sqfp,iatnr(jatom),llmm, &
                    dndpu,dndpd,dndpt,LM(1,1,jatom),LMMAX(JATOM))
!
        do 1008 lm1=1,lmmax(jatom)
         qyy(lm1)=dydp2(k,lm1)
         if(iex.gt.60) qyy(lm1)=0.0d0
 1008   continue
        call srolyl(qrhou,qrhod,qyy,sqfp,iatnr(jatom),llmm, &
                    d2ndpu,d2ndpd,d2ndpt,LM(1,1,jatom),LMMAX(JATOM))
!
        do 1009 lm1=1,lmmax(jatom)
         qyy(lm1)=dydtp(k,lm1)
         if(iex.gt.60) qyy(lm1)=0.0d0
 1009   continue
        call srolyl(qrhou,qrhod,qyy,sqfp,iatnr(jatom),llmm, &
                    dndtpu,dndtpd,dndtpt,LM(1,1,jatom),LMMAX(JATOM))
!
        do 1010 lm1=1,lmmax(jatom)
         qrhou(lm1)=dcsp(ir,lm1,1)
         qrhod(lm1)=dcsp(ir,lm1,2)
         qyy(lm1)=dydt(k,lm1)
         if(iex.gt.60) qyy(lm1)=0.0d0
 1010   continue
        call srolyl(qrhou,qrhod,qyy,sqfp,iatnr(jatom),llmm, &
                    dndrtu,dndrtd,dndrtt,LM(1,1,jatom),LMMAX(JATOM))
!
        do 1011 lm1=1,lmmax(jatom)
         qyy(lm1)=dydp(k,lm1)
         if(iex.gt.60) qyy(lm1)=0.0d0
 1011   continue
        call srolyl(qrhou,qrhod,qyy,sqfp,iatnr(jatom),llmm, &
                    dndrpu,dndrpd,dndrpt,LM(1,1,jatom),LMMAX(JATOM))
!
!--------------------------------         
        endif                                                          
!_vresp
        if (iex.eq.45) then
           do 1022 lm1=1,lmmax(jatom)
              qrhou(lm1)=vclmsp(ir,lm1,jatom,1)
              qrhod(lm1)=vclmsp(ir,lm1,jatom,2)
              qyy(lm1)=yy(k,lm1)
 1022      continue
           call srolyl(qrhou,qrhod,qyy,sqfp,iatnr(jatom),llmm, &
                       vfu,vfd,vft,LM(1,1,jatom),LMMAX(JATOM))
           vfu=vfu*rir
           vfd=vfd*rir
        endif
!_vresp
!                                                                       
!******************************************************
       if(iex.lt.10) then
!
           vxcu=xcpot (ft,fu,pi,iex)
           vxcd=xcpot (ft,fd,pi,iex)
           excu=xcener(ft,fu,pi,iex)
           excd=xcener(ft,fd,pi,iex)
       if(iex.eq.5) then
           RS=(0.75D0/(PI*ft))**THRD
           ZET=(FU-FD)/ft
           CALL CORLSD(RS,ZET,ECLSD,VLSDU,VLSDD,ECRS,ECZET,ALFC)
           vxcu=vxcu+2.d0*VLSDU
           vxcd=vxcd+2.d0*VLSDD
           excu=excu+2.d0*eclsd
           excd=excd+2.d0*eclsd
        endif
!bext
        if(imfield.gt.0) then
         vxcu=vxcu - Bextene
         vxcd=vxcd + Bextene
        endif            
!bext                                  
!!!        if(abs(iefeld).gt.0) then
!!!        z=spt(3,k)*r(ir)+a(3)*pos(3,index)
!!!!.......E-field only for  mult=1 and z-direction (hex or tetrag)
!!!         if(iefeld.lt.0) then
!!!!          set potential to zero if iefeld negative
!!!           vxcu=0.0d0
!!!           vxcd=0.0d0
!!!         endif
!!!!.......odd i
!!!        if(abs(iefeld).lt.999) then
!!!         iefeld1=abs(iefeld)
!!!         do i=1,iefeld1-1,2
!!!          vxcu=vxcu+refeld*8.d0/pisq/i**2*cos(i*z*tpi/a(3))
!!!          vxcd=vxcd+refeld*8.d0/pisq/i**2*cos(i*z*tpi/a(3))
!!!         enddo
!!!!         if(k.lt.4) write(6,*) ir,z,spt(3,k),vxcu
!!!        else
!!!!  efeld 7/8:1/8 geteilt
!!!         iefeld1=mod(abs(iefeld),1000)
!!!         do i=-iefeld1,iefeld1
!!!         if(i.eq.0) cycle
!!!          vxcu=vxcu+refeld*((8.d0*cos(i*pi4)+i*pi*sin(i*pi4)-8.d0) &
!!!                      /2.d0/pisq/i**2 *cos(i*z*tpi/a(3))+ &
!!!               (8.d0*cos(7.d0*pi4*i)+7.d0*i*pi*sin(7.d0*pi4*i)-8.d0) &
!!!                    /14.d0/i**2/pisq *cos(i*z*tpi/a(3)))  +  &
!!!               refeld*((-8.d0*sin(pi4*i)+i*pi*cos(pi4*i)+i*pi) &
!!!                      /2.d0/i**2/pisq *sin(i*z*tpi/a(3)) &
!!!                    -(-8.d0*sin(7.d0*pi4*i)+7.d0*i*pi*cos(7.d0*pi4*i)&
!!!                    +7.d0*i*pi)/14.d0/i**2/pisq *sin(i*z*tpi/a(3)))
!!!          vxcd=vxcd+refeld*((8.d0*cos(i*pi4)+i*pi*sin(i*pi4)-8.d0) &
!!!                      /2.d0/pisq/i**2 *cos(i*z*tpi/a(3))+ &
!!!               (8.d0*cos(7.d0*pi4*i)+7.d0*i*pi*sin(7.d0*pi4*i)-8.d0) &
!!!                    /14.d0/i**2/pisq *cos(i*z*tpi/a(3)))  +  &
!!!               refeld*((-8.d0*sin(pi4*i)+i*pi*cos(pi4*i)+i*pi) &
!!!                      /2.d0/i**2/pisq *sin(i*z*tpi/a(3)) &
!!!                    -(-8.d0*sin(7.d0*pi4*i)+7.d0*i*pi*cos(7.d0*pi4*i)&
!!!                    +7.d0*i*pi)/14.d0/i**2/pisq *sin(i*z*tpi/a(3)))
!!!         enddo
!!!       endif
!!!        endif
!
        if(ir.eq.jri(jatom)) then
         vxcsu(k)=vxcu
         vxcsd(k)=vxcd
        endif
!                     
        tmp1=vxcu*weight(k)
        tmp2=vxcd*weight(k)
        tmp3=excu*weight(k)
        tmp4=excd*weight(k)
        DO 9 LM1=1,LLMM                                                   
         bbu (lm1)=bbu (lm1)+yy(k,lm1)*tmp1
         bbd (lm1)=bbd (lm1)+yy(k,lm1)*tmp2
         bbbu(lm1)=bbbu(lm1)+yy(k,lm1)*tmp3
         bbbd(lm1)=bbbd(lm1)+yy(k,lm1)*tmp4
 9      continue
       endif
!-----------------------
       if(iex.gt.10) then
!
        call grans(rp,qsin(k),qcos(k),dndrt,dndtt,dndpt, &
                   d2ndrt,d2ndtt,d2ndpt,dndrtt,dndrpt,dndtpt, &
                   gx,gy,gz,gmag,g2,ggx,ggy,ggz,gdggt,ir)
        call grans(rp,qsin(k),qcos(k),dndru,dndtu,dndpu, &
                   d2ndru,d2ndtu,d2ndpu,dndrtu,dndrpu,dndtpu, &
                   gxu,gyu,gzu,gmagu,g2u,ggxu,ggyu,ggzu,gdggu,ir)
        call grans(rp,qsin(k),qcos(k),dndrd,dndtd,dndpd, &
                   d2ndrd,d2ndtd,d2ndpd,dndrtd,dndrpd,dndtpd, &
                   gxd,gyd,gzd,gmagd,g2d,ggxd,ggyd,ggzd,gdggd,ir)
!_vresp
        if (iex.eq.45) then
           vxu(k)=vfu
           vxd(k)=vfd
        endif
!_vresp
        CALL VXCLM2(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
                    FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD, &
                    VXU(K),VXD(K),exu(k),exd(k),GDGGT,tauup,taudn, &
                    mod(iex,50),ir)
!
        vxcu=vxu(k)
        vxcd=vxd(k)
        excu=exu(k)
        excd=exd(k)
!           if(k.eq.1) write(6,*) r(ir),excu,excd
!
        if(imfield.gt.0) then
         vxcu=vxcu - Bextene
         vxcd=vxcd + Bextene
        endif   
!
!!!        if(abs(iefeld).gt.0) then
!!!        z=spt(3,k)*r(ir)+a(3)*pos(3,index)
!!!!.......E-field only for  mult=1 and z-direction (hex or tetrag)
!!!         if(iefeld.lt.0) then
!!!!set potential to zero if iefeld negative, should calculate just efield without VXC
!!!           vxcu=0.0d0
!!!           vxcd=0.0d0
!!!         endif
!!!        if(abs(iefeld).lt.999) then
!!!         iefeld1=abs(iefeld)
!!!!.......odd i
!!!         do i=1,iefeld1-1,2
!!!          vxcu=vxcu+refeld*8.d0/pisq/i**2*cos(i*z*tpi/a(3))
!!!          vxcd=vxcd+refeld*8.d0/pisq/i**2*cos(i*z*tpi/a(3))
!!!         enddo
!!!!         if(k.lt.4) write(*,*) ir,z,spt(3,k),vxcu
!!!        else
!!!!  efeld 7/8:1/8 geteilt
!!!         iefeld1=mod(abs(iefeld),1000)
!!!         do i=-iefeld1,iefeld1
!!!         if(i.eq.0) cycle
!!!          vxcu=vxcu+refeld*((8.d0*cos(i*pi4)+i*pi*sin(i*pi4)-8.d0) &
!!!                      /2.d0/pisq/i**2 *cos(i*z*tpi/a(3))+ &
!!!               (8.d0*cos(7.d0*pi4*i)+7.d0*i*pi*sin(7.d0*pi4*i)-8.d0) &
!!!                    /14.d0/i**2/pisq *cos(i*z*tpi/a(3)))  +  &
!!!               refeld*((-8.d0*sin(pi4*i)+i*pi*cos(pi4*i)+i*pi) &
!!!                      /2.d0/i**2/pisq *sin(i*z*tpi/a(3)) &
!!!                    -(-8.d0*sin(7.d0*pi4*i)+7.d0*i*pi*cos(7.d0*pi4*i)&
!!!                    +7.d0*i*pi)/14.d0/i**2/pisq *sin(i*z*tpi/a(3)))
!!!          vxcd=vxcd+refeld*((8.d0*cos(i*pi4)+i*pi*sin(i*pi4)-8.d0) &
!!!                      /2.d0/pisq/i**2 *cos(i*z*tpi/a(3))+ &
!!!               (8.d0*cos(7.d0*pi4*i)+7.d0*i*pi*sin(7.d0*pi4*i)-8.d0) &
!!!                    /14.d0/i**2/pisq *cos(i*z*tpi/a(3)))  +  &
!!!               refeld*((-8.d0*sin(pi4*i)+i*pi*cos(pi4*i)+i*pi) &
!!!                      /2.d0/i**2/pisq *sin(i*z*tpi/a(3)) &
!!!                    -(-8.d0*sin(7.d0*pi4*i)+7.d0*i*pi*cos(7.d0*pi4*i)&
!!!                    +7.d0*i*pi)/14.d0/i**2/pisq *sin(i*z*tpi/a(3)))
!!!         enddo
!!!        endif
!!!        endif
!
        if(ir.eq.jri(jatom)) then
         vxcsu(k)=vxcu
         vxcsd(k)=vxcd
        endif
!                                                                       
        tmp1=vxcu*weight(k)
        tmp2=vxcd*weight(k)
        tmp3=excu*weight(k)
        tmp4=excd*weight(k)                                                  
        DO 3012 LM1=1,LLMM                                                   
         bbu (lm1)=bbu (lm1)+yy(k,lm1)*tmp1
         bbd (lm1)=bbd (lm1)+yy(k,lm1)*tmp2
         bbbu(lm1)=bbbu(lm1)+yy(k,lm1)*tmp3
         bbbd(lm1)=bbbd(lm1)+yy(k,lm1)*tmp4
 3012   continue

       endif
!
 77    CONTINUE
!
!================= VXC, UP ==================
!
      IF(IATNR(JATOM).LT.0) THEN                                        
      IF(ir.EQ.JRI(JATOM)) THEN                                          
          WRITE(6,1540)                                                 
          SIGMA=0.D0                                                      
          DO 12 k=1,npt1                                                
          TEST=0.D0                                                       
          DO 13 LM1=1,LLMM                                              
  13      TEST=TEST+bbu(LM1)*YY(k,LM1)                                   
          IF(MOD(k,2).EQ.0) THEN                                       
            DIFVXC=vxcsu(k)-TEST                                           
            WRITE(6,1550)k,vxcsu(k),TEST,DIFVXC                            
          ENDIF                                                         
  12      SIGMA=SIGMA+weight(k)*(vxcsu(k)-TEST)**2                                 
          SIGMA=SQRT(SIGMA/wnorm)                                        
          WRITE(6,*) 'ATOM',JATOM,'  AT RMT: SIGMA OF V-XC FIT: '        &
         ,SIGMA                                                         
          sigmaj(1)=sigma
          END IF                                                         
          DO 11 LM1=1,LLMM                                                 
   11     VXC(ir,LM1,JATOM,1)=bbu(LM1)                                   
      ELSE                                                              
       IF(ir.EQ.JRI(JATOM)) THEN                                          
          WRITE(6,*) 'k   VXC    FIT'                                   
          SIGMA=0.D0                                                      
          DO 112 k=1,npt1                                               
          TEST=0.D0                                                       
          DO 113 LM1=1,LLMM                                             
 113      TEST=TEST+bbu(LM1)*YY(k,LM1)                                   
          IF(MOD(k,2).EQ.0)WRITE(6,*) k,vxcsu(k),TEST                   
 112      SIGMA=SIGMA+weight(k)*(vxcsu(k)-TEST)**2                                
          SIGMA=SQRT(SIGMA/wnorm)                                        
          WRITE(6,*) 'ATOM ',JATOM,'  AT RMT: SIGMA OF V-XC FIT: '       &
         ,SIGMA                                                         
          sigmaj(1)=sigma
       END IF                                                         
       CALL cub_xc_back(vxc(1,1,jatom,1),ir,llmm, &
                        lm(1,1,jatom),bbu,lmmax(jatom))
      ENDIF                                                         
!
!================= VXC, DN ==================
!
      IF(IATNR(JATOM).LT.0) THEN                                        
      IF(ir.EQ.JRI(JATOM)) THEN                                          
          WRITE(6,1540)                                                 
          SIGMA=0.D0                                                      
          DO 5012 k=1,npt1                                           
          TEST=0.D0                                                       
          DO 5013 LM1=1,LLMM                                              
 5013     TEST=TEST+bbd(LM1)*YY(k,LM1)                                   
          IF(MOD(k,2).EQ.0) THEN                                       
            DIFVXC=vxcsd(k)-TEST                                           
            WRITE(6,1550)k,vxcsd(k),TEST,DIFVXC                            
          ENDIF                                                         
 5012     SIGMA=SIGMA+weight(k)*(vxcsd(k)-TEST)**2                                 
          SIGMA=SQRT(SIGMA/wnorm)                                        
          WRITE(6,*) 'ATOM',JATOM,'  AT RMT: SIGMA OF V-XC FIT: '        &
         ,SIGMA                                                         
          sigmaj(2)=sigma

         END IF                                                         
          DO 5011 LM1=1,LLMM                                                 
 5011     VXC(ir,LM1,JATOM,2)=bbd(LM1)                                   
      ELSE                                                              
       IF(ir.EQ.JRI(JATOM)) THEN                                          
          WRITE(6,*) 'k   VXC    FIT'                                   
          SIGMA=0.D0                                                      
          DO 5112 k=1,npt1                                          
          TEST=0.D0                                                       
          DO 5113 LM1=1,LLMM                                             
 5113     TEST=TEST+bbd(LM1)*YY(k,LM1)                                   
          IF(MOD(k,2).EQ.0)WRITE(6,*) k,vxcsd(k),TEST                   
 5112     SIGMA=SIGMA+weight(k)*(vxcsd(k)-TEST)**2                                
          SIGMA=SQRT(SIGMA/wnorm)                                        
          WRITE(6,*) 'ATOM ',JATOM,'  AT RMT: SIGMA OF V-XC FIT: '       &
         ,SIGMA                                                         
          sigmaj(2)=sigma
       END IF                                                         
         CALL cub_xc_back(vxc(1,1,jatom,2),ir,llmm, &
                        lm(1,1,jatom),bbd,lmmax(jatom))
      ENDIF                                                         
!                                                                       
!================= EXC, UP ==================
!
      IF(IATNR(JATOM).LT.0) THEN                                        
       DO 19 LM1=1,LLMM                                                  
  19   EXC(ir,LM1,JATOM,1)=bbbu(LM1)                                   
      ELSE                                                              
         CALL cub_xc_back(exc(1,1,jatom,1),ir,llmm, &
                        lm(1,1,jatom),bbbu,lmmax(jatom))
      ENDIF                                                         
!                                                                       
!================= EXC, DN ==================
!
      IF(IATNR(JATOM).LT.0) THEN
       DO 5019 LM1=1,LLMM                                                  
 5019  EXC(ir,LM1,JATOM,2)=bbbd(LM1)                                   
      ELSE                                                              
         CALL cub_xc_back(exc(1,1,jatom,2),ir,llmm, &
                        lm(1,1,jatom),bbbd,lmmax(jatom))
      ENDIF                                                         
!                                                                       
      rp=rp*exp(dx(jatom))
 34   CONTINUE                                                          
!********************************************************
!      LOOP OVER ALL RADIAL MESH COMPLETED
!                                                                       

 1530 FORMAT(3X,'CONDITION PARAMETER ATOM',I3,2X,A10,2X,'RCOND= ',       &
             F8.5)                                                      
 1540 FORMAT(3X,'INDEX',6X,'V-XC',7X,'FIT',5X,'DIFFERENCE')             
 1550 FORMAT(I8,2F10.5,F15.5)                                           
!c      iex=iex1
      RETURN                                                            
      END                                                               
