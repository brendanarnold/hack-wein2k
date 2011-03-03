      SUBROUTINE TH1(ROA,ROAX,ROAY,ROAZ,GROA,GDGGU,XLAPA, &
      ROB,ROBX,ROBY,ROBZ,GROB,GDGGD,XLAPB,GRO,GDGGT,EXC,VXCA,VXCB)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(CON43=4.D0/3.D0,CON83=8.D0/3.D0,CON13=1.D0/3.D0)
      PARAMETER(CON76=7.D0/6.D0,CON86=8.D0/6.D0,CON96=9.D0/6.D0, &
                CON106=10.D0/6.D0,CON116=11.D0/6.D0,CON126=2.D0)
      DIMENSION R(7:12),DRDA(7:12),DRDB(7:12),CONR(7:12)
      DIMENSION X(2),DXDA(2),DXDB(2)
      REAL CONST(21)
      INTEGER ICONR(21),ipri
      DATA CONST/-0.728255D0,0.331699D0,-0.102946D1,0.235703D0, &
                 -0.876221D-1,0.140854D0,0.336982D-1,-0.353615D-1, &
                 0.497930D-2,-0.645900D-1,0.461795D-1,-0.757191D-2, &
                 -0.242717D-2,0.428140D-1,-0.744891D-1,0.386577D-1, &
                 -0.352519D0,0.219805D1,-0.372927D1,0.194441D1, &
                 0.128877D0/
      DATA CONR/CON76,CON86,CON96,CON106,CON116,CON126/
      DATA ICONR/7,8,9,10,8,9,10,11,9,10,11,12,9,10,11,12,7,8,9,10,0/
      data ipri/0/
      save ipri
      ipri=0
      if(groa.lt.0.000001d0) then
        write(*,*) 'groa',groa
        ipri=ipri+1
      endif
      EXC=0.D0
      VXCA=0.D0
      VXCB=0.D0
      RO=ROA+ROB
      ROX=ROAX+ROBX
      ROY=ROAY+ROBY
      ROZ=ROAZ+ROBZ
      GRO2A=GROA**2
      GRO2B=GROB**2
      GFGFU=ROX*ROAX+ROY*ROAY+ROZ*ROAZ
      GFGFD=ROX*ROBX+ROY*ROBY+ROZ*ROBZ
      GFUGFD=ROAX*ROBX+ROAY*ROBY+ROAZ*ROBZ
      F43MU=2.D0*RO**CON43*GROA
      F43MD=2.D0*RO**CON43*GROB
!     R PRIME
      DO I=7,12
       R(I)=ROA**CONR(I)+ROB**CONR(I)
       DRDA(I)=CONR(I)*ROA**(CONR(I)-1.D0)
       DRDB(I)=CONR(I)*ROB**(CONR(I)-1.D0)
      ENDDO
!     S PRIME
      S=((ROA-ROB)/RO)**2
      DSDA=4.D0*(ROA-ROB)*ROB/RO**3
      DSDB=-4.D0*(ROA-ROB)*ROA/RO**3
!     X PRIME
      X(1)=(GROA+GROB)/(2.D0*RO**CON43)
      DXDA(1)=-CON43*X(1)/RO
      DXDB(1)=DXDA(1)
      X(2)=(GRO2A+GRO2B)/(2.D0*RO**CON83)
      DXDA(2)=-CON83*X(2)/RO
      DXDB(2)=DXDA(2)
!     Y PRIME
      Y=(GRO2A+GRO2B-2.D0*GFUGFD)/RO**CON83
      DYDA=-CON83*Y/RO
      DYDB=-CON83*Y/RO
      DO I=1,4
       JR=ICONR(I)
       EXC=EXC+CONST(I)*R(JR)
       VXCA=VXCA+CONST(I)*DRDA(JR)
       VXCB=VXCB+CONST(I)*DRDB(JR)
      if(ipri.eq.1) write(6,*) vxca,vxcb
      ENDDO
      DO I=5,8
       JR=ICONR(I)
       JX=1
       EXC=EXC+CONST(I)*R(JR)*X(JX)
      if(ipri.eq.0) then
       VXCA=VXCA+CONST(I)*(DRDA(JR)*X(JX)+R(JR)*DXDA(JX)- &
       (DRDA(JR)*GRO2A+DRDB(JR)*GFUGFD+R(JR)*(XLAPA- &
        CON43*GFGFU/RO-GDGGU/GROA))/F43MU)
       VXCB=VXCB+CONST(I)*(DRDB(JR)*X(JX)+R(JR)*DXDB(JX)- &
       (DRDB(JR)*GRO2B+DRDA(JR)*GFUGFD+R(JR)*(XLAPB- &
        CON43*GFGFD/RO-GDGGD/GROB))/F43MD)
      else
        VXCA=VXCA+CONST(I)*(DRDA(JR)*X(JX)+R(JR)*DXDA(JX))
        VXCB=VXCB+CONST(I)*(DRDB(JR)*X(JX)+R(JR)*DXDB(JX))
        write(6,*) vxca,vxcb,DRDA(JR),X(JX),R(JR),DXDA(JX), &
        DRDA(JR),GRO2A,DRDB(JR),GFUGFD,R(JR),XLAPA, &
        CON43,GFGFU,RO,GDGGU,GROA,F43MU
      endif
      ENDDO
      DO I=9,12
       JR=ICONR(I)
       JX=2
       EXC=EXC+CONST(I)*R(JR)*X(JX)
       VXCA=VXCA+CONST(I)*(DRDA(JR)*X(JX)+R(JR)*DXDA(JX)- &
            (DRDA(JR)*GRO2A+DRDB(JR)*GFUGFD+ &
            R(JR)*(XLAPA-CON83*GFGFU/RO))/RO**CON83)
       VXCB=VXCB+CONST(I)*(DRDB(JR)*X(JX)+R(JR)*DXDB(JX)- &
            (DRDB(JR)*GRO2B+DRDA(JR)*GFUGFD+ &
            R(JR)*(XLAPB-CON83*GFGFD/RO))/RO**CON83)
!     * (DRDB(JR)*GRO2B+DRDA(JR)*GFUGFD+R(JR)*(XLAPB-
!     *  CON43*GFGFD/RO-GDGGD/GROB))/F43MD)
! The line above can't be correct!!! It's from case 5-8
      if(ipri.eq.1) write(6,*) vxca,vxcb
      ENDDO
      DO I=9,12
       JR=ICONR(I)
       JX=2
       EXC=EXC+CONST(I)*R(JR)*X(JX)
       VXCA=VXCA+CONST(I)*(DRDA(JR)*X(JX)+R(JR)*DXDA(JX)- &
            (DRDA(JR)*GRO2A+DRDB(JR)*GFUGFD+ &
            R(JR)*(XLAPA-CON83*GFGFU/RO))/RO**CON83)
       VXCB=VXCB+CONST(I)*(DRDB(JR)*X(JX)+R(JR)*DXDB(JX)- &
            (DRDB(JR)*GRO2B+DRDA(JR)*GFUGFD+ &
            R(JR)*(XLAPB-CON83*GFGFD/RO))/RO**CON83)
      if(ipri.eq.1) write(6,*) vxca,vxcb
      ENDDO
      DO I=13,16
       JR=ICONR(I)
       EXC=EXC+CONST(I)*R(JR)*Y
       VXCA=VXCA+CONST(I)*(DRDA(JR)*Y+R(JR)*DYDA- &
       (2.D0*(DRDA(JR)*(GRO2A-GFUGFD)+DRDB(JR)*(GFUGFD-GRO2B))- &
        R(JR)*(XLAPA-XLAPB-CON83*(GFGFU-GFGFD)))/RO**CON83)
       VXCA=VXCA+CONST(I)*(DRDB(JR)*Y+R(JR)*DYDB- &
       (2.D0*(DRDB(JR)*(GRO2B-GFUGFD)+DRDA(JR)*(GFUGFD-GRO2A))- &
        R(JR)*(XLAPB-XLAPA-CON83*(GFGFD-GFGFU)))/RO**CON83)
      if(ipri.eq.1) write(6,*) vxca,vxcb
      ENDDO
      DO I=17,20
       JR=ICONR(I)
       EXC=EXC+CONST(I)*R(JR)*S
       VXCA=VXCA+CONST(I)*(DRDA(JR)*S+R(JR)*DSDA)
       VXCB=VXCB+CONST(I)*(DRDB(JR)*S+R(JR)*DSDB)
      if(ipri.eq.1) write(6,*) vxca,vxcb
      ENDDO
      EXC=EXC/RO+CONST(21)
      VXCA=VXCA+CONST(21)
      VXCB=VXCB+CONST(21)
      RETURN
      END
