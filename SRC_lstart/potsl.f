      SUBROUTINE POTSL(DV,D,DP,DQ,DR,DPAS,IEX,DEXV,Z,NP,ION,ICUT,JSPIN) 
!                                                                       
! POTENTIEL INTEGRATION PAR UNE METHODE A 4 POINTS                      
! DV POTENTIEL   D DENSITE   DP BLOC DE TRAVAIL  DR POINTS DE TABULATION
! DPAS PAS EXPONENTIEL   DEXV COEFFICIENT MULTIPLICATIF POUR L ECHANGE  
! Z NUMERO ATOMIQUE   NP NOMBRE DE POINTS   ION=Z-NOMBRE D ELECTRONS    
! SI ICUT EST NUL ON CORRIGE EVENTUELLEMENT LE POTENTIEL EN -(ION+1)/R  
! **********************************************************************
!      PARAMETER (NPT=881)                                               
      IMPLICIT REAL*8 (A-H,O-Z)
!     IMPLICIT REAL (D)                                                 
      INCLUDE 'param.inc'
      DIMENSION DV(NPT,2),D(NPT,2),DCSP(NPT,2),DDCSP(NPT,2)
      DIMENSION VLS(2)
      DIMENSION DP(*),DR(*),DQ(*)                    
      DATA THRD/0.3333333333333333D0/
      pi=acos(-1.d0)
      DAS=DPAS/24.                                                      
      DO 1 I=1,NP                                                       
      DQ(I)=D(I,1)+D(I,2)                                               
 1    DV(I,1)=DQ(I)*DR(I)                                               
      DLO=EXP(DPAS)                                                     
      DLO2=DLO*DLO                                                      
      DP(2)=DR(1)*(DQ(2)-DQ(1)*DLO2)/(12.*(DLO-1.))                     
      DP(1)=DV(1,1)/3.-DP(2)/DLO2                                       
      DP(2)=DV(2,1)/3.-DP(2)*DLO2                                       
      J=NP-1                                                            
      DO 3 I=3,J                                                        
 3    DP(I)=DP(I-1)+DAS*(13.*(DV(I,1)+DV(I-1,1))-(DV(I-2,1)+DV(I+1,1))) 
      DP(NP)=DP(J)                                                      
      DV(J,1)=DP(J)                                                     
      DV(NP,1)=DP(J)                                                    
      DO 5 I=3,J                                                        
      K=NP+1-I                                                          
 5    DV(K,1)=DV(K+1,1)/DLO+DAS*(13.*(DP(K+1)/DLO+DP(K))-(DP(K+2)/DLO2+  &
      DP(K-1)*DLO))                                                     
      DV(1,1)=DV(3,1)/DLO2+DPAS*(DP(1)+4.*DP(2)/DLO+DP(3)/DLO2)/3.      
      DLO=-(ION+1) 

      IF (IEX.LT.10) THEN
        DO 7 I=1,NP                                                       
          DV(I,2)=DV(I,1)                                                   
          DO 7 ISPIN=1,JSPIN                                                
!     IF(IEX.LT.2)DEXC=-3.*DEXV*((DR(I)*D(I)/105.27578)**(1./3.))       
!     IF(IEX.GT.1)DEXC=DR(I)*HL1(D(I),DR(I),IEX)                        
            FPRHO=DQ(I)/DR(I)/DR(I)                                           
            FPRHOS=D(I,ISPIN)/DR(I)/DR(I)                                     
            IF(IEX.EQ.5) THEN
	       RS=(3.d0/FPRHO)**THRD
	       ZET=(D(I,1)-D(I,2))/DQ(I)
               CALL CORLSD(RS,ZET,ECLSD,VLS(1),VLS(2),ECRS,ECZET,ALFC)
               DEXC=DR(I)*(VXCSP(FPRHO,FPRHOS,iex)+VLS(ISPIN)) 
	    ELSE
	       IF(IEX.GT.1)DEXC=DR(I)*VXCSP(FPRHO,FPRHOS,iex)
	    ENDIF
!     IF(I.LT.30) WRITE(6,*) I,DR(I),D(I,ISPIN),DEXC                    
            DV(I,ISPIN)=DV(I,ISPIN)-Z+DEXC                                    
            IF (ICUT.NE.0) GO TO 7                                            
            IF (DV(I,ISPIN).GT.DLO) DV(I,ISPIN)=DLO                           
 7      DV(I,ISPIN)=DV(I,ISPIN)/DR(I)
      ELSE
        call drho(np,d,dcsp,ddcsp,dr)
        DO 107  I=1,NP                                                    
        DV(I,2)=DV(I,1)                                                   
!        if (i.eq.5) write (6,*)'dv',dv(i,1),dv(i,2)                        
        DP(I)=0.                                                          
        qrhou=d(i,1)
        qrhod=d(i,2)
        fu=qrhou/4.d0/pi
        fd=qrhod/4.d0/pi
        fu=fu/dr(i)/dr(i)
        fd=fd/dr(i)/dr(i)
        ft=fu+fd
!
        qrhou=dcsp(i,1)
        qrhod=dcsp(i,2)
        dndru=qrhou/4.d0/pi
        dndrd=qrhod/4.d0/pi
        dndrt=dndru+dndrd
!       
        qrhou=ddcsp(i,1)
        qrhod=ddcsp(i,2)
        d2ndru=qrhou/4.d0/pi
        d2ndrd=qrhod/4.d0/pi
        d2ndrt=d2ndru+d2ndrd
!
        dndtt=0.
        dndpt=0.
        d2ndtt=0.
        d2ndpt=0.
        dndrtt=0.
        dndrpt=0.
        dndtpt=0.
!
        dndtu=0.
        dndpu=0.
        d2ndtu=0.
        d2ndpu=0.
        dndrtu=0.
        dndrpu=0.
        dndtpu=0.
!
        dndtd=0.
        dndpd=0.
        d2ndtd=0.
        d2ndpd=0.
        dndrtd=0.
        dndrpd=0.
        dndtpd=0.
        qsin=0.
        qcos=0.
        call grans(dr(i),qsin,qcos,dndrt,dndtt,dndpt, &
                   d2ndrt,d2ndtt,d2ndpt,dndrtt,dndrpt,dndtpt, &
                   gx,gy,gz,gmag,g2,ggx,ggy,ggz,gdggt)
        call grans(dr(i),qsin,qcos,dndru,dndtu,dndpu, &
                   d2ndru,d2ndtu,d2ndpu,dndrtu,dndrpu,dndtpu, &
                   gxu,gyu,gzu,gmagu,g2u,ggxu,ggyu,ggzu,gdggu)
        call grans(dr(i),qsin,qcos,dndrd,dndtd,dndpd, &
                   d2ndrd,d2ndtd,d2ndpd,dndrtd,dndrpd,dndtpd, &
                   gxd,gyd,gzd,gmagd,g2d,ggxd,ggyd,ggzd,gdggd)
!
        CALL VXCLM2(FU,GXU,GYU,GZU,GMAGU,G2U,GGXU,GGYU,GGZU, &
                    FD,GXD,GYD,GZD,GMAGD,G2D,GGXD,GGYD,GGZD, &
                    vxcup,vxcdn,dexc1,dexc2,GDGGT,iex)
        VXCup=DR(I)*VXCup/2.0d0                        
        VXCdn=DR(I)*VXCdn/2.0d0                        
!        if (i.eq.5) write (6,*)'DEXCup,dexcdn,z,icut',dexcup,dexcdn
!     *  , z,icut                        
!        if (i.eq.5) write (6,*)'dv1',dv(i,1),dv(i,2)                        
        DV(I,1)=DV(I,1)-Z+VXCup                                    
        DV(I,2)=DV(I,2)-Z+VXCdn                                    
!        if (i.eq.5) write (6,*)'dv2',dv(i,1),dv(i,2)                        
        IF (ICUT.NE.0) GO TO 109                                            
        IF (DV(I,1).GT.DLO) DV(I,1)=DLO                           
        IF (DV(I,2).GT.DLO) DV(I,2)=DLO                           
 109    DV(I,1)=DV(I,1)/DR(I)
 107    DV(I,2)=DV(I,2)/DR(I)
      endif
      RETURN                                                            
      END                                                               




