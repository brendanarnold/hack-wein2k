      SUBROUTINE  Horb(NV,NAT,NNSL,POS,VI,R,MULT, &
                     LQIND,LQNS,PEI,ROTLOC,ROTIJ,H)
!
      use matrices, only: HSROWS, KZZ, XK, YK, ZK
!      use orb, only   : VORB, IATOM, NLORB, LORB
      use loabc, only : ALO, BLO, CLO, PI12LO, PE12LO
      use lolog, only : NLO, ILO
      use orb, only   : BEXT, VORB, NMOD, NSP, NATORB, IATOM, NLORB, LORB
      use struk, only : NDF
      use albl, only : AL, BL
!      use lapw_timer, only :  time_orb1, time_orb2, time_orb3, time_orb4, &
!         START_TIMER, STOP_TIMER, READ_CPU_TIME
      IMPLICIT NONE
      INCLUDE 'param.inc'
!
!        Arguments
!
      INTEGER            NNSL, NV, NAT
      INTEGER            MULT(NAT), LQIND(NAT)
      INTEGER            LQNS(6,NGAU,NAT)
!      DOUBLE PRECISION   OTIME1, OTIME2, OTIME3, OTIME4
      DOUBLE PRECISION   VI
      DOUBLE PRECISION   R(NAT)
      DOUBLE PRECISION   ROTLOC(3,3,NAT), ROTIJ(3,3,NDF)
      DOUBLE PRECISION   POS(3,NDF)
!             , XK(HSROWS+1), YK(HSROWS+1), ZK(HSROWS+1)
!_REAL      DOUBLE PRECISION   H(HSROWS,HSROWS)
!_COMPLEX      COMPLEX*16         H(HSROWS,HSROWS)
!!_REAL      DOUBLE PRECISION   pot(-3:3,-3:3)
!!_COMPLEX      COMPLEX*16         pot(-3:3,-3:3)
      COMPLEX*16         pot(-3:3,-3:3)
      DOUBLE PRECISION   PEI(LMAX,NAT)
!
!        AL(n,l,i) - coefficient al(kn) for atom i
!                    (including equivalent atoms)
!        BL(n,l,i) - coefficient bl(kn) for atom i
!                    (including equivalent atoms)
!
!!      DOUBLE PRECISION, allocatable :: AL(:,:,:)
!!      DOUBLE PRECISION, allocatable :: BL(:,:,:)
!
      COMPLEX*16, allocatable :: alm(:,:),blm(:,:),clm(:,:)
      COMPLEX*16, allocatable :: almt(:,:),blmt(:,:)
      COMPLEX*16  TMP1,TMP2,TMP3
!_REAL      DOUBLE PRECISION TMP4,TMP5,TMP6,TMP7
!_COMPLEX   COMPLEX*16 TMP4,TMP5,TMP6,TMP7
!
!        YL(lm,:,i) - temporary storage of spherical harmonics for
!                     local orbitals
!
      COMPLEX*16, allocatable :: YL(:,:)
!
!        Local Parameters
!
      COMPLEX*16         CIMAG
      PARAMETER          (CIMAG  = (0.0D+0,1.0D+0))
!
!        Locals
!
      DOUBLE PRECISION   ARGX, ARGY, ARGZ, RHOATM, PI, TPI, ARG
      DOUBLE PRECISION   VEC(3), ROTV1(3), ROTV2(3)
      INTEGER            i1, i2, jlo, i12, i22, j, jlo2, nvorb, jeq, nlv
      INTEGER            iat, indatm, jneq, lvor, l0m0, l0, m0, m1, m2, i
      DOUBLE COMPLEX     czero
      PARAMETER          (CZERO=(0.0D0,0.0D0) )

!
      COMPLEX*16         CIL
      COMPLEX*16, allocatable :: PHS(:)
!      COMPLEX*16         PHS(HSROWS)
!
!        External Subroutines
!
      EXTERNAL           ROTATE, YLM
!
!        Intrinsic Functions
!
      INTRINSIC          DCONJG, EXP, SQRT, ATAN, DBLE, DIMAG
!
!!      allocate( AL(HSROWS, 0:LMAX-1,NDF) )
!!      allocate( BL(HSROWS, 0:LMAX-1,NDF) )
      allocate( ALM(-3:3,HSROWS) )
      allocate( BLM(-3:3,HSROWS) )
      allocate( CLM(-3:3,HSROWS) )
      allocate( ALMt(-3:3,HSROWS) )
      allocate( BLMt(-3:3,HSROWS) )
      allocate( YL(0:NSLMAX*NSLMAX-1, HSROWS) )
      allocate( PHS(HSROWS) )

!      OTIME1=0
!      OTIME2=0
!      OTIME3=0
!      OTIME4=0
      PI = 4.0D+0*ATAN(1.0D+0)
      TPI= 2.0D0*PI
      iat=0
      INDATM = 0
      DO 110 JNEQ = 1, NAT
        nvorb=0
         if(iatom(iat+1).eq.jneq)then
           nvorb=1  
           RHOATM = 4.0D+0*PI*SQRT(VI)*R(JNEQ)*R(JNEQ)
           iat=iat+1
         endif
         DO 120 JEQ = 1, MULT(JNEQ)
          INDATM = INDATM + 1
           if(nvorb.eq.1)then
            do nlv=1,nlorb(iat)
            lvor=lorb(nlv,iat)
             do m1=-lvor,lvor
              do m2=-lvor,lvor
!               pot(m1,m2)=vorb(iat,lvor,m1,m2)
                pot(m2,m1)=vorb(iat,lvor,m1,m2)
              enddo
             enddo
!              CALL START_TIMER(time_orb1)
              DO 170 I = 1, NV+NLO
!
!        calculate phase-factor for Alm, Blm
!
!               ARGX = POS(1,INDATM)*KZZ(1,I)
!               ARGY = POS(2,INDATM)*KZZ(2,I)
!               ARGZ = POS(3,INDATM)*KZZ(3,I)
!               PHS(I) = EXP(CIMAG*2.0D+0*PI*(ARGX+ARGY+ARGZ))
                ARG = ( POS(1,INDATM)*KZZ(1,I) + &
                        POS(2,INDATM)*KZZ(2,I) + &
                        POS(3,INDATM)*KZZ(3,I))*TPI
                PHS(I) = DCMPLX( COS(ARG), SIN(ARG))
!
!        calculate spherical-harmonics for Alm, Blm
!
               VEC(1) = XK(I)
               VEC(2) = YK(I)
               VEC(3) = ZK(I)
               CALL ROTATE (VEC,ROTIJ(1,1,INDATM),ROTV1)
               CALL ROTATE (ROTV1,ROTLOC(1,1,JNEQ),ROTV2)
               CALL YLM (ROTV2,lvor,YL(0,I))
170      continue
!        CALL STOP_TIMER(time_orb1)
!        OTIME1=OTIME1+READ_CPU_TIME(time_orb1)
!
!        Alm, Blm
!
!            CALL START_TIMER(time_orb2)
            CIL = (1.0D+0,0.0D+0)
            l0m0=0
            DO 140 L0 = 0, lvor
!               do 165 i=1,nv+nlo
!               DO 165 M0 = -L0, L0
!               DO 150 M0 = -L0, L0
                do 165 i=1,nv+nlo
                do 165 M0 = -L0,L0
                alm(m0,i)=czero
                blm(m0,i)=czero
                clm(m0,i)=czero
165             continue
               DO 150 M0 = -L0, L0
                if(L0.eq.lvor)then
                do 160 i=1,nv
                TMP1=rhoatm*cil*phs(i)*conjg(yl(l0m0,i))
!             ALM(I,M0)=rhoatm*cil*phs(i)*conjg(yl(l0m0,i))*al(i,l0,jneq)
!             BLM(I,M0)=rhoatm*cil*phs(i)*conjg(yl(l0m0,i))*bl(i,l0,jneq)
                ALM(M0,I)=TMP1*al(i,l0,jneq)
                BLM(M0,I)=TMP1*bl(i,l0,jneq)
160             continue
             do jlo=1,ilo(lvor,jneq)
              call lohns(jneq,mult,i1,i2,l0,jlo)
              do 166 i=nv+i1,nv+i2
               TMP1=rhoatm*cil*phs(i)*conjg(yl(l0m0,i))
!               ALM(I,M0)=rhoatm*cil*phs(i)*conjg(yl(l0m0,i))
               blm(m0,i)=blo(l0,jlo,jneq)*TMP1
               clm(m0,i)=clo(l0,jlo,jneq)*TMP1
               alm(m0,i)=alo(l0,jlo,jneq)*TMP1
166           continue
             enddo
            endif
            l0m0=l0m0+1
150         continue
            cil=cil*cimag
140         continue
!        CALL STOP_TIMER(time_orb2)
!        OTIME2=OTIME2+READ_CPU_TIME(time_orb2)
!   end of ALM, BLM and CLM calculation
!   add orbital term to the hamiltonian
!   <APW||APW>
!        CALL START_TIMER(time_orb3)
                DO J=1,NV
                        DO M0=-LVOR,LVOR
                        TMP1=czero
                        TMP2=czero
                        do M1=-LVOR,LVOR
                                TMP1=TMP1+pot(m1,m0)*alm(m1,j)
                                TMP2=TMP2+pot(m1,m0)*blm(m1,j)
                        ENDDO
                        ALMT(M0,J)=TMP1
                        BLMT(M0,J)=TMP2*pei(lvor+1,jneq)
                        ENDDO
                ENDDO                                
!                DO 260 I = 1, NV
!                 do 261 j=1,i
                 DO 260 J=1,NV
                  do 261 I=J,NV
                   TMP4=0.0D0
                   TMP5=0.0D0
                   do m0=-lvor,lvor
!                   TMP1=czero
!                   TMP2=czero
!                   do m1=-lvor,lvor
!!                    H(i,j) = H(i,J)+pot(m0,m1)*    &
!!     &                 (dconjg(alm(m0,i))*alm(j,m1)+ &
!!     &                 (dconjg(blm(m0,i))*blm(j,m1))*pei(lvor+1,jneq))
!                        TMP1=TMP1+pot(m1,m0)*alm(m1,j)
!                        TMP2=TMP2+pot(m1,m0)*blm(m1,j)
!                   enddo
                   TMP4=TMP4+ALMT(M0,J)*dconjg(alm(m0,i))
                   TMP5=TMP5+BLMT(M0,J)*dconjg(blm(m0,i))
                   enddo
!                   H(i,j)=H(i,j)+TMP4+TMP5*pei(lvor+1,jneq)
                   H(i,j)=H(i,j)+TMP4+TMP5
  261            CONTINUE
  260           CONTINUE
!      CALL STOP_TIMER(time_orb3)
!      OTIME3=OTIME3+ READ_CPU_TIME(time_orb3)
!      CALL START_TIMER(time_orb4)
!
              do jlo=1,ilo(lvor,jneq)
              call lohns(jneq,mult,i1,i2,lvor,jlo)
!  <APW||LO>
              do i=nv+i1,nv+i2
                 do j=1,nv
                   TMP4=0.0D0
                   TMP5=0.0D0
                   TMP6=0.0D0
                   TMP7=0.0D0
                   do m0=-lvor,lvor
                   TMP1=czero
                   TMP2=czero
                   do m1=-lvor,lvor
!                    H(I,J) = H(I,J)+pot(m0,m1)*    &
!     &                 (dconjg(alm(m0,i))*alm(j,m1)+ &
!     &                  dconjg(blm(m0,i))*blm(j,m1)*pei(lvor+1,jneq) &
!     &                 +dconjg(clm(m0,i))*alm(j,m1)*pi12lo(lvor,jneq) &
!     &                 +dconjg(clm(m0,i))*blm(j,m1)*pe12lo(lvor,jneq))
                        TMP1=TMP1+alm(m1,j)*pot(m1,m0)
                        TMP2=TMP2+blm(m1,j)*pot(m1,m0)
                   enddo
                   TMP4=TMP4+TMP1*dconjg(alm(m0,i))
                   TMP5=TMP5+TMP2*dconjg(blm(m0,i))
                   TMP6=TMP6+TMP1*dconjg(clm(m0,i))
                   TMP7=TMP7+TMP2*dconjg(clm(m0,i))
                   enddo
                   H(I,J)=H(I,J)+TMP4+TMP5*pei(lvor+1,jneq)+ &
                        TMP6*pi12lo(lvor,jneq)+TMP7*pe12lo(lvor,jneq)
                  enddo
               enddo
!  <LO||LO>
                      do jlo2=1,jlo
              call lohns(jneq,mult,i12,i22,lvor,jlo2)
!                 write(*,*) jlo,jlo2,i1,i2,i12,i22
              do i=nv+i1,nv+i2
                 do j=nv+i12,min(i,nv+i22)
!                   write(*,*) i,j,jlo,jlo2
                   TMP4=0.0D0
                   TMP5=0.0D0
                   TMP6=0.0D0
                   TMP7=0.0D0
                   do m0=-lvor,lvor
                   TMP1=czero
                   TMP2=czero
                   TMP3=czero
                   do m1=-lvor,lvor
!                    H(I,J) = H(I,J)+pot(m0,m1)*    &
!     &                 ( &
!     &                  dconjg(alm(m0,i))*alm(j,m1)+ &
!     &                  dconjg(blm(m0,i))*blm(j,m1)*pei(lvor+1,jneq) &
!     &                 +(
!     &                   dconjg(alm(m0,i))*clm(j,m1)+ dconjg(clm(m0,i))*alm(j,m1) &
!     &                  )*pi12lo(lvor,jneq) &
!     &                +( &
!     &                   dconjg(blm(m0,i))*clm(j,m1)+ dconjg(clm(m0,i))*blm(j,m1) &
!     &                  )*pe12lo(lvor,jneq) &
!     &                + dconjg(clm(m0,i))*clm(j,m1) &
!                        )
                        TMP1=TMP1+pot(m1,m0)*alm(m1,j)
                        TMP2=TMP2+pot(m1,m0)*blm(m1,j)
                        TMP3=TMP3+pot(m1,m0)*clm(m1,j)
                   enddo
                   TMP4=TMP4+TMP1*dconjg(alm(m0,i))+ dconjg(clm(m0,i))*TMP3
                   TMP5=TMP5+dconjg(blm(m0,i))*TMP2
                   TMP6=TMP6+dconjg(clm(m0,i))*TMP1+ dconjg(alm(m0,i))*TMP3
                   TMP7=TMP7+dconjg(blm(m0,i))*TMP3+ dconjg(clm(m0,i))*TMP2
                   enddo
                   H(I,J)=H(I,J)+TMP5*pei(lvor+1,jneq)+TMP6*pi12lo(lvor,jneq)+&
                        TMP7*pe12lo(lvor,jneq)+TMP4
                  enddo
               enddo
                     ENDdo
                     ENDdo
!      CALL STOP_TIMER(time_orb4)
!      OTIME4=OTIME4+ READ_CPU_TIME(time_orb4)

           enddo
         endif

  120    CONTINUE
  110 CONTINUE
!      deallocate( AL )
!      deallocate( BL )
      deallocate( ALM )
      deallocate( BLM )
      deallocate( CLM )
      deallocate( YL )
      deallocate( PHS )
      deallocate( ALMT )
      deallocate( BLMT )
!      write(6,*)'timer orb1 ',OTIME1
!      write(6,*)'timer orb2 ',OTIME2
!      write(6,*)'timer orb3 ',OTIME3   
!      write(6,*)'timer orb4 ',OTIME4
      RETURN
      END
