      SUBROUTINE MMatrix(JATOM,ivmax,NN_,N_,is)
      use mxyz
      use ablm, a=>alm,b=>blm,c=>clm
      use intu
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
!......COMPUTING : < PSI'| GRAD | PSI >...............................
!................. | PSI > := (ALM(NB) UL + BLM(NB) UPL) YLM .........
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
       IMPLICIT REAL*8 (A-H,O-Z)
       INCLUDE 'param.inc'
       INTEGER JATOM,J,N_(ivmax),NN_(ivmax)
!       COMPLEX*16  A,B,C,MX_,MY_,MZ_,SX_,SY_,SZ_
       COMPLEX*16  CZERO,IMAG
       LOGICAL     lapw(0:lmax2),loor(0:lomax),lloor(0:lmax2)
!ad
!      COMMON /INTU/ Duu1(LMAX1,NATO,2),Duu2(LMAX1,NATO,2), &
!                    Duup1(LMAX1,NATO,2),Duup2(LMAX1,NATO,2), &
!                    Dupu1(LMAX1,NATO,2),Dupu2(LMAX1,NATO,2), &
!                    Dupup1(LMAX1,NATO,2),Dupup2(LMAX1,NATO,2), &
!                    Ruu(LMAX1,NATO,2),Ruup(LMAX1,NATO,2), &
!                    Rupu(LMAX1,NATO,2),Rupup(LMAX1,NATO,2)
!      COMMON /INTUL/Duul1(LOMAX1,NATO,2), Duul2(LOMAX1,NATO,2), &
!                    Dulup1(LOMAX1,NATO,2),Dulup2(LOMAX1,NATO,2), &
!                    Dupul1(LOMAX1,NATO,2),Dupul2(LOMAX1,NATO,2), &
!                    Dulu1(LOMAX1,NATO,2), Dulu2(LOMAX1,NATO,2), &
!                    Dulul1(LOMAX1,NATO,2),Dulul2(LOMAX1,NATO,2), &
!                    Ruul(LOMAX1,NATO,2),  Rulu(LOMAX1,NATO,2), &
!                    Rupul(LOMAX1,NATO,2), Rulup(LOMAX1,NATO,2), &
!                    Rulul(LOMAX1,NATO,2)
!      COMMON /ABLM/ A(NUME,(LMAX+1)**2),B(NUME,(LMAX+1)**2), &
!                    C(NUME,(LMAX+1)**2)
!      COMMON /MXYZ/ MX_(NUMEO),MY_(NUMEO),MZ_(NUMEO), &
!                    SX_(NUMEO),SY_(NUMEO),SZ_(NUMEO)
      common /lolog/   nlo,nlov,nlon,lapw,ilo(0:lomax),loor,lloor
!ad 
       DATA    CZERO/(0.0D0,0.0D0)/IMAG/(0.0D+0,1.0D+0)/ 
       DATA    ONE/1.0D+0/TWO/2.0D+0/THREE/3.0D+0/

!........INLINE FUNCTIONS FOR YLM OTHOGONALITY ....................
!      
       F0(L)   = (TWO*L+ONE)*(TWO*L+THREE)
       F1(L,M) = - SQRT(( L + M + ONE ) * ( L + M + TWO ) / F0(L))
       F3(L,M) =   SQRT(( L - M + ONE ) * ( L - M + TWO ) / F0(L))
       F5(L,M) =   SQRT(( L - M + ONE ) * ( L + M + ONE ) / F0(L))
!..................................................................
         J=JATOM
         IDEX=0 
         do iv=1,ivmax
         SX_(iv)=CZERO
         SY_(iv)=CZERO
         SZ_(iv)=CZERO
         enddo
!ad                                          
!ad......................................................L=0,LMAX-1
!ad
      DO  L=0,LMAX-1 
!ad
        L1=L+1
        L2=L+2
        H1  = (Duu1  (L1,J,is) - L * Ruu  (L1,J,is))
        H2  = (Duup1 (L1,J,is) - L * Ruup (L1,J,is))
        H3  = (Dupu1 (L1,J,is) - L * Rupu (L1,J,is))
        H4  = (Dupup1(L1,J,is) - L * Rupup(L1,J,is))
        H5  = (Duu2  (L1,J,is) + L2* Ruu  (L1,J,is))
        H6  = (Duup2 (L1,J,is) + L2* Rupu (L1,J,is))
        H7  = (Dupu2 (L1,J,is) + L2* Ruup (L1,J,is))
        H8  = (Dupup2(L1,J,is) + L2* Rupup(L1,J,is))
!ad
!ad   for LO's
!ad
      if (L1.LE.LOMAX+1) THEN
        HL1 = (Duul1 (L1,J,is) - L * Ruul(L1,J,is))
        HL2 = (Dupul1(L1,J,is) - L * Rupul(L1,J,is))
        HL3 = (Dulu2 (L1,J,is) + L2* Ruul(L1,J,is))
        HL4 = (Dulup2(L1,J,is) + L2* Rupul(L1,J,is))
        HL5 = (Dulu1 (L1,J,is) - L * Rulu(L1,J,is))
        HL6 = (Dulup1(L1,J,is) - L * Rulup(L1,J,is))
        HL7 = (Duul2 (L1,J,is) + L2* Rulu(L1,J,is))
        HL8 = (Dupul2(L1,J,is) + L2* Rulup(L1,J,is))
        HL9 = (Dulul1(L1,J,is) - L * Rulul(L1,J,is))
        HL0 = (Dulul2(L1,J,is) + L2* Rulul(L1,J,is))
!adr
!       write(*,345) l,Hl1,hl2,hl3,hl4,hl5,hl6,hl7,hl8
!adr 
      end if
!ad
!ad..........................................................M=-L,L
!ad
      DO M=-L,L 
!ad
!ad......IDEX  = (l,m)
!ad......IDEX1 = (l+1,m)
!ad......IDEX11=(l+1,m+1)
!ad......IDEX0 =(l+1,m-1)
!ad
!...IDEX->l,m..IDEX11->l+1,m+1..IDEX10->l+1,m-1..IDEX1->l+1,m
       IDEX  = IDEX+1
       IDEX11= IDEX+2*L+3 
       IDEX10= IDEX+2*L+1 
       IDEX1 = IDEX+2*L+2
!ad
!ad   for LAPW's
!ad
       do iv=1,ivmax
       n=n_(iv)
       nn=nn_(iv)
!ad
!ad   x+iy component
!ad
       SX_(iv)=SX_(iv) &
         +(H1*CONJG(A(nn,IDEX11))*A(n,IDEX)  &
         + H2*CONJG(A(nn,IDEX11))*B(n,IDEX) &
         + H3*CONJG(B(nn,IDEX11))*A(n,IDEX) &
         + H4*CONJG(B(nn,IDEX11))*B(n,IDEX) &
                                           )*F1(L,M)   &
         +(H5*CONJG(A(nn,IDEX))  *A(n,IDEX10) &
         + H6*CONJG(A(nn,IDEX))  *B(n,IDEX10) &
         + H7*CONJG(B(nn,IDEX))  *A(n,IDEX10) &
         + H8*CONJG(B(nn,IDEX))  *B(n,IDEX10) &
                                           )*F3(L,M)
!ad
!ad   x-iy component
!ad
       SY_(iv)=SY_(iv) &
         +(H1*CONJG(A(nn,IDEX10))*A(n,IDEX)  &
         + H2*CONJG(A(nn,IDEX10))*B(n,IDEX)  &
         + H3*CONJG(B(nn,IDEX10))*A(n,IDEX) &
         + H4*CONJG(B(nn,IDEX10))*B(n,IDEX) &
                                           )*F3(L,M)   &
         +(H5*CONJG(A(nn,IDEX))  *A(n,IDEX11) &
         + H6*CONJG(A(nn,IDEX))  *B(n,IDEX11) &
         + H7*CONJG(B(nn,IDEX))  *A(n,IDEX11) &
         + H8*CONJG(B(nn,IDEX))  *B(n,IDEX11) &
                                           )*F1(L,M)
!ad
!ad   z component
!ad
       SZ_(iv)=SZ_(iv) &
         +(H1*CONJG(A(nn,IDEX1)) *A(n,IDEX)  &
         + H2*CONJG(A(nn,IDEX1)) *B(n,IDEX) &
         + H3*CONJG(B(nn,IDEX1)) *A(n,IDEX) &
         + H4*CONJG(B(nn,IDEX1)) *B(n,IDEX) &
         + H5*CONJG(A(nn,IDEX))  *A(n,IDEX1) &
         + H6*CONJG(A(nn,IDEX))  *B(n,IDEX1) &
         + H7*CONJG(B(nn,IDEX))  *A(n,IDEX1) &
         + H8*CONJG(B(nn,IDEX))  *B(n,IDEX1) &
                                          )*F5(L,M)
      enddo
!ad
!ad   for LO's
!ad
!ad   IF (lloor(l)) THEN
      IF (lloor(l)) THEN
!ad
      do iv=1,ivmax
      n=n_(iv)
      nn=nn_(iv)
!ad
!ad   x+iy component
!ad
       SX_(iv)=SX_(iv) &
         +(HL1*CONJG(A(nn,IDEX11))*C(n,IDEX) &
         + HL2*CONJG(B(nn,IDEX11))*C(n,IDEX) &
                                            )*F1(L,M) &
         +(HL3*CONJG(C(nn,IDEX))  *A(n,IDEX10) &
         + HL4*CONJG(C(nn,IDEX))  *B(n,IDEX10) &
                                            )*F3(L,M)
!ad
!ad   x-iy component
!ad
       SY_(iv)=SY_(iv) &
         +(HL1*CONJG(A(nn,IDEX10))*C(n,IDEX) &
         + HL2*CONJG(B(nn,IDEX10))*C(n,IDEX) &
                                            )*F3(L,M) &
         +(HL3*CONJG(C(nn,IDEX))*  A(n,IDEX11) &
         + HL4*CONJG(C(nn,IDEX))*  B(n,IDEX11) &
                                            )*F1(L,M)
!ad
!ad   z component
!ad
       SZ_(iv)=SZ_(iv) &
         +(HL1*CONJG(A(nn,IDEX1)) *C(n,IDEX) &
         + HL2*CONJG(B(nn,IDEX1)) *C(n,IDEX) &
         + HL3*CONJG(C(nn,IDEX))  *A(n,IDEX1) &
         + HL4*CONJG(C(nn,IDEX))  *B(n,IDEX1) &
                                           )*F5(L,M)
!ad
       enddo
!ad
       END IF
       IF (lloor(l+1)) THEN   
 
       do iv=1,ivmax
       n=n_(iv)
       nn=nn_(iv)
!ad
!ad   x+iy component
!ad
       SX_(iv)=SX_(iv) &
         +(HL5*CONJG(C(nn,IDEX11))*A(n,IDEX) &
         + HL6*CONJG(C(nn,IDEX11))*B(n,IDEX) &
                                            )*F1(L,M) &
         +(HL7*CONJG(A(nn,IDEX))  *C(n,IDEX10) &
         + HL8*CONJG(B(nn,IDEX))  *C(n,IDEX10) &
                                            )*F3(L,M)
!ad
!ad   x-iy component
!ad
       SY_(iv)=SY_(iv) &
         +(HL5*CONJG(C(nn,IDEX10))*A(n,IDEX) &
         + HL6*CONJG(C(nn,IDEX10))*B(n,IDEX) &
                                            )*F3(L,M) &
         +(HL7*CONJG(A(nn,IDEX))  *C(n,IDEX11) &
         + HL8*CONJG(B(nn,IDEX))  *C(n,IDEX11) &
         )*F1(L,M)
!ad
!ad   z component
!ad
       SZ_(iv)=SZ_(iv) &
         +(HL5*CONJG(C(nn,IDEX1)) *A(n,IDEX) &
         + HL6*CONJG(C(nn,IDEX1)) *B(n,IDEX) &
         + HL7*CONJG(A(nn,IDEX))  *C(n,IDEX1) &
         + HL8*CONJG(B(nn,IDEX))  *C(n,IDEX1) &
                                           )*F5(L,M)
!ad
       enddo
!ad
       END IF
!ad
       IF ( (lloor(l)).AND.(lloor(l+1)) ) THEN
!ad
       do iv=1,ivmax
       n=n_(iv)
       nn=nn_(iv)
!ad
!ad   x+iy component
!ad
       SX_(iv)=SX_(iv) &
         +(HL9*CONJG(C(nn,IDEX11))*C(n,IDEX) &
                                            )*F1(L,M) &
         +(HL0*CONJG(C(nn,IDEX))  *C(n,IDEX10) &
                                            )*F3(L,M)
!ad
!ad   x-iy component
!ad
       SY_(iv)=SY_(iv) &
        +(HL9*CONJG(C(nn,IDEX10))*C(n,IDEX) &
                                           )*F3(L,M) &
        +(HL0*CONJG(C(nn,IDEX))  *C(n,IDEX11) &
                                           )*F1(L,M)
!ad
!ad   z component
!ad
       SZ_(iv)=SZ_(iv) &
         +(HL9*CONJG(C(nn,IDEX1)) *C(n,IDEX) &
         + HL0*CONJG(C(nn,IDEX))  *C(n,IDEX1) &
                                           )*F5(L,M)
!ad
       enddo
!ad           
       END IF
!ad
!ad    LO's done
!ad
!ad    END IF

       END DO 
!ad
!ad................................................end.......M=-L,L
!ad
       END DO
!ad                                          
!ad................................................end...L=0,LMAX-1
!ad
!ad
!ad    Transformation to Cartesian coordinates
!ad
       do iv=1,ivmax
         MX_(iv)=(SX_(iv)+SY_(iv)) / TWO
         MY_(iv)=(SX_(iv)-SY_(iv)) /IMAG  / TWO
         MZ_(iv)=SZ_(iv)
       enddo

       RETURN
       END
