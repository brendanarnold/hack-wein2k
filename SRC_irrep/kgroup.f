      SUBROUTINE KGROUP(FL,DB1,IZ,IIZ,IORD,ISK,ISKDEN,LKG,IKG, &
                        KKK,NV)
      use felder
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!      INCLUDE 'param.inc'
!
!      DIMENSION        KV(3,NMAT)
      LOGICAL          FL(FLMAX)
      DIMENSION        IZ(3,3,NSYM),IIZ(3,3,NSYM)
      DIMENSION        ISK(3),IKR(3)
      DIMENSION        LKG(NSYM),DB1(3,3),DKR(3)
      DATA             TOLK/1E-6/
!******************************************************************************
!
!.....determine the group G(k) of the allowed wave vector k,
!.....i.e., finding those Ri having the property k*inv(Ri)=k+Km
!.....db1(3,3) transforms Km to units of primitiv vectors.
      DO 2 I=1,NSYM
    2   LKG(I)=0
      IKG=0
      DO 40 I=1,IORD
      DO 42 I1=1,3
 42      IKR(I1)=ISK(1)*IIZ(1,I1,I) &
                +ISK(2)*IIZ(2,I1,I) &
                +ISK(3)*IIZ(3,I1,I) -ISK(I1) 
!
        DO 44 I1=1,3
44         DKR(I1)=(DB1(I1,1)*DBLE(IKR(1)) &
                  + DB1(I1,2)*DBLE(IKR(2)) &
                  + DB1(I1,3)*DBLE(IKR(3)))/DBLE(ISKDEN)
        DTEST=DABS(NINT(DKR(1))-DKR(1)) &
             +DABS(NINT(DKR(2))-DKR(2)) &
             +DABS(NINT(DKR(3))-DKR(3))
!las1
!       DO 441 IN=1,NV
!          ITEST=ABS(IKR(1)-KV(1,IN))
!     &         +ABS(IKR(2)-KV(2,IN))
!     &         +ABS(IKR(3)-KV(3,IN))
!          IF(ITEST.EQ.0) GOTO 451 
! 441   CONTINUE
! 451   CONTINUE
!       DTEST=DBLE(ITEST)
!las2
        IF(DTEST.LE.TOLK) THEN
          IKG=IKG+1
          LKG(IKG)=I 
        ENDIF
   40   CONTINUE
!
!las
!      IF(IKG.LE.1) THEN
!        IF(LKG(1).NE.1) STOP 'kgroup: Identity operator not found'
!        WRITE(6,*) 'Only the identity operator E:'
!        WRITE(6,*) 'The pointgroup is C1'
!        WRITE(6,*) 'Charater table according to Table 3, page 31,',
!     &             'in G.F. Koster, et al.'
!        IF(.NOT.FL(2)) THEN
!          WRITE(6,*) 'All states have G1-symmetry' 
!        ELSE
!          WRITE(6,*) 'All states have G2-symmetry of the double group' 
!        ENDIF
!      ENDIF
!las
      RETURN
      END  
