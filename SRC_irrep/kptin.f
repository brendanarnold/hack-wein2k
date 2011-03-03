      SUBROUTINE KPTIN(FL,NTYPE,NLO,SK,NE,NV,ISK,ISKDEN,KKK)
      USE FELDER
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      REAL*8, ALLOCATABLE ::  APA(:)          ! NMAT
      COMPLEX*16     CZERO,IMAG
      CHARACTER*10   BNAME
      LOGICAL        FL(FLMAX)
      DIMENSION      SK(3),ISK(3)
      DIMENSION      IPRIME(13), ID(3)
      DATA           IPRIME/2,3,5,7,11,13,17,19,23,29,31,37,0/  
      DATA           CZERO/(0.0D0,0.0D0)/,IMAG/(0.0D0,1.0D0)/
!******************************************************************************
!
      NNE=0
      READ(10,END=999) SK(1),SK(2),SK(3),BNAME,NV,NE,WEIGHT
      IF(FL(2)) READ(9,END=999) SK(1),SK(2),SK(3),BNAME,NV,NE,WEIGHT
!
      ALLOCATE ( A(NV,NE), B(NV,NE), EE(NE), KV(3,NV) )
      ALLOCATE ( PH(NSYM,NV), XM(NSYM,NE), L(NSYM,NV) )
      ALLOCATE ( APA(NV) )
!
      READ(10) ((KV(J,I),J=1,3),I=1,NV)
      IF(FL(2)) READ(9) ((KV(J,I),J=1,3),I=1,NV)
      DO 210 J=1,NE
        READ(10) NUM,EE(J)
        IF(FL(2)) READ(9) NUM,EE(J)
        NNE=NNE+1
        EE(NNE)=EE(J)
        DO 216 I=1,NV
 216      B(I,NNE)=CZERO
        IF(FL(1)) THEN
            READ(10) (A(I,NNE),I=1,NV)
          IF(FL(2)) THEN
            READ(9)  (B(I,NNE),I=1,NV)
          ENDIF
        ELSE
          READ(10) (APA(I),I=1,NV)
          DO 240 I=1,NV
 240         A(I,NNE)=CMPLX(APA(I),0.D0)
          IF(FL(2)) THEN
            READ(9) (APA(I),I=1,NV)
            DO  I=1,NV
              B(I,NNE)=CMPLX(APA(I),0.D0)
            ENDDO  
          ENDIF
        ENDIF
!
!.......the complex conjugate
        if(2.eq.1) then
        DO 250 I=1,NV
          A(I,NNE)=CONJG(A(I,NNE))
          B(I,NNE)=CONJG(B(I,NNE))
 250    CONTINUE
        endif
!
  210 CONTINUE
      NV=NV-NLO
      NE=NNE
!
!.....SK(i)=ISK(i)/ISKDEN, where ISK(i) and ISKDEN are integers
      DO 15 I=1,3
  	J=1
	DO WHILE((J.LT.MSTP).AND. &
                 (ABS(J*SK(I)-NINT(J*SK(I))).GT.(1.D0/MSTP)))
  	  J=J+1
	END DO
 	IF(ABS(SK(I)).LT.(1.D0/MSTP)) THEN
	  ISK(I)=0
	  ID(I)=1
	ELSE
	  ISK(I)=NINT(J*SK(I))
	  ID(I)=J
	ENDIF
        IF(I.EQ.1) THEN	
 	  ISKDEN=ID(I)
  	ELSE
 	  ISKDEN=ISKDEN*ID(I)
	ENDIF
   15 CONTINUE
      DO 20 I=1,3
   20   ISK(I)=ISK(I)*ISKDEN/ID(I)
!
!
      IPM=1
      DO WHILE(IPM.NE.0)
        IMM=IPRIME(IPM)
        DO WHILE(  MOD(ISK(1),IMM).EQ.0 .AND.  &
                   MOD(ISK(2),IMM).EQ.0 .AND. &
                   MOD(ISK(3),IMM).EQ.0 .AND. &
                   MOD(ISKDEN ,IMM).EQ.0      )
          DO 30 I=1,3
   30       ISK(I)=INT(ISK(I)/IMM)
          ISKDEN=ISKDEN/IMM
        END DO
        IPM=IPM+1
        IF(IPRIME(IPM).EQ.0) IPM=0
      END DO    
!
!.....K -> (K+k)*ISKDEN
      DO 50 I=1,NV
      DO 50 J=1,3
   50   KV(J,I)=KV(J,I)*ISKDEN+ISK(J)
!
!.....output
      WRITE(6,500)
      WRITE(6,509) KKK,BNAME
      WRITE(6,510) (SK(I),I=1,3)
!
      DEALLOCATE (APA)
      RETURN
 999  KKK=-KKK
 500  FORMAT(/,80('*'))
 509  FORMAT(/,/,'knum =',I3,4X,'kname= ',A10)
 510  FORMAT('k =',3F9.6,/)
      END 
