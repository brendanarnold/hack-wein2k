      SUBROUTINE KPTIN(flcom,A,SS,NE,NV,EE,KV,bname,weight,kkk,isi)
        USE param
      IMPLICIT REAL*8(A-H,O-Z)
      logical flcom
!
      COMPLEX*16     A(NMAT,NUME,2)
!C_COMPLEX      COMPLEX*16     APA(NMAT)
      REAL*8      APA(NMAT)
      CHARACTER*10   BNAME
      DIMENSION      EE(NUME,2),SS(3),KV(3,NMAT,2)
!**********************************************************************
!
      nfile=8+isi
      NNE=0
      READ(nfile,end=999) SS(1),SS(2),SS(3),BNAME,NV,NE,WEIGHT
      IF(NV.GT.NMAT) write(6,*) 'TOO MANY KJS, NV.GT.NMAT ',nv,nmat
      IF(NE.GT.NUME) write(6,*) 'TOO MANY BANDS, NE.GT.NUME',ne,nume
      IF(NV.GT.NMAT) STOP 'TOO MANY KJS: '
      IF(NE.GT.NUME) STOP 'TOO MANY BANDS: '
      READ(nfile)(KV(1,I,isi),KV(2,I,isi),KV(3,I,isi),I=1,NV)
      DO 210 J=1,NE
        READ(nfile)NUM,EE(J,isi)
          NNE=NNE+1
          EE(NNE,isi)=EE(J,isi)
          IF(FLcom) THEN 
            READ(nfile) (A(I,NNE,isi),I=1,NV)
          ELSE
            READ(nfile)(APA(I),I=1,NV)
            DO 200 I=1,NV
              A(I,NNE,isi)=CMPLX(APA(I),0.D0)
  200       CONTINUE
           ENDIF
  210 CONTINUE
      NE=NNE
      return
 999  kkk=-kkk
      return
      END
