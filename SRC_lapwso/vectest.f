      program vectest
! program vectest reads unformatted case.vector files created by
! LAPW1
        USE param
      IMPLICIT REAL*8(A-H,O-Z)
!
      COMPLEX*16     A(NMAT,NUME,2)
      CHARACTER*10   BNAME
      DIMENSION      APA(NMAT),EE(NUME,2),EMM(3),SS(3),KV(3,NMAT)
      DIMENSION       EMM(3),E(LMX,NATO,2)
      DOUBLE PRECISION   DPELO(0:LOMAX,NATO,2), ELO(0:LOMAX,NATO,2)
      ntap=8
      ity=1
      read(5,*)jspin,icmplx,nk
      do 999 isi=1,jspin
      ntap=ntap+1
      READ(ntap) (E(L,ITY,isi),L=1,LMX)
      write(6,*)' lmax',lmax
      write(6,100)(e(l,ity,isi),l=1,lmx)
      READ(ntap) (ELO(L,ITY,isi),L=0,LOMAX)
      write(6,*)' lomax',lomax
      write(6,100)(ELO(L,ITY,isi),L=0,LOMAX)
100   format(5e12.4)
!
!**********************************************************************
!
      do 210 ikpt=1,nk
      READ(ntap,end=998) SS(1),SS(2),SS(3),BNAME,NV,NE,WEIGHT
104   format(3f8.3,a10,2i4,f8.4)
      write(6,107)ss(1),ss(2),ss(3),nv,ne
107   format(' ss,nv,ne',3f8.3,2i5)
      READ(ntap) (KV(1,I),KV(2,I),KV(3,I),I=1,NV)
      write(6,103) (KV(1,I),KV(2,I),KV(3,I),I=1,NV)
103   format(5(3i4,x))
      DO 210 J=1,NE
        READ(ntap) NUM,EE(J,isi)
        write(6,105) NUM,EE(J,isi)
        write(8,205) NUM,EE(J,isi)
205   format(i4,4x,e24.14)
105   format(i4,e14.5)
          IF(icmplx.ne.0) THEN 
            READ(ntap) (A(I,j,isi),I=1,NV)
         write(6,106)  (A(I,j,isi),I=1,NV)
106   format(6e12.3)
          ELSE
            READ(ntap) (APA(I),I=1,NV)
         write(6,106)  (APA(I),I=1,NV)
            DO 200 I=1,NV
              A(I,j,isi)=CMPLX(APA(I),0.D0)
  200       CONTINUE
           ENDIF
       write(6,108)j,num,ee(j,isi),a(1,j,isi),a(nv,j,isi)
108    format(' j,num,e,a1,anv',2i4,5e12.4)
  210 CONTINUE
999   continue
  998 continue
      stop   
      END 
