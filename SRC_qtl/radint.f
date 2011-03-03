	SUBROUTINE RADINT(JATOM,ISPIN)
        USE param
        USE struct
        IMPLICIT REAL*8 (A-H,O-Z)

      COMMON /RADFU/   RF1(NRAD,0:LMAX2,2,nrf),RF2(NRAD,0:LMAX2,2,nrf)
      COMMON /RINTEG/  RI_MAT(0:lmax2,nrf,nrf,2,2)
      DIMENSION  rx(nrad)
!      
!.....set up radial mesh.
      DO  I=1,jrj(JATOM)
      RX(I)=r0(JATOM)*EXP(DX(JATOM)*(i-1))
      enddo
      
      ri_mat(0:lmax2,1:nrf,1:nrf,1:ispin,1:ispin)=0.d0

      DO 20 L=0,LMAX2
      DO 20 IS1=1,ISPIN
      DO 20 IS2=1,ISPIN
      DO 20 IF1=1,NRF
      DO 20 IF2=1,IF1
      CALL RINT13(rf1(1,l,is1,if1),rf2(1,l,is1,if1), &
           rf1(1,l,is2,if2),rf2(1,l,is2,if2), &
           ri_mat(l,if1,if2,is1,is2),JATOM)
      ri_mat(l,if2,if1,is2,is1)=ri_mat(l,if1,if2,is1,is2)
 20   CONTINUE
      END

