!----------------------------------------------------------------------
!
!     READ diagonal MATRIXELEMENTS
!
!----------------------------------------------------------------------
!
      SUBROUTINE READ_DIAG(NCOL)
!      INCLUDE 'param.inc'
      use felder
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 dum(mg0)
!      REAL*4  OPMAT
      CHARACTER*10 KNAME
      COMMON /SWITCH/ ISWITCH
!      COMMON /OPME/ OPMAT(NKPT,INUME,MG), &
!                    EMINo,EMAXo,OML,OM1,MIMA(NKPT,2),NK,KRA
      COMMON /OPME/  EMINo,EMAXo,OML,OM1,NK,KRA
!
      COMMON /SYM2/   TAU(3,NSYM),IORD,IMAT(3,3,NSYM)
      KKK=0
 1110 READ(9,9011,END=2220) KK,NEMIN,NEMAX,emi,emi,kname
      KKK=KKK+1
      DO  NB1=NEMIN,NEMAX
        READ(9,*)
      enddo
      goto 1110
 2220 continue
      rewind (9)
      write(*,*) 'opmat allocated with kkk,nemax,ncol',kkk,nemax,ncol,'(',kkk*nbindex*ncol/1000000,' MB)'
!      allocate (MIMA(kkk,2),OPMAT(kkk,nemax,MG))
      allocate (MIMA(kkk,2),OPMAT(kkk,nemax,ncol))
       opmat=0.0
      KKK=0
      read(9,*) imist
!cad
 1111 continue
      READ(9,9011,END=2222) KK,NEMIN,NEMAX,emi,emi,kname
!jan
!!      if (ncol.gt.mg) stop 'ncol.gt.MG'
!jan
      KKK=KKK+1
      MIMA(KKK,1)=NEMIN
      MIMA(KKK,2)=NEMAX
!cad
      DO 119 NB1=NEMIN,NEMAX
        READ(9,9920,err=9999) NB1a,NB1a, &
            (dum(I),I=1,NCOL)
      do i=1,ncol
       if (ABS(dum(I)).lt.1.d-20) dum(I)=0.d0
       OPMAT(KKK,NB1,I)=dum(I)
      enddo
  119 CONTINUE 
!.....go for next k-point
      GOTO 1111
!.....end of k-points
 2222 CONTINUE 
      NK=KKK
!     WRITE(6,*) NK, 'k-points read'
      RETURN
 9999 write(*,*) KKK,NB1a,NB1a
      stop ' something wrong in unit 4'
 9011 FORMAT(/,6X,I6,15X,2I5,4X,2f5.2,3X,a10,/)
 9920 FORMAT(3X,2I4,9E13.6)
      END

