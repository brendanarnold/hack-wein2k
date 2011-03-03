!----------------------------------------------------------------------
!
!     READ OPTICAL MATRIXELEMENTS
!
!----------------------------------------------------------------------
      SUBROUTINE READOPMAT(NCOL,NEOCC)
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
      nbindex=0
 1110 READ(3,9011,END=2220) KK,NEMIN,NEMAX,emi,emi,kname
      NEOCC=min(NEOCC,NEMAX-1)      
!fb      nbindex=max(nbindex,(1+nemax-nemin)*(nemax-nemin+2)/2)
!      nbindex=max(nbindex,((1+nemax-(nemin+neocc)/2)*(neocc-nemin+1))
      nbindex=max(nbindex,((1+nemax)*(neocc-nemin+1)-(nemin+neocc)*(neocc-nemin+1)/2))
!print*,nbindex,nemin,neocc,nemax
      KKK=KKK+1
      DO  NB1=NEMIN,NEMAX
        DO  NB2=NB1,NEMAX
        READ(3,*) 
        enddo
      enddo
      goto 1110
 2220 rewind (3)
!fb
      write(*,*) 'opmat allocated with kkk,nbindex,ncol',kkk,nbindex,ncol,'(',kkk*nbindex*ncol/1000000,' MB)'
!      allocate (MIMA(kkk,2),OPMAT(kkk,nbindex,MG))
      allocate (MIMA(kkk,2),OPMAT(kkk,nbindex,ncol))
       opmat=0.0
      KKK=0
      read(3,*) imist
!ad
 1111 continue
      READ(3,9011,END=2222) KK,NEMIN,NEMAX,emi,emi,kname
!         write(*,*) kkk,kk,nemin,NEMAX,emi,emi,kname
!!      if(nemax.gt.NUME) then
!!        write(*,*) ' NUME',NUME,' .lt. nemax=',nemax,' at k:',kkk+1
!!        stop
!!      end if
!jan
!!!      if (ncol.gt.mg) stop 'ncol.gt.MG'
!jan
      KKK=KKK+1
      MIMA(KKK,1)=NEMIN
      MIMA(KKK,2)=NEMAX
      NBINDEX=0          
!fb      DO 119 NB1=NEMIN,NEMAX
      DO 119 NB1=NEMIN,NEOCC
        DO 119 NB2=NB1,NEMAX
          NBINDEX=NBINDEX+1    
!          IF (NBINDEX.GT.INUME) THEN                     
!             WRITE(6,78) 
! 78   format(/,' parameter INUME set too small!' / &
!            'check paramter and band indices in *.injoint ')
!             STOP                                         
!          END IF                                        
        READ(3,9920,err=9999) NB1a,NB2a, &
            (dum(I),I=1,NCOL)
!    &      (OPMAT(KKK,NBINDEX,I),I=1,NCOL)
      do i=1,ncol
       if (ABS(dum(I)).lt.1.d-20) dum(I)=0.d0
       OPMAT(KKK,NBINDEX,I)=dum(I)
      enddo
  119 CONTINUE 
!fb
      DO  NB1=NEOCC+1,NEMAX
        DO  NB2=NB1,NEMAX
        READ(3,*) 
        enddo
      enddo
!.....go for next k-point
      GOTO 1111
!.....end of k-points
 2222 CONTINUE 
      NK=KKK
!     WRITE(6,*) NK, 'k-points read'
      RETURN
 9999 write(*,*) KKK,NB1a,NB2a
      stop ' something wrong in unit 3'
 9011 FORMAT(/,6X,I6,15X,2I5,4X,2f5.2,3X,a10,/)
 9920 FORMAT(3X,2I4,9E13.6)
      END

