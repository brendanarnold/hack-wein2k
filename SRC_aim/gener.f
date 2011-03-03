      SUBROUTINE GENER(BR1)
      use atpos
      use crit
      use sphe
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
      DIMENSION BR1(3,3)                                                
      real*8 A,alpha,beta,gamma,gamma1,tmp,toofar
      real*8, allocatable :: Awork(:,:),work(:)
      integer, allocatable :: iwork(:)
      logical ortho
      COMMON /LATT/ A(3),alpha,beta,gamma,gamma1,ortho
!      COMMON /ATPOS/ ATP(3,nnpos),NPOS,NSA,NSB,NSC
!
      toofar=2.D0*(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      tmp = min(a(1),a(2))
      tmp = min(a(3),tmp)
      toofar=max(toofar,tmp*tmp*4.D0)
!
      if(allocated(atp)) deallocate (atp)
      allocate (atp(3,nnpos1))
      allocate (awork(3,nnpos1),iwork(nnpos1),work(nnpos1))
      NPOS=0                                                            
!      JA=-NSA-1                                                         
! 10   JA=JA+1                                                           
!      IF(JA.GT.NSA) GOTO 1                                              
!      JB=-NSB-1                                                         
! 11   JB=JB+1                                                           
!      IF(JB.GT.NSB) GOTO 10                                             
!      JC=-NSC-1                                                         
! 12   JC=JC+1                                                           
!      IF(JC.GT.NSC) GOTO 11                                             
!      if(ortho)write(6,*)'Orthogonal axes'
      DO JA=-NSA,NSA
      DO JB=-NSB,NSB
      DO JC=-NSC,NSC
      NPOS=NPOS+1                                                       
      IF(NPOS.GT.nnpos1) STOP 'NPOS.GT.nnpos1'                                 
      TMP=0
      DO I=1,3                                                        
        ATP(I,NPOS)=BR1(1,I)*JA+BR1(2,I)*JB+BR1(3,I)*JC                   
        TMP=TMP+ATP(I,NPOS)*ATP(I,NPOS)
        AWORK(I,NPOS)=ATP(I,NPOS)
      enddo
      WORK(NPOS)=TMP
      IWORK(NPOS)=NPOS
      if(TMP .GT. TooFar )NPOS=NPOS-1
!      write(6,21)NPOS,(ATP(I,NPOS),I=1,3),JA,JB,JC
!21    format(' Origin ',i3,3F12.6,3i3)
!      GOTO 12
      enddo
      enddo
      enddo                                                           
! 1    CONTINUE                                                          
      WRITE(6,20) NPOS                                          
 20   FORMAT('0NUMBER OF CELL ORIGINS:',I5)
!     Sort the distances in a sloppy way
      call SORTAG (WORK, NPOS, IWORK)
      do j=1,npos
        jj=iwork(j)
        do i=1,3
           atp(i,j)=awork(i,jj)
        enddo
!        write(6,21)j,(atp(i,j),i=1,3)
! 21     format('   Origin ',i3,3F12.5)
      enddo
      deallocate (awork, iwork, work)
      nnpos=npos                             
      ndif=ndat                                                    
      if(allocated(pc)) then
        deallocate (pc,evpc,zpc,pcrb,icpc)
      endif
      allocate( pc(3,NNPOS*NDIF),evpc(3,NNPOS*NDIF), &
         zpc(3,3,NNPOS*NDIF),pcrb(3,NNPOS*NDIF),icpc(NNPOS*NDIF) )
      RETURN                                                            
      END
                                                               
