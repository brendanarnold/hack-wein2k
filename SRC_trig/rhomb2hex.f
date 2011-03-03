      program r2h
	implicit real*8 (a-h,o-z)
!     double precision R(3,3),H(1:3,1:3),ah(3),aho(3),aro(3),
      dimension R(3,3),H(1:3,1:3),ah(3),aho(3),aro(3), &
      RI(3,3)
        DIMENSION AMINV(9),IWORK(3),IWRK(3)
!     integer i,j,l,ii,jj,kk,ll
!     
!.... Changed by P.Blaha
!     Interchange first and second line !!!! 
      r(2,1)=0.8660254037844d0
      r(2,2)=-0.5d0
      r(2,3)=0.d0
      r(1,1)=0.d0
      r(1,2)=1.d0
      r(1,3)=0.d0
!     H(1,1)=0.8660254037844
!     H(1,2)=-0.5
!     H(1,3)=0.
!     H(2,1)=0.
!     H(2,2)=1.
!     H(2,3)=0.
!.... End of changes by P.Blaha
      r(3,1)=0.d0
      r(3,2)=0.d0
      r(3,3)=1.d0
!      
      h(1,1)=0.2886751345948d0
      h(1,2)=-0.5d0
      h(1,3)=0.3333333333333d0
      h(2,1)=0.2886751345948d0
      h(2,2)=0.5d0
      h(2,3)=0.3333333333333d0
      h(3,1)=-.5773502691896d0
      h(3,2)=0.d0
      h(3,3)=0.3333333333333d0
!       DO 712 I=1,3
!       DO 712 J=1,3
!       INDEX=(I-1)*3+J
! 712     AMINV(INDEX)=R(I,J)
!         CALL MINV(AMINV,3,DETR,IWORK,IWRK)
!         DO 713 I=1,3
!         DO 713 J=1,3
!         INDEX=(I-1)*3+J
! 713     RI(I,J)=AMINV(INDEX)
!      
!      do 20 l=1,10
      Print *, 'Please give rx,ry,rz '
      read (*,*)ah(1),ah(2),ah(3)
!      
      aho(1)=ah(1)*H(1,1)+ah(2)*H(2,1)+ah(3)*H(3,1)
      aho(2)=ah(1)*H(1,2)+ah(2)*H(2,2)+ah(3)*H(3,2)
      aho(3)=ah(1)*H(1,3)+ah(2)*H(2,3)+ah(3)*H(3,3)
        write(*,*) 'input  rhomb'
        write(*,100) (ah(j),j=1,3)
!       write(*,*) 'out put hex'
!       write(*,100) (aho(j),j=1,3)
! **********To perfor inversion of matrix******************
      detR=0d0
      do i=1,3
       do j=1,3
        ii=mod(i,3)+1
        jj=mod(j,3)+1
        kk=mod(i+1,3)+1
        ll=mod(j+1,3)+1
        RI(j,i)=(R(ii,jj)*R(kk,ll)-R(ii,ll)*R(kk,jj))
!       write(*,*)RI(j,i),i,j
       end do
       if (i.eq.1) then
        do j=1,3
         detR=detR+RI(j,i)*R(i,j)*(1)**(j+i)
        end do
       endif
       do j=1,3
        RI(j,i)=RI(j,i)/detR
       end do
      end do
! *************************************************
      aro(1)=aho(1)*RI(1,1)+aho(2)*RI(2,1)+aho(3)*RI(3,1)
      aro(2)=aho(1)*RI(1,2)+aho(2)*RI(2,2)+aho(3)*RI(3,2)
      aro(3)=aho(1)*RI(1,3)+aho(2)*RI(2,3)+aho(3)*RI(3,3)
        do j=1,3
        if (aro(j).gt.1.d0) then
        aro(j)=aro(j)-1.d0
        endif
        if (aro(j).lt.0.d0) then
        aro(j)=aro(j)+1.d0
        endif
        end do
!
        write(*,*) 'output hex'
        write(*,100) (aro(j),j=1,3)
!20    continue
100     format(3f12.8)
      stop
      end


