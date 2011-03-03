      subroutine forcint(pi,br2,brat,alat,alpha,lattic)
!ad
      implicit real*8 (a-h,o-z)
      logical ORTHO
      COMMON  /ORTH/     ORTHO
      CHARACTER*4        LATTIC
      dimension br2(3,3),brat(3,3),alat(3),alpha(3)
      do i=1,3
      do j=1,3
      brat(i,j)=0.0d0
      brat(i,j)=br2(i,j)*alat(j)/2/PI
      enddo
      enddo
!
      test=abs(pi/2.d0-alpha(3))
      if(lattic(1:1).eq.'C'.and.(test.gt.1.d-6)) then
      brat(1,1)=1.d0/sin(alpha(3))
      brat(1,2)=0.d0
      brat(1,3)=0.0d0
      brat(2,1)=-cos(alpha(3))/sin(alpha(3))
      brat(2,2)=1.0d0
      brat(2,3)=0.0d0
      brat(3,1)=0.0d0
      brat(3,2)=0.0d0
      brat(3,3)=1.0d0
      endif
!ad
      if(lattic(1:1).eq.'H') then
      brat(1,1)=2.0d0/sqrt(3.0d0)
      brat(1,2)=1.0d0/sqrt(3.0d0)
      brat(1,3)=0.0D0
      brat(2,1)=0.0D0
      brat(2,2)=1.0d0 
      brat(2,3)=0.0D0
      brat(3,1)=0.0D0
      brat(3,2)=0.0D0
      brat(3,3)=1.0d0 
      endif
      if(lattic(1:1).eq.'R') then
      brat(1,1)=1.D0/SQRT(3.D0)
      brat(1,2)=1.d0/sqrt(3.d0)
      brat(1,3)=-2.d0/sqrt(3.d0)
      brat(2,1)=-1.0d0
      brat(2,2)=1.0d0 
      brat(2,3)=0.0D0
      brat(3,1)=1.0D0
      brat(3,2)=1.0D0
      brat(3,3)=1.0d0 
      endif
!ad
      if(ORTHO) then
      brat(1,1)=1.0d0
      brat(1,2)=0.0d0
      brat(1,3)=0.0D0
      brat(2,1)=0.0D0
      brat(2,2)=1.0d0
      brat(2,3)=0.0D0
      brat(3,1)=0.0D0
      brat(3,2)=0.0D0
      brat(3,3)=1.0d0
      endif

!     write(*,*) 'BRAVAIS MATRIX:'
!     write(*,777) brat
!777  format(3(3f12.6,/))
!ad
      return
      end
