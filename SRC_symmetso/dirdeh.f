      subroutine dirdeh(zvec,s,ich,label)
      IMPLICIT REAL*8 (A-H,O-Z)
      logical s(48)
      dimension zvec(3),yvec(3)
      character*2 ich
!
      if(ich.eq.'  ') then
               zvec(1)=0.d0
               zvec(2)=0.d0
               zvec(3)=1.d0
      else if(ich.eq.'my') then
               if(s(19)) then
                  zvec(1)=0.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
                  l=19
               else if(s(16))then
                  zvec(1)=1.d0
                  zvec(2)=0.d0
                  zvec(3)=0.d0
                  l=16
               else if(s(15))then
                  zvec(1)=1.d0
                  zvec(2)=sqrt(3.d0)
                  zvec(3)=0.d0
                  l=15
               else if(s(17))then
                  zvec(1)=1.d0
                  zvec(2)=-sqrt(3.d0)
                  zvec(3)=0.d0
                  l=17
               else if(s(18))then
                  zvec(1)=sqrt(3.d0)
                  zvec(2)=-1.0
                  zvec(3)=0.d0
                  l=18
               else if(s(20))then
                  zvec(1)=sqrt(3.d0)
                  zvec(2)=1.0
                  zvec(3)=0.d0
                  l=20
               else 
                  stop 'm missing '
               endif
      else if(ich.eq.'2y') then
               if(s(8)) then
                  zvec(1)=0.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
                  l=8
               else if(s(5))then
                  zvec(1)=1.d0
                  zvec(2)=0.d0
                  zvec(3)=0.d0
                  l=5
               else if(s(4))then
                  zvec(1)=1.d0
                  zvec(2)=sqrt(3.d0)
                  zvec(3)=0.d0
                  l=4
               else if(s(6))then
                  zvec(1)=1.d0
                  zvec(2)=-sqrt(3.d0)
                  zvec(3)=0.d0
                  l=6
               else if(s(7))then
                  zvec(1)=sqrt(3.d0)
                  zvec(2)=-1.d0
                  zvec(3)=0.d0
                  l=7
               else if(s(9))then
                  zvec(1)=sqrt(3.d0)
                  zvec(2)=1.d0
                  zvec(3)=0.d0
                  l=9
!               else if(s(3))then
!                  zvec(1)=0.d0
!                  zvec(2)=0.d0
!                  zvec(3)=1.d0
               else 
                  stop '2 missing '
               endif

      else if(ich.eq.'mz_my_') then

      else if(ich.eq.'2z_my ') then

      else if(ich.eq.'2z_2y ') then

      else if(ich.eq.'mz') then
!         write(*,*) ' CAUTION with LOC ROT MAT'
               if(s(14)) then
                  zvec(1)=0.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
               else if(s(16))then
                  zvec(1)=1.d0
                  zvec(2)=0.d0
                  zvec(3)=0.d0
               else if(s(15))then
                  zvec(1)=1.d0
                  zvec(2)=sqrt(3.d0)
                  zvec(3)=0.d0
               else if(s(17))then
                  zvec(1)=1.d0
                  zvec(2)=-sqrt(3.d0)
                  zvec(3)=0.d0
               else if(s(18))then
                  zvec(1)=sqrt(3.d0)
                  zvec(2)=-1.0
                  zvec(3)=0.d0
               else if(s(19))then
                  zvec(1)=0.d0
                  zvec(2)=1.0
                  zvec(3)=0.d0
               else if(s(20))then
                  zvec(1)=sqrt(3.d0)
                  zvec(2)=1.0
                  zvec(3)=0.d0
               else 
                  stop 'm missing '
               endif

      else if(ich.eq.'2z') then
               if(s(3)) then
                  zvec(1)=0.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
               else if(s(5))then
                  zvec(1)=1.d0
                  zvec(2)=0.d0
                  zvec(3)=0.d0
               else if(s(4))then
                  zvec(1)=1.d0
                  zvec(2)=sqrt(3.d0)
                  zvec(3)=0.d0
               else if(s(6))then
                  zvec(1)=1.d0
                  zvec(2)=-sqrt(3.d0)
                  zvec(3)=0.d0
               else if(s(7))then
                  zvec(1)=sqrt(3.d0)
                  zvec(2)=-1.d0
                  zvec(3)=0.d0
               else if(s(8))then
                  zvec(1)=0.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
               else if(s(9))then
                  zvec(1)=sqrt(3.d0)
                  zvec(2)=1.d0
                  zvec(3)=0.d0
               else 
                  stop '2 missing '
               endif
      end if
!
      end
