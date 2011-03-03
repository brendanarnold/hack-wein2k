      subroutine dirdef(zvec,s,ich,label)
      IMPLICIT REAL*8 (A-H,O-Z)
      logical s(48)
      dimension zvec(3),yvec(3)
      character*2 ich
!
!       write(*,*) ich,s
 1    continue
      if(ich.eq.'  ') then
               zvec(1)=0.d0
               zvec(2)=0.d0
               zvec(3)=1.d0
      else if(ich.eq.'my') then
               if(s(7)) then
                  zvec(1)=0.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
                  l=7
               else if(s(8))then
                  zvec(1)=1.d0
                  zvec(2)=0.d0
                  zvec(3)=0.d0
                  l=8
               else if(s(12))then
                  zvec(1)=1.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
                  l=12
               else if(s(13))then
                  zvec(1)=-1.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
                  l=13
               else if(s(14))then
                  zvec(1)=1.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
                  l=14
               else if(s(15))then
                  zvec(1)=0.d0
                  zvec(2)=1.d0
                  zvec(3)=1.d0
                  l=15
               else if(s(16))then
                  zvec(1)=-1.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
                  l=16
               else if(s(17))then
                  zvec(1)=0.d0
                  zvec(2)=-1.d0
                  zvec(3)=1.d0
                  l=17
               else if(s(6))then
                  zvec(1)=0.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
                  l=6
               else 
                  stop 'm missing '
               endif
      else if(ich.eq.'mz') then
               if(s(6)) then
                  zvec(1)=0.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
               else if(s(14))then
                  zvec(1)=1.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
               else if(s(15))then
                  zvec(1)=0.d0
                  zvec(2)=1.d0
                  zvec(3)=1.d0
               else if(s(16))then
                  zvec(1)=-1.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
               else if(s(17))then
                  zvec(1)=0.d0
                  zvec(2)=-1.d0
                  zvec(3)=1.d0
               else if(s(8))then
                  zvec(1)=1.d0
                  zvec(2)=0.d0
                  zvec(3)=0.d0
               else if(s(12))then
                  zvec(1)=1.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
               else if(s(13))then
                  zvec(1)=-1.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
               else if(s(7))then
                  zvec(1)=0.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
               else 
                  stop 'm missing '
               endif
      else if(ich.eq.'-3') then
               if(s(28)) then
                  zvec(1)=1.d0
                  zvec(2)=1.d0
                  zvec(3)=1.d0
               else if(s(29))then
                  zvec(1)=-1.d0
                  zvec(2)=-1.d0
                  zvec(3)=1.d0
               else if(s(30))then
                  zvec(1)=1.d0
                  zvec(2)=-1.d0
                  zvec(3)=1.d0
               else if(s(31))then
                  zvec(1)=1.d0
                  zvec(2)=-1.d0
                  zvec(3)=1.d0
               else 
                  stop '-3 missing '
               endif

      else if(ich.eq.'3z') then
               if(s(24)) then
                  zvec(1)=1.d0
                  zvec(2)=1.d0
                  zvec(3)=1.d0
               else if(s(25))then
                  zvec(1)=1.d0
                  zvec(2)=1.d0
                  zvec(3)=-1.d0
               else if(s(26))then
                  zvec(1)=1.d0
                  zvec(2)=-1.d0
                  zvec(3)=-1.d0
               else if(s(27))then
                  zvec(1)=1.d0
                  zvec(2)=-1.d0
                  zvec(3)=1.d0
               else 
                  stop '-3 missing '
               endif

      else if(ich.eq.'2y') then
               if(s(4)) then
                  zvec(1)=0.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
                  l=4
               else if(s(3))then
                  zvec(1)=1.d0
                  zvec(2)=0.d0
                  zvec(3)=0.d0
                  l=3
               else if(s(18))then
                  zvec(1)=1.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
                  l=18
               else if(s(19))then
                  zvec(1)=-1.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
                  l=19
               else if(s(20))then
                  zvec(1)=1.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
                  l=20
               else if(s(21))then
                  zvec(1)=0.d0
                  zvec(2)=1.d0
                  zvec(3)=1.d0
                  l=21
               else if(s(22))then
                  zvec(1)=-1.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
                  l=22
               else if(s(23))then
                  zvec(1)=0.d0
                  zvec(2)=-1.d0
                  zvec(3)=1.d0
                  l=23
               else if(s(5))then
                  zvec(1)=0.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
                  l=5
               else 
                  stop '2y missing '
               endif
      else if(ich.eq.'2x') then
               if(s(3))then
                  zvec(1)=1.d0
                  zvec(2)=0.d0
                  zvec(3)=0.d0
                  l=3
               else if(s(4)) then
                  zvec(1)=0.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
                  l=4
               else if(s(18))then
                  zvec(1)=1.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
                  l=18
               else if(s(19))then
                  zvec(1)=-1.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
                  l=19
               else if(s(20))then
                  zvec(1)=1.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
                  l=20
               else if(s(21))then
                  zvec(1)=0.d0
                  zvec(2)=1.d0
                  zvec(3)=1.d0
                  l=21
               else if(s(22))then
                  zvec(1)=-1.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
                  l=22
               else if(s(23))then
                  zvec(1)=0.d0
                  zvec(2)=-1.d0
                  zvec(3)=1.d0
                  l=23
               else if(s(5))then
                  zvec(1)=0.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
                  l=5
               else 
                  stop '2x missing '
               endif
      else if(ich.eq.'2z') then
               if(s(5))then
                  zvec(1)=0.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
               else if(s(3))then
                  zvec(1)=1.d0
                  zvec(2)=0.d0
                  zvec(3)=0.d0
               else if(s(4)) then
                  zvec(1)=0.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
               else if(s(18))then
                  zvec(1)=1.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
               else if(s(19))then
                  zvec(1)=-1.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
               else if(s(20))then
                  zvec(1)=1.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
               else if(s(21))then
                  zvec(1)=0.d0
                  zvec(2)=1.d0
                  zvec(3)=1.d0
               else if(s(22))then
                  zvec(1)=-1.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
               else if(s(23))then
                  zvec(1)=0.d0
                  zvec(2)=-1.d0
                  zvec(3)=1.d0
               else 
                  stop '2 missing '
               endif

      else if(ich.eq.'4z') then
               if(s(11)) then
                  zvec(1)=0.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
               else if(s(9))then
                  zvec(1)=1.d0
                  zvec(2)=0.d0
                  zvec(3)=0.d0
               else if(s(10))then
                  zvec(1)=0.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
               else 
                  stop '4 missing '
               endif

      else if(ich.eq.'-4') then
               if(s(34)) then
                  zvec(1)=0.d0
                  zvec(2)=0.d0
                  zvec(3)=1.d0
               else if(s(32))then
                  zvec(1)=1.d0
                  zvec(2)=0.d0
                  zvec(3)=0.d0
               else if(s(33))then
                  zvec(1)=0.d0
                  zvec(2)=1.d0
                  zvec(3)=0.d0
               else 
                  stop '-4 missing '
               endif
      else
               stop 'wrong ich'
      end if
!
      if(label.eq.0) return
      s(l)=.false.
      label=label-1
      goto 1
      end
