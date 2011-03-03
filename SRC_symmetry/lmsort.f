      subroutine lmsort(lm,lmmax)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension lm(2,49),lm1(2,49),lm2(4)
      logical minus,plus
      lm1(1,1)=0
      lm1(2,1)=0
      j=2
      do 10 l=1,6
         do 20 m=0,l
            plus=.false.
            minus=.false.
            do 30 i=2,lmmax
               if (iabs(lm(1,i)).eq.l.and.lm(2,i).eq.m) then
                  if (lm(1,i).eq.l) then
                     plus=.true.
                     lm2(1)=lm(1,i)
                     lm2(2)=lm(2,i)
                  else 
                     minus=.true.
                     lm2(3)=lm(1,i)
                     lm2(4)=lm(2,i)
                  endif
               endif
  30        continue
            if (plus) then
               lm1(1,j)=lm2(1)
               lm1(2,j)=lm2(2)
               j=j+1
            endif
            if (minus) then
               lm1(1,j)=lm2(3)
               lm1(2,j)=lm2(4)
               j=j+1
            endif
  20     continue
  10  continue
      do 40 j=2,lmmax
         lm(1,j)=lm1(1,j)
         lm(2,j)=lm1(2,j)
  40  continue
      return
      end
