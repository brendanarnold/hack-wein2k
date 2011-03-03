      subroutine kurki(la,lp,mu,mp,sign,lm,lmmax)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer sign,lm(2,48)
      logical log
!     if sign=-1 and m=0 do not add LM combination
      log=.true.
      isign=sign
!      if(sign.eq.(-2)) then
      if(sign.lt.0) then
        isign=-1
        log=.false.
      endif
      do 10 l1=0,6
      l=(l1*la+lp)*isign
      if(abs(l).gt.6) goto 10
      do 11 m1=0,6
      m=m1*mu+mp
      if(m.gt.abs(l)) goto 10
      lmmax=lmmax+1
      lm(1,lmmax)=l
      lm(2,lmmax)=m
      if(m.eq.0) then
        if(log) then
           lm(1,lmmax)=abs(lm(1,lmmax))
        else
           lmmax=lmmax-1
        endif
      endif
 11   continue
 10   continue
      return
      end
