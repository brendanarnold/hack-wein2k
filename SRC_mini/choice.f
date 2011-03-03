      subroutine choice(xa,xb,xe,ya,yb,ye,za,zb,ze,ichoic)
      implicit double precision (a-h,o-z)
      tiny=1.d-6
      tinyd=0.05d0
      if (DABS(-xa+xb).lt.tiny) goto 100
      if (DABS(-ya+yb).lt.tiny) goto 200
      if (DABS(-za+zb).lt.tiny) goto 300
      disx=1+(xa-xe)/(-xa+xb)
      disy=1+(ya-ye)/(-ya+yb)
      disz=1+(za-ze)/(-za+zb)
      write(6,*) 'disx,disy,disz',disx,disy,disz
!      write(23,*) 'disx,disy,disz',disx,disy,disz
      if ((dabs(disx-disy).gt.tinyd).and.(dabs(disx-disz).gt.tinyd) &
             .and.(dabs(disy-disz).gt.tinyd)) then
         ichoic=1
         return
      endif
      if ((dabs(disx-disy).lt.tinyd).and.(dabs(disx-disz).lt.tinyd) &
             .and.(dabs(disy-disz).lt.tinyd)) then
         ichoic=2
         return
      endif
      if (dabs(disx-disy).lt.tinyd) then
         ichoic=6
         return
      endif
      if (dabs(disx-disz).lt.tinyd) then
         ichoic=5
         return
      endif
      if (dabs(disy-disz).lt.tinyd) then
         ichoic=6
         return
      endif
      goto 900 
!
  100 if (DABS(-ya+yb).lt.tiny) goto 150
      if (DABS(-za+zb).lt.tiny) then
         write(6,*)'dabs of x,y,z', &
                   DABS(-xa+xb),DABS(-ya+yb),DABS(-za+zb)
         ichoic=3
      else
         disy=1+(ya-ye)/(-ya+yb)
         disz=1+(za-ze)/(-za+zb)
         write(6,*) 'dabs of x,disy,disz',DABS(-xa+xb),disy,disz
         if (dabs(disy-disz).lt.tiny) ichoic=3
         if (dabs(disy-disz).gt.tiny) ichoic=7
      endif
      return
!
  150 if (DABS(-za+zb).lt.tiny) then
         write(6,*)'dabs of x,y,z', &
         DABS(-xa+xb),DABS(-ya+yb),DABS(-za+zb)
         ichoic=0
      else
         write(6,*)'dabs of x,y,z', &
         DABS(-xa+xb),DABS(-ya+yb),DABS(-za+zb)
         ichoic=4
      endif
      return
!
  200 if (DABS(-za+zb).lt.tiny) then
         write(6,*)'dabs of x,y,z', &
                   DABS(-xa+xb),DABS(-ya+yb),DABS(-za+zb)
         ichoic=2
      else
         disx=1+(xa-xe)/(-xa+xb)
         disz=1+(za-ze)/(-za+zb)
         write(6,*) 'dabs of y,disx,disz',DABS(-ya+yb),disx,disz
         if (dabs(disx-disz).lt.tiny) ichoic=2
         if (dabs(disx-disz).gt.tiny) ichoic=6
      endif
      return
!
  300 disx=1+(xa-xe)/(-xa+xb)
      disy=1+(ya-ye)/(-ya+yb)
      write(6,*) 'dabs of z,disx,disy',DABS(-za+zb),disx,disy
      if (dabs(disx-disy).lt.tiny) ichoic=2
      if (dabs(disx-disy).gt.tiny) ichoic=5
      return
!
!        error handling
!
  900 INFO = 1
!
!        error in selecting ichoic
!
      CALL OUTERR('CHOICE','error in selecting ichoic')
      CALL OUTERR('CHOICE','no option found')
      CALL OUTERR('CHOICE','info=1')
      GOTO 999
  999 STOP 'MINI - Error'
      end
