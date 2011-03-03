 subroutine fiteta(qm,fi,teta)

   real*8   qm(3),fi,teta
   real*8   mx,my,mz,pi,azero,teta_old,fi_old
   

   mz=qm(3)
   mx=qm(1)
   my=qm(2)

   pi=acos(-1.d0)   
   azero=1.0d-10

   if ((abs(mx).lt.azero).and.(abs(my).lt.azero)) then
      if (mz.gt.-azero) then
         fi=0.0d0
         teta=0.0d0
      else
         fi=0.0d0
         teta=pi
      endif
   else 
      
      if (abs(mx).gt.azero) then
         fi=atan(abs(my/mx))
         fi_old=fi
         if ((my.ge.0).and.(mx.ge.0)) fi=fi
         if ((my.ge.0).and.(mx.lt.0)) fi=pi-fi
         if ((my.lt.0).and.(mx.lt.0)) fi=pi+fi
         if ((my.lt.0).and.(mx.ge.0)) fi=2.0d0*pi-fi
      else
         fi=pi/2.0d0 
         fi_old=fi
         if (my.gt.0) fi=fi
         if (my.lt.0) fi=pi+fi
      endif
      
      if (abs(mz).gt.azero) then
         teta=atan(sqrt(mx**2+my**2)/abs(mz))
         teta_old=teta
         if (mz.gt.0.0d0) teta=teta
         if (mz.lt.0.0d0) teta=pi-teta
      else
         teta=pi/2.d0
         teta_old=teta
      endif
      
   endif

 end subroutine fiteta
