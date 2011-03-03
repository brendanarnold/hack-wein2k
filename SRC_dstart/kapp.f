

         subroutine kapp(akk,akap,rmt,h1)


         implicit real*8 (a-h,o-z)

         dimension akap(4)
         logical ent

         f1(xx)=xx*xx*exp(-aa*xx)/(eps+exp(-aa*xx))
         f2(xx)=xx*exp(-aa*xx)/(eps+exp(-aa*xx))*sin(akk*xx)


         npkt=40
         
         aa=akap(1)
         eps=akap(3)


         if (akk.lt.1.e-6) then

            dx0=1./aa/npkt
            
            iptx=int(rmt/dx0/2.+50)
	    iptx=2*iptx
            dx=rmt/iptx
            
            da=f1(0.d0)+f1(rmt)
            ent=.true.
            do  ir=1,iptx-1
               x=rmt*ir/iptx   
               w=2
               if (ent) w=4
               ent=.not.ent
               da=da+f1(x)*w
            enddo
            h1=dx*da*akap(2)/3.
            akap(4)=h1
!cccc            akap(4)=h1*2.
!           write(6,100) h1*2.
 100        format(3x,' ladung in der kappe:  ',f10.5)
            return
         else
            aaf=dmax1(akk,aa)
            dx0=1./aaf/npkt
            
            iptx=int(rmt/dx0/2.+50)
	    iptx=2*iptx
            dx=rmt/iptx
            
            da=f2(0.d0)+f2(rmt)
            ent=.true.
            do  ir=1,iptx-1
               x=rmt*ir/iptx   
               w=2
               if (ent) w=4
               ent=.not.ent
               da=da+f2(x)*w
            enddo
            h1=dx*da*akap(2)/3.
	    endif
            return
            end
