      SUBROUTINE SOMM2(r,DP,DPAS,DA,np1,NP)                           
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DP(np),r(np)                                       

      rr1=r(np1)
      g1=dlog(dp(np1)/rr1/rr1)
      da=0.d0
      
      do  ir=np1+1,np
        rr2=r(ir)
        g2=dlog(dp(ir)/rr2/rr2)  
        ac=-(g2-g1)/(rr2-rr1)
        bc=dexp((g1+g2)/2.d0+ac*(rr1+rr2)/2.d0)

        x=rr2
        h=-exp(-ac*x)*(x*x/ac+2.d0*x/ac/ac+2.d0/ac/ac/ac)
        x=rr1
        h=h+exp(-ac*x)*(x*x/ac+2.d0*x/ac/ac+2.d0/ac/ac/ac)
        da=da+h*bc
        rr1=rr2
        g1=g2
      enddo
      return
      end
