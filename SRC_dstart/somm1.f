SUBROUTINE SOMM1(r0,DP,DPAS,DA,np1,NP)                           
  IMPLICIT REAL*8 (A-H,O-Z)
  DIMENSION DP(np)                                       
  logical ent
  da=dp(np1)*exp((np1-1)*dpas)
  da=da+dp(np)*exp((np-1)*dpas)
  ent=.true.
  do i=np1+1,np-1
     fak=2.
     if (ent) fak=4.
     ent=.not.ent
     da=da+fak*dp(i)*exp((i-1)*dpas)
  enddo
  da=da*r0*dpas/3.
  RETURN                                                            
END SUBROUTINE SOMM1
