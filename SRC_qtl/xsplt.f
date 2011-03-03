! Subroutine XSPLT calculates the density matrix of l orbital
! q. number. For k-points coupled to original one by operation
! idet(i)=-1 it includes the effect of time invertion on the
! spherical harminics ( Tinv Y(l,m)=(-1)**m Y(l,-m) )
      SUBROUTINE XSPLT(nemin,nemax,L,isym,iso)    
      USE param
      USE struct
      USE sym2
      USE abc
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16     sum  
      COMMON /RINTEG/  RI_MAT(0:lmax2,nrf,nrf,2,2)
!---------------------------------------------------------------------  
 5     format(6f10.5)
        N=2*L+1
	do 100 num=nemin,nemax
        do 100 is1=1,iso
	do 100 is2=1,iso
        do 100 irf1=1,nrf
        do 100 irf2=1,nrf
        do 100 m1=-l,l
        do 100 m2=-l,l

         ly=L+1+m1
	 lpy=l+1+m2
         ind1=L+1+idet(isym)*m1+N*(is1-1)
         ind2=L+1+idet(isym)*m2+N*(is2-1)
          sum= &
           alm(ly,num,irf1,is1)*conjg(alm(lpy,num,irf2,is2))* &
          ri_mat(L,IRF1,IRF2,IS1,IS2)
       if (idet(isym).lt.0) sum=dconjg(sum)
          xqtl(ind1,ind2,num)=xqtl(ind1,ind2,num)+ &
          idet(isym)**(m1+m2)*sum/IORD
 100   continue
      RETURN                                                            
      END                                                               
