      subroutine gen_change(change,cth,sth,cfi,sfi)
      implicit none
!
!     Returns in change(i,j) the derivatives d(u_i)/d(x_j)
!     where u_1=r, u_2=theta u_3=phi
!     and   x_1=x, x_2=y, x_3=z
!     Takes care of the special case theta=0
!
      real*8 change,cth,sth,cfi,sfi

      dimension change(3,3)
!.....calculate sines and cosines of the polar angles of the vector v.  

      change(1,1)=sth*cfi
      change(1,2)=sth*sfi
      change(1,3)=cth
      change(2,1)=cth*cfi
      change(2,2)=cth*sfi
      change(2,3)=-sth
      change(3,1)=-sfi
      change(3,2)=cfi
      change(3,3)=0.D0
      return                                                            
      end                                                               

