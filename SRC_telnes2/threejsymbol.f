!BOP
! !ROUTINE: ThreeJSymbol
! !INTERFACE:
    real*8 function ThreeJSymbol(L1, L2, L3, M1, M2, M3)
! !INPUT/OUTPUT PARAMETERS:
!   L1,L2,L3,M1,M2,M3   :   integers specifying the 3j-symbol (L1   L2   L3)
!                                                              M1   M2   M3
! !DESCRIPTION:
!   Calculates Wigner 3j symbols using the Racah formula.
! !REVISION HISTORY:
!   Updated November 2004 (Kevin Jorissen)
!EOP

      implicit none
      integer,intent(in) ::  L1, L2, L3, M1, M2, M3
      integer k,kmin,kmax
      real*8 tjs
      real*8, external :: fact

      
      if ((abs(m1).gt.l1).or.(abs(m2).gt.l2).or.(abs(m3).gt.l3)) then
         ThreeJSymbol = dble(0)
         return
      else
         tjs = dble(0)
         kmin = max0(0, -(l3-l2+m1), -(l3-l1-m2))
         kmax = min0(l1+l2-l3, l1-m1, l2+m2)
         do k = kmin, kmax
            tjs = tjs + ((-1)**k) / fact(k) / fact(l1+l2-l3-k)  &
                 / fact(l1-m1-k) / fact(l2+m2-k) &
                 / fact(l3-l2+m1+k) / fact(l3-l1-m2+k)
         enddo
         
         tjs = tjs * (-1)**(l1-l2-m3)
         tjs = tjs * dsqrt (fact(l1+l2-l3) * fact(l1-l2+l3)  &
              * fact(-l1+l2+l3) / fact(l1+l2+l3+1) &
              * fact(l1+m1) * fact(l1-m1) * fact(l2+m2) * fact(l2-m2) &
              * fact(l3+m3) * fact(l3-m3) )
         
         ThreeJSymbol = tjs
      endif

      return
      end
      
