      subroutine con_ev(n_ene,n_kpt, &
        efermi,eigen,nevl)
!
      IMPLICIT REAL*8 (A-H,O-Z)
!      INCLUDE 'param.inc'
!
      dimension eigen(nevl,*),n_ene(*)
      do 10 j=1,n_kpt
      do 10 i=1,n_ene(j)
 10   eigen(i,j)=(eigen(i,j)-efermi)*13.6058d0
      efermi=0.0d0
      return
      end
