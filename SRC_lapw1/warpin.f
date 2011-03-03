      COMPLEX*16 FUNCTION WARPIN(X,Y,Z)
      use out, only      : IORD, NKK, KZZ1, KKK, IMAT, TAU, TAUP, WARP, WARP1, WARP2, WARP3
      IMPLICIT NONE
      INTEGER            X,  Y,  Z
!
      INCLUDE 'param.inc'
!
!        Common blocks
!
!
!        IMAT(1:3,1:3,j) - matrix representation of (space group)
!                          symmetry operation j
!        IORD            - number of symmetry operations of space group
!        KKK(1:3,j)      - star vectors of reciprocal lattice
!                          vector KZZ(:)
!        KZZ1(1:3)       - reciprocal lattice vector
!        NKK             - number of vectors in KKK (size of the star)
!        TAU(1:3,j)      - non-primitive translation vector for symmetry
!                          operation j
!        TAUP(j)         - tauphase factor of vector KKK(:,j)
!        WARP(i,j,k)     - interstitial warpin for K vector (i,j,k)
!
!      INTEGER            IORD, NKK
!      INTEGER            KZZ1(3), KKK(3,NSYM), IMAT(3,3,NSYM)
!      DOUBLE PRECISION   TAU(3,NSYM)
!      COMPLEX*16         TAUP(NSYM)
!      COMPLEX*16         WARP(-KMAX1:KMAX1,-KMAX2:KMAX2,-KMAX3:KMAX3)
!      COMMON  /OUT/      WARP, TAU, TAUP, IORD, IMAT, KZZ1, KKK, NKK
!      SAVE    /OUT/
!
      IF ((ABS(X) .LE. WARP1) .AND. &
          (ABS(Y) .LE. WARP2) .AND. &
          (ABS(Z) .LE. WARP3)) THEN
         WARPIN = WARP(X,Y,Z)
      ELSE
         WARPIN = (0.0D+0,0.0D+0)
      ENDIF
!
      RETURN
!
!        End of 'WARPIN'
!
      END
