!
!     Constant parameter definition
!
      INTEGER   NCOM, NNPOS1, NRAD, NSYM , PARTMAX
!      INTEGER   NPT00, NPT, NGAUSS, LMAX2
      INTEGER   NGAUSS, LMAX2
      INTEGER   MAXFOLL
!      PARAMETER          (NATO=   6)
!      PARAMETER          (NDIF=   50)
      PARAMETER          (NNPOS1= 5000)
      PARAMETER          (NRAD=  881)
!      PARAMETER          (NPT00= 190) 
!      PARAMETER          (NPT=NRAD+NPT00)  
      PARAMETER          (NSYM=   48)
!      PARAMETER          (NWAV=   20000)
!      PARAMETER          (NXDIM= 10)
!      PARAMETER          (NYDIM= 10)
      PARAMETER		 (NGAUSS= 200)
      PARAMETER          (LMAX2=8)
      PARAMETER          (NCOM=   (lmax2+1)*(lmax2+1))
      PARAMETER          (PARTMAX= 8)
      PARAMETER          (MAXFOLL=12000)
!
!     Convergence terms
      real*8 dmin, drrmin
!     Smallest acceptable step to take
      parameter          ( dmin=1.D-3)
!     Convergence for location of surface in a.u.
      parameter          ( drrmin=5.D-4)
!
!     Trust region parameters
      integer ismall,ilarge
      real*8 ftol
!     Smallest fraction of trust radius to use
      parameter (ismall = 4)
!     Largest multiple of trust radius to use
      parameter (ilarge = 1)
!     Smallest change in predicted f to be viable
      parameter (ftol = 1.D-6)

