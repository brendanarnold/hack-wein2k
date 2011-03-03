!
!     Constant parameter definition
!
      INTEGER        MG, MF
!     mg: usually 21, but eq. ((lxdos+1)**2)*((Lxdos+1)**2 + 1)
!      PARAMETER          (MG=      273)
      PARAMETER          (MG=      41)
      PARAMETER          (MF=      MG)
!     For "cross-DOS" set lxdos to 3
      parameter (lxdos=1)
      parameter (lxdos2=(lxdos+1)**2)
