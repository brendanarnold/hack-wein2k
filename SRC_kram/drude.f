!ad
!ad _____________________________________________________________________
!ad
!ad    written by Robert Abt
!ad    updated by Claudia Ambrosch-Draxl November 1998
!ad _____________________________________________________________________
!ad
!ad  calculates the dielectric tensor
!ad  for given plasma frequency and damping within the Drude model
!ad
!ad      eps(w) = 1 - pl^2 / ( w^2 + i w gamma )
!ad
!ad      pl2 = pl ^ 2  .... squared plasmafrequency
!ad      w             .... frequency
!ad      gamma         .... damping energy
!ad
!ad  Since eps(w) is divergent at w -> 0 it is better to deals with 
!ad  the optical conductivty rather than the dielectric function
!
!            gamma in this notation is not 
!                  the full width at half hight
!                  -> thats gamma/2 !
!
!    eps2 = Im(eps) = pl2 gamma / w / (w^2 + gamma^2)
!
!    eps1 = Re(eps) = 1 - pl2 /(w^2 + gamma^2)
!
!ad
!ad __________________________ IMAGINARY PART __________________________
!ad
      real*8 function eps2 (pl2,gamma,w)

!ad
      real*8 w,gamma,pl2,h
      h= w * ( w*w + gamma*gamma )
      if (h.gt.1.e-9) then
        eps2 = gamma * pl2 / h
      else
        eps2 = 999.999
      endif
      if(w.lt.1) then
!     write(6,99) pl2,gamma,w,eps2
 99   format(' wpl2, gamma,energy,eps2: ',4f12.5)
      endif
      return
      end
!ad
!ad ____________________________________________________________________
!ad
!ad
!ad _____________________________ REAL PART ____________________________
!ad
!ad
      real*8 function eps1 (pl2,gamma,w)
      real*8 w,gamma,pl2
        eps1 = 1.d0 - pl2 / ( w*w + gamma*gamma )
      return
      end
!ad
!ad ____________________________________________________________________
!ad
!ad
!ad _________________________ OPTICAL CONDUCTIVITY _____________________
!ad
!ad
      real*8 function sig1 (pl2,gamma,w)
      real*8 w,gamma,pl2
        sig1 = gamma * pl2 / ( w*w + gamma*gamma )
      return
      end
!ad
!ad ____________________________________________________________________
!ad

