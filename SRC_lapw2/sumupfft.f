      SUBROUTINE sumupfft(ifft1,ifft2,ifft3,w,fft,sumfft)
!
! Accumulates fft fields
! GM 6/6-00
      IMPLICIT NONE
      INTEGER    :: i1,i2,i3
      INTEGER    :: ifft1,ifft2,ifft3
      REAL*8     :: w
      COMPLEX*16 :: fft(ifft1,ifft2,ifft3),sumfft(ifft1,ifft2,ifft3)
      INTRINSIC DCONJG
      do i3=1,ifft3
         do i2=1,ifft2
            do i1=1,ifft1
               sumfft(i1,i2,i3)=sumfft(i1,i2,i3)+ &
                    w*DBLE(fft(i1,i2,i3)*DCONJG(FFT(I1,I2,I3)))
!                    w*fft(i1,i2,i3)*DCONJG(FFT(I1,I2,I3))
            enddo
         enddo
      enddo
      return
      end


