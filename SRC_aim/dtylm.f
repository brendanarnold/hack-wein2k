      subroutine dtylm(cth,sth,cfi,sfi,xy,lomax,y,dy)
      implicit none

!.....given the ylm in y, returns the derivative  
!      with respect to theta                                         
!     Small changes by L. D. Marks, January 2006
      complex*16 y(*),dy(*),efi

      real*8 tsfpi,xy,cth,sth,cfi,sfi
      real*8 pi,tpi,fpi,sqpi,sqtpi,tmp1
     
      integer l,m,i,lomax

      logical deb

      common /DEBUG/ deb
      COMMON /CTES/ pi,tpi,fpi,sqpi,sqtpi
      include 'param.inc'
      real*8 RTLM(0:LMAX2,-LMAX2:LMAX2)
      real*8 RTLL(0:LMAX2,-LMAX2:LMAX2)
      logical GetRTLM
!      COMMON /RTPARS/RTLM,RTLL,GetRTLM
      DATA GetRTLM/.true./
      save rtlm,rtll,getrtlm
!.....calculate sines and cosines of the polar angles of the vector v.  
      tsfpi=4.D0*sqpi
!
!     Setup terms if needed
      if(GetRTLM)then
         do L=0,LMAX2
            DO M=-L,L
               IF (  m.eq.1 .or. m.eq.-1 ) then
                  RTLL(L,M)=-m*sqrt(L*(L+1.D0)*(2.D0*L+1.D0))/tsfpi
               else
                  RTLL(L,M)=0.D0
               endif
               RTLM(L,M)=sqrt(L*(L+1.D0)-m*(m+1.D0))
            ENDDO
         ENDDO
         GetRTLM=.false.
      ENDIF
!
      if (xy.eq.0.d0) then
!........v points in the z direction (special formula for dy 
!........to avoid divergence)
        if(deb) write(6,*) ':DYTLM xy.lt.snull xy=',xy
        do l=0,lomax
          do m=-l,l
             i=l*(l+1)+1+m
!            if ( m.eq.1 .or. m.eq.-1 ) then
!              dy(i)=-m*sqrt(l*(l+1.)*(2.*l+1.))/tsfpi
!            else
!              dy(i)=0.
!            endif
             DY(I)=RTLL(L,M)
          enddo
        enddo
      else
!........v does not point in the z direction
        efi=dcmplx(cfi,-sfi)
        tmp1=cth/sth
        do l=0,lomax
          do m=-l,l-1
            i=l*(l+1)+1+m
!            dy(i)=m*cth/sth*y(i)+efi*RTLM(L,M)*y(i+1)
            dy(i)=m*tmp1*y(i)+efi*RTLM(L,M)*y(i+1)
          enddo
          i=l*(l+1)+1+l
!          dy(i)=l*cth/sth*y(i)
          dy(i)=l*tmp1*y(i)
        enddo
      endif
      return                                                            
      end                                                               

