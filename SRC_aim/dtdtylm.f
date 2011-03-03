      subroutine dtdtylm(cth,sth,cfi,sfi,xy,lomax,y,dy,ddy)
      implicit none

!.....given the dtylm in y,  returns the second derivative  
!      with respect to theta 
!     Small changes by L. D. Marks, January 2006
      complex*16 y(*),dy(*),ddy(*),efi


      real*8 tsfpi,xy,cth,sth,cfi,sfi
      real*8 pi,tpi,fpi,sqpi,sqtpi,tmp1,tmp2
     
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

!      write(6,*) ':DYTLM xy=',xy
      if (xy.eq.0.d0) then                                         
!........v points in the z direction (special formula for dy 
!........to avoid divergence)
        if (deb) write(6,*) ':DYTLM xy.lt.snull xy=',xy
        do l=0,lomax
          do m=-l,l
             i=l*(l+1)+1+m
!            if ( m.eq.1 .or. m.eq.-1 ) then
!              ddy(i)=-m*sqrt(l*(l+1.)*(2.*l+1.))/tsfpi
!            else
!              ddy(i)=0.
!            endif
             ddy(i)=RTLL(L,M)
          enddo
        enddo
      else
!........v does not point in the z direction
        efi=dcmplx(cfi,-sfi)
        tmp1=cth/sth
        tmp2=1.d0/(sth*sth)
        do l=0,lomax
          do m=-l,l-1
            i = l*(l+1)+1+m
            ddy(i) = m*(tmp1*dy(i)-y(i)*tmp2)+ &
               efi*RTLM(L,M)*dy(i+1)
!            ddy(i) = m*cth/sth*dy(i)-m*y(i)/(sth*sth)+ &
!               efi*RTLM(L,M)*dy(i+1)
!            ddy(i) = m*cth/sth*dy(i)-m*y(i)/(sth*sth)+ &
!               efi*sqrt(l*(l+1.)-m*(m+1.))*dy(i+1)
          enddo
          i = l*(l+1)+1+l
!          ddy(i) = l*(cth/sth*dy(i)-y(i)/(sth*sth))
          ddy(i) = l*(tmp1*dy(i)-y(i)*tmp2)
        enddo
      endif
      return                                                            
      end                                                               

