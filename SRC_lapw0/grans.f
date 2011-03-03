      subroutine grans(r,qsin,qcos,dndr,dndt,dndp, &
                       d2ndr,d2ndt,d2ndp,dndrt,dndrp,dndtp, &
                       gx,gy,gz,gmag,g2,ggx,ggy,ggz,gdgg,ir)                
!
      implicit real*8 (a-h,o-z)
!
      r2=r*r
      qsin2=qsin*qsin
!
      gx=dndr
      gy=dndt/r
      gz=dndp/r/qsin
!
      gmag=dsqrt(gx*gx+gy*gy+gz*gz)
!
      g2x=d2ndr+2.d0*dndr/r
      g2y=(d2ndt+qcos/qsin*dndt)/r2
      g2z=d2ndp/r2/qsin2
      g2=g2x+g2y+g2z
!      if(ir.eq.4)write(6,'(8e15.5)') g2,g2x,g2y,g2z,d2ndr,2.d0*dndr/r
!     *, dndr,r
!
      ggx=dndr*d2ndr+ &
          dndt*(dndrt-dndt/r)/r2+ &
          dndp*(dndrp-dndp/r)/r2/qsin2
      ggy=(dndr*dndrt+ &
           dndt*d2ndt/r2+ &
           dndp*(dndtp-qcos/qsin*dndp)/r2/qsin2)/r
      ggz=(dndr*dndrp+ &
           dndt*dndtp/r2+ &
           dndp*d2ndp/r2/qsin2)/r/qsin
      ggx=ggx/gmag  
      ggy=ggy/gmag  
      ggz=ggz/gmag  
!
      gdgg=gx*ggx+gy*ggy+gz*ggz
!
      return
      end
