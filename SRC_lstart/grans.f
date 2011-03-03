      subroutine grans(r,qsin,qcos,dndr,dndt,dndp, &
                       d2ndr,d2ndt,d2ndp,dndrt,dndrp,dndtp, &
                       gx,gy,gz,gmag,g2,ggx,ggy,ggz,gdgg)                
!
      implicit real*8 (a-h,o-z)
!
!
      gx=dndr
      gy=0.
      gz=0.
!
      gmag=ABS(gx)
!
      g2x=d2ndr+2.d0*dndr/r
      g2y=0.
      g2z=0.
      g2=g2x
!
      ggx=dndr*d2ndr
      ggy=0.
      ggz=0.
      ggx=ggx/gmag  
      ggy=ggy/gmag  
      ggz=ggz/gmag  
!
      gdgg=gx*ggx
!
      return
      end
