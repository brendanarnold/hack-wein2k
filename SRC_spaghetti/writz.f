      subroutine writz(ymin1,sizez,yval1,lyflag)
      IMPLICIT REAL*8 (A-H,O-Z)
      call movet (-1.5d0,yval1-sizez*0.2d0)
      if (lyflag.eq.1) then
	 write(11,100) sizez*20.d0,ymin1
      else
	 write(11,110) sizez*20.d0,ymin1
      endif
      return
 110  format('/Times-Roman findfont',f10.5,' scalefont setfont',/, &
      '(',f5.1,') show')
 100  format('/Times-Roman findfont',f10.5,' scalefont setfont',/, &
      '(',f6.3,') show')
      end

