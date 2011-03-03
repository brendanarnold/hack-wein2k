      SUBROUTINE QTL(NEMIN,NEMAX,icase,iso)
        USE param
        USE case
        USE abc
        implicit real*8 (a-h,o-z)

        complex*16 czero,sum

      data czero /(0d0,0d0)/
!----------------------------------
	nn=(2*lcase(icase)+1)*iso
  	do num=nemin,nemax
!       write(85,*) 'NUM=',num
!	write(85,5)((real(xqtl(k,l,num)),l=1,6),k=1,6)
!       write(85,*)
!       write(85,5)((imag(xqtl(k,l,num)),l=1,6),k=1,6)
!5      format(6f8.4,/)
	do i=1,nn
	sum=czero
	do k=1,nn
 	do l=1,nn
 	sum=sum+dconjg(cf(i,k,icase))*xqtl(k,l,num)*cf(i,l,icase)
 	end do
	end do
	xden(i,num)=dble(sum)
	end do
	end do

	end
	
	
