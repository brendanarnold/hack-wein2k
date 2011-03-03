      subroutine finish( &
       xwert,fwert,ewert,mult,nato,alat,nparam, &
       title,bestr,stopfil)
!
      implicit real*8 (a-h,o-z)
!
      dimension  alat(3)
      dimension  xwert(nato*3), fwert(nato*3)
      dimension  mult(nato)
!
      character*(*) title,bestr
      character*(*) stopfil
!
!.......output of final result:
        write(6,1050)title,bestr
1050    format(///20('=')/5x,a/20('=')// &
        10x,a)
        write(6,'(//6h  #   ,6h  atom,6h   xyz,6h  mult,10h    coord.,10h    (rel.),10h    F*mult,10h    F/atom)')
        write(6,1060)(j,(j-1)/3+1,mod(j-1,3)+1,mult((j-1)/3+1), &
        xwert(j),xwert(j)/alat(mod(j-1,3)+1), &
        fwert(j),fwert(j)/dble(mult((j-1)/3+1)), &
        j=1,nparam)
1060    format(4i6,4f10.5)
!
       write(6,1070)ewert
1070   format(/10x,'minimum energy = ',f20.10,' Ryd ')
!
        open( unit=99,file=stopfil)
!
        return
        end
