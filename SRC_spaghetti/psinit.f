      subroutine psinit(xoffs,yoffs,xcm,ycm,dlwid)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension      bccol(50,3)
!
      write(11,120) 
      write(11,110) (xoffs+xcm+3)*28.3465,(yoffs+ycm+2)*28.3465
      write(11,130)
      write(11,*) '/CMM { 0.283465 mul} def'
      write(11,*) '/IK { CMM exch CMM exch lineto} def '
      write(11,*) '/HK { CMM exch CMM exch moveto} def '
      write(11,1)
  1   format('/CM { 28.3465 mul} def',/, &
             '/T { CM exch CM exch translate} def',/, &
             '/RL { CM exch CM exch rlineto} def',/, &
             '/RM { CM exch CM exch rmoveto} def',/, &
             '/L { CM exch CM exch lineto} def',/, &
             '/M { CM exch CM exch moveto} def',/, &
             '/BOX {stroke newpath  M ',/, &
             'length 0 rlineto 0 hight rlineto length neg 0 rlineto ',/, &
             'closepath 1 setgray fill 0 setgray } def')
      write(11,100) xoffs,yoffs
!
      call defins(bccol)
      write(11,140) dlwid
 140  format(f5.2,' setlinewidth')
!
      write(11,'(//,a2)',advance='NO') '% '
      call wrtdate(11)
      write(11,2)
  2   format('%Color definitions   ',/, &
             '/bdef {bind def} bind def',/, &
             '/ldef {load def} bind def',/, &
             '/sr /setrgbcolor ldef')
      do i=1,4
        write(11,3) i,bccol( i,1),bccol( i,2),bccol( i,3)
      end do
      do i=10,21
        write(11,4) i,bccol( i,1),bccol( i,2),bccol( i,3)
      end do
      write(11,5)
  3   format('/c0',i1,' {',3f5.2,' sr} bdef')
  4   format('/c', i2,' {',3f5.2,' sr} bdef')
  5   format('%End Color definitions',/)
!
      return

 130  FORMAT('%%Pages: 1',/,'%%EndComments')
 120  FORMAT('%!PS-Adobe-2.0',/,'%%Orientation: Portrait')
 110  FORMAT('%%BoundingBox: 0 0 ',f10.5,f10.5)
 100  FORMAT('newpath',f10.5,' CM',f10.5,' CM translate')
      end
