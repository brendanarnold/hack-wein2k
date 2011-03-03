      subroutine defins(bccol)
      use irr_param
      IMPLICIT REAL*8 (A-H,O-Z)
!      INCLUDE 'param.inc'
!
      dimension        bccol(50,3)
      character*40     sfont(7)
      common /prtft/   sfont
!
!-------------------------------------------------------------------
!.....line width definitions
!     dlwid=line width, default
!     dawid=line width of axis
!     dfwid=line width of Fermi level
!
      if(dlwid.le.0d0) dlwid=1.d0
      dawid=dlwid*1.0
      dfwid=dlwid*0.3     
!
!
!.....font definitions:
!     sfont(1): strnam = text    
!           2          = symbols
!           3          = fermi 
!           4          = labels 
!           5          = header      
!
      if(ifont.eq.0) then
        do i=1,5
        sfont(i) = '                                        ' 
        end do
      elseif(ifont.eq.1) then
        sfont(1) = '/Times-Roman                            '
        sfont(2) = '/Symbol                                 '
        sfont(3) = '/Times-Roman                            '
        sfont(4) = '/Times-Roman                            '
        sfont(5) = '/Times-Roman                            '
      elseif(ifont.eq.2) then
        sfont(1) = '/Times-Roman                            '
        sfont(2) = '/Symbol                                 '
        sfont(3) = '/Times-Roman-Italic                     '
        sfont(4) = '/Times-Roman                            '
        sfont(5) = '/Times-Roman                            '
      elseif(ifont.eq.3) then
        sfont(1) = '/Helvetica /ISOLatin1Encoding           '
        sfont(2) = '/Symbol                                 '
        sfont(3) = '/Helvetica-Italic                       '
        sfont(4) = '/Helvetica                              '
        sfont(5) = '/Helvetica                              ' 
      else
         stop 'error in input: iprtf shall be <= 3'
      endif
!
!
!
!.....color definition 
!     bccol(    1,rbg) = used for one-color   plot
!     bccol( 2- 4,rbg) = used for three-color plot
!     bccol(10-22,rgb) = used for multi-color plot
!                        and for separating irr. representations     
!
!.....1-color plot
      bccol( 1,1) = 0.0
      bccol( 1,2) = 0.0
      bccol( 1,3) = 1.0
!
!.....3-color plot
      bccol( 2,1) = 0.0
      bccol( 2,2) = 0.0
      bccol( 2,3) = 1.0
      bccol( 3,1) = 1.0
      bccol( 3,2) = 0.0
      bccol( 3,3) = 0.0
      bccol( 4,1) = 0.0
      bccol( 4,2) = 1.0
      bccol( 4,3) = 0.0
!
!.....multi-color plot
      bccol(10,1) = 0.2
      bccol(10,2) = 0.2
      bccol(10,3) = 0.2
      bccol(11,1) = 0.0
      bccol(11,2) = 0.0
      bccol(11,3) = 1.0
      bccol(12,1) = 1.0
      bccol(12,2) = 0.0
      bccol(12,3) = 0.0
!
      bccol(13,1) = 0.0
      bccol(13,2) = 0.9
      bccol(13,3) = 0.0
      bccol(14,1) = 0.6
      bccol(14,2) = 0.6
      bccol(14,3) = 0.0
      bccol(15,1) = 0.7
      bccol(15,2) = 0.7
      bccol(15,3) = 0.7
!
      bccol(16,1) = 1.0
      bccol(16,2) = 0.0
      bccol(16,3) = 1.0
      bccol(17,1) = 0.9
      bccol(17,2) = 0.9
      bccol(17,3) = 0.0
      bccol(18,1) = 0.0
      bccol(18,2) = 0.9
      bccol(18,3) = 0.9
!
      bccol(19,1) = 0.0
      bccol(19,2) = 0.0
      bccol(19,3) = 0.0
      bccol(20,1) = 1.0
      bccol(20,2) = 0.5
      bccol(20,3) = 0.0
      bccol(21,1) = 0.0
      bccol(21,2) = 1.0
      bccol(21,3) = 0.0
!
!
      end   
