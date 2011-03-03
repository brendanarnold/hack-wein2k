      subroutine pgrpnr(grpnam,igrp)
      IMPLICIT REAL*8 (A-H,O-Z)
!
      character*3    grpnam
!--------------------------------------------------------------
       
      igrp=0
      if(grpnam.eq.'C1 ') igrp= 1
      if(grpnam.eq.'Ci ') igrp= 2
      if(grpnam.eq.'C2 ') igrp= 3
      if(grpnam.eq.'Cs ') igrp= 4
!
      if(grpnam.eq.'C2h') igrp= 5
      if(grpnam.eq.'D2 ') igrp= 6
      if(grpnam.eq.'C2v') igrp= 7
      if(grpnam.eq.'D2h') igrp= 8
!
      if(grpnam.eq.'C4 ') igrp= 9
      if(grpnam.eq.'S4 ') igrp=10
      if(grpnam.eq.'C4h') igrp=11
      if(grpnam.eq.'D4 ') igrp=12
!
      if(grpnam.eq.'C4v') igrp=13
      if(grpnam.eq.'D2d') igrp=14
      if(grpnam.eq.'D4h') igrp=15
      if(grpnam.eq.'C3 ') igrp=16
!
      if(grpnam.eq.'C3i') igrp=17
      if(grpnam.eq.'D3 ') igrp=18
      if(grpnam.eq.'C3v') igrp=19
      if(grpnam.eq.'D3d') igrp=20
!
      if(grpnam.eq.'C6 ') igrp=21
      if(grpnam.eq.'C3h') igrp=22
      if(grpnam.eq.'C6h') igrp=23
      if(grpnam.eq.'D6 ') igrp=24
!
      if(grpnam.eq.'C6v') igrp=25
      if(grpnam.eq.'D3h') igrp=26
      if(grpnam.eq.'D6h') igrp=27
      if(grpnam.eq.'T  ') igrp=28
!
      if(grpnam.eq.'Th ') igrp=29
      if(grpnam.eq.'O  ') igrp=30
      if(grpnam.eq.'Td ') igrp=31
      if(grpnam.eq.'Oh ') igrp=32
!
      if(igrp.eq.0) stop 'pgrpnr: cannot find point group nr'
      return
      end
