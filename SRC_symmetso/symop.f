      subroutine symop(i,x,y,z,x1,y1,z1,oname,nz,lattic) 
      IMPLICIT REAL*8 (A-H,O-Z)
      character*20 oname 
      character*4 lattic
!      
      if(lattic(1:1).eq.'H') goto 1000
      if(lattic(1:1).eq.'R') goto 1000
!      
      if(i.eq.1) goto 10
      if(i.eq.2) goto 20
      if(i.eq.3) goto 30
      if(i.eq.4) goto 40
      if(i.eq.5) goto 50
      if(i.eq.6) goto 60
      if(i.eq.7) goto 70
      if(i.eq.8) goto 80
      if(i.eq.9) goto 90
      if(i.eq.10) goto 100
      if(i.eq.11) goto 110
      if(i.eq.12) goto 120
      if(i.eq.13) goto 130
      if(i.eq.14) goto 140
      if(i.eq.15) goto 150
      if(i.eq.16) goto 160
      if(i.eq.17) goto 170
      if(i.eq.18) goto 180
      if(i.eq.19) goto 190
      if(i.eq.20) goto 200
      if(i.eq.21) goto 210
      if(i.eq.22) goto 220
      if(i.eq.23) goto 230
      if(i.eq.24) goto 240
      if(i.eq.25) goto 250
      if(i.eq.26) goto 260
      if(i.eq.27) goto 270
      if(i.eq.28) goto 280
      if(i.eq.29) goto 290
      if(i.eq.30) goto 300
      if(i.eq.31) goto 310
      if(i.eq.32) goto 320
      if(i.eq.33) goto 330
      if(i.eq.34) goto 340
      if(i.eq.35) goto 350
      if(i.eq.36) goto 360
      if(i.eq.37) goto 370
      if(i.eq.38) goto 380
      if(i.eq.39) goto 390
      if(i.eq.40) goto 400
      if(i.eq.41) goto 410
      if(i.eq.42) goto 420
      if(i.eq.43) goto 430
      if(i.eq.44) goto 440
      if(i.eq.45) goto 450
      if(i.eq.46) goto 460
      if(i.eq.47) goto 470
      if(i.eq.48) goto 480
!
!     Identitaet
!
 10   oname='1'
      x1= x
      y1= y
      z1= z
      return 
!
!     Inversion 
!
 20   oname='-1'
      x1=-x
      y1=-y
      z1=-z
      return 
!
!     2-zaehlige Drehachsen  
!
  30  oname='2 || x'
      x1= x
      y1=-y
      z1=-z 
      return
  40  oname='2 || y'
      x1=-x
      y1= y
      z1=-z
      return
  50  oname='2 || z'
      x1=-x
      y1=-y
      z1= z                      
      return
 180  oname='2 || 110'
      x1= y
      y1= x
      z1= -z
      return
 190  oname='2 || -110'
      x1=-y
      y1=-x
      z1= -z
      return
 200  oname='2 || 101'
      x1= z
      y1=-y
      z1= x 
      return
 210  oname='2 || 011'
      x1=-x
      y1= z
      z1= y 
      return
 220  oname='2 || -101'
      x1=-z
      y1=-y
      z1=-x 
      return
 230  oname='2 || 0-11'
      x1=-x
      y1=-z
      z1=-y 
      return
!
!     4-zaehlige Drehachsen
!
  90  oname='4 || x'
      x1= x
      y1= z
      z1=-y 
      nz=3 
      return 
 100  oname='4 || y'
      x1=-z
      y1= y
      z1= x
      nz=3
      return 
 110  oname='4 || z'
      x1= y
      y1=-x
      z1= z
      nz=3
      return 
 350  oname='4 || x'
      x1= x
      y1=-z
      z1= y 
      nz=3 
      return 
 360  oname='4 || y'
      x1= z
      y1= y
      z1=-x
      nz=3
      return 
 370  oname='4 || z'
      x1=-y
      y1= x
      z1= z
      nz=3
      return 
!
!
!     Spiegelebenen
!
  60  oname='m n z'
      x1= x
      y1= y
      z1=-z
      return
  70  oname='m n y'
      x1= x
      y1=-y
      z1= z
      return
  80  oname='m n x'
      x1=-x
      y1= y
      z1= z
      return
 120  oname='m n 110'
      x1=-y
      y1=-x
      z1= z
      return
 130  oname='m n -110'
      x1= y
      y1= x
      z1= z
      return
 140  oname='m n 101'
      x1=-z
      y1= y
      z1=-x
      return
 150  oname='m n 011'
      x1= x
      y1=-z
      z1=-y
      return
 160  oname='m n -101'
      x1= z
      y1= y
      z1= x
      return
 170  oname='m n 0-11'
      x1= x
      y1= z
      z1= y
      return
!
!    3 zaehlige Achsen
!
 240  oname='3 || 111'
      x1= z
      y1= x
      z1= y 
      nz=2
      return
 250  oname='3 || -1-1-1'
      x1=-z
      y1= x
      z1=-y 
      nz=2
      return
 260  oname='3 || 1-1-1'
      x1=-z
      y1=-x
      z1= y 
      nz=2
      return
 270  oname='3 || 1-11'
      x1= z
      y1=-x
      z1=-y 
      nz=2
      return
 380  oname='3 || 111'
      x1= y
      y1= z
      z1= x 
      nz=2
      return
 390  oname='3 || -1-1-1'
      x1= y
      y1=-z
      z1=-x 
      nz=2
      return
 400  oname='3 || 1-1-1'
      x1=-y
      y1= z
      z1=-x 
      nz=2
      return
 410  oname='3 || 11-1'
      x1=-y
      y1=-z
      z1= x 
      nz=2
      return
!
!    S6 zaehlige Achsen
!
 280  oname='S6 || 111'
      x1=-z
      y1=-x
      z1=-y 
      nz=5
      return
 290  oname='S6 || -1-11'
      x1= z
      y1=-x
      z1= y 
      nz=5
      return
 300  oname='S6 || 1-1-1'
      x1= z
      y1= x
      z1=-y 
      nz=5
      return
 310  oname='S6 || 1-11'
      x1=-z
      y1= x
      z1= y 
      nz=5
      return
 420  oname='S6 || 111'
      x1=-y
      y1=-z
      z1=-x 
      nz=5
      return
 430  oname='S6 || -1-11'
      x1=-y
      y1= z
      z1= x 
      nz=5
      return
 440  oname='S6 || 1-1-1'
      x1= y
      y1=-z
      z1= x 
      nz=5
      return
 450  oname='S6 || 1-11'
      x1= y
      y1= z
      z1=-x 
      nz=5
      return
!
!    S4-zaehlige Drehachsen
!
 320  oname='S4 || x'
      x1=-x
      y1= z
      z1=-y 
      nz=3 
      return 
 330  oname='S4 || y'
      x1=-z
      y1=-y
      z1= x
      nz=3
      return 
 340  oname='S4 || z'
      x1= y
      y1=-x
      z1=-z
      nz=3
      return 
 460  oname='S4 || x'
      x1=-x
      y1=-z
      z1= y 
      nz=3 
      return 
 470  oname='S4 || y'
      x1= z
      y1=-y
      z1=-x
      nz=3
      return 
 480  oname='S4 || z'
      x1=-y
      y1= x
      z1=-z
      nz=3
      return 
!...........................hex part .........................
 1000 continue
      if(i.eq.1) goto 11
      if(i.eq.2) goto 21
      if(i.eq.3) goto 31
      if(i.eq.4) goto 41
      if(i.eq.5) goto 51
      if(i.eq.6) goto 61
      if(i.eq.7) goto 71
      if(i.eq.8) goto 81
      if(i.eq.9) goto 91
      if(i.eq.10) goto 101
      if(i.eq.11) goto 111
      if(i.eq.12) goto 121
      if(i.eq.13) goto 131
      if(i.eq.14) goto 141
      if(i.eq.15) goto 151
      if(i.eq.16) goto 161
      if(i.eq.17) goto 171
      if(i.eq.18) goto 181
      if(i.eq.19) goto 191
      if(i.eq.20) goto 201
      if(i.eq.21) goto 211
      if(i.eq.22) goto 221
      if(i.eq.23) goto 231
      if(i.eq.24) goto 241
      return
!
!     Identitaet
!
 11   oname='1'
      x1= x
      y1= y
      z1= z
      return 
!
!     Inversion 
!
 21   oname='-1'
      x1=-x
      y1=-y
      z1=-z
      return 
!
!     2-zaehlige Drehachsen  
!
  31  oname='2 || z'
      x1=-x
      y1=-y
      z1= z 
      return
  41  oname='2 || 120'
      x1=-x+y
      y1= y
      z1=-z
      return
  51  oname='2 || 210'
      x1= x
      y1= x-y
      z1=-z                      
      return
  61  oname='2 || 1-10'
      x1=-y
      y1=-x
      z1=-z
      return
  71  oname='2 || 100'
      x1= x-y
      y1=-y
      z1=-z
      return
  81  oname='2 || 010'
      x1=-x
      y1=-x+y
      z1=-z 
      return
  91  oname='2 || 110'
      x1= y
      y1= x
      z1=-z 
      return
!
!     6-zaehlige Drehachsen
!
 101  oname='6 || z'
      x1= x-y
      y1= x
      z1= z 
      nz=5
      return 
 111  oname='6 || -z'
      x1= y
      y1=-x+y
      z1= z
      nz=5
      return 
 121  oname='S6 (-3) || z'
      x1= x-y
      y1= x
      z1=-z
      nz=5
      return 
 131  oname='S6 (-3) || -z'
      x1= y
      y1=-x+y
      z1=-z
      nz=5
      return 
!
!     Spiegelebenen
!
 141  oname='m (s-h) n  z'
      x1= x
      y1= y
      z1=-z
      return
 151  oname='m (s-d)  n  120'
      x1= x-y
      y1=-y
      z1= z
      return
 161  oname='m (s-d)  n  210'
      x1=-x
      y1=-x+y
      z1= z
      return
 171  oname='m (s-d)  n  1-10'
      x1= y
      y1= x
      z1= z
      return
 181  oname='m (s-v)  n  100'
      x1=-x+y
      y1= y
      z1= z
      return
 191  oname='m (s-v)  n  010'
      x1= x
      y1= x-y
      z1= z
      return
 201  oname='m (s-v)  n  110'
      x1=-y
      y1=-x
      z1= z
      return
!
!    3 zaehlige Achsen
!
 211  oname='3 || z'
      x1=-y
      y1= x-y
      z1= z 
      nz=2
      return
 221  oname='3 || -z'
      x1=-x+y
      y1=-x
      z1= z 
      nz=2
      return
 231  oname='S3 (-6) || z'
      x1=-y
      y1= x-y
      z1=-z 
      nz=2
      return
 241  oname='S3 (-6) || -z'
      x1=-x+y
      y1=-x
      z1=-z 
      nz=2
      return
!...........................rhomb part .........................
 2000 continue
      if(i.eq.1) goto 12
      if(i.eq.2) goto 22
      if(i.eq.3) goto 32
      if(i.eq.4) goto 42
      if(i.eq.5) goto 52
      if(i.eq.6) goto 62
      if(i.eq.7) goto 72
      if(i.eq.8) goto 82
      if(i.eq.9) goto 92
      if(i.eq.10) goto 102
      if(i.eq.11) goto 112
      if(i.eq.12) goto 122
!
!     Identitaet
!
 12   oname='1'
      x1= x
      y1= y
      z1= z
      return 
!
!     Inversion 
!
 22   oname='-1'
      x1=-x
      y1=-y
      z1=-z
      return 
!
!     2-zaehlige Drehachsen  
!
  32  oname='2 || y'
      x1=-x
      y1= y
      z1=-z 
      return
  42  oname='2 || 2-10'
      x1=0.5*x-sqrt(3.d0)/2.d0*y
      y1=-0.5*y-sqrt(3.d0)/2.d0*x
      z1=-z
      return
  52  oname='2 || 210'
      x1= 0.5*x+sqrt(3.d0)/2.d0*y
      y1= -0.5*y+sqrt(3.d0)/2.d0*x
      z1=-z                      
      return
!
!     6-zaehlige Drehachsen
!
  62  oname='S6 (-3) || z'
      x1= -0.5d0*x+sqrt(3.d0)/2.d0*y
      y1= -0.5d0*y-sqrt(3.d0)/2.d0*x
      z1=-z 
      nz=5
      return 
 72   oname='S6 (-3) || -z'
      x1= -0.5d0*x+sqrt(3.d0)/2.d0*y
      y1= -0.5d0*y-sqrt(3.d0)/2.d0*x
      z1=-z
      nz=5
      return 
!
!     Spiegelebenen
!
 82   oname='m (s-h) n  y'
      x1= x
      y1=-y
      z1= z
      return
 92   oname='m (s-d)  n  2-10'
      x1= -0.5d0*x+sqrt(3.d0)/2.d0*y
      y1=  sqrt(3.d0)/2.d0*x-0.5d0*y
      z1= z
      return
 102  oname='m (s-d)  n  210'
      x1= -0.5d0*x-sqrt(3.d0)/2.d0*y
      y1= -sqrt(3.d0)/2.d0*x+0.5d0*y
      z1= z
      return
!
!    3 zaehlige Achsen
!
 112  oname='3 || z'
      x1=-0.5d0*x-sqrt(3.d0)/2.d0*x
      y1= sqrt(3.d0)/2.d0*x-0.5d0*y
      z1= z 
      nz=2
      return
 122  oname='3 || -z'
      x1=-0.5d0*x+sqrt(3.d0)/2.d0*x
      y1=-sqrt(3.d0)/2.d0*x-0.5d0*y
      z1= z 
      nz=2
      return
      end
