      Program xyz2struct
!          (spos,mult,isplit,name,rotloc,iatnr,gt)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (NDIF=1000)
      dimension x(ndif),y(ndif),z(ndif),mult(ndif)
      character*10 label(NDIF), name
      character*2 gt
!c
!c
      read(5,'(a10)') label(1)
      i=1
      ieq=1
      mult(ieq)=1
!c.....the format may need adaption for each case, would need lineparser
 10   continue
!      read(5,'(a10,f7.4,2x,f7.4,1x,f7.4)',err=99,end=99) &
      read(5,*,err=99,end=99) &
                                        label(i),x(i),y(i),z(i)
      write(*,*)  label(i),x(i),y(i),z(i)
      if(i.gt.1) then
         if(label(i).eq.label(i-1)) then
             mult(ieq)=mult(ieq)+1
         else
             ieq=ieq+1
             mult(ieq)=1
         endif
      endif
      i=i+1
      goto 10

 99   nat=ieq
      write(21,'(a)') 'xyz2struct'
      write(21,'(a1,26x,i3)') 'P',nat
      write(21,'(a)') 'MODE OF CALC=RELA'
      write(21,'(6f10.5)') 5.,5.,5.,90.,90.,90.
!c
      i=0
      do iatnr=1,ieq
      i=i+1
      write(21,1010) IATNR,x(i),y(i),z(i)
!c      write(21,1010) IATNR,x(i),z(i),y(i)
      write(21,1012) mult(iatnr),8
 1010 FORMAT('ATOM',I4,': X=',F10.8,' Y=',F10.8,' Z=',F10.8)            
 1011 FORMAT(4X,I4,': X=',F10.8,' Y=',F10.8,' Z=',F10.8)                
 1012 FORMAT('          MULT=',i2,'          ISPLIT=',i2)
      DO 55 M=2,MULT(iatnr)              
      i=i+1                       
 55   write(21,1011) IATNR,x(i),y(i),z(i)
      write(21,1510) label(i)
      write(21,1013) 1.,0.,0.,0.,1.,0.,0.,0.,1.          
 1013 format('LOCAL ROT MATRIX:   ',3f10.7,/,20x,3f10.7,/,20x,3f10.7)
 1510 FORMAT(A10,' NPT=  781  R0=0.00001000 RMT=    2.0000   Z: 00.0')  
      enddo
!c
      write(21,*) 0                                          
!c
      end
