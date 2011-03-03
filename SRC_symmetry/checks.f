      SUBROUTINE  checks(a,trans)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /SYM2/ IZ(3,3,48),TAU(3,48),IORD                       
      dimension a(3,3),trans(3)
!
!.....compare symmetry operations with structfile
      do 10 i=1,iord
        do 15 j1=1,3
!           if(abs(trans(j1)-tau(j1,i)).gt.0.01) goto 10
        do 15 j2=1,3
 15        if(nint(a(j1,j2)).ne.iz(j1,j2,i)) goto 10 
!ccc 15        if(nint(a(j2,j1)).ne.iz(j1,j2,i)) goto 10 
      write(6,*) '          this is operation',i,'  in struct file'
        do 16 j1=1,3
           if(abs(trans(j1)-tau(j1,i)).gt.0.01) then
           write(6,20) trans,tau(1,i),tau(2,i),tau(3,i)
 20        format(//'          WARNING !!!!',/, &
           '    transl. differs:',3f10.4,3x,3f10.4,/)
           endif
 16           continue
      return  
 10   continue
      write(6,100) 
 100  format(//,'           WARNING !!!!!',/, &
      '   this symmetry operation was not found in struct!!',/)
      return
      end
