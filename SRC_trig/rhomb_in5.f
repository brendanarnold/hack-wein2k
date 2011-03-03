      PROGRAM trigo                                                        
!
! Calculates the vertices of a plane for lapw5 
! Specify 3 atoms (points) in rhombohedral coordinates within the desired planeC                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
!
      DIMENSION DIF(3),XX(3),PP1(3),Pp2(3),orig(3),xend(3),yend(3)  
      DIMENSION help(3),pnn(3,3),d1(3),d2(3)                     
      COMMON /GENER  /  BR2(3,3),br1(3,3)                                     
      COMMON /STRUK  /  POS(3,3),A(3)
!-----------------------------------------------------------------------
!                                                                       
      write(*,*) ' Input for lapw5 for rhombohedral structures'
      write(*,*) ' specify hex. lattice parameters a and c'
      read(*,*) a(1),a(3)
!.....SET UP LATTICE, AND IF REQUIRED ROTATION MATRICES                 
!.....RHOMBOHEDRAL CASE
 80   BR2(1,1)=1/2.d0/sqrt(3.d0)*a(1)
      BR2(1,2)=-1/2.d0*a(1)                                                    
      BR2(1,3)=1/3.d0*a(3)                                                     
      BR2(2,1)=1/2.d0/SQRT(3.d0)*a(1)                                          
      BR2(2,2)=1*0.5d0*a(1)                                                
      BR2(2,3)=1/3.d0*a(3)                                                     
      BR2(3,1)=-1/SQRT(3.d0)*a(1)                                         
      BR2(3,2)=0.d0                                                
      BR2(3,3)=1/3.d0*a(3)                                                     
!                                                                       
      write(*,*) 'Bravais Matrix:'
      write(*,999) br2
 999  format(3f15.5)
!      CALL DIRLAT (nat,ortho,alpha,beta,gamma)                                
!                                           
      pi=4.d0*atan(1.d0)
       write(*,*) 'give 3 trigonal coordinates'
       read(*,*) (pos(i,1),i=1,3)
       read(*,*) (pos(i,2),i=1,3)
       read(*,*) (pos(i,3),i=1,3)
!......convert to orthogonal
        do i=1,3
        help(1)=pos(1,i)
        help(2)=pos(2,i)
        help(3)=pos(3,i)
        pnn(1,i)=help(1)*BR2(1,1)+help(2)*BR2(2,1)+help(3)*BR2(3,1)            
        pnn(2,i)=help(1)*BR2(1,2)+help(2)*BR2(2,2)+help(3)*BR2(3,2)            
        pnn(3,i)=help(1)*BR2(1,3)+help(2)*BR2(2,3)+help(3)*BR2(3,3)           
!c        write(*,*) 'orth:',pnn(1,i),pnn(2,i),pnn(3,i)
        enddo
!.....calc 2 dif vectors
       do i=1,3
       d1(i)=pnn(i,2)-pnn(i,1)
       d2(i)=pnn(i,3)-pnn(i,1)
       enddo
        write(*,*) 'd1:',(d1(i),i=1,3)
        write(*,*) 'd2:',(d2(i),i=1,3)
!.....normal vector
      pp1(1)= d1(2)*d2(3)-d1(3)*d2(2)
      pp1(2)=-d1(1)*d2(3)+d1(3)*d2(1)
      pp1(3)= d1(1)*d2(2)-d1(2)*d2(1)
!.....2nd normal vector:
      pp2(1)= d1(2)*pp1(3)-d1(3)*pp1(2)
      pp2(2)=-d1(1)*pp1(3)+d1(3)*pp1(1)
      pp2(3)= d1(1)*pp1(2)-d1(2)*pp1(1)
!        write(*,*) 'n1:',(pp1(i),i=1,3)
!        write(*,*) 'n2:',(pp2(i),i=1,3)
      xlen=sqrt(d1(1)**2+d1(2)**2+d1(3)**2)
      ylen=sqrt(pp2(1)**2+pp2(2)**2+pp2(3)**2)
!.....normalization to x-lenght
        pp2(1)=pp2(1)/ylen*xlen
        pp2(2)=pp2(2)/ylen*xlen
        pp2(3)=pp2(3)/ylen*xlen
      ylen=3.d0*sqrt(pp2(1)**2+pp2(2)**2+pp2(3)**2)
      write(*,*) 'distance between atom 1 and 2',xlen
!   
      write(*,*) 'give x and y-offset scaling factors (in units of d12)'      
      read(*,*) alam,amu
!     origin of plot, x-end, y-end
      do i=1,3
      orig(i)=pnn(i,1)-alam*d1(i)-amu*pp2(i)
      xend(i)=pnn(i,2)+alam*d1(i)-amu*pp2(i)
      yend(i)=orig(i)+pp2(i)+2.d0*amu*pp2(i)
      enddo
!c        write(*,*) 'orig_orth:',(orig(i),i=1,3)
!c        write(*,*) 'xend_orth:',(xend(i),i=1,3)
!c        write(*,*) 'yend_orth:',(yend(i),i=1,3)
      TEST1=d1(1)*pp2(1)+d1(2)*pp2(2)+d1(3)*pp2(3)                         
      IF(abs(TEST1).GT.0.001) then
       write(*,*) ' directions not orthogonal ',test1
      endif
!.....back to trigonal coordinates
      call reclat(br2,br1,0)
      hh=br1(1,2)
      br1(1,2)= br1(2,1)
      br1(2,1)=hh
      hh=br1(1,3)
      br1(1,3)= br1(3,1)
      br1(3,1)=hh
      hh=br1(2,3)
      br1(2,3)= br1(3,2)
      br1(3,2)=hh

      write(*,*) 'inv.Bravais Matrix:'
      write(*,999) br1
        help(1)=orig(1)
        help(2)=orig(2)
        help(3)=orig(3)
        orig(1)=help(1)*BR1(1,1)+help(2)*BR1(2,1)+help(3)*BR1(3,1)            
        orig(2)=help(1)*BR1(1,2)+help(2)*BR1(2,2)+help(3)*BR1(3,2)            
        orig(3)=help(1)*BR1(1,3)+help(2)*BR1(2,3)+help(3)*BR1(3,3)           
        help(1)=xend(1)
        help(2)=xend(2)
        help(3)=xend(3)
        xend(1)=help(1)*BR1(1,1)+help(2)*BR1(2,1)+help(3)*BR1(3,1)            
        xend(2)=help(1)*BR1(1,2)+help(2)*BR1(2,2)+help(3)*BR1(3,2)            
        xend(3)=help(1)*BR1(1,3)+help(2)*BR1(2,3)+help(3)*BR1(3,3)           
        help(1)=yend(1)
        help(2)=yend(2)
        help(3)=yend(3)
        yend(1)=help(1)*BR1(1,1)+help(2)*BR1(2,1)+help(3)*BR1(3,1)            
        yend(2)=help(1)*BR1(1,2)+help(2)*BR1(2,2)+help(3)*BR1(3,2)            
        yend(3)=help(1)*BR1(1,3)+help(2)*BR1(2,3)+help(3)*BR1(3,3)           
!cc        write(*,*) 'orig_trig:',(orig(i),i=1,3)
!cc        write(*,*) 'xend_trig:',(xend(i),i=1,3)
!c        write(*,*) 'yend_trig:',(yend(i),i=1,3)
        write(*,*) 'trigonal coordinates for lapw5'
        write(*,*) (int(orig(i)*1000000),i=1,3),1000000
        write(*,*) (int(xend(i)*1000000),i=1,3),1000000
        write(*,*) (int(yend(i)*1000000),i=1,3),1000000
                                                                      
      STOP                                                        
!                                                                       
 1011 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1012 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,7X,I2,8X,I2)     
 1020 FORMAT(///,3X,'ERROR IN EBND  : MULT(JATOM)=0 ...',  &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)           
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5)                               
 1051 FORMAT(20X,3F10.8)                                                
 1510 FORMAT(A80)                                                       
 1511 FORMAT(A4,24X,I2,/,13X,A4)                                        
      END                                                               
      SUBROUTINE RECLAT (A,B,IOPT)

!  CALCULATES RECIPROCAL LATTICE VECTORS. THEIR PRODUCT WITH DIRECT
!  LATTICE VECTORS IS 1 IF IOPT=0 OR 2*PI IF IOPT=1

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(3,3),B(3,3)
      PI=ACOS(-1.D0)
      B(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
      B(2,1)=A(3,2)*A(1,3)-A(1,2)*A(3,3)
      B(3,1)=A(1,2)*A(2,3)-A(2,2)*A(1,3)
      B(1,2)=A(2,3)*A(3,1)-A(3,3)*A(2,1)
      B(2,2)=A(3,3)*A(1,1)-A(1,3)*A(3,1)
      B(3,2)=A(1,3)*A(2,1)-A(2,3)*A(1,1)
      B(1,3)=A(2,1)*A(3,2)-A(3,1)*A(2,2)
      B(2,3)=A(3,1)*A(1,2)-A(1,1)*A(3,2)
      B(3,3)=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      C=1.D0
      IF (IOPT.EQ.1) C=2.D0*PI
      DO 20 I=1,3
         CI=C/(A(1,I)*B(1,I)+A(2,I)*B(2,I)+A(3,I)*B(3,I))
         B(1,I)=B(1,I)*CI
         B(2,I)=B(2,I)*CI
         B(3,I)=B(3,I)*CI
  20  CONTINUE
      END
