      parameter(mesh=120)
      dimension a(mesh,mesh),a1(mesh,mesh),a2(mesh,mesh),ah(4)
      dimension rold(mesh)
!.....nx,ny  interpotated mesh (eg. 2*meshx-1,....
!.....invert:  0  no,  1 data von (-nx+1) bis (nx-1),....
!.....flip:    0  no,  1 interchange x and y
      read(*,*,err=999) meshx,meshy,x0,y0,nx,ny, invert,iflip
       write(*,*) meshx,meshy,x0,y0,nx,ny, invert,iflip
      if(nx.gt.mesh) stop 'mesh parameter'
      if(ny.gt.mesh) stop 'mesh parameter'
      iunit=10
 1    read(*,*,end=99)
      if(iflip.eq.0) then
      do j=1,meshy
      do i=1,meshx
      read(*,*,err=10) x,y,z,r,a(i,j)
      enddo
      enddo   
      else if(iflip.eq.1) then
      do j=1,meshy
      do i=j,meshx
      read(*,*,err=10) x,y,z,r,a(i,j)
      write(*,*) j,i,a(i,j)
       a(j,i)=a(i,j)
      enddo
      enddo   
      endif
!   
!  interpolation along x
      do i=1,meshx
        rold(i)=real((i-1))/(meshx-1)
      enddo
      do j=1,meshy
      do i=1,nx
        iold=2
 55     iold=iold+1
        if(real((iold-1))/(meshx-1).lt.real((i-1))/(nx-1)) goto 55
        rnew=real((i-1))/(nx-1)
!        if(iold.ge.meshx-3) iold=meshx-3
!           write(*,*) j,i,iold,rold(iold-2),rnew
                call INTER(rold(iold-2),a(iold-2,j),4,rnew &
                          ,a1(i,j),dps,.false.) 
      enddo
      enddo
!  interpolation along y
      do i=1,meshy
        rold(i)=real((i-1))/(meshy-1)
      enddo
      do j=1,nx
      do i=1,ny
        iold=2
 56     iold=iold+1
        if(real((iold-1))/(meshy-1).lt.real((i-1))/(ny-1)) goto 56
           rnew=real((i-1))/(ny-1)
!           if(iold.ge.meshy-3) iold=meshy-3
!           write(*,*) j,i,iold,rold(iold-2),rnew
           do n=0,3
           ah(n+1)=a1(j,iold-2+n) 
           enddo  
           call INTER(rold(iold-2),ah(1),4, &
                          rnew,a2(j,i),dps,.false.) 
      enddo
      enddo
           
      iunit=iunit+1
      if(invert.eq.0) then
      write(iunit,*)  nx,ny,x0,y0
      write(iunit,11) ((a2(i,j),j=1,ny),i=1,nx)
      else
      write(iunit,*)  2*nx-1,2*ny-1,x0,y0
      write(iunit,11)((a2(abs(i)+1,abs(j)+1),j=-ny+1,ny-1),i=-nx+1,nx-1)
      endif
 11   format(5e16.8)
      goto 1
 99   stop
 999  write(*,*) 'please insert    meshx,meshy,xlen,ylen, nxint,nyint,invert'
      stop
 10   write(*,*) 'error reading file. Propably wrong mesh'
      end

      SUBROUTINE INTER(R,P,N,RS,PS,DPS,DERIV)                          
!                                                                       
!     ******************************************************************
!                                                                       
!     INTERPOLATE VIA LAGRANGE                                          
!                                                                       
!     ******************************************************************
!                                                                       
!      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL DERIV,NODRIV                                              
      DIMENSION R(N),P(N)                                               
      NODRIV=.NOT.DERIV                                                 
      PS=0.E0                                                           
      DPS=0.E0                                                          
      DO 1 J=1,N                                                        
      TERM=1.E0                                                         
      DENOM=1.E0                                                        
      DTERM=0.E0                                                        
      DO 2 I=1,N                                                        
      IF(I.EQ.J) GO TO 2                                                
      DENOM=DENOM*(R(J)-R(I))                                           
      TERM=TERM*(RS-R(I))                                               
      IF(NODRIV) GO TO 2                                                
      DTERM1=1.E0                                                       
      DO 3 K=1,N                                                        
      IF(K.EQ.J.OR.K.EQ.I) GO TO 3                                      
      DTERM1=DTERM1*(RS-R(K))                                           
    3 CONTINUE                                                          
      DTERM=DTERM+DTERM1                                                
    2 CONTINUE                                                          
      IF(NODRIV) GO TO 1                                                
      DPS=DPS+DTERM*P(J)/DENOM                                          
    1 PS=PS+TERM *P(J)/DENOM                                            
      RETURN                                                            
      END                                                               
