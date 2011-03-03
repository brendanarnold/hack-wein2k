      program kreis
!     calc hex k vectors with fixed magnitude and varying angle from 0 to 30 
      dimension a(2,2),b(2,2),r(2),rk(2),ik(2,0:30)
!
      pi=acos(-1.0)
      a(1,1)=5.637825
      a(2,1)=-3.255000 
      a(2,2)=6.510000 
      a(1,2)=0.0
!
      b(1,1)=    .177373
      b(2,1)=0.0
      b(1,2)=.088687
      b(2,2)=.153610
!
      rk(1)=1./6.
      rk(2)=1./6.
!
      r(1)=rk(1)*b(1,1)+rk(2)*b(1,2)
      r(2)=rk(1)*b(2,1)+rk(2)*b(2,2)
      rmag=sqrt(r(1)**2+r(2)**2)
      write(6,*) r,rmag
!
      do i=0,30
      r(1)=rmag*cos(i*pi/180.0)
      r(2)=rmag*sin(i*pi/180.0)
      rk(1)=r(1)*a(1,1)+r(2)*a(2,1)
      rk(2)=r(1)*a(1,2)+r(2)*a(2,2)
      ik(1,i)=nint(rk(1)*99998)
      ik(2,i)=nint(rk(2)*99998)
      write(*,'(i10,4i5,f5.2)') i,nint(rk(1)*99998),nint(rk(2)*99998) &
                                ,0,99998,1.0
      enddo
      do i=0,30
      write(*,'(i10,4i5,f5.2)') i+31,ik(1,i),ik(2,i),49999,99998,1.0
      enddo
!
      end
