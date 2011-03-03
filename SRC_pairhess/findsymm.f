!     Find symmetry operations that relate reduced to larger set
      subroutine findsymm(pred,pos,N1,N2,br)
      implicit real*8 (a-h,o-z)
      include 'params.inc'
      dimension pos(3,48*natmax),pred(3,natmax)
      dimension br(3,3),deb(3)
      parameter (tol=1.D-4)
!
      index=0
      do j=1,nat
!
!        Skip search for the first atom
         index=index+1
         ind=index
         icode(1,index)=j
         icode(2,index)=0
         do k=2,mult(j)
            index=index+1
!           Search over operations
            do i=1,iord
!             Positions after applying symmetry operations
              x = iz(1,1,i)*pred(1,j) + iz(2,1,i)*pred(2,j) + iz(3,1,i)*pred(3,j) +  &
              tau(1,i) 
              y = iz(1,2,i)*pred(1,j) + iz(2,2,i)*pred(2,j) + iz(3,2,i)*pred(3,j) +  &
              tau(2,i)
              z = iz(1,3,i)*pred(1,j) + iz(2,3,i)*pred(2,j) + iz(3,3,i)*pred(3,j) + &
              tau(3,i)
!             Distances
              dx=x-pos(1,index)
              dy=y-pos(2,index)
              dz=z-pos(3,index)
!             Remove integer elements
              DX0=MOD(DX,1.0D0)
              DY0=MOD(DY,1.0D0)
              DZ0=MOD(DZ,1.0D0)
!
              if(dx0.gt.0.50001D0)DX0=DX0-1.D0
              if(dy0.gt.0.50001D0)DY0=DY0-1.D0
              if(dZ0.gt.0.50001D0)DZ0=DZ0-1.D0
!
              if(dx0.lt.-0.50001D0)DX0=DX0+1.D0
              if(dy0.lt.-0.50001D0)DY0=DY0+1.D0
              if(dZ0.lt.-0.50001D0)DZ0=DZ0+1.D0

!             Within tol a.u. on all axes
              if( (alat(1)*abs(dx0) .lt. tol) .and. &
                  (alat(2)*abs(dy0) .lt. tol) .and. &
                  (alat(3)*abs(dz0) .lt. tol) )then
!                 This one works, so store information
                     icode(1,index)=j
                     icode(2,index)=i
                     icode(3,index)=nint(dx0-dx)
                     icode(4,index)=nint(dy0-dy)
                     icode(5,index)=nint(dz0-dz)
!                     write(6,*)'Got ',(icode(ii,index),ii=1,5)
                  goto 100
               endif
            enddo
!       Aaargh, did not find the equivalent site
!       Die, gracefully
            write(6,*)'Cannot find operation for atom ',j,' Equiv ',k
            write(6,*)'Please check the symmetry operations with x sgroup'
            write(6,66)'Base atom ',(pred(L,j),L=1,3)
            write(6,66)'Equiv atom',(pos(L,index),L=1,3)
66          format(a,3f14.8)
            do i=1,iord
                write(6,1001)i
!               Positions after applying symmetry operations
                deb(1) = iz(1,1,i)*pred(1,j) + iz(2,1,i)*pred(2,j) + iz(3,1,i)*pred(3,j) +  &
                        tau(1,i) 
                deb(2) = iz(1,2,i)*pred(1,j) + iz(2,2,i)*pred(2,j) + iz(3,2,i)*pred(3,j) +  &
                        tau(2,i) 
                deb(3) = iz(1,3,i)*pred(1,j) + iz(2,3,i)*pred(2,j) + iz(3,3,i)*pred(3,j) + &
                        tau(3,i) 
                do L=1,3
                        write(6,1002)(iz(kk,L,i),kk=1,3),tau(L,i),deb(L)
                enddo
            enddo
1001        format('Operation ',i3)
1002        format(3I5,2F12.6)
            stop 'Symmetry problem, please check with x sgroup'
100      continue
         enddo
      enddo
      return
      end

