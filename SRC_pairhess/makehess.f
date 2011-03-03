      subroutine makehess(translat,nuse,distances,pred,hess,br1,N1,N2,steps)
!     Pairwise Hessian initializer
!     Version 1.0, Feb 2006
!     J. Rondinelli, B. Deng and L. D. Marks
!
      implicit real*8 (a-h,o-z)
      include 'params.inc'
      dimension pred(3,natmax), g1(3,natmax),g2(3,natmax),steps(3)
      dimension hess(3*natmax,3*natmax),xred(3,natmax)
      dimension br1(3,3)
      integer nuse(natmax),translat(4,neigmax,48*natmax)
      dimension distances(neigmax,natmax*48)
!     Central difference hessian using gradients
!
!      write(6,*)'Steps ',steps
!      write(6,*)'Latt ',alat
      do j=1,nat*3
         do k=1,nat*3
            hess(j,k)=0.0D0
         enddo
      enddo
      ind1=0
      do j=1,nat
         do i=1,3
            ind1=ind1+1
            if( .not. isfree(i,j) )goto 100
!           Copy: sometimes might be needed
            do n=1,nat
               do ii=1,3
                  xred(ii,n)=pred(ii,n)
               enddo
            enddo
!           Shift
            xred(i,j)=pred(i,j)+steps(i)
!           Gradient
            call grad(translat,nuse,distances,xred,g1,br1,N1,N2,steps)
!
!           Copy: sometimes might be needed
            do n=1,nat
               do ii=1,3
                  xred(ii,n)=pred(ii,n)
               enddo
            enddo
!           Shift
            xred(i,j)=pred(i,j)-steps(i)
!
            call grad(translat,nuse,distances,xred,g2,br1,N1,N2,steps)
!
!           Hessian terms
            ind2=0
            do k=1,nat
               do L=1,3
                  ind2=ind2+1
!       The mini code uses atomic units, not fractional units, so rescale here
!       Also, multiply up by 2, kind of better for TiO2
                  if( isfree(L,K) ) &
                  hess(ind1,ind2)=(g1(L,K)-g2(L,K))/(steps(i)*alat(i)*alat(L))
!                  write(78,78)ind1,ind2,hess(ind1,ind2)
78                format('Hess ',2i3,D15.8)
               enddo
            enddo
100         continue
          enddo
      enddo
      return
      end
!
      subroutine grad(translat,nuse,distances,pred,gred,br1,N1,N2,steps)
      implicit real*8 (a-h,o-z)
      include 'params.inc'
      dimension pred(3,natmax), gred(3,natmax),steps(3)
      dimension xred(3,natmax)
      dimension br1(3,3)
      integer nuse(natmax),translat(4,neigmax,48*natmax)
      dimension distances(neigmax,natmax*48)
!
!     Central different gradients using energies
!
!     For each point, shift by step and calculate gradient
      gred(1:3,1:nat)=0.D0
      do j=1,nat
         do i=1,3
            if( .not. isfree(i,j) ) goto 100
!           Copy: sometimes might be needed
            do n=1,nat
               do ii=1,3
                  xred(ii,n)=pred(ii,n)
               enddo
            enddo
!           Shift
            xred(i,j)=pred(i,j)+steps(i)

!
            E1=energy(translat,distances,br1,nuse,xred,N2)
!
!           Copy: sometimes might be needed
            do n=1,nat
               do ii=1,3
                  xred(ii,n)=pred(ii,n)
               enddo
            enddo
!           Shift
            xred(i,j)=pred(i,j)-steps(i)

!
            E2=energy(translat,distances,br1,nuse,xred,N2)
            gred(i,j)=(E1-E2)/(2.D0*steps(i))
!            write(78,*)'Shift -',i,j,gred(i,j)
!            if(abs(gred(i,j)) .gt.1.D0)&
!                  write(78,78)(xred(LL,i),LL=1,3),(xred(LL,j),LL=1,3)
!78                format('Panic ',(3F15.7))
100         continue
            enddo
      enddo
      return
      end
