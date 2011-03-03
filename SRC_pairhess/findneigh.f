        subroutine findneigh(translat,nuse,distances,pred,br,N1,N2,rmax)
!     Pairwise Hessian initializer
!     Version 1.0, Feb 2006
!     J. Rondinelli, B. Deng and L. D. Marks
!
!     Find the neighbors
        implicit real*8 (a-h,o-z)
        include 'params.inc'
        dimension pred(3,natmax),pall(3,48*natmax)
        dimension idist(4,neigmax*natmax*48),iwork(neigmax*48)
        dimension dist(neigmax),distances(neigmax,natmax*48)
        dimension trans(3,neigmax)
        integer translat(4,neigmax,natmax*48),nuse(natmax)
        dimension br(3,3)
        character*67 ErrMsg
!
        rmax0=rmax
        rmax=rmax*rmax
!
        write(6,1000)' Nearest Neighbour listing'
        write(6,1000)' *************************'
!
!       Fixup positions to numerical accuracy
        call fixup3(pred)
!
!       Expand the set to all equivalent atoms
        call expandset(pred,pall)
!
!       Loop for first atom
        do j=1,nat
           N=0
           WMULT(J)=DBLE(MULT(J))/DBLE(IORD)
!          Unit cell translations
           do 10 i1=-3,3
           do 10 i2=-3,3
           do 10 i3=-3,3
!          Convert to a.u. in cartesians
           x11=(pred(1,j)+i1)*BR(1,1)+(pred(2,j)+i2)*BR(2,1)+(pred(3,j)+i1)*BR(3,1)
           x12=(pred(1,j)+i1)*BR(1,2)+(pred(2,j)+i2)*BR(2,2)+(pred(3,j)+i2)*BR(3,2)
           x13=(pred(1,j)+i1)*BR(1,3)+(pred(2,j)+i2)*BR(2,3)+(pred(3,j)+i3)*BR(3,3)
!          Loop over second atom
           do 10 k=1,N2
               x21=pall(1,k)*BR(1,1)+pall(2,k)*BR(2,1)+pall(3,k)*BR(3,1)
               x22=pall(1,k)*BR(1,2)+pall(2,k)*BR(2,2)+pall(3,k)*BR(3,2)
               x23=pall(1,k)*BR(1,3)+pall(2,k)*BR(2,3)+pall(3,k)*BR(3,3)
!
!              Distance
               p11=x21-x11
               p12=x22-x12
               p13=x23-x13
               d=p11*p11+p12*p12+p13*p13
!              In target range ?
               if((d .gt. 0.25).and. (d .lt. rmax))then
                     n=n+1
                     if(n .gt. neigmax)then
!                       Oops, too many
                        WRITE (ERRMSG,9001)j
9001                    format('Error, too many neighbors for atom ',i3)
                        CALL OUTERR('PairHess',ERRMSG)
                        write(6,9001)j
                        write(*,9001)j
                        write(*,9002)'Current value of rmax ',rmax0
                        STOP 'PairHess error, too many neighbors, reduce rmax'
                      endif
!                    Store information for later
                     iwork(n)=n
                     dist(n)=sqrt(d)
                     idist(1,n)=i1
                     idist(2,n)=i2
                     idist(3,n)=i3
                     idist(4,n)=k
                     trans(1,n)=p11
                     trans(2,n)=p12
                     trans(3,n)=p13
              endif
10         continue
!
!          Die if we have no neighbors
           if(n .eq. 0)then
             WRITE (ERRMSG,9000)j
9000         format('Error, no neighbors for atom ',i3)
             CALL OUTERR('PairHess',ERRMSG)
             write(6,9000)j
             write(*,9000)j
             write(*,9002)'Current value of rmax ',rmax0
9002         format(a,F10.6)
             STOP 'PairHess error, no neighbors, increase rmax'
           endif
!          Sort so they are increasing
           call sortag(dist, N, iwork)
!          Save
           nuse(j)=n
!
!        Info that the user might want to check
!        This should be consistent with x nn output
         write(6,*)'Atom ',j,' Neighbours ',n
         write(6,6)' Position ',(pall(k,j),k=1,3)
6        format(a,3F14.8)
           d0=dist(1)
           au=0.52918D0
           do i=1,n
              d=dist(i)
              distances(i,j)=d
              k=iwork(i)
              translat(1,i,j)=idist(1,k)
              translat(2,i,j)=idist(2,k)
              translat(3,i,j)=idist(3,k)
              translat(4,i,j)=idist(4,k)
              atrans(1,i,j)=trans(1,k)
              atrans(2,i,j)=trans(2,k)
              atrans(3,i,j)=trans(3,k)
              wt= exp(-(d-d0)/d0*scl)
              atrans(4,i,j)=wt
              if(wt.gt.cutoff)then
                  write(6,1000)' Distance to atom ',idist(4,k),d*au,wt
                  nuse(j)=i
              else
                  goto 200
              endif
           enddo
200        write(6,1000)
      enddo
      return
1000  format(a,i3,F12.6,' Angs, weight ',F10.5)
      end
