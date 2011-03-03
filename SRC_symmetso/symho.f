!..... program sym is used for calculation of transformation parameters
!..... of so wavefunctions. It transforms the symmetry operations
!..... to frame of reference z||M. Variable det=1 ... M not reversed
!..... det=-1 ... M reversed ( must by later coupled to time  inversion ).
!..... det= 0 ... M changes direction                                     
       
      SUBROUTINE SYMHO(lattic,alat,alpha,theta,phi,iord,imat,br1,isigk)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer isigk(48),imat(3,3,48)
      character*4 lattic
      LOGICAL print
      DIMENSION rot(3,3),trans(3,3),rotinv(3,3), rot1(3,3)
      DIMENSION alat(3),alpha(3),br1(3,3),br1inv(3,3),det(48)
!      parameter (pi=3.141592654d0) 
      pi=acos(-1.d0)
      print=.false.	
     print=.true.	
      ct=cos(theta)
      cf=cos(phi)
      st=sin(theta)
      sf=sin(phi)
      rot(1,1)=cf*ct
      rot(1,2)=-sf
      rot(1,3)=cf*st
      rot(2,1)=sf*ct
      rot(2,2)=cf
      rot(2,3)=sf*st
      rot(3,1)=-st
      rot(3,2)=0
      rot(3,3)=ct
      call INVERSS(rot,rotinv)
      if(print)then
      write(6,*)'SPIN COORD SYSTEM:'
      write(6,5)((rot(i,j),j=1,3),i=1,3)
      write(6,5)((rotinv(i,j),j=1,3),i=1,3)
      endif
!cccc
        alpha(1)=alpha(1)*180.d0/pi           
        alpha(2)=alpha(2)*180.d0/pi           
        alpha(3)=alpha(3)*180.d0/pi           
      call br1dm(br1,lattic,alat,alpha)
      call INVERSS(br1,br1inv)
      if(print)then
      write(6,*)' BR1'
      write(6,5)((br1(i,j),j=1,3),i=1,3)
      write(6,*)' BR1inv'
      write(6,5)((br1inv(i,j),j=1,3),i=1,3)
      endif
      do 100 i=1,IORD
        do j=1,3
         do k=1,3
         rot1(j,k)=0.            
           do l=1,3
           do m=1,3
            rot1(j,k)=rot1(j,k)+   &
                      br1(j,l)*imat(l,m,i)*br1inv(m,k)
!!changed P.B. 27.6.05                  br1(j,l)*imat(m,l,i)*br1inv(m,k)
            enddo
           enddo 
         enddo
        enddo
         call transform(trans,rotinv,rot1,rot)
 5       format(3f12.6)
 	if (print) then
	WRITE(6,*)'****SYMMETRY OPERATION ',I
	write(6,*)'in global cartesian coordinate system:'
	write(6,5)((rot1(ii,j),j=1,3),ii=1,3)
        write(6,*)
	write(6,*)'in spin coordinate system:'
        write(6,5)((trans(ii,jj),jj=1,3),ii=1,3)
        write(6,*)
	end if
         det(i)=trans(1,1)*trans(2,2)-trans(1,2)*trans(2,1)
	 if (abs(1.-abs(det(i))).gt.1d-2) then
         isigk(i)=0
         else if (det(i).gt.0.999) then
         isigk(i)=1
         else
         isigk(i)=-1
	 end if
 100  continue
      end
      
      subroutine transform(T,Pinv,A,P)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION T(3,3),P(3,3),A(3,3),Pinv(3,3)
      
      do i=1,3
         do j=1,3
            sum=0
            do k=1,3
               do l=1,3
                  sum=sum+Pinv(i,k)*A(k,l)*P(l,j)
               end do
            end do
            T(i,j)=sum
         end do
!       write(6,*)'Transf. matrix:',(T(i,k),k=1,3)
      end do
      return
      end

        real*8 function ACOSS(x)
        implicit real*8 (A-H,O-Z)

        if ((abs(x)-1.).gt.1d-4) then
        write(6,*)'x=',x
        stop 'ACOSS ERROR'
        end if
        if (x.ge.1) then
        acoss=0.0
        else if (x.le.-1.) then
        acoss=acos(-1.0)
        else
        acoss=acos(x)
        end if
        end

















