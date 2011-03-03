!..... progaram sym is used for calculation of transformation parameters
!..... of so wavefunctions. It transforms the symmetry operations
!..... to frame of reference z||M. Variable det=1 ... M not reversed
!..... det=-1 ... M reversed ( must by later coupled to time  inversion ).
!.... Angle phase(I) is used for spin part of transformation matrix 
!..... acting of spinor function.
       
      SUBROUTINE SYM(THETA,FI)
      USE param
      USE struct
      USE sym2
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GENER/  BR1(3,3),BR2(3,3)

      DIMENSION rot(3,3),trans(3,3),rotinv(3,3)
      
      pi=acos(-1.d0)
      CALL symoper

      ct=cos(theta)
      cf=cos(fi)
      st=sin(theta)
      sf=sin(fi)
      rot(1,1)=cf*ct
      rot(1,2)=-sf
      rot(1,3)=cf*st
      rot(2,1)=sf*ct
      rot(2,2)=cf
      rot(2,3)=sf*st
      rot(3,1)=-st
      rot(3,2)=0
      rot(3,3)=ct
      write(6,*)'SPIN COORD SYSTEM:'
      write(6,5)((rot(i,j),j=1,3),i=1,3)
	write(6,*)
      call INVERSSYMDEF(rot,rotinv)
      do 100 i=1,IORD
	write(6,*)'in global cartezian coordinates:'
	 write(6,5)((opimat(ii,jj,i),jj=1,3),ii=1,3)
 5	 format(3f12.6)
          write(6,*)'____________________________'
         call transform(trans,rotinv,opimat(1,1,i),rot)
          write(6,*)'in spin coordinates:'
          write(6,5)((trans(ii,jj),jj=1,3),ii=1,3)
         det=trans(1,1)*trans(2,2)-trans(1,2)*trans(2,1)
!.... calculation of the phase shift of so - functions under
!...  sym operations: for operations det=1 is the spin matrix 
!... diagonal: | exp(-i*phase/2),             0|
!.......       |               0,exp(i*phase/2)|
!.... for operations det=-1 is the spin matrix of the pure rotation
!.... off-diagonal (i.e. reverts spinor up->dn, dn->up) and must
!.... be combined with time inversion ( which reverts the spinor again )
!.... in this case the phase(i) defined in SYM includes also
!.... part of the time inversion ( the other part of time inversion
!.... - complex conjugation of wave function is included in XSPLT )

	 DD=det*trans(3,3)
	 if (abs(trans(1,1)).gt.1d0) then
	 a=trans(1,1)/abs(trans(1,1))
         else
	 a=trans(1,1)
         end if

	 b=trans(1,2)

	 if (abs(b).gt.1d-8) then
	 phase(i)=-DD*b/(abs(-DD*b))*acoss(DD*a)
	 else
	 phase(i)=acoss(DD*a)
	 end if
 	 if (det.lt.0.) then
 	 phase(i)=phase(i)+pi
 	 end if

	 write(6,*)i,phase(i),det
	 write(6,*)
	 if (abs(1.-abs(det)).gt.1d-4) then
	 write(6,*)'symm. operation ',i,' so-det=',det
	 STOP
	 end if
	 if (det.gt.0.) then
	 idet(i)=1
	 else
	 idet(i)=-1
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














