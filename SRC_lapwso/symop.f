      SUBROUTINE symop(theta,fi)
!     for .orhto. systems   opimat == imat
!     for not.ortho. systems opimat == BR1 . imat . BR1^-1 
      USE rotmat
      USE struct
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 BR1in(3,3),opimat(3,3),sglo(3,3),rotinv(3,3),trans(3,3)

      REAL*8 rloc(3,3),sloc(3,3),th(10),ph(10)
      LOGICAL          ORTHO
      COMMON /ORTH/   ORTHO
      COMMON /GENER/  BR1(3,3),BR2(3,3)

!.........inverssymdef is a subroutine in sph-UP.frc and
!.........calculates the inverse of an 3*3 matrix.......
      pi=acos(-1.d0)
      CALL INVERSSYMDEF(BR1,BR1in)
      ct=dcos(theta)
      cf=dcos(fi)
      st=dsin(theta)
      sf=dsin(fi)
      sglo(1,1)=cf*ct
      sglo(1,2)=-sf
      sglo(1,3)=cf*st
      sglo(2,1)=sf*ct
      sglo(2,2)=cf
      sglo(2,3)=sf*st
      sglo( 3,1)=-st
      sglo(3,2)=0
      sglo(3,3)=ct
 5	format(3f12.7)
      call INVERSSYMDEF(sglo,rotinv)
	
	indj=0
        do 100 jatom=1,nat
	do 100 mu=1,mult(jatom)
	indj=indj+1
          do j=1,3
          do k=1,3
          opimat(j,k)=0. 
          end do
          end do
 
          do j=1,3
          do k=1,3
            do l=1,3
            do m=1,3
          opimat(j,k)=opimat(j,k)+                   &
               BR1(j,l)*ROTIJ(l,m,indj)*BR1in(m,k)
            end do
            end do
          end do
          end do
         call transform(trans,rotinv,opimat,sglo)

	 do 200 i=1,3
	 do 200 j=1,3
	 rloc(i,j)=0.d0
	 do  200 k=1,3
	 rloc(i,j)=rloc(i,j)+opimat(k,j)*rotloc(i,k,jatom)
 200     continue

         tt=determinant(rloc)
         do 300 i=1,3
         do 300 j=1,3
         sloc(i,j)=0.d0
         do  300 k=1,3
         sloc(i,j)=sloc(i,j)+tt*rloc(i,k)*sglo(k,j)
 300     continue

!        write(6,*)'ATOM:',indj
! write(6,*)'local coord. system:'
! do ii=1,3
! write(6,5)(rloc(ii,jj),jj=1,3)
! end do
! write(6,*)'rotij:'
! do ii=1,3
! write(6,5)(rotij(ii,jj,indj),jj=1,3)
! end do
! write(6,*)'opimat:'
! do ii=1,3
! write(6,5)(opimat(ii,jj),jj=1,3)
! end do
! write(6,*)
! write(6,*)'spin with respect to local coord.:'
! do ii=1,3
! write(6,5)(sloc(ii,jj),jj=1,3)
! end do
! write(6,*)
	
	 call euler(sloc,a,b,c)
	 call couple(indj,a,b,c)


         det(indj)=trans(1,1)*trans(2,2)-trans(1,2)*trans(2,1)	
         if (abs(1.-abs(det(indj))).gt.1.d-2) then
         write(6,*)'WRONG SYMMETRY'
         write(6,555)indj,det(indj)
555     format(' atom',i3,' det',f10.5)
!       stop
        endif
	
	DD=determinant(trans)
        if (DD.lt.-0.5) then
         do ii=1,3
         do jj=1,3
         trans(ii,jj)=-trans(ii,jj)
         end do
         end do
        end if
         call euler(trans,a,b,c)
!       write(6,3)indj,a*180/pi,b*180/pi,c*180/pi
  3     format(i2,1x,'Euler angles: a,b,c:  ',3f7.1)
  4     format('det=',f3.0,' phase=',f8.4)
         if (det(indj).gt.0.5) then
         phase(indj)=a+c
         else
         phase(indj)=a-c
         end if
!       write(6,4)det(indj),phase(indj)
 100	continue
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
      end do
      return
      end

      SUBROUTINE INVERSSYMDEF(A,AINV)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(3,3),AINV(3,3)
        det= determinant(a)
      AINV(1,1) =(   A(2,2) * A(3,3) - A(2,3) * A(3,2) ) / det
      AINV(2,1) =( - A(2,1) * A(3,3) + A(2,3) * A(3,1) ) / det
      AINV(3,1) =(   A(2,1) * A(3,2) - A(2,2) * A(3,1) ) / det
      AINV(1,2) =( - A(1,2) * A(3,3) + A(1,3) * A(3,2) ) / det
      AINV(2,2) =(   A(1,1) * A(3,3) - A(1,3) * A(3,1) ) / det
      AINV(3,2) =( - A(1,1) * A(3,2) + A(1,2) * A(3,1) ) / det
      AINV(1,3) =(   A(1,2) * A(2,3) - A(1,3) * A(2,2) ) / det
      AINV(2,3) =( - A(1,1) * A(2,3) + A(1,3) * A(2,1) ) / det
      AINV(3,3) =(   A(1,1) * A(2,2) - A(1,2) * A(2,1) ) / det
      RETURN
      END

      REAL*8 FUNCTION determinant(A)
      REAL* 8 A(3,3)

      determinant=a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)  &
            +a(1,3)*a(2,1)*a(3,2)-a(3,1)*a(2,2)*a(1,3)  &
            -a(1,1)*a(3,2)*a(2,3)-a(2,1)*a(1,2)*a(3,3)
      END
