      subroutine constrained_forces(numeq,nat,A,b,F,pos,mult)
!---  calculate constrained forces 
!     constraint is linear of the form  xA=b; here
!     nat - number of notequivalent atoms
!     numeq - number of constraints (columns of A(:,:))
!     x - coordinates of notequivalent atoms, array posneq,  dim=(3*nat) 
!     A - constant matrix specified by user, on exit contains 
!         orthonormalized columns,  dim=(3*nat,numeq)
!     b - constant vector  specified by user, unchanged on exit, dim=(numeq)
!         used only to check that atomic positions satisfy constraint   
!     F - forces, constrained on exit
      
      integer numeq,nat,mult(nat)
      integer i,j,jdp,info
      integer,allocatable :: ldp(:)
      real*8 A(3*nat,numeq), b(numeq), pos(3,*), F(3,nat)
      real*8,allocatable :: tau(:),posneq(:,:),wrk(:)
      real*8 tmp
      real*8 dnrm2
      external dgemv,dnrm2,dgeqrf,dorgqr
      
      if(numeq.lt.0) return
      allocate(posneq(3,nat),wrk(3*nat),tau(numeq),ldp(numeq))

!     check if atomic positions satisfy constraint
      j=1
      do i=1,nat
         posneq(1,i)=pos(1,j)
         posneq(2,i)=pos(2,j)
         posneq(3,i)=pos(3,j)
         j=j+mult(i)
      enddo
      do i=1,numeq
         wrk(i)=b(i)
      enddo
!       write(6,*) (wrk_constr(i), i=1,3*nat)
      call dgemv('T',3*nat,numeq,1d0,A,3*nat,posneq,1,-1d0,wrk,1)
      tmp=dnrm2(numeq,wrk,1)
      if(sqrt(tmp/nat) .gt. 1d-5) then
         write(6,*) "WARNING! atom positions don't satisfy constraint", sqrt(tmp/nat)
         print*, "WARNING! atom positions don't satisfy constraint", sqrt(tmp/nat)
      endif  
!     orthonormalize columns of A using QR decomposition
      call dgeqrf(3*nat,numeq,A,3*nat,tau,wrk,3*nat,info)
!     check for linear dependent equations, nullify these if exist
      ldp(1:numeq)=0
      jdp=0
      do i=1,numeq
         if(abs(A(i,i)).lt.1d-10) then
           ldp(i)=1
           jdp=1
         endif
      enddo   
      if(jdp.eq.1) then
         write(6,*) 'WARNING!', &
           ' Linear dependent constraints have been found:'
         do i=1,numeq
            if(ldp(i).eq.1) write(6,'(i6)') i
         enddo   
      endif
!     generate Q      
      call dorgqr(3*nat,numeq,numeq,A,3*nat,tau,wrk,3*nat,info)
      if(jdp.eq.1) then ! nullify linear dependent vectors 
         do i=1,numeq
            if( ldp(i).eq.1 ) A(1:3*nat,i)=0d0
         enddo
      endif   
!     project F on Q, result in wrk_constr
      call dgemv('T',3*nat,numeq,1d0,A,3*nat,F,1,0d0, wrk,1)
!     remove from F projection on Q
      call dgemv('N',3*nat,numeq,-1d0,A,3*nat,wrk,1,1d0,F,1)
      deallocate(posneq,wrk,tau,ldp)
      return
      end
