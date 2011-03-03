       subroutine pdsyr2m (uplo,transa,transb,m,n, &
                          alpha,a,ia,ja,desca,b,ib,jb,descb,&
			  beta,c,ic,jc,descc)
!
! Author: Dieter Kvasnicka
!         Institute of Materials Chemistry
!         Vienna University of Technology, Austria
!         2001-2002
!
!  Purpose
!  =======
!
!  pdsyr2m performs one of the matrix-matrix operations
!
!     op1(C) := alpha*op2( A )*op2( B ) + beta*op1(C),
!
!  where  op1( X ) means that only the UPPER or only the LOWER
!  part of X is accessed
!
!  and  op2( X ) is one of
!
!     op2( X ) = X   or   op2( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with 
!  op2( A ) an m by n matrix, 
!  op2( B ) an n by m matrix, and 
!  C an m by m matrix.
!
!     uplo: 'U'pper or 'L'ower
!     transa,transb: op2( X ) = X   or   op2( X ) = X', with X=A,B
!     m: number of rows of op2(A) and op2(B) and of
!        rows and columns of C
!     n: number of columns of op2(A) and op2(B)
!     alpha
!     op2(A): mXn matrix
!     op2(B): mXn matrix
!     beta
!     c: symmetric mXm matrix, 
!        stored either above (uplo='U') or below (uplo='L') the
!        diagonal
!     ia,ja,ib,jb,ic,jc: first and second index of A, B, C, respectively
!     desca,descb,descc: matrix descriptor of A, B, C, respectively
!                        see scalapack users guide at 
!                        www.netlib.org/scalapack/slug/index.html
!
       implicit none
       character*1 uplo,transa,transb
       integer m,n,ia,ja,ib,jb,ic,jc
       integer desca(9), descb(9), descc(9)
       double precision alpha,beta
       double precision a(*), b(*), c(*)
#ifdef Parallel

       integer nblock, nprow, npcol, myrow, mycol
       double precision half
       parameter (half=0.5d0)

       integer i,j,k,ii,jj,kk,imax,jmax,kmax
       logical upper,atrans,btrans
       double precision t

       logical lsame
       external lsame

       upper=lsame(uplo, 'U')
       atrans=lsame(transa, 'T')
       btrans=lsame(transb, 'T')

       nblock=desca(5)
       if(desca(6).ne.nblock .or. descb(5).ne.nblock .or. &
	  descb(6).ne.nblock .or. descc(5).ne.nblock .or. &
	  descc(6).ne.nblock) &
	 stop 'wrong parameter in matrix descriptors in pdsyr2m'
       call blacs_gridinfo( desca(2), nprow, npcol, myrow, mycol)
       nblock=nblock*npcol

       if(upper .and. (atrans .and. .not.btrans)) then
         do jj=1,m,nblock
           jmax=min(m,jj+nblock-1)
	   call pdsyr2k(uplo,transa,jmax-jj+1,n,half*alpha,a,1,jj,desca, &
			b,1,jj,descb,beta,c,jj,jj,descc)
           if(jj.ne.1) &
             call pdgemm(transa,transb,jj-1,jmax-jj+1,n, &
                    alpha,a,1,1,desca,b,1,jj,descb, &
                    beta,c,1,jj,descc)
         end do
       else if(upper .and. (.not.atrans .and. btrans)) then
         do jj=1,m,nblock
           jmax=min(m,jj+nblock-1)
	   call pdsyr2k(uplo,transa,jmax-jj+1,n,half*alpha,a,jj,1,desca, &
			b,jj,1,descb,beta,c,jj,jj,descc)
           if(jj.ne.1) &
             call pdgemm(transa,transb,jj-1,jmax-jj+1,n, &
                    alpha,a,1,1,desca,b,jj,1,descb, &
                    beta,c,1,jj,descc)
         end do
       else
         stop 'wrong parameter in pdsyr2m'
       end if
#endif
       end

