      subroutine dsymm2(side,uplo,n,m,alpha,a,lda,b,ldb,beta,c,ldc)
!
!     Version of mm.f adapted for 
!     Symmetric Matrix
!     
!     This Version of a symmetric matrix multiplication is
!     on some computers faster than the Level 3 Blas.
!
!     Parameters:
!     n: number of rows of b and c
!     m: number of columns of b and c
!     a: Matrix A, nXn 
!     lda: leading dimension of a
!     b: Matrix B, nXm
!     ldb: leading dimension of b
!     c: target Matrix C, nXm
!     ldc: leading dimension of c
!
!     ####################################
      implicit none
      character*1 uplo,side
      integer n,m,lda,ldb,ldc
      double precision a(lda,n),b(ldb,m),c(ldb,m)
      double precision alpha,beta

      integer nblock
      parameter(nblock=64)
      double precision zero,one
      parameter (zero=0.0d0,one=1.0d0)

      double precision aa(nblock,nblock)
      double precision factor
      integer nn,i,j,k,ii,jj,kk
      integer imax,kmax
      logical upper

      logical lsame
      external lsame

      if(.NOT.lsame(side, 'L'))  &
         stop 'wrong parameter side in dsymm2'
      if(beta.NE.zero) stop 'wrong parameter beta in dsymm2'
      if(alpha.NE.one) stop 'wrong parameter alpha in dsymm2'
      upper=lsame(uplo, 'U')

      do ii = 1,n,nblock
         imax = min(n,ii+nblock-1)
         factor = zero
	 do kk = 1,n,nblock
            kmax = min(n,kk+nblock-1)

            do k = kk,kmax
               do i = ii,imax

        	  if( (upper.AND.(k.GE.i)).OR. &
                      ((.NOT.upper).AND.(i.GE.k))) THEN
                      aa(i-ii+1,k-kk+1) = a(i,k)
		    else
                      aa(i-ii+1,k-kk+1) = a(k,i)
		  end if

               end do
            end do

	    CALL DGEMM('N','N',imax-ii+1,m,kmax-kk+1, &
                       one,aa,nblock,b(kk,1),ldb, &
                       factor,c(ii,1),ldc)
	    factor=one

         end do
      end do
      end

