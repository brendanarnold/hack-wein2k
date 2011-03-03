      subroutine zhemm2(side,uplo,n,m,alpha,a,lda,b,ldb,beta,c,ldc)
!
!     Version of mm.f adapted for 
!     Complex Hermitian Matrix
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
      complex*16 a(lda,n),b(ldb,m),c(ldb,m)
      complex*16 alpha,beta

      integer nblock
      parameter(nblock=64)
      complex*16 czero,cone
      parameter (czero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0))

      complex*16 aa(nblock,nblock)
      complex*16 factor
      integer i,k,ii,kk
      integer imax,kmax
      logical upper

      logical lsame
      external lsame

      if(.NOT.lsame(side, 'L'))  &
         stop 'wrong parameter side in zhemm2'
      if(beta.NE.czero) stop 'wrong parameter beta in zhemm2'
      if(alpha.NE.cone) stop 'wrong parameter alpha in zhemm2'
      upper=lsame(uplo, 'U')

!      do j=1,m
!	do i=1,n
!	  c(i,j)=czero
!	end do
!      end do
      do ii = 1,n,nblock
         imax = min(n,ii+nblock-1)
         factor = czero
	 do kk = 1,n,nblock
            kmax = min(n,kk+nblock-1)

            do k = kk,kmax
               do i = ii,imax

        	  if( (upper.AND.(k.GE.i)).OR. &
                      ((.NOT.upper).AND.(i.GE.k))) THEN
                      aa(i-ii+1,k-kk+1) = a(i,k)
		    else
                      aa(i-ii+1,k-kk+1) = dconjg(a(k,i))
		  end if

               end do
            end do

	    CALL ZGEMM('N','N',imax-ii+1,m,kmax-kk+1, &
                       cone,aa,nblock,b(kk,1),ldb, &
                       factor,c(ii,1),ldc)
	    factor=cone

         end do
      end do
      end

