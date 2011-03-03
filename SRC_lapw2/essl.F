       subroutine c_transpose(a,lda,m,n0,n,b,ldb)

       integer lda,m,n,ldb
       complex *16 a(lda,n),b(ldb,m)
!
! for IBM essl use:
!       call zgetmo(a,lda,m,n,b,ldb)
! otherwise:
        do ie=n0,n
        do jl=1,m
        b(ie,jl)=a(jl,ie)
        enddo
        enddo
!       
       return
       end
