SUBROUTINE read_vec(more_kpoints,nemin,nemax,kkk,n,jatom,ne,etot,so,nspin1)
  USE param
  USE parallel
  USE com; USE char
  USE lo
  USE xa; USE xa3
  IMPLICIT NONE
#ifdef Parallel
  INCLUDE 'mpif.h'
#endif

  LOGICAL               :: more_kpoints,so
  INTEGER               :: ne,nemin,nemax,kkk,n,jatom,nspin1,n_pw
  INTEGER               :: iblock_tot,ios,i,num
  REAL*8                :: etot,s,t,z,vnorm1
  REAL*8                :: vnorm(nume)

!_REAL  REAL*8,ALLOCATABLE      ::    a_buf(:)
!_COMPLEX  COMPLEX*16,ALLOCATABLE  ::    a_buf(:)
  INTEGER,ALLOCATABLE   :: kx_buf(:),ky_buf(:),kz_buf(:)

  CHARACTER *10    BNAME                                

#ifdef Parallel
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

  IF (myid.EQ.0) READ(10,IOSTAT=ios) s,t,z,bname,n,ne
  more_kpoints=.FALSE.

#ifdef Parallel
  CALL MPI_BCAST(ios,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(ne,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(s,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(t,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(z,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
#endif
  IF (ios /= 0) RETURN
  more_kpoints=.TRUE.
!  IF(myid.EQ.0) WRITE(6,205) s,t,z,n,ne,bname         

  iblock_tot=(n-(nlo+nlon+nlov)-1)/iblock+1
  ibpp=(n-(nlo+nlon+nlov)-1)/(iblock*npe)+1

  ALLOCATE(kx_buf(ibpp*iblock*npe),ky_buf(ibpp*iblock*npe),kz_buf(ibpp*iblock*npe))
  ALLOCATE(a_buf(ibpp*iblock*npe))
  kx_buf=0 ; ky_buf=0 ; kz_buf=0 



  kkk=kkk+1                                                         
  IF(jatom.EQ.1) nk=nk+1                                            
  IF(myid.EQ.0) THEN
     WRITE(31,205) s,t,z,n,ne,bname                                  
     READ(10) (kx_buf(i),ky_buf(i),kz_buf(i),i=1,n-(nlo+nlon+nlov)),(kxlo(i),kylo(i),kzlo(i),i=1,(nlo+nlon+nlov))
  ENDIF
#ifdef Parallel
  CALL MPI_SCATTER(kx_buf, ibpp*iblock, MPI_INTEGER, kx, ibpp*iblock, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_SCATTER(ky_buf, ibpp*iblock, MPI_INTEGER, ky, ibpp*iblock, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_SCATTER(kz_buf, ibpp*iblock, MPI_INTEGER, kz, ibpp*iblock, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#else
  DO i=1,n-(nlo+nlon+nlov)
     kx(i)=kx_buf(i)
     ky(i)=ky_buf(i)
     kz(i)=kz_buf(i)
  ENDDO
#endif
  DO i=1,ibpp*iblock
     bkx(i)=(s+kx(i))
     bky(i)=(t+ky(i))
     bkz(i)=(z+kz(i))
  ENDDO
  IF(myid.EQ.0) THEN
     DO i=1,nlo+nlon+nlov
        bkxlo(i)=(s+kxlo(i))
        bkylo(i)=(t+kylo(i))
        bkzlo(i)=(z+kzlo(i))
     ENDDO
  ENDIF

!ccc mod for so
  vnorm1=1.d0
  IF(so.AND.(nspin1.EQ.1)) vnorm1=0.5d0
  READ(12,202,iostat=ios) (vnorm(i),i=1,ne)
  IF (ios /= 0) THEN
     DO i=1,ne
        vnorm(i)=1.d0
     ENDDO
  ENDIF
  nemin=1                                                           
  nemax=0                                                           
  DO
     IF(myid.EQ.0) THEN
        READ(10) num,e(num)                                               
        READ(10) (a_buf(i),i=1,n-(nlo+nlon+nlov)),(a_lo(i,num),i=1,(nlo+nlon+nlov))
        IF(e(num).LT.emin) nemin=nemin+1                                  
        IF(e(num).LE.ef) nemax=nemax+1                              
        IF(e(num).GE.emin.AND.e(num).LE.ef) THEN
           weight(num)=weigh(kkk,num)
           IF(ABS(weight(num)).GT.1.d-8.OR.modus.EQ.'QTL  ') THEN
              !ccc modifications for so
              IF(jatom.EQ.1) THEN
                 etot=etot+e(num)*weight(num)*vnorm(num)*vnorm1            
                 xwt=xwt+weight(num)*vnorm(num)*vnorm1
                 NB(NK)=NB(NK)+1
              ENDIF
           ELSE
              nemax=nemax-1
           ENDIF
        ENDIF
     ENDIF
#ifdef Parallel
  CALL MPI_BCAST(num,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!_REAL  CALL MPI_SCATTER(a_buf,ibpp*iblock,MPI_DOUBLE_PRECISION,a(1,num),ibpp*iblock,MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)
!_COMPLEX  CALL MPI_SCATTER(a_buf,ibpp*iblock,MPI_DOUBLE_COMPLEX,a(1,num),ibpp*iblock,MPI_DOUBLE_COMPLEX,0, MPI_COMM_WORLD, ierr)
#else
     a(1:ibpp*iblock,num)=a_buf(1:ibpp*iblock)
#endif
     IF(num.EQ.ne) EXIT
  ENDDO
#ifdef Parallel
  CALL MPI_BCAST(nemin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(nemax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

  DEALLOCATE(kx_buf,ky_buf,kz_buf,a_buf)
  RETURN
202   FORMAT(4e20.12)
205 FORMAT(/,1X,' K-POINT:',3F8.4,1X,I5,I4,2X,A10)                      
END SUBROUTINE read_vec


