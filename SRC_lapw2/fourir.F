SUBROUTINE fourir(nnlo,nwave)
  !                                  
  ! 1. Evaluates fourir expansion of charge density and kinetic energy density
  !    in interstitial regions by FFT
  ! 2. Calculates surface part of forces from interstitial kinetic energy 
  !    density
  ! 3. Eq.(B9) MBSSN-01 (subroutine fsumai1)
  ! 8/8-2000 GKHM
  USE param
  USE defs
  USE parallel
  USE char; USE com
  USE reclat
  USE struk
  IMPLICIT NONE
#ifdef Parallel
  INCLUDE 'mpif.h'
#endif


  INTEGER                  :: iff1,iff2,iff3,iff1t,iff2t,iff3t
  INTEGER                  :: ios,nwave,kmax(3),i,j,jj,jx,n,ne,kkk,jmin,je
  INTEGER                  :: iproc,num,isig,je1,ii,i1,i2,i3,ia1,ia2,nnlo
  INTEGER                  :: jsub, jnum, j1, j2
  INTEGER,ALLOCATABLE      :: keigen(:,:), jl(:), ju(:)
#if defined (SCSL)
  INTEGER                  :: isys(0:1)
#endif
  REAL*8                   :: e(nume),emist,s,t,z,e_buf,fftin,volin
  CHARACTER*10     kname
  LOGICAL force
!_REAL      REAL*8,allocatable ::           A(:,:)
!_COMPLEX      COMPLEX*16,allocatable ::           A(:,:)
  COMPLEX*16, ALLOCATABLE  :: fft(:,:,:),tp(:,:,:),sumfft(:,:,:),vsumfft(:,:,:),tloc(:,:,:)
  DOUBLE PRECISION,ALLOCATABLE :: dwork(:)
  COMPLEX*16, ALLOCATABLE  :: cwork(:)
  
  COMPLEX*16,ALLOCATABLE   :: rho1(:),vrho1(:),ekin1(:)
  COMPLEX*16,ALLOCATABLE   :: rhok_buf(:),vrhok_buf(:),ekink_buf(:)
  COMPLEX*16,ALLOCATABLE   :: rhok(:),vrhok(:),ekink(:)
!_REAL  REAL*8,allocatable       :: a_buf(:)
!_COMPLEX  COMPLEX*16,allocatable        :: a_buf(:)


  ALLOCATE (a(nmat,nume),a_buf(nmat))
  ALLOCATE (rho1(nwave*nsym),vrho1(nwave*nsym),ekin1(nwave*nsym))
  ALLOCATE (rhok(nwave),vrhok(nwave),ekink(nwave))
  ALLOCATE (rhok_buf(nwave),vrhok_buf(nwave),ekink_buf(nwave))
#ifdef Parallel
  ALLOCATE (jl(0:npe-1), ju(0:npe-1))
#endif

#if defined (SCSL)
  isys(0) = 1
#endif
  
  vrhok(1:nwave)=zeroc
  rhok(1:nwave)=zeroc
  vrhok_buf(1:nwave)=zeroc
  rhok_buf(1:nwave)=zeroc
  vrho1(1:nwave*nsym)=zeroc
  rho1(1:nwave*nsym)=zeroc
  ekink_buf(1:nwave)=zeroc
  ekin1(1:nwave*nsym)=zeroc
  
  ALLOCATE(keigen(3,nmat))

  kmax(1:3)=0
          
!.....FIND MAX RECIPROCAL LATTICE VECTORS
  IF(myid.EQ.0) THEN
     REWIND 10
     DO i=1,nat
        READ(10) emist
        READ(10) emist
     ENDDO
     READ(10,end=998) s,t,z,kname,n,ne
     n=n-nnlo
     READ(10) (keigen(1,i),keigen(2,i),keigen(3,i),i=1,n)
     REWIND 10
     DO i=n,1,-1
        kmax(1)=MAX(keigen(1,i),kmax(1))
        kmax(2)=MAX(keigen(2,i),kmax(2))
        kmax(3)=MAX(keigen(3,i),kmax(3))
     ENDDO
     WRITE(6,*) 'max k indices:',kmax(1),kmax(2),kmax(3)
     kmax(1)=kmax(1)+1
     kmax(2)=kmax(2)+1
     kmax(3)=kmax(3)+1
  ENDIF
#ifdef Parallel
  CALL MPI_BCAST(kmax, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif  
!     set iff values for FFT calculation
  iff1=2*kmax(1)
  iff2=2*kmax(2)
  iff3=2*kmax(3)
  iff1=(iff1+1)*2
  iff2=(iff2+1)*2
  iff3=(iff3+1)*2
  iff1t=iff1
  iff2t=iff2
  iff3t=iff3
  CALL IFFLIM(iff1t,iff1)
  CALL IFFLIM(iff2t,iff2)
  CALL IFFLIM(iff3t,iff3)
  IF(myid.eq.0) WRITE(6,*) 'n,iff1,iff2,iff3',n,iff1,iff2,iff3
  
  ALLOCATE (fft(iff1,iff2,iff3),tp(iff1,iff2,iff3), &
       sumfft(iff1,iff2,iff3),vsumfft(iff1,iff2,iff3), &
       tloc(iff1,iff2,iff3))
#if defined (SCSL)
  ALLOCATE (dwork(2 * MAX(iff1, iff2, iff3)),                        &
            cwork((2*iff1+256) + (2*iff2+256) + (2*iff3+256)))
#else
  ALLOCATE (dwork(4*(iff1+iff2+iff3)+15),cwork(iff1+iff2+iff3))
#endif
  fft=zeroc
  tp=zeroc
  sumfft=zeroc
  vsumfft=zeroc
  tloc=zeroc
  
  force=.false.
  IF(modus.EQ.'FOR ')   force=.TRUE.
  IF(myid.EQ.0) WRITE(6,*) '    FORCE-CALCULATION:',force
  
  DO i=1,nat
     READ(10) emist
     READ(10) emist
  ENDDO
  kkk=0

!.....SUM OVER DIFFERENT K-POINTS AND THEIR WEIGHTS
  DO 
     if(myid.EQ.0) READ(10,IOSTAT=ios) S,T,Z,KNAME,N,NE
#ifdef Parallel
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
     CALL MPI_BCAST(ios,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
     CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
     CALL MPI_BCAST(ne,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
     CALL MPI_BCAST(s,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
     CALL MPI_BCAST(t,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
     CALL MPI_BCAST(z,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
#endif
     IF (ios /= 0) EXIT

!.....READ IN WAVEFUNCTIONS
     n=n-nnlo
     KKK=KKK+1
     IF(NE.GT.NUME) GOTO 900
     IF(myid.EQ.0) READ(10) (KEIGEN(1,I),KEIGEN(2,I),KEIGEN(3,I),I=1,N)
#ifdef Parallel
     CALL MPI_BCAST(keigen, 3*n, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
     JMIN=1
     JE=0
! Einlesen der Eigenwerte und Eigenvektoren
     if (myid == 0) then
        DO
           READ(10) num,e_buf
           READ(10,IOSTAT=ios) (a_buf(i),i=1,n)
           if (ios /= 0) stop 'end.'
           e(je+1)=e_buf
           a(1:n,je+1)=a_buf(1:n)
           IF(e_buf.LT.emin) jmin=jmin+1
           IF(e_buf.LE.ef) je=je+1
           IF(e_buf.GE.emin.AND.e_buf.LE.ef) THEN
              IF(ABS(weigh(kkk,num)).LT.1.d-8) je=je-1
           ENDIF
           IF(num.EQ.ne) EXIT
        ENDDO         ! loop over eigenvalues
     endif
#ifdef Parallel
     CALL MPI_BCAST(jmin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     CALL MPI_BCAST(je, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
! Determine j region processor iproc is responsible for
     jsub = (je-jmin+npe)/npe
     do iproc=0,npe-1
	jl(iproc) = iproc*jsub + 1
	ju(iproc) = min((iproc+1)*jsub,je)
     end do
     IF(myid.eq.0) THEN 
        DO iproc=1,npe-1
           j1 = jl(iproc)
           j2 = ju(iproc)
           jnum = j2-j1+1
	   if (jnum < 1) cycle
           CALL MPI_SEND(e(j1), jnum, MPI_DOUBLE_PRECISION, iproc, 0, MPI_COMM_WORLD, ierr) 
!_REAL           CALL MPI_SEND(a(1,j1), nmat*jnum, MPI_DOUBLE_PRECISION, iproc, 0, MPI_COMM_WORLD, ierr)
!_COMPLEX           CALL MPI_SEND(a(1,j1), nmat*jnum, MPI_DOUBLE_COMPLEX, iproc, 0, MPI_COMM_WORLD, ierr)
        ENDDO
     ELSE
	j1 = jl(myid)
	j2 = ju(myid)
	jnum = j2-j1+1
	if (jnum < 1) cycle
        CALL MPI_RECV( e(j1), jnum, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD,statusmpi, ierr)
!_REAL        CALL MPI_RECV( a(1,j1), nmat*jnum, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD,statusmpi, ierr)
!_COMPLEX        CALL MPI_RECV( a(1,j1), nmat*jnum, MPI_DOUBLE_COMPLEX, 0, MPI_ANY_TAG, MPI_COMM_WORLD,statusmpi, ierr)
     ENDIF
     jmin = jl(myid)
     je = ju(myid)
#endif
!     write (60+myid,*) 'jmin, je = ', jmin, je
     IF(je.EQ.0) CYCLE

! FFT eigenvectors and accumulate weighted square in sumfft
     isig=-1
     DO je1=jmin,je
        CALL setfft1(n,iff1,iff2,iff3,a(1,je1),fft,keigen)
#if defined (SCSL)
        CALL zzfft3d(0,iff1,iff2,iff3,0.0D0,0,1,1,0,1,1,cwork,dwork,    &
                     isys)
        CALL zzfft3d(isig,iff1,iff2,iff3,1.0D0,fft,iff1,iff2,fft,iff1,  &
                     iff2,cwork,dwork,isys)
#else
        CALL c3fft(iff1,iff2,iff3,fft,iff1,iff2,isig, &
             cwork,dwork,ierr)
#endif
        CALL sumupfft(iff1,iff2,iff3,weigh(kkk,je1),fft,sumfft)
        CALL sumupfft(iff1,iff2,iff3,-weigh(kkk,je1)*e(je1),fft, &
             vsumfft)
        IF (force) THEN
           CALL sumupfft(iff1,iff2,iff3,-weigh(kkk,je1)*e(je1),fft, &
                tloc)
           DO ii=1,3
              CALL setfft2(n,iff1,iff2,iff3,a(1,je1),tp,ii,s,t,z,keigen)
#if defined (SCSL)
              CALL zzfft3d(isig,iff1,iff2,iff3,1.0D0,tp,iff1,iff2,      &
                           tp,iff1,iff2,cwork,dwork,isys)
#else
              CALL c3fft(iff1,iff2,iff3,tp,iff1,iff2,isig, &
                   cwork,dwork,ierr)
#endif
              CALL sumupfft(iff1,iff2,iff3,weigh(kkk,je1),tp,tloc)
           ENDDO
        ENDIF
     ENDDO
  ENDDO
!.....END LOOP OVER ALL K-POINTS
!      do i3=1,iff3
!         do i2=1,iff2
!            do i1=1,iff1
!               tloc(i1,i2,i3)=tloc(i1,i2,i3)+vsumfft(i1,i2,i3)
!            enddo
!         enddo
!      enddo

! Back transform cumulative grids
  isig=1
#if defined (SCSL)
  CALL zzfft3d(isig,iff1,iff2,iff3,1.0D0,sumfft,iff1,iff2,sumfft,iff1,  &
               iff2,cwork,dwork,isys)
  CALL zzfft3d(isig,iff1,iff2,iff3,1.0D0,vsumfft,iff1,iff2,vsumfft,iff1,&
               iff2,cwork,dwork,isys)
#else
  CALL c3fft(iff1,iff2,iff3,sumfft,iff1,iff2,isig, &
       cwork,dwork,ierr)
  CALL c3fft(iff1,iff2,iff3,vsumfft,iff1,iff2,isig, &
       cwork,dwork,ierr)
#endif
  IF(force) THEN
#if defined (SCSL)
     CALL zzfft3d(isig,iff1,iff2,iff3,1.0D0,tloc,iff1,iff2,tloc,iff1,   &
                  iff2,cwork,dwork,isys)
#else
     CALL c3fft(iff1,iff2,iff3,tloc,iff1,iff2,isig, &
          cwork,dwork,ierr)
#endif
  ENDIF
  fftin = 1.d0/dble(iff1*iff2*iff3)
  DO i3=1,iff3
     DO i2=1,iff2
        DO i1=1,iff1
           sumfft(i1,i2,i3)=sumfft(i1,i2,i3)*fftin
           vsumfft(i1,i2,i3)=vsumfft(i1,i2,i3)*fftin
        ENDDO
     ENDDO
  ENDDO
  IF (force) THEN
     DO i3=1,iff3
        DO i2=1,iff2
           DO i1=1,iff1
              tloc(i1,i2,i3)=tloc(i1,i2,i3)*fftin
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  call getfft(nwave,iff1,iff2,iff3,rho1,sumfft,kmax)
  call getfft(nwave,iff1,iff2,iff3,vrho1,vsumfft,kmax)
  if(force) then
     call getfft(nwave,iff1,iff2,iff3,ekin1,tloc,kmax)
  endif

!.....SUM OVER ALL RECPR. LATTIC VECTORS OF ONE STAR
  ia1=1
  volin = 1.0d0/vol
  DO j=1,nwave
     ia2=ia1+inst(j)-1
     DO jj=ia1,ia2
        rhok_buf(j)=rhok_buf(j)+rho1(jj)*CONJG(tauk(jj))
        vrhok_buf(j)=vrhok_buf(j)+vrho1(jj)*CONJG(tauk(jj))
        ekink_buf(j)=ekink_buf(j)+ekin1(jj)*CONJG(tauk(jj))
     ENDDO
     ia1=ia2+1
     vrhok_buf(j)=vrhok_buf(j) * volin
     rhok_buf(j)=rhok_buf(j) * volin
     ekink_buf(j)=ekink_buf(j) * volin
  ENDDO
!cccc
#ifdef Parallel
  CALL MPI_REDUCE(rhok_buf, rhok, nwave, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr ) 
  CALL MPI_REDUCE(vrhok_buf, vrhok, nwave, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr ) 
  CALL MPI_REDUCE(ekink_buf, ekink, nwave, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr ) 
  IF(myid.NE.0) GOTO 998
#else
  DO j=1,nwave
     vrhok(j)=vrhok_buf(j)
     rhok(j)=rhok_buf(j)
     ekink(j)=ekink_buf(j)
  ENDDO
#endif
  WRITE(8,*) '   VALENCE CHARGE DENSITY IN INTERSTITIAL '
  WRITE(8,*)
  WRITE(8,2061)  NWAVE
  WRITE(8,2071)  ( (KZZ(JX,J),JX=1,3),RHOK(J), J=1,NWAVE)
  WRITE(28,*) '   VALENCE CHARGE DENSITY IN INTERSTITIAL '
  WRITE(28,*)
  WRITE(28,2061)  NWAVE
  WRITE(28,2071)  ( (KZZ(JX,J),JX=1,3),VRHOK(J), J=1,NWAVE)
  WRITE(6,204) NWAVE,kmax

     close(8)

  if(force) call fsumai1(ekink,nwave)
  
998 CONTINUE
#ifdef Parallel
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
#endif
  DEALLOCATE(a,a_buf)
  DEALLOCATE(rhok,vrhok,ekink,rhok_buf,vrhok_buf,ekink_buf)
  DEALLOCATE(fft,tp,sumfft,vsumfft,tloc,dwork,cwork)
  DEALLOCATE(rho1,vrho1,ekin1)
  DEALLOCATE(keigen)
  RETURN
  !
  !        Error messages
  !
  !
900 CALL OUTERR('FOURIR','VALUE OF PARAMETER NUME TOO SMALL FOR NE.')
  STOP 'FOURIR2 - Error'
  !
204 FORMAT(I10,' FOURIER COEFFICIENTS CALCULATED',/ &
         ' LARGEST COMPONENTS:',3I6)
2061 FORMAT(1X,'        ',I10,1X,'NUMBER OF PW')
2071 FORMAT(3X,3I5,2E19.12)
  
END SUBROUTINE fourir
