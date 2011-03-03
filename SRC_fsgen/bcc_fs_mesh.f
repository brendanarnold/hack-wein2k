!     calc hex k vectors for fermi surfaces in 3 planes 
      dimension a(2,2),b(2,2),r(2),rk(2),ik(2,0:30)
!
      write(*,*) 'generates fermi surface mesh in bcc case:'
      write(*,*) 'plane Gamma-H-N:                   1'
      write(*,*) 'plane Gamma-H-N + normal-offset    2'
      write(*,*) 'plane Gamma-N-H-N:                 3'
      write(*,*) 'select 1-3:'
      read(*,*) ind
      if(ind.eq.1) then
        write(*,*)'select mesh divisions: (30,36,42,...)'
        write(*,*)'set nx,ny to nmesh+1; invers=1,flip=0'
        read(*,*) mesh
        i=0
        do iy=0,mesh
        do ix=0,mesh
        i=i+1
        if(i.eq.1)then
          write(2,'(i10,4i5,3f5.2,5x,2i4)') i,0,0,0,mesh*2,1.0,-7.,1.5, &
                                     mesh/2+1,mesh/3+1
        else
          write(2,'(i10,4i5,f5.2)') i,iy,iy,2*ix,mesh*2,1.0
        endif
        enddo
        enddo
        write(2,'(a)') 'END'
      else if(ind.eq.2) then
        write(*,*)'select mesh divisions: (30,36,42,...); offset (0-mesh)'
        write(*,*)'set nx,ny to nmesh+1; invers=1,flip=0'
        read(*,*) mesh,moff
        i=0
        do iy=0,mesh
        do ix=0,mesh
        i=i+1
        if(i.eq.1)then
          write(2,'(i10,4i5,3f5.2,5x,2i4)') i,moff,-moff,0,mesh*2,1.0,-7., &
                                            1.5,mesh/2+1,mesh/3+1
        else
        write(2,'(i10,4i5,f5.2)') i,iy+moff,iy-moff,2*ix,mesh*2,1.0
        endif
        enddo
        enddo
        write(2,'(a)') 'END'
      else if(ind.eq.3) then
        write(*,*)'select mesh divisions: (nx (even))'
        write(*,*)'set nx,ny to nmesh+1; invers=1,flip=1'
       read(*,*) mesh
        i=0
        do iz=0,mesh
        do ix=iz,mesh
        i=i+1
        if(i.eq.1)then
          write(2,'(i10,4i5,3f5.2,5x,2i4)') i,0,0,0,mesh*2,1.0,-7., &
                                            1.5,mesh/2+1,mesh/2+1
        else
          write(2,'(i10,4i5,f5.2)') i,ix+iz,ix-iz,0, &
                                    mesh*2,1.0
        endif
        enddo
        enddo
        write(2,'(a)') 'END'
      endif

      end
