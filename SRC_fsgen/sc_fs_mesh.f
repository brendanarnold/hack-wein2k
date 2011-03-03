!     calc hex k vectors for fermi surfaces in 3 planes 
      dimension a(2,2),b(2,2),r(2),rk(2),ik(2,0:30)
!
      write(*,*) 'generates fermi surface mesh in sc case:'
      write(*,*) 'plane Gamma-X-M:     1'
      write(*,*) 'x-z plane:           2'
      write(*,*) 'plane Gamma-K-A:     3'
      write(*,*) 'monoclinic plane Gamma-Z-E-A:     4'
      write(*,*) 'orthorhombic xy-plane (up to 1):     5'
      write(*,*) 'select 1-5:'
      read(*,*) ind
      if(ind.eq.1) then
        write(*,*)'select mesh divisions: (30,36,42,...)'
        read(*,*) mesh
        i=0
        do iy=0,mesh
        do ix=0,mesh
        i=i+1
        if(i.eq.1)then
          write(2,'(i10,4i5,3f5.2,5x,2i4)') i,0,0,0,mesh*2,1.0,-7.,1.5, &
                                     mesh/2+1,mesh/3+1
        else
          write(2,'(i10,4i5,f5.2)') i,iy,iy,ix,mesh*2,1.0
        endif
        enddo
        enddo
        write(2,'(a)') 'END'
      else if(ind.eq.2) then
        write(*,*)'select mesh divisions: (even nx,nz)'
        read(*,*) meshx,meshz
        i=0
        do iz=0,meshz/2
        do ix=0,meshx/2
        i=i+1
        if(i.eq.1)then
          write(2,'(i10,4i5,3f5.2,5x,2i4)') i,0,0,0,meshx*meshz,1.0,-7., &
                                            1.5,meshx/2+1,meshz/2+1
        else
        write(2,'(i10,4i5,f5.2)') i,ix*meshz,0,iz*meshx,meshx*meshz,1.0
        endif
        enddo
        enddo
        write(2,'(a)') 'END'
      else if(ind.eq.3) then
        write(*,*)'select mesh divisions: (nx (mult of 3), nz (even))'
        read(*,*) meshx,meshz
        i=0
        do iz=0,meshz/2
        do ix=0,meshx/3
        i=i+1
        if(i.eq.1)then
          write(2,'(i10,4i5,3f5.2,5x,2i4)') i,0,0,0,meshx*meshz,1.0,-7., &
                                            1.5,meshx/2+1,meshz/2+1
        else
          write(2,'(i10,4i5,f5.2)') i,ix*meshz,ix*meshz,iz*meshx, &
                                    meshx*meshz,1.0
        endif
        enddo
        enddo
        write(2,'(a)') 'END'
      else if(ind.eq.4) then
        write(*,*)'select mesh divisions: (nx (even), nz (even))'
        read(*,*) meshx,meshz
        i=0
        do iz=0,meshz/2
        do ix=0,meshx/2
        i=i+1
        if(i.eq.1)then
          write(2,'(i10,4i5,3f5.2,5x,2i4)') i,0,0,0,meshx*meshz,1.0,-7., &
                                            1.5,meshx/2+1,meshz/2+1
        else
          write(2,'(i10,4i5,f5.2)') i,ix*meshz,ix*meshz,iz*meshx, &
                                    meshx*meshz,1.0
        endif
        enddo
        enddo
        write(2,'(a)') 'END'
      else if(ind.eq.5) then
        write(*,*)'select mesh divisions: (even nx,ny)'
        read(*,*) meshx,meshy
        i=0
        do iy=0,meshy
        do ix=0,meshx
        i=i+1
        if(i.eq.1)then
          write(2,'(i10,4i5,3f5.2,5x,2i4)') i,0,0,0,meshx*meshy,1.0,-7., &
                                            1.5,meshx+1,meshy+1
        else
        write(2,'(i10,4i5,f5.2)') i,ix*meshy,iy*meshx,0,meshx*meshy,1.0
        endif
        enddo
        enddo
        write(2,'(a)') 'END'
      endif

      end
