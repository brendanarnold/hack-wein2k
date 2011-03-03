   subroutine read_octave(strname,domult)

     use nstruct     

     character*50 strname
     character*70 line
     logical domult

      open(unit=10,file='struct.write',status='unknown')


       call jumpto('nat')
       read(10,*) n_nat

       call init_nstruct(n_nat)

       call jumpto('a')
       read(10,*) n_aa,n_bb,n_cc
       call jumpto('alpha')
       read(10,*) n_alpha
       call jumpto('aname')
       do i=1,n_nat
          read(10,*)
          read(10,'(a)') n_aname(i)
       enddo
       call jumpto('lattic')
       read(10,*)
       read(10,'(a)') n_lattic
       call jumpto('pos')
       do i=1,n_nat
          read(10,*) n_pos(1:3,i)
       enddo
       call jumpto('jrj')
       do i=1,n_nat
          read(10,*) n_jrj(i)
       enddo
       call jumpto('r0')
       do i=1,n_nat
          read(10,*) n_r0(i)
       enddo
       call jumpto('rmt')
       do i=1,n_nat
          read(10,*) n_rmt(i)
       enddo
       call jumpto('zz')
       do i=1,n_nat
          read(10,*) n_zz(i)
       enddo

      n_ndif=n_nat        

      do i=1,n_nat
         n_iatnr(i)=-i
      enddo

      n_isplit(1:n_nat)=8

      do i=1,n_nat
         n_rotloc(1:3,1:3,i)=reshape((/1.0d0,0.0d0,0.0d0,&
                                     0.0d0,1.0d0,0.0d0,&
                                     0.0d0,0.0d0,1.0d0/),&
                                     (/3,3/))
      enddo
      n_nsym=0
      n_iord=0
      n_mult(1:n_nat)=1
      n_iz(1:3,1:3,1)=reshape((/1.0d0,0.0d0,0.0d0,&
                                0.0d0,1.0d0,0.0d0,&
                                0.0d0,0.0d0,1.0d0/),&
                                 (/3,3/))
      n_inum(1)=1
      n_tau(1:3,1)=0.0d0

      n_title='octavetmp'

      if (domult) call makemult

   end subroutine read_octave

    
    subroutine write_octave(strname)

     use nstruct     
     use struct     

     character*50 strname

      open(unit=10,file='struct.read',status='unknown')

      write(10,'(a)') '# Created rwoctave'
      write(10,'(2a)') '# name: ', trim(strname)
      write(10,'(a)') '# type: struct'
      write(10,'(a)') '# length: 12'
      write(10,'(a)') '# name: a'
      write(10,'(a)') '# type: matrix'
      write(10,'(a)') '# rows: 1'
      write(10,'(a)') '# columns: 3'
      write(10,'(3f16.9)') n_aa,n_bb,n_cc
      write(10,'(a)') '# name: alpha'
      write(10,'(a)') '# type: matrix'
      write(10,'(a)') '# rows: 1'
      write(10,'(a)') '# columns: 3'
      write(10,'(3f16.9)') n_alpha
      write(10,'(a)') '# name: aname'
      write(10,'(a)') '# type: string'
      write(10,'(a,i6)') '# elements:',n_nat 
      do i=1,n_nat
         write(10,'(a,i3)') '# length:',10
         write(10,'(a)') n_aname(i) 
      enddo
      write(10,'(a)') '# name: jrj'
      write(10,'(a)') '# type: matrix'
      write(10,'(a,i6)') '# rows:', n_nat
      write(10,'(a)') '# columns: 1'
      write(10,'(i6)') n_jrj
      write(10,'(a)') '# name: lat2car'
      write(10,'(a)') '# type: matrix'
      write(10,'(a)') '# rows: 3'
      write(10,'(a)') '# columns: 3'
      write(10,'(3f16.9))') transpose(br1_dir)
      write(10,'(a)') '# name: brlat'
      write(10,'(a)') '# type: matrix'
      write(10,'(a)') '# rows: 3'
      write(10,'(a)') '# columns: 3'
      write(10,'(3f16.9))') transpose(br2_dir)
      write(10,'(a)') '# name: lattic'
      write(10,'(a)') '# type: string'
      write(10,'(a)') '# elements: 1'
      write(10,'(a,i3)') '# length:',len(trim(lattic))
      write(10,'(a)') trim(lattic)
      write(10,'(a)') '# name: nat'
      write(10,'(a)') '# type: scalar'
      write(10,'(i6)') n_nat     
      write(10,'(a)') '# name: pos'
      write(10,'(a)') '# type: matrix'
      write(10,'(a,i6)') '# rows:',n_nat
      write(10,'(a)') '# columns: 3'
      do i=1,n_nat
         write(10,'(3f16.9)') n_pos(1:3,i)
      enddo
      write(10,'(a)') '# name: r0'
      write(10,'(a)') '# type: matrix'
      write(10,'(a,i6)') '# rows:',n_nat
      write(10,'(a)') '# columns: 1'
      do i=1,n_nat
         write(10,'(f16.9)') n_r0(i)
      enddo
      write(10,'(a)') '# name: rmt'
      write(10,'(a)') '# type: matrix'
      write(10,'(a,i6)') '# rows:',n_nat
      write(10,'(a)') '# columns: 1'
      do i=1,n_nat
         write(10,'(f16.9)') n_rmt(i)
      enddo
      write(10,'(a)') '# name: zz'
      write(10,'(a)') '# type: matrix'
      write(10,'(a,i6)') '# rows:',n_nat 
      write(10,'(a)') '# columns: 1'
      do i=1,n_nat
         write(10,'(f16.9)') n_zz(i)
      enddo

      close(10)

    end subroutine write_octave

    subroutine write_octave_sym(strname)

     use nstruct     
     use struct     

     character*50 strname
      open(unit=10,file='struct.read.sym',status='unknown')

      write(10,'(a)') '# Created by rwoctave'
      write(10,'(2a)') '# name: ', trim(strname)
      write(10,'(a)') '# type: struct'
      write(10,'(a)') '# length: 3'
      write(10,'(a)') '# name: nsym'
      write(10,'(a)') '# type: scalar'
      write(10,'(i6)') nsym     
      write(10,'(a)') '# name: rot'
      write(10,'(a)') '# type: matrix'
      write(10,'(a)') '# ndims: 3'
      write(10,'(a,i5)') ' 3 3 ',nsym
      do i=1,nsym
      write(10,'(i5)') transpose(iz(1:3,1:3,i))
      enddo
      write(10,'(a)') '# name: tr'
      write(10,'(a)') '# type: matrix'
      write(10,'(a)') '# ndims: 2'
      write(10,'(a,i5)') ' 3 ',nsym
      write(10,'(f10.6)') tau(1:3,1:nsym)

      close(10)

    end subroutine write_octave_sym

    subroutine jumpto(var)

       character var*(*)
       integer iunit

       character h*1,line*100,var1*50 

       iunit=10 
       rewind(iunit)

1      continue

       read(iunit,'(a)',end=2) line 
!write(*,*) line
       if (line(1:7).eq.'# name:') then
          var1=adjustl(line(8:100))
!          write(*,'(a,a)') 'var1=',var1
!          write(*,'(a,a)') ' var=',var
          if (trim(var1).eq.trim(var)) then
3            continue
             read(iunit,'(a)') line
!             write(*,'(a)') line(1:7) 
             if (line(1:7).eq.'# type:') then
                var1=adjustl(line(8:100))
!                write(*,'(a,a)') 'var1=',var1 
                if (trim(var1).eq.'cell') then
                   read(iunit,*)
                   read(iunit,*)
                   read(iunit,*)
                   goto 3
                elseif (trim(var1).eq.'matrix') then
                   read(iunit,*)
                   read(iunit,*)
                   goto 2
                elseif (trim(var1).eq.'scalar') then
                   goto 2
                elseif (trim(var1).eq.'string') then
                   read(iunit,*)
                   goto 2
                endif
             endif 
          endif
       endif

       goto 1

2     continue

     end subroutine jumpto


     subroutine makemult

       use nstruct

       character*10, allocatable:: aname(:)
       integer, allocatable:: mult(:),jrj(:)  
       real*8, allocatable :: r0(:),rmt(:),pos(:,:),zz(:)

       allocate(aname(n_ndif),mult(n_ndif),jrj(n_ndif))
       allocate(r0(n_ndif),rmt(n_ndif),pos(3,n_ndif),zz(n_ndif))

       nmult=0 
       mult=0
       jj=0

       n_ndif=n_nat

       do i=1,n_ndif
          do j=1,nmult
             if (trim(aname(j)).eq.trim(n_aname(i))) then
                jj=jj+1
                mult(j)=mult(j)+1
                pos(1:3,jj)=n_pos(1:3,i)
                goto 1
             endif
          enddo
          nmult=nmult+1
          mult(nmult)=mult(nmult)+1
          aname(nmult)=n_aname(i)
          r0(nmult)=n_r0(i)
          rmt(nmult)=n_rmt(i)
          zz(nmult)=n_zz(i)
          jrj(nmult)=n_jrj(i)
          jj=jj+1
          pos(1:3,jj)=n_pos(1:3,i)
1         continue

       end do

       n_nat=nmult       
       n_pos(1:3,1:n_ndif)=pos(1:3,1:n_ndif)
       n_r0(1:n_nat)=r0(1:n_nat)
       n_rmt(1:n_nat)=rmt(1:n_nat)
       n_jrj(1:n_nat)=jrj(1:n_nat)
       n_zz(1:n_nat)=zz(1:n_nat)
       n_mult(1:n_nat)=mult(1:n_nat)
       n_aname(1:n_nat)=aname(1:n_nat)
       

     end subroutine makemult
