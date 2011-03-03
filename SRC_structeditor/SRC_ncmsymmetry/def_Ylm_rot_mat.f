 subroutine def_Ylm_rot_mat

  use Ylm_rot
  use rotations, only: lmax,ideb
  use struct, only: rotloc, br1_dir,nat

  integer  mod,i,l,j,jj,ipgo
  real*8   sym(3,3),det
  real*8   a,b,c
  real*8   a1,b1,c1
  real*8   a2,b2,c2
  real*8   a3,b3,c3
  complex*16 mat(-lmax:lmax,-lmax:lmax)
  complex*16 mat1(-lmax:lmax,-lmax:lmax)
  complex*16 mat2(-lmax:lmax,-lmax:lmax)
  complex*16 mat3(-lmax:lmax,-lmax:lmax)
  complex*16 matti(-lmax:lmax,-lmax:lmax)
  real*8  rot1(3,3),rot2(3,3),rot1_inv(3,3)
  real*8  br1_dir_inv(3,3)

  call inversa(br1_dir,br1_dir_inv)

  if (ideb.gt.0) then
     write(6,'(a)') '**** def Ylm rot mat ****'
     write(6,*)
     write(6,'(a)') 'br1_dir:'
     write(6,'(3f15.8)') transpose(br1_dir)
     write(6,*)         
     write(6,'(a)') 'br1_dir_inv:'
     write(6,'(3f15.8)') transpose(br1_dir_inv)
     write(6,*)         
  endif

  do i=1,nat
 
     if (ideb.gt.0) then    
        write(6,*) '----------------------------------------'
        write(6,*) '     atom=',i
        write(6,*) '----------------------------------------'
     endif

     rot1(1:3,1:3)=rotloc(1:3,1:3,i)
     rot2(1:3,1:3)=rotloc(1:3,1:3,i)
!rotloc read as transposed, see module.f
     rot1=transpose(rot1)
     rot2=transpose(rot2)
! now we have local->globar coordinates transf.

     call euler(1,rot1,a1,b1,c1)
     call euler(1,transpose(rot2),a3,b3,c3)

     if (ideb.gt.0) then
        write(6,'(a)') 'rot1:'
        write(6,'(3f15.8)') transpose(rot1)
        write(6,'(a,3f15.8)') 'rot1 eulers:', a1,b1,c1
        write(6,*)         
        write(6,'(a)') 'rot2:'
        write(6,'(3f15.8)') transpose(rot2)
        write(6,'(a,3f15.8)') 'rot2 eulers:', a3,b3,c3
        write(6,*)         
     endif

     do ipgo=1,npgopat(i)

        sym(1:3,1:3)=pg_sym_oper(i,ipgo,1:3,1:3)
        
        if (ideb.gt.0) then
           write(6,'(a)') 'sym :'
           write(6,'(3f15.8)') transpose(sym)
           if (pgtinv(i,ipgo).gt.0.5d0) write(6,'(a)') 'timeinversion present' 
           write(6,*)      
        endif

        call determinant(sym,det)

        if (det.lt.0d0) then
           sym(1:3,1:3)=-sym(1:3,1:3)
!!$           write(6,'(a)') 'sym is  reflection, so multiply by inversion:'
!!$           write(6,'(3f15.8)') transpose(sym)
!!$           write(6,*)      
        endif
        
        sym=matmul(matmul(br1_dir,sym),br1_dir_inv)
        
        call euler(1,sym,a2,b2,c2)

        if (ideb.gt.0) then
           write(6,'(a,3f15.8)') 'sym eulers:', a2,b2,c2
           write(6,*)           
        endif

        do l=0,lmax

           call find_rot_mat(l,a1,b1,c1,mat1,lmax)
           call find_rot_mat(l,a2,b2,c2,mat2,lmax)
           if (det.lt.0) call apply_inversion_Ylm(l,lmax,mat2)
           call find_rot_mat(l,a3,b3,c3,mat3,lmax)
           
           mat=matmul(mat2,mat1)
           mat=matmul(mat3,mat)        
           Ylm_rot_mat(i,ipgo,l,-l:l,-l:l)=mat(-l:l,-l:l)

           if (ideb.gt.0) then
              write(6,'(a,i3)') 'Ylm rotation matrix for l=',l
              do j=-l,l
                 write(6,'(100(2f8.4,3x))') (mat(j,jj),jj=-l,l)
              enddo
              write(6,*)
           endif

           if (pgtinv(i,ipgo).gt.0.5d0) then
              call def_timeinv_Ylm(l,lmax,matti)
              mat3=matmul(matti,mat3)
           endif

           mat=matmul(mat2,mat1)
           mat=matmul(mat3,mat)        

           Ylm_rot_mat_dmat(i,ipgo,l,-l:l,-l:l)=mat(-l:l,-l:l)
!!$        write(6,'(a,i3)') 'Ylm for dmat l=',l
!!$        do j=-l,l
!!$           write(6,'(100(2f8.4,3x))') (mat(j,jj),jj=-l,l)
!!$        enddo
!!$        write(6,*)

        enddo
     enddo
  enddo

 end subroutine def_Ylm_rot_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

 subroutine find_rot_mat(l,a,b,c,mat,lmax)

  integer l,lmax
  real*8  a,b,c
  complex*16 mat(-lmax:lmax,-lmax:lmax)
  integer n,m
  complex*16 imag
  real*8  dlmat 

  imag=(0.0d0,1.0d0)
  mat=(0.0d0,0.0d0)
  do n=-l,l
     do m=-l,l
        call find_dlmat(l,n,m,b,dlmat) 
        mat(n,m)=exp(-imag*n*c)*dlmat*exp(-imag*m*a)
     enddo
  enddo

 end subroutine find_rot_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine find_dlmat(l,n,m,b,dlmat)

  integer l,n,m
  real*8  dlmat,b,fact,b2
  integer k

  b2=b/2.0d0
  dlmat=0.0d0
  do k=0,2*l
 
     if (((l-n-k).ge.0).and.((l+m-k).ge.0).and.((k-m+n).ge.0).and.&
         ((l+n).ge.0).and.((l+m).ge.0).and.((l-n).ge.0).and.((l-m).ge.0)) then

         dlmat=dlmat+((-1.0d0)**(k-m+n))*&
             (cos(b2)**(2*l+m-n-2*k))*(sin(b2)**(2*k+n-m))*&
             sqrt(fact(l+n)*fact(l+m)*fact(l-n)*fact(l-m))/&
             (fact(l-n-k)*fact(l+m-k)*fact(k)*fact(k-m+n))     

     endif

  enddo

 end subroutine find_dlmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 real*8 function fact(n)
   
   integer n,j

   if (n.eq.0) then
      fact=1
   else
      fact=1
      do j=1,n
         fact=fact*j
      end do
   end if

 end Function fact
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine def_timeinv_Ylm(l,lmax,mat)

   integer l,i,j,jj
   complex*16 mat(-lmax:lmax,-lmax:lmax)
   complex*16 tmp(-l:l,-l:l)
   complex*16      ti(-l:l,-l:l)

   mat=(0.0d0,0.0d0)
    do i=-l,l
       do j=-l,l
          ti(i,j)=0.0d0
          if (i.eq.-j) ti(i,j)=(-1.0d0)**i
       enddo
    enddo

!!$        write(6,'(a,i3)') 'time inversion ti for l=',l
!!$        do j=-l,l
!!$           write(6,'(100(2f8.4,3x))') (ti(j,jj),jj=-l,l)
!!$        enddo
!!$        write(6,*)

    mat(-l:l,-l:l)=ti(-l:l,-l:l)

  end subroutine def_timeinv_Ylm

 subroutine apply_inversion_Ylm(l,lmax,mat)

   integer l
   complex*16 mat(-lmax:lmax,-lmax:lmax)

   mat=((-1.0d0)**l)*mat
   
 end subroutine apply_inversion_Ylm
  
