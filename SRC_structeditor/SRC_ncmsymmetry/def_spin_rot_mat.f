 subroutine def_spin_rot_mat

  use Ylm_rot
  use struct, only: br1_dir,br1_rec,nat
  use rotations
  use defs

  integer  i,ipgo
  real*8   sym(3,3),det,a,b,c
  complex*16   cimag
  real*8   sprot1(3,3),sprot1_inv(3,3),sprot2(3,3),symrot(3,3),sprots(3,3)
  real*8   sprot(3,3)
  real*8  br1_dir_inv(3,3)
  real*8  br1_rec_inv(3,3)
  complex*16 sigmay(2,2),spin_rot(2,2)
   
  cimag=(0.0d0,1.0d0)

  if (ideb.gt.0) then
     write(6,*)
     write(6,*)
     write(6,'(a)') '**** def spin rot mat ****'
     write(6,*)
  endif

  iat=1
  do i=1,nat
     
     if (i.gt.1) iat=iat+mult(i-1)

     do ipgo=1,npgopat(i)

        sym(1:3,1:3)=pg_sym_oper(i,ipgo,1:3,1:3)
        call determinant(sym,det)

        if (ideb.gt.0) then
           write(6,'(a)') 'sym :'
           write(6,'(3f15.8)') transpose(sym)
           if (pgtinv(i,ipgo).gt.0.5d0) write(6,'(a)') 'timeinversion present' 
           write(6,*)      
        endif

        if (det.lt.0d0) then
           sym(1:3,1:3)=-sym(1:3,1:3)
           if (ideb.gt.0) then
              write(6,'(a)') 'sym is  reflection, so multiply by inversion:'
              write(6,'(3f15.8)') transpose(sym)
              write(6,*)      
           endif
        endif
        
        a=fi(iat)!/pi2a 
        b=teta(iat)!/pi2a
        sprot1(1,1)=cos(a)*cos(b)
        sprot1(1,2)=-sin(a)
        sprot1(1,3)=cos(a)*sin(b)
        sprot1(2,1)=sin(a)*cos(b)
        sprot1(2,2)=cos(a)
        sprot1(2,3)=sin(a)*sin(b)
        sprot1(3,1)=-sin(b)
        sprot1(3,2)=0.0d0
        sprot1(3,3)=cos(b)
        
        call euler(0,sprot1,a,b,c)

        if (ideb.gt.0) then
           write(6,'(a,3f15.8)') 'sprot1 eulers:', a*pi2a,b*pi2a,c*pi2a
           write(6,'(a)') 'sprot1:'
           write(6,'(3f15.8)') transpose(sprot1)
           write(6,*)
        endif
        spin_rot(1,1)=exp(-cimag*(a+c)*0.5d0)*cos(b*0.5d0)
        spin_rot(2,1)=exp(cimag*(a-c)*0.5d0)*sin(b*0.5d0)
        spin_rot(1,2)=-exp(-cimag*(a-c)*0.5d0)*sin(b*0.5d0)
        spin_rot(2,2)=exp(cimag*(a+c)*0.5d0)*cos(b*0.5d0)
        spin_rot_mat(i,ipgo,1,1:2,1:2)=spin_rot(1:2,1:2)
        
        
        a=fi(iat)!/pi2a 
        b=teta(iat)!/pi2a
        sprot2(1,1)=cos(a)*cos(b)
        sprot2(1,2)=-sin(a)
        sprot2(1,3)=cos(a)*sin(b)
        sprot2(2,1)=sin(a)*cos(b)
        sprot2(2,2)=cos(a)
        sprot2(2,3)=sin(a)*sin(b)
        sprot2(3,1)=-sin(b)
        sprot2(3,2)=0.0d0
        sprot2(3,3)=cos(b)
        
        call euler(0,sprot2,a,b,c)
        if (ideb.gt.0) then
           write(6,'(a,3f15.8)') 'sprot2 eulers:', a*pi2a,b*pi2a,c*pi2a
           write(6,'(a)') 'sprot2:'  
           write(6,'(3f15.8)') transpose(sprot2)
           write(6,*)
        endif
        spin_rot(1,1)=exp(-cimag*(a+c)*0.5d0)*cos(b*0.5d0)
        spin_rot(2,1)=exp(cimag*(a-c)*0.5d0)*sin(b*0.5d0)
        spin_rot(1,2)=-exp(-cimag*(a-c)*0.5d0)*sin(b*0.5d0)
        spin_rot(2,2)=exp(cimag*(a+c)*0.5d0)*cos(b*0.5d0)
        spin_rot_mat(i,ipgo,2,1:2,1:2)=spin_rot(1:2,1:2)
        
        
        call inversa(br1_dir,br1_dir_inv)
        call inversa(sprot1,sprot1_inv)
        
        sym=matmul(matmul(br1_dir,sym),br1_dir_inv)
        call euler(0,sym,a,b,c)
        if (ideb.gt.0) then
           write(6,'(a,3f15.8)') 'sym in cartesian eulers:', a*pi2a,b*pi2a,c*pi2a
           write(6,'(a)') 'sym :'  
           write(6,'(3f15.8)') transpose(sym)
           write(6,*)
           write(6,*)
        endif
        spin_rot(1,1)=exp(-cimag*(a+c)*0.5d0)*cos(b*0.5d0)
        spin_rot(2,1)=exp(cimag*(a-c)*0.5d0)*sin(b*0.5d0)
        spin_rot(1,2)=-exp(-cimag*(a-c)*0.5d0)*sin(b*0.5d0)
        spin_rot(2,2)=exp(cimag*(a+c)*0.5d0)*cos(b*0.5d0)
        spin_rot_mat(i,ipgo,3,1:2,1:2)=spin_rot(1:2,1:2)
        
     enddo
   
  enddo

end subroutine def_spin_rot_mat
   
 subroutine apply_rotsym(sdm,i,ipgo)

   use Ylm_rot

   integer i,ipgo
   complex*16 sdm(2,2)
   complex*16 rot(2,2)
   complex*16 roti(2,2)

   rot(1:2,1:2)=spin_rot_mat(i,ipgo,1,1:2,1:2)
   roti=transpose(conjg(rot))
   sdm=matmul(matmul(rot,sdm),roti)

   rot(1:2,1:2)=spin_rot_mat(i,ipgo,3,1:2,1:2)
   roti=transpose(conjg(rot))
   sdm=matmul(matmul(rot,sdm),roti)

   rot(1:2,1:2)=spin_rot_mat(i,ipgo,2,1:2,1:2)
   roti=transpose(conjg(rot))
   sdm=matmul(matmul(roti,sdm),rot)

   if (pgtinv(i,ipgo).gt.0.5d0) call apply_timeinv(sdm)
   
 end subroutine apply_rotsym
 

 subroutine apply_rotsym_rec(sdm,i,ipgo)

   use Ylm_rot

   integer i,ipgo
   complex*16 sdm(2,2)
   complex*16 rot(2,2)
   complex*16 roti(2,2)

   rot(1:2,1:2)=spin_rot_mat(i,ipgo,3,1:2,1:2)
   roti=transpose(conjg(rot))
   sdm=matmul(matmul(rot,sdm),roti)

   if (pgtinv(i,ipgo).gt.0.5d0) call apply_timeinv(sdm)
   
 end subroutine apply_rotsym_rec

 subroutine apply_timeinv(sdm)
   
   complex*16 sdm(2,2)
   complex*16 sigmay(2,2),sdm_t(2,2)
   

   sdm=conjg(sdm)
   sigmay(1,1)=(0.0d0,0.0d0)
   sigmay(1,2)=(0.0d0,-1.0d0)
   sigmay(2,1)=(0.0d0,1.0d0)
   sigmay(2,2)=(0.0d0,0.0d0)
   sdm=matmul(matmul(sigmay,sdm),sigmay)

 end subroutine apply_timeinv
