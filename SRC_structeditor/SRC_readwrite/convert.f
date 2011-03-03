 subroutine convert

  use struct
  use nstruct     

    integer        icp,iri
    character*50   text
    integer        i,j,ii,kk,l,k
    character      tnum*10
    real*8         tr(3,4),tra(3,3),trainv(3,3)
    integer, allocatable :: iio(:,:)
    real*8 pi
    real*8   tmpv(3)

      pi=acos(-1.0d0) 

      call latgen_struct

      n_title=title
      n_lattic=lattic
      n_irel=irel
      n_aa=aa
      n_bb=bb
      n_cc=cc
      n_alpha(1:3)=alpha(1:3)      
      n_nat=ndif
      ntr=1
      tr(1:3,1)=(/0.0d0,0.0d0,0.0d0/)

      call init_nstruct(n_nat)
      ii=0
      kk=0
      allocate(iio(nat,maxval(mult)))
      do i=1,nat
         do j=1,mult(i)
            ii=ii+1
            iio(i,j)=ii
         enddo
      enddo
      ii=0
      do i=1,nat
         do j=1,mult(i)
            ii=ii+1
            n_pos(1:3,ii)=pos(1:3,iio(i,j))
            do l=1,3
               if (n_pos(l,ii).gt.1.0d0) n_pos(l,ii)=n_pos(l,ii)-1.0d0
               if (n_pos(l,ii).lt.0.0d0) n_pos(l,ii)=n_pos(l,ii)+1.0d0
            enddo
            call num2text(ii,tnum) 
            n_aname(ii)=trim(aname(i))
            n_jrj(ii)=jrj(i)
            n_iatnr(ii)=-ii
            n_isplit(ii)=isplit(i)
            n_r0(ii)=r0(i)
            n_dx(ii)=dx(i)
            n_rmt(ii)=rmt(i)
            n_zz(ii)=zz(i)             
            
            n_rotloc(1:3,1:3,ii)=reshape((/1.0d0,0.0d0,0.0d0,&
                 0.0d0,1.0d0,0.0d0,&
                 0.0d0,0.0d0,1.0d0/),&
                 (/3,3/))
         enddo
      enddo

      n_nsym=1
      n_iord=1
      n_mult(1:n_nat)=1
      n_iz(1:3,1:3,1)=reshape((/1.0d0,0.0d0,0.0d0,&
                                0.0d0,1.0d0,0.0d0,&
                                0.0d0,0.0d0,1.0d0/),&
                                 (/3,3/))
      n_inum(1)=1
      n_tau(1:3,1)=0.0d0

 end subroutine convert

 subroutine num2text(ii,tnum)

   integer ii
   character tnum*(*),tnum1*100
   
   write(tnum1,*) ii
   tnum=adjustl(tnum1)
   
 end subroutine num2text
