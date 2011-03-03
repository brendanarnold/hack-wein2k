
  subroutine test_symetry

    use struct
    use rotations

    integer i
    real*8   symr(3,3),syms(3,3),t(3)
    logical ok
    integer, allocatable:: eqiv(:),eqivall(:,:)
    integer n_nsym,n_nat
    integer, allocatable:: oldid(:),oldiat(:)

    real*8,allocatable       :: n_r0(:),n_dx(:),n_rmt(:),n_zz(:)
    real*8,allocatable       :: n_rotloc(:,:,:)
    real*8,allocatable       :: n_tau(:,:)
    character*10,allocatable :: n_aname(:)
    integer,allocatable      :: n_mult(:),n_jrj(:),n_iatnr(:),n_isplit(:)
    integer,allocatable      :: n_iz(:,:,:),n_inum(:)
    real*8,allocatable       :: n_pos(:,:)
    real*8, allocatable     :: n_fi(:),n_teta(:)
    integer, allocatable    :: n_ncmopt(:)
    integer, allocatable    :: n_tmatom(:)
    logical aone

    allocate(eqiv(ndif),eqivall(ndif,nsym))
    allocate(oldid(ndif),oldiat(ndif))

    allocate(n_aname(ndif),n_mult(ndif),n_jrj(ndif),&
         n_r0(ndif),n_dx(ndif),n_rmt(ndif),&
         n_zz(ndif),n_rotloc(3,3,ndif),n_iatnr(ndif),&
         n_isplit(ndif),n_pos(3,ndif),&
         n_iz(3,3,nsym),n_tau(3,nsym),n_inum(nsym))
    allocate(n_fi(ndif),n_teta(ndif),n_ncmopt(ndif),n_tmatom(ndif))

    write(6,'(//,a)') '----------------------------------'
    write(6,'(a)')    '  Filtering symmetry operations'
    write(6,'(a,/)')  '----------------------------------'

    write(6,'(a,i5)') 'number of symetry operations:',nsym
    write(6,*)
      
      n_nsym=0
      
      do iti=0,1
         
         if (iti.eq.1.and.nmagmod) cycle

         do i=1,nsym
            symr(1:3,1:3)=iz(1:3,1:3,i)
            t(1:3)=tau(1:3,i)
            symr=transpose(symr)
            syms=symr

            call test_operation(syms,symr,t,iti,ok,eqiv,i)
            
            if (ok) then               
               n_nsym=n_nsym+1
               write(6,'(/,a)') '*****************************'
               write(6,'(a,i4)') 'operation accepted:',n_nsym
               write(6,'(3f10.5)') transpose(symr)
               write(6,'(3f10.5)') t
               write(6,'(a,i3)') 'time inversion:',iti
               write(6,'(a)') 'equivalency (atom,symr(atom)):'
               write(6,'(2i5)') (id,eqiv(id),id=1,ndif)
               write(6,'(a,/)') '*****************************'
               eqivall(1:ndif,n_nsym)=eqiv(1:ndif)
               n_iz(1:3,1:3,n_nsym)=transpose(symr)
               n_tau(1:3,n_nsym)=t(1:3)
               tinv(n_nsym)=iti
               n_inum(n_nsym)=n_nsym
               
            endif

         enddo
      enddo

      call find_newmult(eqivall,n_nsym,n_mult,n_nat,oldid,oldiat)

      id=0
      do i=1,n_nat
         n_aname(i)=aname(oldiat(i))
         n_jrj(i)=jrj(oldiat(i))
         n_r0(i)=r0(oldiat(i))
         n_dx(i)=dx(oldiat(i))
         n_rmt(i)=rmt(oldiat(i))
         n_zz(i)=zz(oldiat(i))
         n_rotloc(1:3,1:3,i)=rotloc(1:3,1:3,oldiat(i))
!!! assume non cubic site
         n_iatnr(i)=-i         
         n_isplit(i)=isplit(oldiat(i))
         do j=1,n_mult(i)
            id=id+1
            n_pos(1:3,id)=pos(1:3,oldid(id))
            n_fi(id)=fi(oldid(id))
            n_teta(id)=teta(oldid(id))
            n_ncmopt(id)=ncmopt(oldid(id))
            do im=1,nmag
               if (oldid(id).eq.tmatom(im)) n_tmatom(im)=id
            enddo
         enddo
      enddo

      deallocate(aname,mult,jrj,r0,dx,rmt,zz,rotloc,iatnr,isplit)

      nat=n_nat
      nsym=n_nsym
      iord=nsym
      
      allocate(aname(nat),mult(nat),jrj(nat),r0(nat),dx(nat),rmt(nat))
      allocate(zz(nat),rotloc(3,3,nat),iatnr(nat),isplit(nat))
      deallocate (iz,tau)
      allocate (iz(3,3,nsym),tau(3,nsym))

      alpha=alpha*pi2a

      aname=n_aname
      mult=n_mult
      jrj=n_jrj
      r0=n_r0
      dx=n_dx
      rmt=n_rmt
      zz=n_zz
      rotloc=n_rotloc
      iatnr=n_iatnr
      isplit=n_isplit
      pos(1:3,1:ndif)=n_pos(1:3,1:ndif)
      tmatom=n_tmatom
      iz(:,:,1:nsym)=n_iz(:,:,1:nsym)
      tau(:,1:nsym)=n_tau(:,1:nsym)
      fi=n_fi
      teta=n_teta
      ncmopt=n_ncmopt

      sgroup='  nie mam pojêcia'

      call write_inncmsym(eqivall)


      deallocate(eqiv)

  end subroutine test_symetry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write_inncmsym(eqivall)

    use struct

    integer eqivall(ndif,nsym)
    integer i,i1,i2

    lmax=6
    nop=nsym*ndif
    write(8,'(i6,10x,a)') lmax,'# lmax'
    write(8,'(i6,10x,a)') nop,'# navat'
    do i=1,nsym
       do j=1,ndif
          write(8,'(2i5)') j,eqivall(j,i)
          write(8,'(3i5)') ((iz(i1,i2,i),i1=1,3),i2=1,3)
          write(8,'(3f10.5)') (tau(i1,i),i1=1,3)
          write(8,'(i5)') tinv(i) 
       enddo
    enddo
    
  end subroutine write_inncmsym

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine find_newmult(eqivall,n_nsym,n_mult,n_nat,oldid,oldiat)

    use struct
    
    integer eqivall(ndif,nsym),n_mult(ndif)
    integer oldid(ndif),oldiat(ndif)
    integer n_nsym,n_nat
    integer eqiv(ndif,ndif),indat(ndif)
    integer i,j,id,is,im,ia
    logical id_c_is_eqiv, id_is_eqiv
    integer ia_eqiv 
    integer idtmp
    logical nextcycle

    id=0
    do i=1,nat
       do j=1,mult(i)
          id=id+1
          indat(id)=i
       enddo
    enddo

    n_nat=1
    n_mult(n_nat)=1
    eqiv(n_nat,1)=1
    do id=1,ndif
       id_is_eqiv=.false.
       do is=1,n_nsym
          id_c=eqivall(id,is)
          id_c_is_eqiv=.false.
          do ia=1,n_nat
             do im=1,n_mult(ia)
                if (id.eq.eqiv(ia,im)) then
                   id_is_eqiv=.true.
                   ia_eqiv=ia
                endif
                if (id_c.eq.eqiv(ia,im)) then
                   id_c_is_eqiv=.true.
                   ia_eqiv=ia
                endif
             enddo
          enddo
          if (id_is_eqiv.and.(.not.id_c_is_eqiv)) then
             n_mult(ia_eqiv)=n_mult(ia_eqiv)+1
             eqiv(ia_eqiv,n_mult(ia_eqiv))=id_c
             id_c_is_eqiv=.true.
          endif
          if (id_c_is_eqiv.and.(.not.id_is_eqiv)) then
             n_mult(ia_eqiv)=n_mult(ia_eqiv)+1
             eqiv(ia_eqiv,n_mult(ia_eqiv))=id
             id_is_eqiv=.true.
          endif
       enddo
       if (.not.id_is_eqiv) then
          n_nat=n_nat+1
          n_mult(n_nat)=1
          eqiv(n_nat,1)=id
       endif
    enddo

    do i=1,n_nat
1      continue
       nextcycle=.false.
       do j=1,n_mult(i)-1
          if (eqiv(i,j).gt.eqiv(i,j+1)) then
             idtmp=eqiv(i,j)
             eqiv(i,j)=eqiv(i,j+1)
             eqiv(i,j+1)=idtmp
             nextcycle=.true.
          endif
       enddo
       if (nextcycle) goto 1
    enddo

    write(6,'(//,a,i5)') 'new nat=',n_nat
    id=0
    do i=1,n_nat
       oldiat(i)=indat(eqiv(i,1))
       do j=1,n_mult(i)
          id=id+1
          oldid(id)=eqiv(i,j)
          write(6,'(4(a,i3,5x))') &
               'atom=',i,&
               'index=',id,&
               'old atom=',indat(eqiv(i,j)),&
               'old index=',eqiv(i,j)
       enddo
    enddo

    
  end subroutine find_newmult

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_operation(syms,symr,t,iti,ok,eqiv,i)

    use struct
    use rotations

    real*8   symr(3,3),syms(3,3),t(3)
    logical ok
    integer iat,ind,im,ind1,ind2
    real*8  tol,xyz(3),xyz2(3)
    integer eqiv(ndif)
    real*8, allocatable  :: ct(:,:)
    real*8        car2lat(3,3)
    real*8        mm1(3),mm2(3)
    logical vectors_equal
    real*8 det
    integer iti
    real*8 vectmod,ncmqt(3),symq(3,3)

    call determinant(syms,det)
    ! magnetic moment is pseudovector
    if (det.lt.0) syms=(-1.0d0)*syms
    
    ! time inversion inverts mm
    if (iti.eq.1) syms=syms*(-1.0d0)

    symq=symr
    if (iti.eq.1) symq=symr*(-1.0d0)

    if (ideb.gt.0) then
       write(6,'(/,a,i8)') 'testing operation (real):',i
       write(6,100) transpose(symr)
       write(6,101) t
       write(6,'(a,i3)') 'time inversion:',iti       
       write(6,'(/,a)') 'testing operation (spin):'
       write(6,100) transpose(syms)
       write(6,'(a,f10.5)') 'determinant   :',det
100    format('rotation    ',3f10.5)
101    format('translation ',3f10.5)
    endif

    call inversa(br1_dir,car2lat)

    eqiv(1:ndif)=0
    tol=5.0d-4
    ok=.true.

    if (lattic(1:1).eq.'P') then
       nct=1
       allocate(ct(1:3,nct))
       ct(1:3,1)=(/0.0d0,0.0d0,0.0d0/)
    else if (lattic(1:1).eq.'F') then
       nct=4
       allocate(ct(1:3,nct))
       ct(1:3,1)=(/0.0d0,0.0d0,0.0d0/)
       ct(1:3,2)=(/0.5d0,0.5d0,0.0d0/)
       ct(1:3,3)=(/0.5d0,0.0d0,0.5d0/)
       ct(1:3,4)=(/0.0d0,0.5d0,0.5d0/)
    else if (lattic(1:1).eq.'B') then
       nct=2
       allocate(ct(1:3,nct))
       ct(1:3,1)=(/0.0d0,0.0d0,0.0d0/)
       ct(1:3,2)=(/0.5d0,0.5d0,0.5d0/)
    else if (lattic(1:1).eq.'H') then
       nct=1
       allocate(ct(1:3,nct))
       ct(1:3,1)=(/0.0d0,0.0d0,0.0d0/)
    else if (lattic(1:1).eq.'R') then
       nct=1
       allocate(ct(1:3,nct))
       ct(1:3,1)=(/0.0d0,0.0d0,0.0d0/)
    else if (lattic(1:3).eq.'CXY') then
       nct=2
       allocate(ct(1:3,nct))
       ct(1:3,1)=(/0.0d0,0.0d0,0.0d0/)
       ct(1:3,2)=(/0.5d0,0.5d0,0.0d0/)
    else if (lattic(1:3).eq.'CYZ') then
       nct=2
       allocate(ct(1:3,nct))
       ct(1:3,1)=(/0.0d0,0.0d0,0.0d0/)
       ct(1:3,2)=(/0.0d0,0.5d0,0.5d0/)
    else if (lattic(1:3).eq.'CXZ') then
       nct=2
       allocate(ct(1:3,nct))
       ct(1:3,1)=(/0.0d0,0.0d0,0.0d0/)
       ct(1:3,2)=(/0.5d0,0.0d0,0.5d0/)
    endif

    ind1=0
    ind0=0
    do iat=1,nat
       if (iat.gt.1) ind0=ind0+mult(iat-1)
       do im1=1,mult(iat)
          ind1=ind1+1
          xyz(1:3)=pos(1:3,ind1)           
          xyz=matmul(symr,xyz)
          xyz=xyz+t+10.0d0
          xyz(1:3)=mod(xyz(1:3),1.0d0)
          if (ideb.gt.1) write(6,'(a,3f10.6)') 'xyx=',xyz
          ind2=ind0
          do im2=1,mult(iat)
             ind2=ind2+1
             do ict=1,nct
                xyz2(1:3)=(pos(1:3,ind2)+ct(1:3,ict))
                if (ideb.gt.1) write(6,'(a,3f10.6)') 'xyx2=',xyz2
                if (vectors_equal(xyz,xyz2,tol)) then
                   eqiv(ind1)=ind2
                   goto 1
                endif
             enddo
          enddo
          ok=.false.
          return
!!$          write(6,*)
!!$          write(6,*)
!!$          write(6,'(a)') 'cos nie tak z symr w test_operation:'
!!$          write(6,*)
!!$          write(6,'(a)') 'symr:'
!!$          write(6,'(3f10.5)') transpose(symr)
!!$          write(6,'(a)') 'translaction:'
!!$          write(6,'(3f10.5)') t
!!$          write(6,'(a)') 'xyz:'
!!$          write(6,'(3f10.5)') xyz          
!!$          stop
1         continue

          if (.not.matom(ind1)) cycle
          if (.not.matom(ind2)) cycle

          write(6,'(a,2i5)') 'atoms (ind1,ind2):',ind1,ind2
          mm1(1)=sin(teta(ind1))*cos(fi(ind1))
          mm1(2)=sin(teta(ind1))*sin(fi(ind1))
          mm1(3)=cos(teta(ind1))
          if (ideb.gt.0) then
             write(6,'(a)') 'polar angles'
             write(6,'(i5,3f10.5)') ind1,fi(ind1)*pi2a,teta(ind1)*pi2a
             write(6,'(i5,3f10.5)') ind2,fi(ind2)*pi2a,teta(ind2)*pi2a
             write(6,'(a)') 'cartesian cell coordinate'
             write(6,'(i5,3f10.5)') ind1,mm1
          endif   
          mm2(1)=sin(teta(ind2))*cos(fi(ind2))
          mm2(2)=sin(teta(ind2))*sin(fi(ind2))
          mm2(3)=cos(teta(ind2))
          if (ideb.gt.0) then
             write(6,'(i5,3f10.5)') ind2,mm2
          endif   
          mm1=matmul(car2lat,mm1)
          mm2=matmul(car2lat,mm2)
          if (ideb.gt.0) then
             write(6,'(a)') 'unit cell coordinate'
             write(6,'(i5,3f10.5)') ind1,mm1
             write(6,'(i5,3f10.5)') ind2,mm2
          endif   
          mm1=matmul(syms,mm1)
          if (ideb.gt.0) then
             write(6,'(a)') 'afrer symetry:'
             write(6,'(i5,3f10.5)') ind1,mm1
             write(6,'(i5,3f10.5)') ind2,mm2
             write(6,*)vectors_equal(mm1,mm2,tol)
          endif   
          if (.not.vectors_equal(mm1,mm2,tol)) then
             ok=.false.
             return
          endif

          ncmqt(1:3)=matmul(symq,ncmq)
          if (ideb.gt.0) then
             write(6,'(a,3f8.4)') 'input q=',ncmq
             write(6,'(a,3f8.4)') 'symet q=',ncmqt
             write(6,*) vectors_equal(ncmqt,ncmq,tol)
          endif   
          if (.not.vectors_equal(ncmqt,ncmq,tol)) then
             ok=.false.
             return
          endif

       enddo
    enddo
    
  deallocate(ct)

  end subroutine test_operation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  logical function vectors_equal(xyz,xyz2,tol)

    real*8 xyz(3), xyz2(3),tol
    real*8 dxyz(3)
    integer k

    dxyz(1:3)=abs(xyz(1:3)-xyz2(1:3))
    do k=1,3
       if (abs(mod(dxyz(k),1.0d0)).lt.tol) dxyz(k)=0.0d0
!       if (abs(dxyz(k)-1.0d0).lt.tol) dxyz(k)=0.0d0
    enddo
    vectors_equal=dxyz(1).lt.tol.and.dxyz(2).lt.tol.and.dxyz(3).lt.tol

  end function vectors_equal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function vectmod(v)

    real*8 v(3)
    
    vectmod=sqrt(v(1)**2+v(2)**2+v(3)**3)

  end function vectmod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function aone(a)

    real*8   a(3,3)
    real*8   b(3,3)
    real*8   one(3,3)
    real*8  zero
    logical ok
    integer i,j,n


    zero=1.0d-8

    one(1,1:3)= (/1.0d0,0.0d0,0.d0/)
    one(2,1:3)= (/0.0d0,1.0d0,0.d0/)
    one(3,1:3)= (/0.0d0,0.0d0,1.d0/)
    
    b=a-one
    
    ok=.false.
    do i=1,3
       do j=1,3
          if (abs(b(i,j)).gt.zero) ok=.true.
       enddo
    enddo

    aone=.not.ok
    
  end function aone
