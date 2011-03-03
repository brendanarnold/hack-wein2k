  subroutine fix_nonmagatom_1

    use struct
    use Ylm_rot, only: npgop,pg_sym_oper,pgtinv,npgopat
    use rotations

    complex*16   a(3,3),w(3),vl(3,3),vr(3,3),work(9)
    real*8       rwork(6)
    integer  ipiv(3),info
    logical  anotone
    integer  ione
    complex*16    vect(npgop,3,3),cvect(3,3),detc
    integer       nvect(npgop),ncvect
    real*8    small,tol,dir(3)
    integer i,j,iat,iop,ii
    real*8  br(3,3)
    integer, allocatable ::  itranop(:)
    integer ntranop
    real*8   sym(3,3),syms(3,3),qm(3),qmold(3),det
    logical vectors_equal
    logical  :: warn
    logical magi,magi1,magia
    real*8  :: car2lat(3,3),p2a

    write(6,'(//,a)')'-----------------------------'    
    write(6,'(a)')'  Fixing nonmagnetic atoms (some of equiv magnetic)'
    write(6,'(a,/)')'-----------------------------'    

    allocate (itranop(nsym))

    p2a=45.0d0/atan(1.0d0)

    br=br1_dir

    call inversa(br,car2lat)

    small=1.0d-8
    tol=5.d-4
    warn=.false.

    iat=0
    do i=1,nat
       magi1=.false.
       magia=.true.
       do j=1,mult(i)
          iat=iat+1
          magi1=magi1.or.matom(iat)
          magia=magia.and.matom(iat)
       enddo
       magi=magi1.and.(.not.magia)
       if (magi) then
          iat=iat-mult(i)
          do j=1,mult(i)
             iat=iat+1
             if (matom(iat)) iatm=iat
          enddo
          iat=iat-mult(i)
          do j=1,mult(i)
             iat=iat+1
             if (.not.matom(iat)) then
                ind1=iatm
                ind2=iat
                dir(1)=sin(teta(ind1))*cos(fi(ind1))
                dir(2)=sin(teta(ind1))*sin(fi(ind1))
                dir(3)=cos(teta(ind1))
                dir=matmul(car2lat,dir)                
                call find_tranop(ind1,ind2,itranop,ntranop)                
                write(6,'(2(a,i4),3x,a,i3)') 'magnetic atom',ind1,' nonmagnetic atom', ind2,&
                     '  number of trnsfor.',ntranop
                do jj=1,ntranop                   
                   sym(1:3,1:3)=iz(1:3,1:3,itranop(jj))
                   syms=transpose(sym)
                   call determinant(syms,det)
                   ! magnetic moment is pseudovector
                   if (det.lt.0) syms=syms*(-1.0d0)             
                   ! time inversion inverts mm
                   if (tinv(itranop(jj)).gt.0.5d0) syms=syms*(-1.0d0)                   
                   qm=matmul(syms,dir)                   
!                   write(6,'(a,3f10.5)')    '             magnetic moment:',qm                   
                   if (jj.eq.1) then
                      qmold=qm
                   else
                      if (.not.vectors_equal(qm,qmold,tol)) then
                       !  write(6,'(a)') 'warning:  vectors not equal'
                         warn=.true.
                      endif
                   endif                   
                   qm=matmul(br,qm)
                   call  fiteta(qm,fi(ind2),teta(ind2))
                   write(6,'(a,i4,2f10.5)') 'fi, teta of nonmagnetic atom:',ind2,fi(ind2)*p2a,teta(ind2)*p2a
                enddo                
                matom(ind2)=.true.
             endif
          enddo
       endif
    enddo

    deallocate (itranop)

    if (warn) stop 'check output, more then one possible choices for some mm'

  end subroutine fix_nonmagatom_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fix_nonmagatom_2

    use struct
    use Ylm_rot, only: npgop,pg_sym_oper,pgtinv,npgopat
    use rotations

    complex*16   a(3,3),w(3),vl(3,3),vr(3,3),work(9)
    real*8       rwork(6)
    integer  ipiv(3),info
    logical  anotone
    integer  ione
    complex*16    vect(npgop,3,3),detc
    real*8 cvect(3,3)
    integer       nvect(npgop),ncvect
    real*8    small,tol,dir(3)
    integer i,j,iat,iop,ii
    real*8  br(3,3)
    integer, allocatable ::  itranop(:)
    integer ntranop
    real*8   sym(3,3),syms(3,3),qm(3),qmold(3),det
    logical vectors_equal
    logical  :: warn
    logical magi,magi1,magia
    real*8  :: car2lat(3,3),p2a


    write(6,'(//,a)')'-----------------------------'    
    write(6,'(a)')'  Fixing nonmagnetic atoms (all equiv nonmagnetic)'
    write(6,'(a,/)')'-----------------------------'    

    allocate (itranop(nsym))


    p2a=45.0d0/atan(1.0d0)

    br=br1_dir
    call inversa(br,car2lat)

    small=1.0d-8
    tol=5.d-4
    warn=.false.

    iat=1
    do i=1,nat

       if (i.gt.1) iat=iat+mult(i-1)

       if (.not.matom(iat)) then
          
          iop=0
          do j=1,npgopat(i)
             a(1:3,1:3)=pg_sym_oper(i,j,1:3,1:3)

             call determinant_c(a,detc)
             det=dble(detc)
             if (det.lt.0.0d0) a=-a
             if (pgtinv(i,j).gt.0.5d0) a=-a

             !if (anotone(a,npgopat(i))) then
 
                write(6,'(a)') 'spin op mat'
                write(6,'(3(2f10.5,2x))') a
                write(6,'(a,f6.3)') 'det=',det
                if (pgtinv(i,j).gt.0.5d0)  write(6,'(a)') 'time inv present'

                call ZGEEV( 'N', 'V', 3, A, 3, W, VL, 3, VR, 3,&
                     WORK, 9, RWORK, INFO )
          
                if (info.eq.0) then
                   ione=0
                   iop=iop+1
                   do ii=1,3
                      if (abs(w(ii)-1.0d0).lt.small) then
                         ione=ione+1
                         vect(iop,ione,1:3)=vr(1:3,ii)  
                      endif
                      write(6,'(a,i3,a,i3,a,2f10.5,a,3(2f10.5,3x))') &
                           'atom=',i,&
                           '  pgoper=',j,&
                           '  eigval=',w(ii),&
                           '  eigvec=',vr(1:3,ii)  
                   enddo
                   write(6,*)
                   if (ione.eq.0) then
                      write(6,'(//,a)') 'error in fix_nonmagatom: vect not found, same of the pg'
                      write(6,'(   a)') '              operation does not preserve any direction'
                      write(6,'(   a)') 'Try to fix direction of one mm'
                      write(0,'(//,a)') 'error in fix_nonmagatom: vect not found, same of the pg'
                      write(0,'(   a)') '              operation does not preserve any direction'
                      write(0,'(a,//)') 'Try to fix direction of one mm'
                      stop 'error in fix_nonmagatom: vect not found'
                   endif
                   nvect(iop)=ione
                else
                   write(6,'(//,a,i5)') 'error in fix_nonmagatom',info
                   stop 'error in fix_nonmagatom'
                endif
             !endif
          enddo

          call find_common_vect(vect,nvect,iop,cvect,ncvect)

          write(6,*)
          write(6,'(2(a,i3))') 'atom=',i,'   ncvect=',ncvect
          do ii=1,ncvect
             write(6,'(a,3f10.5)') '       cvect=',cvect(ii,1:3)
          enddo
          write(6,*)

          dir(1:3)=cvect(1,1:3)

          ind1=iat
          do ii=1,mult(i)
             ind2=iat+ii-1

             call find_tranop(ind1,ind2,itranop,ntranop)

             write(6,'(a,2i4,a,i3)') 'equiv atoms:',ind1,ind2,'  number of trnsfor.',ntranop

             do jj=1,ntranop

                sym(1:3,1:3)=iz(1:3,1:3,itranop(jj))
                syms=transpose(sym)
                call determinant(syms,det)
                ! magnetic moment is pseudovector
                if (det.lt.0) syms=syms*(-1.0d0)             
                ! time inversion inverts mm
                if (tinv(itranop(jj)).gt.0.5d0) syms=syms*(-1.0d0)
                
                write(6,'(a)') 'spin op mat'
                write(6,'(3f10.5)') transpose(syms) 
                
                qm=matmul(syms,dir)

                write(6,'(a,3f10.5)') '       magnetic moment (unit cell):',qm

                if (jj.eq.1) then
                   qmold=qm
                else
                   if (.not.vectors_equal(qm,qmold,tol)) then
                      write(6,'(a)') 'warning:  vectors not equal'
                      warn=.true.
                   endif
                endif
                   
             enddo

             qm=matmul(br,qm)
             write(6,'(a,3f10.5)') '       magnetic moment (cartesian):',qm
             call  fiteta(qm,fi(ind2),teta(ind2))
             write(6,'(a,i4,2f10.5)') 'fi, teta of nonmagnetic atom',ind2,fi(ind2)*p2a,teta(ind2)*p2a
          enddo
       endif
    enddo

    deallocate (itranop)
    if (warn) stop 'check output'

  end subroutine fix_nonmagatom_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_tranop(ind1,ind2,itranop,ntranop)

    use struct
    use rotations

    integer ind1,ind2
    integer itranop(*),ntranop
    
    real*8   symr(3,3),t(3)
    logical ok
    real*8  tol,xyz(3),xyz2(3)
    real*8, allocatable  :: ct(:,:)
    real*8        car2lat(3,3)
    logical vectors_equal

    call inversa(br1_dir,car2lat)

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

    tol=5.d-4
    ntranop=0
    do isym=1,nsym

       symr(1:3,1:3)=iz(1:3,1:3,isym)
       t(1:3)=tau(1:3,isym)
       symr=transpose(symr)
       
       xyz(1:3)=pos(1:3,ind1)           
       xyz=matmul(symr,xyz)
       xyz=xyz+t+2.0d0
       xyz(1:3)=mod(xyz(1:3),1.0d0)
       
       ok=.false.

       do ict=1,nct
          xyz2(1:3)=(pos(1:3,ind2)+ct(1:3,ict))
          xyz2=xyz2+2.0d0
          xyz2(1:3)=mod(xyz2(1:3),1.0d0)
          
          if (vectors_equal(xyz,xyz2,tol)) ok=.true.
       enddo
          
       if (ok) then 
          ntranop=ntranop+1         
          itranop(ntranop)=isym
       endif
   
    enddo
    
  deallocate(ct)

end subroutine find_tranop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_common_vect(vect,nvect,iop,cvect,ncvect)

    use Ylm_rot, only: npgop

    complex*16    vect(npgop,3,3)
    real*8   cvect(3,3)
    integer       nvect(npgop),iop,ncvect
    real*8    small
    real*8  v1(3),v2(3)
    logical vecteq,ok1,ok
    integer i,j,ii
    real*8, allocatable ::   a(:,:),s(:),u(:,:),vt(:,:),work(:),r(:),rp(:),x(:),xp(:)
    integer mr3,nrow,lwork,info
    real*8 wn(3),w1(3),w2(3)
    logical warn,jest

    small=1.0d-8
    warn=.false.

    write(6,'(//,a)')'-----------------------------'    
    write(6,'(a)')'  Searching for common vectors'
    write(6,'(a,/)')'-----------------------------'    

    nrow=0

    do i=1,iop
       if (nvect(i).eq.3) then
          write(6,'(a,i4)')     'operation=',i
          write(6,'(a,3f10.5)') 'space: ',real(vect(i,1,1:3)) 
          write(6,'(a,3f10.5)') '       ',real(vect(i,2,1:3)) 
          write(6,'(a,3f10.5)') '       ',real(vect(i,3,1:3))           
       endif
       if (nvect(i).eq.2) then 
          write(6,'(a,i4)')     'operation=',i
          write(6,'(a,3f10.5)') 'plane: ',real(vect(i,1,1:3)) 
          write(6,'(a,3f10.5)') '       ',real(vect(i,2,1:3)) 
          nrow=nrow+1
       endif       
       if (nvect(i).eq.1) then 
          write(6,'(a,i4)')     'operation=',i
          write(6,'(a,3f10.5)') 'vecto: ',real(vect(i,1,1:3)) 
          nrow=nrow+1
       endif
    enddo   

    do i=1,iop
       do j=1,nvect(i)
          do ii=1,3
             if (imag(vect(i,j,ii)).gt.small) then
                write(6,*) 'vect are not real in find_common_vect'
                stop 'vect are not real in find_common_vect'                
             endif
          enddo
       enddo
    enddo

    if (iop.eq.1) then
       ncvect=3
       cvect(1:3,1:3)=vect(1,1:3,1:3)
       return
    endif   

    if (nrow.eq.0) then
       write(6,*) 'nrow.eq.0 in find_common_vect'
       stop 'nrow.eq.0 in find_common_vect'
    endif

    mr3=min(nrow,3)
    lwork=2*MAX(3*MIN(nrow,3)+MAX(nrow,3),5*MIN(nrow,3))
    allocate(a(nrow,3))
    allocate(s(mr3))
    allocate(u(nrow,nrow))
    allocate(vt(3,3))
    allocate(work(lwork))
    allocate(r(nrow))
    allocate(rp(nrow))
    allocate(x(3))
    allocate(xp(3))

    nrow=0
    do i=1,iop
       if (nvect(i).eq.2) then
          w1=vect(i,1,1:3)
          w2=vect(i,2,1:3)
          call vecprod(w1,w2,wn)
          nrow=nrow+1
          a(nrow,1:3)=wn(1:3)
          r(nrow)=0.0d0
       endif
       if (nvect(i).eq.1) then
          nrow=nrow+1
          a(nrow,1:3)=vect(i,1,1:3)
          r(nrow)=1.0d0
       endif
    enddo

    write(6,'(a,i4)') 'nrow=',nrow
    do i=1,nrow
       write(6,'(a,3f10.5)') 'A=',A(i,1:3)
    enddo
    do i=1,nrow
       write(6,'(a,f10.5)') 'R=',r(i)
    enddo

    call DGESVD('A','A',nrow,3,A,nrow,S,U,nrow,VT,3,WORK,LWORK,INFO)

    write(6,'(a)') '---------------------------'
    write(6,'(a,100f10.5:)') 'S=',s
    do i=1,nrow
       write(6,'(a,100f10.5:)') 'U=',u(i,1:nrow)
    enddo
    do i=1,3
       write(6,'(a,3f10.5:)') 'VT=',vt(i,1:3)
    enddo

    rp=matmul(transpose(u),r)
    write(6,'(a,100f10.5:)') 'Rp=',rp(1:nrow)

    do i=1,mr3
       if ((s(i).lt.small).and.(rp(i).gt.small)) then
          write(6,*) 'no solution SVD in find_common_vect'
          stop 'no solution SVD in find_common_vect'
       endif
    enddo
    do i=mr3+1,nrow
       if (rp(i).gt.small) then
          write(6,*) 'no solution SVD in find_common_vect'
          stop 'no solution SVD in find_common_vect'
       endif
    enddo

    xp(1:3)=0.0d0
    do i=1,mr3
       if (s(i).gt.small) xp(i)=rp(i)/s(i)       
    enddo

    ok=.false.
    do i=1,3
       if (abs(xp(i)).gt.small) ok=.true.
    enddo

    if (.not.ok) then
       jest=.false.
       do i=1,mr3
          if (s(i).lt.small) then
             x(1:3)=vt(i,1:3)       
             jest=.true.
          endif   
       enddo
       if (.not.jest) then
          if (mr3+1.gt.3) then 
             write(6,*) 'no solution SVD in find_common_vect'
             stop 'no solution SVD in find_common_vect'
          endif
       endif
       x(1:3)=vt(mr3+1,1:3)
    else
       x(1:3)=matmul(transpose(vt),xp)
    endif

    write(6,'(a,3f10.5)') 'solution=',x

    ncvect=1
    cvect(1,1:3)=x(1:3)

    write(6,'(//,a)')'------------------------------------'    
    write(6,'(a)')'  Congratulations, found something!!'
    write(6,'(a,/)') '------------------------------------'    


  end subroutine find_common_vect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function vecteq(v1,v2)

    complex*16  v1(3),v2(3)
    real*8    small
    logical ok
    integer i

    small=1.0d-8

    ok=.true.
    do i=1,3
       if(abs(v1(i)-v2(i)).gt.small) ok=.false.
    enddo
    vecteq=ok

  end function vecteq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function anotone(a,n)

    complex*16   a(3,3)
    complex*16   b(3,3)
    complex*16   one(3,3)
    real*8  zero
    logical ok
    integer i,j,n


    if (n.eq.1) then
       anotone=.true.
       return
    endif

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

    anotone=ok
    
  end function anotone
