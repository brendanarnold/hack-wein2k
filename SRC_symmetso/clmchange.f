       subroutine  clmchange(ifile,ifw,nat,ieq,jrj,iord,iz,tau,ksym)
       implicit real*8(a-h,o-z)
      INCLUDE 'param.inc'
       complex*16 imag,ro2
       integer iz(3,3,48),ksym(48),ksym2(48),kv(3),kvtest(3,48)
       dimension TAU(3,noper)
! this program serves to read clm/potential files in FLAPW WIEN program
! then to construct the new files in the enlarged part of BZ
! when s-o coupling and magnetization are present
! jri(nat) is number of radial points for atom nat
        CHARACTER*19     nkktext
        character a(111)
        dimension jrj(*),ieq(*)
        character, allocatable :: al(:,:)
        real*8,allocatable ::rho(:,:),ro1(:,:)
        integer,allocatable :: lmmax(:),kv1(:,:)
        allocate(lmmax(nat))
!
        imag=(0.d0,1.d0)
        tpi=2.d0*acos(-1.d0)
!        read(5,*)nat,(jri(j),j=1,nat)
!        read(5,*)cclm1
333     format(' number of atoms',i3,'number of points',2i7)
444     format(' file 1 vsn*',f5.1)
        read(ifile,1980,end=999,err=999)(a(i),i=1,79)
        a(78)='s'
        a(79)='o'
      write(6,*) 'Generating clm file',ifile
        write(ifw,1980)(a(i),i=1,79)
        read(ifile,1980)(a(i),i=1,79)
        write(ifw,1980)(a(i),i=1,79)
        read(ifile,1980)a(1)
        write(ifw,1980)a(1)
        jatso=0
        do 79 jatom=1,nat
        do 80 jeq=1,ieq(jatom)
        if(jeq.gt.1) write(6,*) 'duplicating atom',jatom
        jatso=jatso+1
        if(jeq.eq.1) read(ifile,2000)(a(i),i=1,15),j
       write(ifw,1990)jatso
 1990 FORMAT(3X,'ATOM NUMBER=',I3,5X,10A4)                             
        if(jeq.eq.1)then
          read(ifile,2000)(a(i),i=1,15),lmmax(jatom)
          allocate (al(lmmax(jatom),111),rho(lmmax(jatom),nrad))
        endif
       write(ifw  ,2111)lmmax(jatom)
 2111 FORMAT(3X,'NUMBER OF LM',I3)                                   
        if(jeq.eq.1)then
        do 35 lm1=1,lmmax(jatom)
       read(ifile,2000)a(1)
       read(ifile,2000)a(1)
        read(ifile,1980)(al(lm1,i),i=1,25)
        read(ifile,1980)a(1)
        read(ifile,2020)(rho(lm1,j),j=1,jrj(jatom))
 35     continue
       endif
       do 36 lm1=1,lmmax(jatom)
       write(ifw,2000)a(1)
       write(ifw,2000)a(1)
       write(ifw    ,1980)(al(lm1,i),i=1,25)
       write(ifw,    1980)a(1)
        write(ifw,2020)(rho(lm1,j),j=1,jrj(jatom))
36     continue
       do 88 i=1,6
       if(jeq.eq.1)read(ifile,1980)a(1)
        write(ifw,1980)a(1)
 88    continue
 80     continue
        deallocate (al,rho)
 79     continue
        if(ifile.le.26) then
          deallocate(lmmax)
          return
        endif
        read(ifile,1980)(a(i),i=1,79)
        write(ifw,1980)(a(i),i=1,79)
        read(ifile,1980)a(1)
        write(ifw,1980)a(1)
!
!        read(ifile,2030)(a(i),i=1,13),nwav1
  read(ifile,'(a19)') nkktext
  read(nkktext,'(9x,i10)',err=6767) nwav1
  goto 6768
 6767 read(nkktext,'(13x,i6)') nwav1
 6768 continue
        allocate (ro1(2,nwav1*iord),kv1(3,nwav1*iord))

        do 140 j=1,nwav1
 140    read(ifile,2070)(kv1(jx,j),jx=1,3),(ro1(k,j),k=1,2)
 345    format(i5,2e15.7)
 2070  format(3x,3i5,2e19.12)
 2030   format(13a1,i6)
 1980   format(80a1)
 2000   format(15a1,i3)
 2020   format(3x,4e19.12)
!  **** clm file read *****
        nw=nwav1
      write(6,*)' nw     k     j         knew      rhonew              kold'
!      print *, nwav1,iord,(ksym(i),i=1,iord)
        do 20 j=1,nwav1
!        do 20 j=1,10
        do 22 k=1,iord
        ksym2(k)=ksym(k)
        do 22 i1=1,3
        kvtest(i1,k)=0
        do 22 i2=1,3
        kvtest(i1,k)=kvtest(i1,k)+iz(i1,i2,k)*kv1(i2,j)
22      continue
          do 23 k1=1,iord
          if(ksym(k1).eq.0)then
           do k2=1,iord
             if(ksym2(k2).eq.1) then
             if(kvtest(1,k1).eq.kvtest(1,k2).and. &
              kvtest(2,k1).eq.kvtest(2,k2).and. &
              kvtest(3,k1).eq.kvtest(3,k2)) goto 23
             endif
           enddo
! vector kv is new, generate its little-star and block equal in ksym2
           do k=1,iord
             if(ksym(k).eq.1) then
               do  i1=1,3
               kv(i1)=0
               do  i2=1,3
               kv(i1)=kv(i1)+iz(i1,i2,k)*kvtest(i2,k1)
               enddo
               enddo
                 do k2=1,iord
                  if(kv(1).eq.kvtest(1,k2).and. &
                     kv(2).eq.kvtest(2,k2).and. &
                     kv(3).eq.kvtest(3,k2)) ksym2(k2)=1
                 enddo
             endif
            enddo
      
        nw=nw+1

         TK = 0.0D+0
         DO l = 1, 3
            TK = TK + TAU(l,k1)*kv1(l,J)*TPI
         enddo
        do 25 l=1,3
        kv1(l,nw)=kvtest(l,k1)
25      continue
        ro2=cmplx(ro1(1,j),ro1(2,j))* EXP(DCMPLX(TK)*(-IMAG))
        ro1(1,nw)=dble(ro2)
        ro1(2,nw)=aimag(ro2)
        if(ifw.eq.55) write(6,500)nw,j,k,(kv1(l,nw),l=1,3), &
        (ro1(l,nw),l=1,2),(kv1(l,j),l=1,3)
500     format(3i6,3x,3i4,2f10.6,3i4)

          endif
23      continue

21      continue
20      continue
        write(ifw,2031) nw
 2031   format('         ',i10,' NUMBER OF PW')
        do  j=1,nw
        write(ifw,2070)(kv1(jx,j),jx=1,3),(ro1(k,j),k=1,2)
        enddo
        deallocate(lmmax)
        deallocate (ro1,kv1)
        return
!
 999    write(6,*) ' Warning: clmfile',ifile,' skipped'
        deallocate(lmmax)
        return
        end
