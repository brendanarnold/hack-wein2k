      SUBROUTINE CRITICS(index,stwo,sthree,sfour,dstmax)
!
!     AIM (1999-2000) by Javier D. Fuhr and Jorge O. Sofo
!
!     Instituto Balseiro and Centro Atomico Bariloche
!     S. C. de Bariloche - Rio Negro
!     Argentina
!     e-mail: fuhr@cab.cnea.gov.ar
!             sofo@cab.cnea.gov.ar
!
!
!     Search for critical points starting between
!     atom INDEX and its neighbors.
!     
!     STWO, STHREE, SFOUR: logical switches (TRUE/FALSE) indicating
!                          to search CP starting in the middle point
!                          of two, three or four atoms. One of these
!                          atoms is the INDEX.
!     DSTMAX: maximum distance to search for neighbors
!
      use crit
      use sphe
      use atpos
      implicit none
      INCLUDE 'param.inc'

      real*8 dstmax,dlimit,br1,br2,br3,br4
      real*8 a,alpha,beta,gamma,gamma1
      real*8 v1,v2,vt,dif,v3,v4,vi
      real*8 olddist,dmax,dist

      integer index,nnpos,ibest,ipbest
      integer nb,nshell,n1,n2,n3
      integer l,jj,nc,j,k,ii,ndif

      logical stwo,sthree,sfour,srho,sgrho,shrho
      logical deb,ortho

      parameter (dlimit=1.d-4)
      COMMON /DEBUG/ deb
      COMMON /BRAV/ br1(3,3),br2(3,3),br3(3,3),br4(3,3)
!      COMMON /SPHE/ POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO), &
!         RNOT(NATO),ndat,IATNR(NDIF),IOP(NDIF),MULT(NATO)
!      COMMON /ATPOS/ ATP(3,nnpos),NPOS,NSA,NSB,NSC
      COMMON /LATT/ A(3),alpha,beta,gamma,gamma1,ortho
!      COMMON /CRIT/ pc(3,NNPOS*NDIF),evpc(3,NNPOS*NDIF), &
!         zpc(3,3,NNPOS*NDIF),pcrb(3,NNPOS*NDIF),npc,icpc(NNPOS*NDIF), &
!         npcm3


      dimension v1(3),v2(3),vt(3),dif(3),v3(3),v4(3),vi(3)
      real*8,allocatable :: dists(:)
      integer,allocatable :: nnat(:),inat(:)
      integer,allocatable :: nr(:)
      integer,allocatable :: ibat(:),ipibat(:)
      ndif=ndat
      nnpos=npos
      allocate ( dists(NNPOS*NDIF))
      allocate ( nnat(NNPOS*NDIF),inat(NNPOS*NDIF))
      allocate ( nr(NNPOS*NDIF))
      allocate ( ibat(NNPOS*NDIF),ipibat(NNPOS*NDIF))
      srho=.true.
      sgrho=.false.
      shrho=.false.

      do jj=1,3
        vi(jj)=br4(1,jj)*pos(1,index)+ &
           br4(2,jj)*pos(2,index)+br4(3,jj)*pos(3,index)
      enddo
      
      nc=0
      do j=1,npos
        do k=1,ndat
          dist=0.d0
          do ii=1,3
            dif(ii)=pos(ii,index)-pos(ii,k)-atp(ii,j)
            dist=dist+dif(ii)*dif(ii)
          end do
          dist=sqrt(dist)
          if(.not.((dist.gt.dstmax).or.(dist.lt.0.001))) then
            nc=nc+1
            dists(nc)=dist
            nnat(nc)=k
            inat(nc)=j
          endif
        enddo
      enddo
      call ord(dists,nr,nc)
      nb=0
      olddist=0.d0
      nshell=0
      do n1=1,nc
        n2=nr(n1)
        n3=nnat(n2)
        if(dists(n1).lt.(2.d0*dists(1))) then
          if((dists(n1)-olddist).gt.dlimit) then
            nshell=nshell+1
            olddist=dists(n1)
            if(nshell.eq.5) goto 31
          endif
          nb=nb+1
          ibat(nb)=n3
          ipibat(nb)=inat(n2)
          write(6,*) ':NEIG ',index,n3,inat(n2),dists(n1)
        else
          goto 31
        endif
      enddo
 31   continue

      npc=0
      npcm3=0
      dmax=1.d-2
      
!
!.....SEARCH BETWEEN EACH PAIR OF ATOMS
!
      if(stwo) then
        do j=1,nb
          do ii=1,3
            v1(ii)=pos(ii,index)
            v2(ii)=pos(ii,ibat(j))+atp(ii,ipibat(j))
            vt(ii)=(v1(ii)+v2(ii))*0.5D0
          enddo
          call newcrit(vt,vi)
        enddo
      endif
!
!.....SEARCH BETWEEN EACH THREE ATOMS
!
      if(sthree) then
        do j=1,nb
          do k=j+1,nb
            do ii=1,3
              v1(ii)=pos(ii,index)
              v2(ii)=pos(ii,ibat(j))+atp(ii,ipibat(j))
              v3(ii)=pos(ii,ibat(k))+atp(ii,ipibat(k))
              vt(ii)=(v1(ii)+v2(ii)+v3(ii))/3.d0
            enddo
            call newcrit(vt,vi)
          enddo
        enddo
      endif

!
!.....SEARCH BETWEEN EACH FOUR ATOMS
!
      if(sfour) then
        do j=1,nb
          do k=j+1,nb
            do l=j+1,nb
              do ii=1,3
                v1(ii)=pos(ii,index)
                v2(ii)=pos(ii,ibat(j))+atp(ii,ipibat(j))
                v3(ii)=pos(ii,ibat(k))+atp(ii,ipibat(k))
                v4(ii)=pos(ii,ibat(l))+atp(ii,ipibat(l))
                vt(ii)=(v1(ii)+v2(ii)+v3(ii)+v4(ii))*0.25D0 !/4.d0
              enddo
              call newcrit(vt,vi)
            enddo
          enddo
        enddo
      endif

      return

 120  format('0NSHELL PARAMETERS IN X, Y AND Z-DIRECTION:',3I5)
 121  format(a4)
 202  format(':PC2',3F10.6,3E12.4,I4,2E12.4)
 203  format(':PC3',3F10.6,3E12.4,I4,2E12.4)
 204  format(':PC4',3F10.6,3E12.4,I4,2E12.4)
 205  FORMAT('AUTOVAL: ',3F16.8)
 206  FORMAT('AUTOVEC: ',/,3F16.8,/,3F16.8,/,3F16.8)
 207  FORMAT('POS: ',3F16.8)
 208  FORMAT('POS in base: ',3F16.8)
      end


      SUBROUTINE critaround(index)
      use crit
      use sphe
      implicit none
      INCLUDE 'param.inc'

      real*8 pi,tpi,fpi,sqpi,sqtpi
      real*8 dcthe,dphi,cthe0,phi0,cthe,phi,v(3),rmin
      real*8 sthe,r,dmax,dir,ngrho0,ngrho1,ngrho2
      real*8 stpos,stgrho,stngrho,stlen
      real*8 themin,themax,phimin,phimax,h0,frmin,r0,dr0
      real*8 tmin,tmax,pmin,pmax
      real*8 cthemin,cthemax

      integer index,nthe,nphi,npmax,ibest,ipbest
      integer i,j,iat,ipos,iatinit,iposinit,nstep,k,l
      integer index0,nsimax

      logical deb,srch,inter,ssave,ldeb

      COMMON /DEBUG/ deb
!      COMMON /SPHE/ POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO), &
!         RNOT(NATO),ndat,IATNR(NDIF),IOP(NDIF),MULT(NATO)
      COMMON /CTES/ pi,tpi,fpi,sqpi,sqtpi
      COMMON /FOLL/ stpos(3,MAXFOLL),stgrho(3,MAXFOLL), &
           stngrho(MAXFOLL),stlen(MAXFOLL),ldeb
      COMMON /SRF/ themin,themax,phimin,phimax,h0,frmin,r0,dr0, &
         index0,nsimax
!      COMMON /CRIT/ pc(3,NNPOS*NDIF),evpc(3,NNPOS*NDIF), &
!         zpc(3,3,NNPOS*NDIF),pcrb(3,NNPOS*NDIF),npc,icpc(NNPOS*NDIF), &
!         npcm3

      read(5,*) nthe,tmin,tmax
      read(5,*) nphi,pmin,pmax
!     Eliminate rounding errors
      do j=-8,8
        if(j .ne. 0)then
                if(abs(tmin-tpi/j) .lt. 1d-5)tmin=tpi/j
                if(abs(tmax-tpi/j) .lt. 1d-5)tmax=tpi/j
                if(abs(pmin-tpi/j) .lt. 1d-5)pmin=tpi/j
                if(abs(pmax-tpi/j) .lt. 1d-5)pmax=tpi/j
        endif
      enddo
      npc=0
      npcm3=0
      index0=index
      h0=0.08D0
      frmin=0.01D0
      rmin=0.1D0
      npmax=17
      srch=.false.
      call search_at(pos(1,index),r,inter,iatinit,iposinit,ibest,ipbest)
      dmax=0.1D0
      dir=-1.0D0
      ssave=.true.

      cthemin=cos(tmin)
      cthemax=cos(tmax)
      dcthe=(cthemax-cthemin)/nthe
      dphi=(pmax-pmin)/nphi
      cthe0=cthemin+dcthe*0.5D0
      phi0=pmin+dphi*0.5D0

      do i=1,nthe
         do j=1,nphi
            cthe=cthe0+(i-1)*dcthe
            sthe=sqrt(1.0-cthe*cthe)
            phi=phi0+(j-1)*dphi
            v(1)=pos(1,index)+rmin*sthe*cos(phi)
            v(2)=pos(2,index)+rmin*sthe*sin(phi)
            v(3)=pos(3,index)+rmin*cthe
            write(6,*) 'Starting follow from: ',v
            call follow(v,npmax,srch,iatinit,iposinit,iat,ipos, &
                 dmax,nstep,dir,ssave)
            write(6,*) 'nstep=',nstep
            if(nstep.lt.3) goto 13
            ngrho0=stngrho(1)
            ngrho1=stngrho(2)
            do k=3,nstep-1
               ngrho2=stngrho(k)
               if((ngrho0.gt.ngrho1).and.(ngrho2.gt.ngrho1)) then
                  do l=1,3
                     v(l)=stpos(l,k)
                  enddo
                  write(6,*) 'Starting search from: ',v
                  call newcrit(v,pos(1,index))
               endif
               ngrho0=ngrho1
               ngrho1=ngrho2
            enddo
            do l=1,3
               v(l)=stpos(l,nstep)
            enddo
            write(6,*) 'Starting search from: ',v
            call newcrit(v,pos(1,index))
!      stop 'END'
 13      enddo
      enddo
      return
      end

      SUBROUTINE NEWCRIT(vt,vi)
      use sphe
      use crit
      use atpos
      implicit none
      INCLUDE 'param.inc'

      real*8 br1,br2,br3,br4
      real*8 a,alpha,beta,gamma,gamma1
      real*8 vt,dif,vi,ev,z,chg,grho,hrho
      real*8 dmax,r,diff,dist,distmin,distmin2

      integer iat,ipos,ires,kk,i1,i2,i3,ibest,ipbest
      integer jj,k,ii,ist,ndif,nnpos,IATNR_MIN,IATNR_MIN2

      logical srho,sgrho,shrho
      logical deb,inter,found,ortho

      COMMON /DEBUG/ deb
      COMMON /BRAV/ br1(3,3),br2(3,3),br3(3,3),br4(3,3)
!      COMMON /SPHE/ POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO), &
!         RNOT(NATO),ndat,IATNR(NDIF),IOP(NDIF),MULT(NATO)
!      COMMON /ATPOS/ ATP(3,nnpos),NPOS,NSA,NSB,NSC
      COMMON /LATT/ A(3),alpha,beta,gamma,gamma1,ortho
!      COMMON /CRIT/ pc(3,NNPOS*NDIF),evpc(3,NNPOS*NDIF), &
!         zpc(3,3,NNPOS*NDIF),pcrb(3,NNPOS*NDIF),npc,icpc(NNPOS*NDIF), &
!         npcm3

      dimension vt(3),dif(3),vi(3)
      dimension ev(3),z(3,3),grho(3),hrho(3,3)
      real*8,allocatable :: dist1(:,:)
      ndif=ndat
      nnpos=npos
      allocate ( dist1(nnpos,ndif))
      dmax=1.d-2
      srho=.true.
      sgrho=.false.
      shrho=.false.

!      write(6,*) 'in newcrit vt=',vt

      inter=.true.
      do jj=1,3
         diff=br4(1,jj)*vt(1)+br4(2,jj)*vt(2)+ &
              br4(3,jj)*vt(3)
         diff=diff-vi(jj)
         if(ortho) diff=diff/a(jj)
         if(abs(diff).gt.(1.5d0)) then
!ccc-test            inter=.false.
            write(6,*) 'Initial point greater than 1.5 a'
         endif
      enddo
      if(inter) then
         call critic(vt,ev,z,dmax,ires)
         if(ires.eq.0) then
            found=.false.
            if(npc.gt.0) then
               do jj=1,npc
                  dist=0.d0
                  do kk=1,3
                     dist=dist+(vt(kk)-pc(kk,jj))*(vt(kk)-pc(kk,jj))
                  enddo
                  dist=sqrt(dist)
                  if(dist.lt.(1.d-4)) then 
                     found=.true.
                     write(6,*) 'CP already found'
                  endif
               enddo
            endif
            if(.not.found) then
               do jj=1,3
                  pcrb(jj,npc+1)=br4(1,jj)*vt(1)+ &
                       br4(2,jj)*vt(2)+br4(3,jj)*vt(3)
                  diff=pcrb(jj,npc+1)-vi(jj)
                  if(ortho) diff=diff/a(jj)
                  if(abs(diff).gt.(1.5d0)) then
                     found=.true.
                     write(6,*) 'CP distance greater than 1.5 a'
                  endif
               enddo
            endif
            if(.not.found) then
               npc=npc+1
               do jj=1,3
                  pc(jj,npc)=vt(jj)
                  evpc(jj,npc)=ev(jj)
                  do kk=1,3
                     zpc(kk,jj,npc)=z(kk,jj)
                  enddo
               enddo
               i1=ev(1)/abs(ev(1))
               i2=ev(2)/abs(ev(2))
               i3=ev(3)/abs(ev(3))
               icpc(npc)=i1+i2+i3
               if (icpc(npc).eq.-3) then
                  npcm3=npcm3+1
               endif
               write(6,*) 'New critical point found'
               write(6,207) (pc(ii,npc),ii=1,3)
               write(6,208) (pcrb(ii,npc),ii=1,3)
               write(6,205) ev
               write(6,206) ((zpc(ii,jj,npc),ii=1,3),jj=1,3)
               distmin=999.     
               distmin2=998. 
               iatnr_min=0    
               do jj=1,npos
                  do k=1,ndat
                     dist=0.d0
                     do ii=1,3
                        dif(ii)=pc(ii,npc)-pos(ii,k)-atp(ii,jj)
                        dist=dist+dif(ii)*dif(ii)
                     end do
                     dist1(jj,k)=sqrt(dist)
                     if(distmin.gt.dist1(jj,k)) then
                       iatnr_min2=iatnr_min
                       distmin2=distmin
                       iatnr_min=iatnr(k)
                       distmin=dist1(jj,k)
                     else if(distmin2.gt.dist1(jj,k)) then
                       distmin2=dist1(jj,k)
                       iatnr_min2=iatnr(k)
                     endif
                  enddo
               enddo
               do jj=1,npos
                  do k=1,ndat
                     if(dist1(jj,k).lt.distmin* 1.5d0) then
                        write(6,'(''distance to atom'',i4,f10.4,'' au'',f10.4,'' ang'')') &
                             iatnr(k),dist1(jj,k),dist1(jj,k)*.529177d0       
                     endif
                  enddo
               enddo
               call vgh_rho(vt,chg,grho,hrho,srho,sgrho,shrho,r, &
                    iat,ipos,ist,ibest,ipbest)
               write(22,202) (pcrb(jj,npc),jj=1,3), &
                    (ev(jj),jj=1,3),icpc(npc),ev(1)+ev(2)+ev(3),chg 
                 if(icpc(npc).eq.-1) then
                   write(6,202) (pcrb(jj,npc),jj=1,3), &
                    (ev(jj),jj=1,3),icpc(npc),ev(1)+ev(2)+ev(3),chg &
                    ,iatnr_min,distmin,iatnr_min2,distmin2
                 else
                   write(6,202) (pcrb(jj,npc),jj=1,3), &
                    (ev(jj),jj=1,3),icpc(npc),ev(1)+ev(2)+ev(3),chg 
                 endif
!$$$            else
!$$$               write(6,*) 'Critical point not found'
            endif
         endif
      endif
      
      return

 202  format(':PC ',3F10.6,3E12.4,I4,2E12.4,i4,f6.3,i4,f6.3)
 205  FORMAT('AUTOVAL: ',3F16.8)
 206  FORMAT('AUTOVEC: ',/,3F16.8,/,3F16.8,/,3F16.8)
 207  FORMAT('POS: ',3F16.8)
 208  FORMAT('POS in base: ',3F16.8)
      end


      SUBROUTINE CRITIC(v,ev,z,dmax,ires)

!     
!     Search for a critical point starting from point v
!     
!     It returns:
!       v:  the final point
!       ev: eigenvalues of the Hessian in the final point
!       z:  eigenvectors of the Hessian in the final point
!       ires: if ires==0 => CP found
!             if ires==1 => CP not found within the maximum steps
!
      use sphe
      implicit none
      INCLUDE 'param.inc'

      real*8 v,ev,z,dmax,br1,br2,br3,br4
      real*8 vt,chg,grho,hrho,vold,y,dc
      real*8 dv,dg,dltcmax,vnorm,r,rhomin

      integer ires,istep,iat,ipos,ibest,ipbest
      integer ipiv,lwork,jj,i,j,id,nrot,info,MAXSTEP,ist

      logical deb,srho,sgrho,shrho

      COMMON /BRAV/ br1(3,3),br2(3,3),br3(3,3),br4(3,3)
      COMMON /DEBUG/ deb
!      COMMON /SPHE/ POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO), &
!         RNOT(NATO),ndat,IATNR(NDIF),IOP(NDIF),MULT(NATO)
      dimension v(3),vt(3),grho(3),hrho(3,3),vold(3),y(3,3)
      dimension ipiv(3),dc(3)
      dimension z(3,3),ev(3)
      PARAMETER (MAXSTEP=1500)
!      data rhomin/1.d-6/
      parameter(rhomin=1.D-6)      
      lwork=30
      istep=0
!      write(6,308) v
!      do jj=1,3
!        vt(jj)=br4(1,jj)*v(1)+
!     $     br4(2,jj)*v(2)+br4(3,jj)*v(3)
!      enddo
!      write (6,309) vt
      srho=.true.
      sgrho=.true.
      shrho=.true.
      ires=0
      
      call vgh_rho(v,chg,grho,hrho,srho,sgrho,shrho,r,iat,ipos,ist,ibest,ipbest)
      if(ist.eq.1.or.chg.lt.rhomin) then
        ires=1
        goto 44
      endif
      
!      write(6,200) grho
!      write(6,201) ((hrho(i,j),j=1,3),i=1,3)

      dg=1.0d0
      dv=1.0d0
      dg = vnorm(grho)

      do while((dv.gt.(1.d-8)).and.(dg.gt.(1.d-7)) &
         .and.(istep.lt.MAXSTEP))
!      do while((dg.gt.(1.d-7))
!     $   .and.(istep.lt.MAXSTEP))
        istep=istep+1
        do i=1,3
          vold(i)=v(i)
        enddo

        call ludcmp(hrho,3,3,ipiv,id,info)
        if(info.ne.0) then
          write(6,*) 'Error inverting hrho:'
          do i=1,3
            write(6,*) (hrho(i,j),j=1,3)
          enddo
          ires=1
          goto 44
!          stop 'ERROR INVERTING HESSIAN'
        endif
        do  i=1,3
          do j=1,3
            y(i,j)=0.
          enddo 
          y(i,i)=1.
        enddo 
        do j=1,3 
          call lubksb(hrho,3,3,ipiv,y(1,j))
        enddo         

        dltcmax = 0.d0
        do i = 1, 3
          dc(i) = 0.d0
          do j = 1, 3
            dc(i) = dc(i) + grho(j) * y(j,i) 
          end do
          dltcmax = dltcmax + dc(i)**2
        end do
        dltcmax = sqrt(dltcmax)
        if(dltcmax.gt.dmax) then
          do i=1,3
            dc(i)=dc(i)*dmax/dltcmax
          end do
        endif
        do i=1,3
          v(i) = v(i) - dc(i)
        end do

        call vgh_rho(v,chg,grho,hrho,srho,sgrho,shrho,r,iat,ipos,ist,ibest,ipbest)
        if(ist.eq.1.or.chg.lt.rhomin) then
          ires=1
          goto 44
        endif
        dg = vnorm(grho)
        
        if (deb) then
          write(6,202) ((y(i,j),j=1,3),i=1,3)
          write(6,206) dc
          write(6,*) 'STEP ',istep
          write(6,205) v
          write(6,200) grho
          write(6,*) ':DGRAD ',dg
          write(6,201) ((hrho(i,j),j=1,3),i=1,3)
        end if
        dv=0.d0
        do i=1,3
          dv=dv+(v(i)-vold(i))*(v(i)-vold(i))
        enddo
      end do

      if(istep.ge.MAXSTEP) then
!        write(6,*) 'istep.gt.MAXSTEP !!!!'
        ires=1
      endif

      do jj=1,3
        vt(jj)=br4(1,jj)*v(1)+ &
           br4(2,jj)*v(2)+br4(3,jj)*v(3)
      enddo
!      write(6,208) v
!      write(6,209) vt

      call vgh_rho(v,chg,grho,hrho,srho,sgrho,shrho,r,iat,ipos,ist,ibest,ipbest)
      if(ist.eq.1.or.chg.lt.rhomin) then
        ires=1
        goto 44
      endif
!      write(6,200) grho
!      write(6,201) ((hrho(i,j),j=1,3),i=1,3)
      call ludcmp(hrho,3,3,ipiv,id,info)
      if(info.ne.0) then
        write(6,*) 'Error inverting hrho:'
        do i=1,3
          write(6,*) (hrho(i,j),j=1,3)
        enddo
        ires=1
        goto 44
!        stop 'ERROR INVERTING HESSIAN'
      endif
      do  i=1,3
        do j=1,3
          y(i,j)=0.
        enddo 
        y(i,i)=1.
      enddo 
      do j=1,3 
        call lubksb(hrho,3,3,ipiv,y(1,j))
      enddo         
!      write(6,202) ((y(i,j),j=1,3),i=1,3)

      call vgh_rho(v,chg,grho,hrho,srho,sgrho,shrho,r,iat,ipos,ist,ibest,ipbest)
      if(ist.eq.1.or.chg.lt.rhomin) then
        ires=1
        goto 44
      endif
      call jacobi(hrho,3,3,ev,y,nrot)
      do i=1,3
        do j=1,3
          z(i,j)=y(j,i)
        enddo
      enddo

!      write(6,203) (ev(i),i=1,3)
!      write(6,204) ((z(i,j),i=1,3),j=1,3)

 44   return

 200  FORMAT(':GRAD ',3F16.8)
 201  FORMAT(':HESSIAN ',/,3F16.8,/,3F16.8,/,3F16.8)
 202  FORMAT(':HESSIAN^(-1) ',/,3F16.8,/,3F16.8,/,3F16.8)
 203  FORMAT(':AUTOVAL ',3F16.8)
 204  FORMAT(':AUTOVEC ',/,3F16.8,/,3F16.8,/,3F16.8)
 205  FORMAT(':POS ',3F16.8)
 206  FORMAT(':DC ',3F16.8)
 208  FORMAT(':POSOUT ',3F16.8)
 209  FORMAT(':RBPOSOUT ',3F16.8)
 308  FORMAT(':POSIN ',3F16.8)
 309  FORMAT(':RBPOSIN ',3F16.8)
      end



      subroutine ord(a,nr,imax)             
!     ORDERS ELEMENTS IN ARRAY A INCREASING IN SIZE          
!       REORDERS CORRESPONDING INDICES (NR)              
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL cont                                             
      DIMENSION a(*),nr(*)                                   
      do i=1,imax                                                    
        nr(i)=i    
      end do
  100 i=1                                                 
      cont=.false.                                        
  110 i=i+1                                                    
      if(a(i).lt.a(i-1)) then                                 
!       INTERCHANGE I AND (I-1) ELEMENT IN ARRAYS A AND NR  
        hlp=a(i)                                           
        a(i)=a(i-1)                                           
        a(i-1)=hlp                                             
        nhlp=nr(I)                                                
        nr(i)=nr(i-1)                                      
        nr(i-1)=nhlp                                           
        cont=.true.                                          
      endif                                          
      if(i.lt.imax) goto 110                                     
      if(cont) goto 100                                  
      return                                                 
      end                                                

