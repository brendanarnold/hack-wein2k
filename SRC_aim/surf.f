      SUBROUTINE SURF(iatinit,readsurf)

!
!     AIM (1999-2000) by Javier D. Fuhr and Jorge O. Sofo
!
!     Instituto Balseiro and Centro Atomico Bariloche
!     S. C. de Bariloche - Rio Negro
!     Argentina
!     e-mail: fuhr@cab.cnea.gov.ar
!             sofo@cab.cnea.gov.ar
!
!     Calculates the surface corresponding to the atom IATINIT
!
      use crit
      use sphe
      implicit none

      real*8 cth,th,ph,wcth,wph,rs,themin,themax,phimin,phimax
      real*8 dr0,h,frmin,r0
      real*8 ct1,ct2,t1,t2,grho,v,vr,cutoff,rsmax,rsmin,rthe0,xy,xyz
      real*8 vcth,vph,vth,rthe,theta,phi,rr

      integer iatinit,nth,nph,index,nsimax
      integer i,j,ii,k,npmax,ith,iph,ith2,iph2


      logical deb,srch,stemp,readsurf
!     real*4 t1,t2

      INCLUDE 'param.inc'
      COMMON /debug/ deb
      COMMON /INTEG/ cth(NGAUSS),th(NGAUSS),ph(NGAUSS),wcth(NGAUSS), &
         wph(NGAUSS),rs(NGAUSS,NGAUSS),nth,nph
      COMMON /SRF/ themin,themax,phimin,phimax,h,frmin,r0,dr0, &
         index,nsimax
!      COMMON /CRIT/ pc(3,NNPOS*NDIF),evpc(3,NNPOS*NDIF), &
!         zpc(3,3,NNPOS*NDIF),pcrb(3,NNPOS*NDIF),npc,icpc(NNPOS*NDIF), &
!         npcm3
!      COMMON /SPHE/ POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO), &
!         RNOT(NATO),ndat,IATNR(NDIF),IOP(NDIF),MULT(NATO)

      dimension grho(3),v(3),vr(3)
      integer igmode
      logical trust
      common /ldmopts/igmode,trust

      rewind(22)

      ct1=cos(themin)
      ct2=cos(themax)
      call gauleg(ct1,ct2,cth,wcth,nth)
!      call gauleg(phimin,phimax,ph,wph,nph)
      call gaulep(phimin,phimax,ph,wph,nph)
!      write(*,*)'Quadrature ',nth,nph
!      do i=1,nth
!        write(6,*)'Theta quadrature ',i,cth(i),wcth(i)
!      enddo
!      do i=1,nph
!        write(6,*)'Phi quadrature ',i,ph(i),wph(i)
!      enddo
      do i=1,nth
        th(i)=acos(cth(i))
        if(.not.readsurf) then
          do j=1,nph
            rs(i,j)=0.d0
          enddo
        endif
      enddo
      
      npmax=17
      cutoff=0.0D0
      rsmax=0.0D0
      rsmin=100.0D0
      rthe0=r0
      srch=.false.
      
      do i=1,3
        v(i)=pos(i,index)
      end do

!     write(6,*) 'npmax in surf= ',npmax

      ith=0
      iph=0

      if(trust)goto 500
      if(.not.readsurf) then
!       srch=.true.
        do i=1,npc
          if((icpc(i).eq.-1)) then
            do j=1,3
              vr(j)=pc(j,i)-v(j)
            enddo
            xy=vr(1)*vr(1)+vr(2)*vr(2)
            xyz=xy+vr(3)*vr(3)
            xyz=sqrt(xyz)
            if (xy.lt.1.D-10) then
              vcth=1.d0
              if(vr(3).lt.0.d0) vcth=-vcth
              vph=0.d0
            else
              vcth=vr(3)/xyz
              vph=atan2(vr(2),vr(1))
            endif
            vth=acos(vcth)
            if(vth.lt.th(1)) then
              ith=0
            else
              if(vth.gt.th(nth)) then
                ith=nth
              else
                do ii=2,nth
                  if(vth.lt.th(ii)) then
                    ith=ii-1
                    goto 13
                  endif
                enddo
              endif
 13         endif
            if(vph.lt.ph(1)) then
              iph=0
            else
              if(vph.gt.ph(nph)) then
                iph=nph
              else
                do ii=2,nph
                  if(vph.lt.ph(ii)) then
                    iph=ii-1
                    goto 14
                  endif
                enddo
              endif
 14         endif
            do j=-1,2
              do k=-1,2
                ith2=ith+j
                iph2=iph+k
                stemp=(iph2.gt.0).and.(iph2.lt.nph+1)
                stemp=stemp.and.(ith2.gt.0).and.(ith2.lt.nth+1)
                if(stemp)then
                  if((rs(ith2,iph2).eq.0.d0)) then
                  rthe=rthe0
                  theta=th(ith2)
                  phi=ph(iph2)
                  if(deb) write(6,*) ':CALCULATING NP',theta,phi, &
                     rthe,npmax
                  call cputim(t1)
                  call rsurf(rr,grho,theta,phi,rthe,iatinit,npmax,srch)
                  call cputim(t2)
                  t2=t2-t1
                  rs(ith2,iph2)=rr
                  write(6,203) i,j,k,theta,phi,rr, &
                     wcth(ith2)*wph(iph2),t2
                  rthe0=rr
                  endif
                endif
              enddo
            enddo
          endif
        enddo

        rthe0=r0
!       srch=.true.

        if((nth.gt.8).and.(nph.gt.8)) then
          rthe=r0
          do i=1,2
            theta=th(i)
            do j=1,nph
              phi=ph(j)
              call cputim(t1)
              if(rs(i,j).eq.0.d0) then
                if(deb) write(6,*) ':CALC NP',theta,phi,rthe,npmax
                call rsurf(rr,grho,theta,phi,rthe,iatinit,npmax,srch)
                rs(i,j)=rr
              end if
              call cputim(t2)
              t2=t2-t1
              write(6,204) theta,phi,rs(i,j),wcth(i)*wph(j),t2
              if (j.eq.1) rthe0=rs(i,j)
              if (rr.lt.rsmin) rsmin=rs(i,j)
              if (rr.gt.rsmax) rsmax=rs(i,j)
              rthe=rs(i,j)
            end do
          end do
          rthe=r0
          do i=nth-1,nth
            theta=th(i)
            rthe=r0
            do j=1,nph
              phi=ph(j)
              call cputim(t1)
              if(rs(i,j).eq.0.d0) then
                if(deb) write(6,*) ':CALC NP',theta,phi,rthe,npmax
                call rsurf(rr,grho,theta,phi,rthe,iatinit,npmax,srch)
                rs(i,j)=rr
              end if
              call cputim(t2)
              t2=t2-t1
              write(6,204) theta,phi,rs(i,j),wcth(i)*wph(j),t2
              if (j.eq.1) rthe0=rs(i,j)
              if (rr.lt.rsmin) rsmin=rs(i,j)
              if (rr.gt.rsmax) rsmax=rs(i,j)
              rthe=rs(i,j)
            end do
          end do
          rthe=rs(2,1)
          do j=1,2
            phi=ph(j)
            do i=3,nth-2
              theta=th(i)
              t2=0.0
              call cputim(t1)
              if(rs(i,j).eq.0.d0) then
                if(deb) write(6,*) ':CALC NP',theta,phi,rthe,npmax
                call rsurf(rr,grho,theta,phi,rthe,iatinit,npmax,srch)
                rs(i,j)=rr
              end if
              call cputim(t2)
              t2=t2-t1
              write(6,204) theta,phi,rs(i,j),wcth(i)*wph(j),t2
              if (j.eq.1) rthe0=rs(i,j)
              if (rr.lt.rsmin) rsmin=rs(i,j)
              if (rr.gt.rsmax) rsmax=rs(i,j)
              rthe=rs(i,j)
            end do
          end do
          rthe=rs(2,nph-1)
          do j=nph-1,nph
            phi=ph(j)
            do i=3,nth-2
              theta=th(i)
              t2=0.0
              call cputim(t1)
              if(rs(i,j).eq.0.d0) then
                if(deb) write(6,*) ':CALC NP',theta,phi,rthe,npmax
                call rsurf(rr,grho,theta,phi,rthe,iatinit,npmax,srch)
                rs(i,j)=rr
              end if
              call cputim(t2)
              t2=t2-t1
              write(6,204) theta,phi,rs(i,j),wcth(i)*wph(j),t2
              if (j.eq.1) rthe0=rs(i,j)
              if (rr.lt.rsmin) rsmin=rs(i,j)
              if (rr.gt.rsmax) rsmax=rs(i,j)
              rthe=rs(i,j)
            end do
          end do
        endif
      endif     
      rthe0=r0
!     srch=.true.

500   do i=1,nth
        theta=th(i)
        rthe=rthe0
        do j=1,nph
          phi=ph(j)
          t2=0.0
          call cputim(t1)
          if(rs(i,j).eq.0.d0) then
            if (deb) write(6,*) ':CALCULATING NP',i,j, &
               theta,phi,rthe,npmax,rs(i,j)
            call rsurf(rr,grho,theta,phi,rthe,iatinit,npmax,srch)
            rs(i,j)=rr
          end if
          call cputim(t2)
          t2=t2-t1
          write(21,200) theta,phi,rs(i,j),wcth(i)*wph(j)
          write(6,204) theta,phi,rs(i,j),wcth(i)*wph(j),t2
          if (j.eq.1) rthe0=rs(i,j)
          if (rr.lt.rsmin) rsmin=rs(i,j)
          if (rr.gt.rsmax) rsmax=rs(i,j)
          rthe=rs(i,j)
        end do
      end do
      write(21,202) rsmin,rsmax

      return

 200  format(2F16.12,2D20.12)
 201  format(':RSUR ',2F12.8,2E16.8)
 202  format(2F15.10)
 203  format(':RSUR PC ',3i3,4E16.8,F10.4)
 204  format(':RSUR ',2F12.8,2E16.8,F10.4)
      end
      



      subroutine init_surf(iat,readcrit)
      use crit
      use sphe
      use atpos
      implicit none
      INCLUDE 'param.inc'
      
      real*8 cth,th,ph,wcth,wph,rs,br1,br2,br3,br4
      real*8 themin,themax,phimin,phimax,h,frmin,r0,dr0
      real*8 a,alpha,beta,gamma,gamma1
      real*8 v,dstmax,tpi,pi

      integer iat,nth,nph,index
      integer i,nsimax

      logical deb,stwo,sthree,sfour,ortho,readcrit

      COMMON /debug/ deb
      COMMON /INTEG/ cth(NGAUSS),th(NGAUSS),ph(NGAUSS),wcth(NGAUSS), &
         wph(NGAUSS),rs(NGAUSS,NGAUSS),nth,nph
      COMMON /BRAV/ BR1(3,3),BR2(3,3),br3(3,3),br4(3,3) 
      COMMON /SRF/ themin,themax,phimin,phimax,h,frmin,r0,dr0, &
         index,nsimax
!      COMMON /SPHE/ POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO), &
!         RNOT(NATO),ndat,IATNR(NDIF),IOP(NDIF),MULT(NATO)
!      COMMON /CRIT/ pc(3,NNPOS*NDIF),evpc(3,NNPOS*NDIF), &
!         zpc(3,3,NNPOS*NDIF),pcrb(3,NNPOS*NDIF),npc,icpc(NNPOS*NDIF), &
!         npcm3
!      COMMON /ATPOS/ ATP(3,nnpos),NPOS,NSA,NSB,NSC
      COMMON /LATT/ A(3),alpha,beta,gamma,gamma1,ortho

      DIMENSION v(3)
      
      read(5,*) index
      read(5,*) nth,themin,themax
      read(5,*) nph,phimin,phimax
      pi=2.D0*acos(-1.D0)
      tpi=2.D0*pi
!     Eliminate rounding errors
      do i=-8,8
        if(i .ne. 0)then
                if(abs(themin-tpi/i) .lt. 1d-5)themin=tpi/dble(i)
                if(abs(themax-tpi/i) .lt. 1d-5)themax=tpi/dble(i)
                if(abs(phimin-tpi/i) .lt. 1d-5)phimin=tpi/dble(i)
                if(abs(phimax-tpi/i) .lt. 1d-5)phimax=tpi/dble(i)
        else
                if(abs(themin) .lt. 1d-5)themin=0.D0
                if(abs(themax) .lt. 1d-5)themax=0.D0
                if(abs(phimin) .lt. 1d-5)phimin=0.D0
                if(abs(phimax) .lt. 1d-5)phimax=0.D0
        endif
      enddo
      write(6,*) 'ATOM ',index
      write(6,*) 'NI THEMAX THEMIN ',nth,themin,themax
      write(6,*) 'NJ PHIMIN PHIMAX ',nph,phimin,phimax
      
      read(5,*) h,frmin,nsimax
      if(h.lt.0.02d0) frmin=0.02d0
      if(frmin.gt.1.d0) frmin=1.d0
      write (6,*) 'H ',h,'  frmin=',frmin
      read(5,*) r0,dr0
      write(6,*) 'R0=',r0,'  dr=',dr0
      read(5,*) NSA,NSB,NSC
      write(6,120) NSA,NSB,NSC
      call gener(BR2)

      do i=1,3
        v(i)=pos(i,index)
      end do

      iat=index
      write(21,300) index,v
      write(21,300) nth,themin,themax,frmin
      write(21,301) nph,phimin,phimax

      if (.not.readcrit) then
!  FIND CRITICAL POINTS
        dstmax=min(30.d0,min(nsa*a(1),nsb*a(2),nsc*a(3)))
        
        stwo=.true.
        sthree=.false.
        sfour=.false.
        call critics(index,stwo,sthree,sfour,dstmax)
      endif
!      stop 'SURF: CRITICAL POINTS FINISHED'

 120  format('0NSHELL PARAMETERS IN X, Y AND Z-DIRECTION:' &
         ,3I5)         
 300  FORMAT(i3,3F18.12)
 301  FORMAT(i3,2F18.12)

      return
      end
