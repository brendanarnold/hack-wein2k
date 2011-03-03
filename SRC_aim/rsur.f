      subroutine rsur()
!
!     AIM (1999-2000) by Javier D. Fuhr and Jorge O. Sofo
!
!     Instituto Balseiro and Centro Atomico Bariloche
!     S. C. de Bariloche - Rio Negro
!     Argentina
!     e-mail: fuhr@cab.cnea.gov.ar
!             sofo@cab.cnea.gov.ar
!
      use sphe
      use atpos
      implicit none
      include 'param.inc'


      real*8 br1,br2,br3,br4,themin,themax
      real*8 phimin,phimax,h,frmin,r0,dr0,a,alpha,beta,gamma,gamma1
      real*8 grho,v,dstmax,cutoff,rr,t1,t2,r,theta,phi

      integer index,i,npmax,ibest,ipbest
      integer iat,ipos,np,nsimax

      logical ortho,inter,deb,srch
      logical stwo,sthree,sfour

      COMMON /DEBUG/ deb
!      COMMON /SPHE/ POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO), &
!         RNOT(NATO),ndat,IATNR(NDIF),IOP(NDIF),MULT(NATO)
      COMMON /BRAV/ BR1(3,3),BR2(3,3),br3(3,3),br4(3,3)
      COMMON /SRF/ themin,themax,phimin,phimax,h,frmin,r0,dr0, &
         index,nsimax
      COMMON /LATT/ A(3),alpha,beta,gamma,gamma1,ortho
!      COMMON /ATPOS/ ATP(3,nnpos),NPOS,NSA,NSB,NSC
      dimension grho(3),v(3)

      srch=.false.
      deb=.false.

      read(5,*) index,np
      read(5,*) h,frmin
      if(h.lt.0.02d0) frmin=0.02d0
      if(frmin.gt.1.d0) frmin=1.d0
      write (6,*) 'H ',h,'  rmin=',frmin
      read(5,*) r0,dr0
      write(6,*) 'R0=',r0,'  dr=',dr0
      read(5,*) NSA,NSB,NSC
      write(6,120) NSA,NSB,NSC
      call gener(BR2)
      do i=1,3
        v(i)=pos(i,index)
      end do
      npmax=9
      cutoff=0.0
      nsimax=0

      dstmax=min(20.d0,min(nsa*a(1),nsb*a(2),nsc*a(3)))

      stwo=.true.
      sthree=.true.
      sfour=.false.
      write (6,*) 'finding critical points'
      call critics(index,stwo,sthree,sfour,dstmax)

      do i=1,np
        read(5,*) theta,phi
        call search_at(v,r,inter,iat,ipos,ibest,ipbest)
        call cputim(t1)
        call rsurf(rr,grho,theta,phi,r0,iat,npmax,srch)
        call cputim(t2)
        t2=t2-t1
        write(21,200) theta,phi,rr
        write(6,204) theta,phi,rr,t2
      enddo

      return

 120  format('0NSHELL PARAMETERS IN X, Y AND Z-DIRECTION:' &
         ,3I5)         
 200  format(2F12.8,2F15.10)
 201  format(':RSUR ',2F14.10,F15.10)
 204  format(':RSUR ',2F14.10,2F15.10,F10.4)
      end

