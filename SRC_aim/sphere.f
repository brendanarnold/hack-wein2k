      SUBROUTINE SPHERE(v,r,chg,grho,hrho,iat,srho,sgrho,shrho,ist)

!
!     AIM (1999-2000) by Javier D. Fuhr and Jorge O. Sofo
!
!     Instituto Balseiro and Centro Atomico Bariloche
!     S. C. de Bariloche - Rio Negro
!     Argentina
!     e-mail: fuhr@cab.cnea.gov.ar
!             sofo@cab.cnea.gov.ar
!
!     This routine calculates the value, the gradient, and the Hessian
!       of the electron density for a point inside a MT.
!
      use sphe
      use rad
      use atpos
      implicit none
!     
!.... POINT IN SPHERE IATNR(IAT)
!     
      complex*16 yl,dtyl,dtdtyl,dfdtyl

      real*8 v,r,chg,grho,hrho,tau,otau,oiz,invoiz,br1,br2,br3,br4,a
      real*8 alpha,beta,gamma,gamma1
      real*8 vt,vt1,vtold,xy
      real*8 cth,sth,cfi,sfi,mat1,mat2

      integer iat,ist,iz,iord,j
      integer jatom,jatom1,index,ir,jj,ipos

      logical ortho,srho,sgrho,shrho,deb,er,transp,notransp

      INCLUDE 'param.inc'
      COMMON /DEBUG/ deb
      COMMON /SYM2/ tau(3,NSYM),otau(3,NSYM),oiz(3,3,NSYM), &
         invoiz(3,3,NSYM),iz(3,3,NSYM),iord
!      COMMON /ATPOS/ ATP(3,nnpos),NPOS,NSA,NSB,NSC
      COMMON /BRAV/ BR1(3,3),BR2(3,3),br3(3,3),br4(3,3)
      COMMON /LATT/ A(3),alpha,beta,gamma,gamma1,ortho
!      COMMON /SPHE/ POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO), &
!         RNOT(NATO),ndat,IATNR(NDIF),IOP(NDIF),MULT(NATO)
      COMMON /YLMS/ yl((lmax2+4)*(lmax2+4)),dtyl((lmax2+3)*(lmax2+3)), &
         dtdtyl((lmax2+2)*(lmax2+2)),dfdtyl((lmax2+2)*(lmax2+2))
!      common /rad/ rm(nrad,nato),jri(nato)

      DIMENSION v(3),vt(3),vt1(3),vtold(3)
      DIMENSION grho(3),hrho(3,3),mat1(3,3),mat2(3,3)

      if(.not.(srho.or.sgrho.or.shrho)) goto 25
      transp=.true.
      notransp=.false.
      ist=0
      JATOM=IABS(IATNR(IAT))
!     
      if(deb) WRITE(6,*) 'JATOM,IAT,IOP,POINT,DIFFVECTOR ', &
         JATOM,IAT,IOP(IAT),V ,VT
      vt(1)=v(1)
      vt(2)=v(2)
      vt(3)=v(3)

!     TRANSFER COORDINATES FROM SPHERE IAT TO JATOM               
      if(ortho) then
        call ROTATE(VT,IZ(1,1,IOP(IAT)),TAU(1,IOP(IAT)),A)
      else
        call rotato(vt,oiz(1,1,iop(iat)),otau(1,iop(iat)))
      end if
!     WRITE(6,*) 'ROTATED VECTOR',VT 

!     REDUCE VECTOR TO SMALLEST POSSIBLE ONE
      index=1
      jatom1=1
 31   if (jatom1.NE.jatom) then 
        index=index+mult(jatom1)
        jatom1=jatom1+1
        goto 31
      end if 
      do j=1,3
        vtold(j)=vt(j)
      enddo
      call REDUC(vt,index,ipos,er)
      if (er) then
!        if (deb) then
!          write(6,61) 'Error in reducing vector',vtold,vt
!          write(6,*) 'Error in reducing vector hex coord',vt1
!          write(6,61) 'Original vector',v
!          write(6,61) 'Output from reduc ',vt
61       format(a,6F12.6)
!        end if
!        if (deb) then
!          do jj=1,3
!            vt1(jj)=br4(1,jj)*v(1)+br4(2,jj)*v(2)+br4(3,jj)*v(3)
!          enddo 
!          write(6,*) 'Original vector hex coord',vt1
!        end if
        ist=1
      goto 25
!     stop 'ERROR !!'
      endif
      if(deb) WRITE(6,*) 'REDUCED VECTOR',VT

!     TRANSFER COORDINATES INTO LOCAL COORDINATE SYSTEM    
      call ROTAT(VT,ROTLOC(1,1,JATOM))
      if(deb) WRITE(6,*) 'LOCAL VECTOR', VT

!     
      if (r.lt.5.0E-5) r=5.0E-5
      ir=1+log(r/rnot(jatom))/dx(jatom)
      if (ir.lt.1) ir=1
      if (ir.gt.jri(jatom)) stop 'ERROR ir.gt.jri(jatom)'

      call cossin(sth,cth,sfi,cfi,xy,vt)
!      xy=vt(1)*vt(1)+vt(2)*vt(2)
!      xyz=xy+vt(3)*vt(3)
!      if (xyz.lt.snull) then
!        xy=0.d0
!        xyz=0.d0
!        cth=one
!        sth=zero                  
!        cfi=one                  
!        sfi=zero
!      else  
!        if (xy.lt.snull) then
!          if (deb) write (6,*) ':SPHERE xy.lt.snull'
!          write (6,*) ':SPHERE xy.lt.snull'
!          xy=0.d0
!          cth=one
!          if (vt(3).lt.zero)  cth=-cth 
!          sth=zero                  
!          cfi=one                  
!          sfi=zero
!        else
!          xy=sqrt(xy)      
!          xyz=sqrt(xyz)   
!          cth=vt(3)/xyz    
!          sth=xy/xyz      
!          cfi=vt(1)/xy    
!          sfi=vt(2)/xy
!        endif
!      endif
      call ylm(cth,sth,cfi,sfi,lmax2+1,yl)
      if(sgrho.or.shrho) then
        call dtylm(cth,sth,cfi,sfi,xy,lmax2,yl,dtyl)
      endif
      if(shrho) then
        call dtdtylm(cth,sth,cfi,sfi,xy,lmax2,yl,dtyl,dtdtyl)
      endif
      if (srho) then
        CALL CHARGE(CHG,IR,R,JATOM,IATNR(IAT))
        if (deb) write(6,*) 'rhosphere',chg
      endif
      if (sgrho) then
        call grhosphe(grho,ir,r,jatom,xy,cth,sth,cfi,sfi,iatnr(iat))
        if (deb) write(6,*) 'sphere grho',grho(1),grho(2),grho(3)

        call rotat_back(grho,rotloc(1,1,jatom))
        if (deb) write(6,*)'after rotat_back',grho(1),grho(2),grho(3)

        call rotate_back(grho,invoiz(1,1,iop(iat)))
        if (deb) write(6,*)'after rotate_back',grho(1),grho(2),grho(3)
      endif
      if (shrho) then
        call hrhosphe(hrho,ir,r,jatom,xy,cth,sth,cfi,sfi,iatnr(iat),ist)
        if (ist.ne.0) then
          goto 25
        endif
        if (deb) write (6,*) 'sphere hrho: ',hrho

        call mat3prod(mat1,oiz(1,1,iop(iat)),rotloc(1,1,jatom), &
           notransp,transp)
        if (deb) write (6,*) 'sphere oiz . Tr[rotloc]: ',mat1

        call mat3prod(mat2,mat1,hrho,notransp,notransp)
        if (deb) write (6,*) 'sphere oiz . Tr[rotloc] . hrho: ',mat2

        call mat3prod(mat1,mat2,rotloc(1,1,jatom),notransp,notransp)
        if (deb) write (6,*) &
           'sphere oiz . Tr[rotloc] . hrho . rotloc: ',mat1

        call mat3prod(hrho,mat1,invoiz(1,1,iop(iat)),notransp,notransp)
        if (deb) write (6,*) &
           'sphere oiz . Tr[rotloc] . hrho . rotloc . invoiz: ',mat1
      endif

 25   return
      end

