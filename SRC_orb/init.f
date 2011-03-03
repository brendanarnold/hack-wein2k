      SUBROUTINE INIT(chmix)
      use opr
      use ldau
      use struct
      use exex
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'param.inc'
!
      CHARACTER*80    test
      character*5   chmix
      character*4     spin(-1:2)
      character*6     test1,test2
      common/print/   ipr,ipr3j,iprv
      common/updn/    nup,spin
      CHARACTER*4     LATTIC,amod(3)
      CHARACTER*80    TITLE
!      CHARACTER*10    ANA
      COMMON/CHAR/     TITLE,LATTIC
!      COMMON/STRUCT/  AA,BB,CC,VOL,RO(NATO),DXM(NATO),RMT(NATO), &
!                      ANA(NATO),JRJ(NATO),NATOM(NATO),NTYPE
      COMMON/REST/    PIA(3),ALPHA(3)
!                                    ,POS(3,NDIF)
!      real*8          U(labc,nato),J(labc,nato)
!      common /ldau/   U,J,nsic
!      complex*16      dmat(nato,3,-3:3,-3:3,2)
!      common/opr/   amix,Bexten,VSP(NRAD,NATO),EOP(LABC,NATO),zz(nato), &
!                    pop(labc,nato),al(3,labc,nato),alm(labc,nato), &
!                    dmat
!      common/opi/     nmod,natorb,iatom(nato),nlorb(nato),nmodop, &
!                      lorb(labc,nato),NCALC(NATO)
!      common/rotmat/  rotloc(3,3,nato)
      common/angles/  theta,phi,thetaloc,philoc
      DIMENSION       xms(3)
      character*34   chmod(5)
      REAL*8,ALLOCATABLE    :: ee(:,:)
      pi=4.d0*atan(1.d0)
!**********************************************************************
      test1='ATOMIC'
! nmod defines type of Vorb 1..LDA+U, 2..OP, 3..Bext
      write(6,199)spin(nup)
      chmod(1)=' LDA+U'
      chmod(2)=' Orbital Polarization'
      chmod(3)=' Interaction with Bext'
      amod(1)='LDAU'
      amod(2)='magn'
      amod(3)='Bext'
! input for all potentials
! nmod defines the type of potential 1...LDA+U, 2...OP, 3...Bext
! natorb is number of atoms to which Vorb is applied
! amix is coefficient for mixing of Vorb
! chmix is PRATT or BROYD
! ipr is printing option, larger ipr, longer output
      ipr3j=0
      iprv=0
      read(5,*)nmod,natorb,ipr
      read(5,702)chmix,amix
 702  format(a5,f8.2)    
      write(6,200)chmod(nmod)
      allocate (iatom(natorb),nlorb(natorb),lorb(labc,natorb),ncalc(natorb))
      allocate (U(labc,natorb),J(labc,natorb),pop(labc,natorb))
      allocate ( dmat(natorb,labc,-3:3,-3:3,2),al(3,labc,natorb),alm(labc,natorb))

      allocate (fxck(42,natorb),fupk(0:6,natorb),nlm1(natorb),lm1(2,42))
      allocate (espxc(natorb))

      do i=1,natorb
! iatom - index of atom in struct file
! nlorb number of orbital moments for which Vorb will be applied
! lorb .. orbital numbers
       read(5,*)iatom(i),nlorb(i),(lorb(k,i),k=1,nlorb(i))
       write(6,202)iatom(i),(lorb(k,i),k=1,nlorb(i))
       if(nlorb(i).gt.labc) stop 'nlorb gt labc'
      enddo
! input for LDA+U
      if(nmod.eq.1)then
       read(5,*)nldau
       if(nldau.lt.4)then
        nexex=0
        if(nldau.eq.1)then
         write(21,*)'       Fully Localized Limit method'
         write(6,*)' Fully Localized Limit method'
        endif
        if(nldau.eq.0)then
         write(21,*)'       Around the mean field method'
         write(6,*)' Around the mean field method'
        endif
        if(nldau.eq.2)then
         write(21,*)'       Mean field Hubbard method'
         write(6,*)' Mean field Hubbard method'
        endif
       do i=1,natorb
          do li=1,nlorb(i)
           read(5,*)U(li,i),J(li,i)
           iat=iatom(i)
           l=lorb(li,i)
           write(6,207)iat,l,U(li,i),J(li,i)
           write(21,207)iat,l,U(li,i),J(li,i)
          enddo
       enddo
       elseif(nldau.eq.4)then       !EECE starts
        nexex=1
        nsp=1           !standard: include nonspherical Vxc to exclude it put nsp=0
        hybr=0.         ! for bybrids mixing in hybr*(Fock exchange)+(1-hybr)*(DFT exchange)
!       read(5,*,err=333)nsp,hybr
        read(5,*,err=333)hybr
333    continue
!      if(nsp.eq.0)then
!       write(6,*)' Only spherical Vxc included'
!       write(21,*)' Only spherical Vxc included'
!      else
!       write(6,*)' Full Vxc potential included'
!       write(21,*)' Full Vxc potential included'
!      endif
       if(hybr.gt.0.000001)then
        write(6,332)hybr
        write(21,332)hybr
332     format('  Hybrid potential with ',f9.4,' of Fock exchange')
       endif
       
       do i=1,natorb
        read(50,*)ia,nlm1(i)
        if(nlm1(i).gt.42) stop 'nlm1 too large'
        write(6,*)' Exact exchange or hybrid metod for atom',ia
        do ll=0,6,2
         read(50,*)ii,k,fupk(ll,i)
         if((i.eq.1).and.(ll.eq.0))then
          write(6,550)i,ia,ll,fupk(ll,i)
          write(21,550)i,ia,ll,fupk(ll,i)
         else
          write(6,560)i,ia,ll,fupk(ll,i)
          write(21,560)i,ia,ll,fupk(ll,i)
         endif
550   format(3i4,f12.6,' iat, atom, k, F^k Slater integral from lapw0, Ry ')
560   format(3i4,f12.6)
        enddo
       enddo
        do i=1,natorb
         espxc(i)=0.
         do k=1,nlm1(i)
          read(50,*)jatom,lmo,(lm1(ii,k),ii=1,2),fxck(k,i)
          if(k.eq.1)then
           write(6,552)i,jatom,lmo,(lm1(ii,k),ii=1,2),fxck(k,i)
           write(21,552)i,jatom,lmo,(lm1(ii,k),ii=1,2),fxck(k,i)
          else
           write(21,562)i,jatom,lmo,(lm1(ii,k),ii=1,2),fxck(k,i)
           write(6,562)i,jatom,lmo,(lm1(ii,k),ii=1,2),fxck(k,i)
          endif
552  format(5i4,f12.6,' iatorb, type, lm1, l, m, <Vxc(l,m)>')
562  format(5i4,2f12.6)
         enddo
         read(50,*)jatom,espxc(i)
         write(6,572)i,jatom,espxc(i)
         write(21,572)i,jatom,espxc(i)
572  format(2i4,f12.6,' iatorb, type, xc energy ')
        enddo
       endif                   !EECE ends
      endif
! end of input for LDA+U
! input for OP
      if(nmod.eq.2)then
       read(5,*)nmodop
        if(nmodop.eq.1)then
          write(6,*)' Different Lz for spin up, spin down'
        else
          write(6,*)' Lz is the same for spin up and spin down'
        endif
       do i=1,natorb
         read(5,*)Ncalc(i)
          if(ncalc(i).eq.1)then
           write(21,*)'       Orb. pol. parameters calculated ab initio'
           write(6,*)' Orb. pol. parameters calculated ab initio'
          else
             do li=1,nlorb(i)
               l=lorb(li,i)
               read(5,*)pop(li,i)
               write(21,203)iatom(i),l,pop(li,i)
               write(6,203)iatom(i),l,pop(li,i)
             enddo
          endif
        enddo
        read(5,*)(xms(i),i=1,3)
       endif
! end of input for OP
       write(6,*)' end of OP input'
! input for Bext
      if(nmod.eq.3)then
       read(5,*)Bex
! corresponding energy in Ry (muB*Bext(T) in Ry):
       Bexten=Bex*4.2543812D-6
       write(6,204)Bex,Bexten
       write(21,204)Bex,Bexten
       read(5,*)(xms(i),i=1,3)
      endif
! ********* input data read **************
! begin to write new orbital potential
      write(12,198)nmod,nup,natorb,Bexten,spin(nup)
!....read struct file
      READ(20,5010) TITLE
      READ(20,541) LATTIC,NTYPE
      nato=ntype
      allocate (VSP(NRAD,NATO),EOP(LABC,NATO),zz(nato),mult1(nato))
      allocate (rotloc(3,3,nato),pos(3,48*nato),IATNR(NATO))
      allocate (RO(NATO),DXM(NATO),RMT(NATO), &
                      ANA(NATO),JRJ(NATO),NATOM(NATO))
      ALLOCATE (ee(0:labc,ntype))
      ee(0:labc,1:ntype)=0.0d0
      READ(20,501) TIJUNK
      READ(20,1020) AA,BB,CC,ALPHA(1),ALPHA(2),ALPHA(3)
      IF(ABS(ALPHA(1)).LT.1.D-5) ALPHA(1)=90.0D0
      IF(ABS(ALPHA(2)).LT.1.D-5) ALPHA(2)=90.0D0
      IF(ABS(ALPHA(3)).LT.1.D-5) ALPHA(3)=90.0D0
      INDEX=0
      DO 10 ITY=1,NTYPE
      INDEX=INDEX+1
      READ(20,543) IATNR(ITY),(POS(I,INDEX),I=1,3)
      READ(20,544) NATOM(ITY)
      mult1(ity)=natom(ity)
      DO 12 iINA=2,NATOM(ITY)
      INDEX=INDEX+1
      READ(20,543) IATNR(ITY),(POS(I,INDEX),I=1,3)
   12 CONTINUE
      READ(20,545) ANA(ITY),JRJ(ITY),RO(ITY),RMT(ITY),ZZ(ITY)
      DXM(ITY)=DLOG(RMT(ITY)/RO(ITY))/(JRJ(ITY)-1)
      READ (20,5060) ((ROTLOC(I1,I2,ITY),I1=1,3),I2=1,3)
   10 CONTINUE
      write(6,*)'STRUCT file read'
! end of struct file reading
!
! calculate polar and azimuthal angles of M or Bext 
      if((nmod.eq.2).or.(nmod.eq.3))then
!.....Generate symmop
      CALL LATGEN2(LATTIC) 
! angles of M or Bext in global orthogonal system
        CALL ANGLE(XMS,THETA,PHI)
        th=theta*180./pi
        ph=phi*180./pi
        write(21,212)amod(nmod),(xms(i),i=1,3),th,ph
        write(6,212)amod(nmod),(xms(i),i=1,3),th,ph
      endif
!
! nothing more required for Bext
      if(nmod.eq.3)return
!
! input for LDA+U and OP
      if((nmod.eq.1).or.(nmod.eq.2))then
!GM edit--------------------------------
      iat=1
      DO jat=1,ntype
         READ(14,'(5(f9.5))',iostat=ios) (ee(jj,jat),jj=0,labc)
         IF(ios /= 0) EXIT
         READ(14,*) 
         IF(jat==iatom(iat)) THEN
            DO k=1,nlorb(iat)
              eop(k,iat)=ee(lorb(k,iat),jat)
              IF(eop(k,iat)>150.) eop(k,iat)=eop(k,iat)-200.d0
            ENDDO
            iat=iat+1
         ENDIF
         IF (iat>natorb) EXIT
      ENDDO
!GM edit--------------------------------
!
      mtap=18
! read spherical potential
      READ(mtap,532)
      nnlo=0
      DO 100 ITY=1,NTYPE
       READ(mtap,531)
       READ(mtap,533) (VSP(I,ITY),I=1,JRJ(ITY))
       READ(mtap,531)
  100 CONTINUE
      endif
      write(6,*)' VSP read'
! end of input for LDA+U and OP
!
! read the density matrices  For each atom the density matrices diagonal
! in L were determined in LAPWDM
        do jat=1,natorb
          do ljat=1,nlorb(jat)
           alm(ljat,jat)=0.
          enddo
        enddo
!
      ndm=2
      if(nup.eq.0)ndm=1
      if(nup.eq.2)ndm=1
      do isp=1,ndm   
       ifile=11-isp
        if(nup.eq.2)ifile=11
        do jat=1,natorb
         read(ifile,*)jatom
          if(iatom(jat).ne.jatom)then
           write(6,*)' Conflict in atom indexes: iatom',iatom(jat), &
                  'ne jatom',jatom
           stop
          endif
          do ljat=1,nlorb(jat)
! alx, aly, alz are components of orbital momentum in global system
!   for atom=jat, l=ljat, spin of potential (isp=1) or opposite (isp=2)
           read(ifile,*)ll,alx,aly,alz   
! alm is projection of total orbital momentum on M or Bext direction
           alm(ljat,jat)=alm(ljat,jat)+ &
             sin(theta)*(alx*cos(phi)+aly*sin(phi))+cos(theta)*alz
!
           l=lorb(ljat,jat)
! al(i,ljat,jat) are lx,ly,lz in global system for atom jat, l=ljat
!                and for the spin of potential
            if(isp.eq.1)then
             al(1,ljat,jat)=alx
             al(2,ljat,jat)=aly
             al(3,ljat,jat)=alz
             write(6,205)jatom,l,(al(i,ljat,jat),i=1,3)
            endif
            if(l.ne.ll)then
             write(6,*)' Conflict in atom orb. number: lorb',l,'ne ll',ll
             stop
            endif
            do m=-l,l
             read(ifile,201)(dmat(jat,ljat,m,mp,isp),mp=-l,l)
            enddo
            if(ipr.gt.1)then
             if(isp.eq.1)then
             if(nup.eq.2)then
              write(6,*)'Atom',jatom,' density matrix UPDN block, L=',l
              else
              write(6,*)'Atom',jatom,' density matrix for spin ', &
                          spin(nup),', L=',l
             endif
             else
              write(6,*)'Atom',jatom,' density matrix for spin ', &
                          spin(-nup),', L=',l
             endif
             write(6,*)' Real part'
              do m=-l,l
               write(6,208)(dble(dmat(jat,ljat,m,mp,isp)),mp=-l,l)
              enddo
             write(6,*)' Imaginary part'
              do m=-l,l
               write(6,208)(aimag(dmat(jat,ljat,m,mp,isp)),mp=-l,l)
              enddo
             write(6,*)
             endif
            if(isp.eq.2)write(6,197)jatom,l,alm(ljat,jat)
            enddo
           enddo
          enddo
! end of density matrices input
! for para-calculations copy dmat(...,1) on dmat(....,2)
      if(nup.eq.0)then
        do jat=1,natorb
          do ljat=1,nlorb(jat)
            do m=-l,l
             do mp=-l,l
              dmat(jat,ljat,m,mp,2)=dmat(jat,ljat,m,mp,1)
             enddo
            enddo
          enddo
        enddo
      endif
      deallocate(ee)
      RETURN
197   format(' atom',i3,' L=',i2,' projection of L on M=',f12.6)
198   format(3i3,e14.6,' nmod, nsp, natorb,', &
             ' muB*Bext (Ry), spin ',a4)
199   format &
      (' Calculation of orbital potential for spin block: ',a4)
200   format(' Type of potential:',a45)
201   format(2(2e16.8,2x))
202   format(' Vorb applied to atom',i4,' orbit. numbers',3i4)
203   format(7x,' Atom',i3,' L=',i2, &
             ' Orb. pol. parameter set equal to ',f8.4,' Ry')
204    format(7x,' Bext=',f10.5,' T; muB*Bext=',e12.5,' Ry')
205   format(' Atom',i3,' L=',i2,' spin of potential; Lx, Ly, Lz=' &
              ,3f10.6)
206   format(' Atom',i3, &
            ' L =',i2,' E=',f9.4,' Ry')
207   format(7x,' Atom',i3,' L=',i3,' U=',f7.3,' J=',f7.3,' Ry')
208   format(7f11.5)
209   format(10x,a6)
210   format(/10x,a10/)
211   format(12x,i2,2x,f10.4)
212   format(/1x,a4,' in global crystal system',3f10.5,/, &
          ' angles in global orthogonal system (M,z)=' &
          ,f8.3,' (M,x)=',f8.3,' deg'/)
  501 FORMAT(A10)
  531 FORMAT(/////)
  532 FORMAT(//)
  533 FORMAT((3X,4E19.12))
  546 FORMAT(4I3)
  547 FORMAT(3F10.4)
  542 FORMAT(3f10.7)
  541 FORMAT(a4,23x,i3)
  543 FORMAT(4X,I4,4X,f10.7,3x,f10.7,3x,f10.7)
  544 FORMAT(15x,i2)
  545 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
  549 FORMAT(A5)
  553 FORMAT(2x,13f8.3)
 1020 FORMAT(6F10.7,10X,F10.7)
 5010 FORMAT(A80)
 5060 FORMAT(20X,3F10.8)
!
      END 
