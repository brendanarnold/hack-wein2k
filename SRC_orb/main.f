PROGRAM main    
  use opr
  use ldau
  use exex
  use struct, only : IATNR
  IMPLICIT REAL*8(A-H,O-Z)
  INCLUDE 'param.inc'
  CHARACTER*80    DEFFN,ERRFN,FNAME
  LOGICAL         FL(5)
  character*4     amod(5),chldau(0:5)
  character*5     chmix
  character*11    status,form
  complex*16      Ldzeta(-3:3,-3:3)
  !      complex*16      dmat(nato,3,-3:3,-3:3,2)
  complex*16      dm(-3:3,-3:3,2),ddm(-3:3,-3:3)
  REAL*8                 :: ddm2(-3:3,-3:3)
  complex*16,allocatable ::      vorb(:,:,:,:),vold(:,:,:,:)
  !      complex*16      vorb(nato,3,-3:3,-3:3),vold(nato,3,-3:3,-3:3)
  complex*16      vu(-3:3,-3:3),vu_amf(-3:3,-3:3),vu_FLL(-3:3,-3:3)
  character*4     spin(-1:2)
  !      dimension       NATOM(NATO)
  common/print/   ipr,ipr3j,iprv
  common/updn/    nup,spin
  complex*16 imag
!  common/exex/nexex,nsp,fxck(42,20),fupk(0:6,20),fxc,vxcl(-3:3,-3:3), &
!    fupka(0:6),hybr,nlm1(20),lm1(42,20),espxc(20),exc
      COMMON/FACT/FCT(0:100),nfac   
  common/angles/  theta,phi,thetaloc,philoc
  real*8 ja
  REAL*8 :: C_KUB(0:10,0:10)
  dimension x(3),xloc(3)
  imag=(0.,1.)
  pi=4.d0*atan(1.d0)

  C_KUB=0.0d0
  C_KUB(0,0)=1.d0
  C_KUB(3,2)=1.d0
  C_KUB(4,0)=.5d0*SQRT(7.d0/3.d0)
  C_KUB(4,4)=.5*SQRT(5.d0/3.d0)
  C_KUB(6,0)=.5d0*SQRT(.5d0)
  C_KUB(6,2)=.25d0*SQRT(11.d0)
  C_KUB(6,4)=-.5d0*SQRT(7.d0/2.d0)
  C_KUB(6,6)=-.25d0*SQRT(5.d0)
  C_KUB(7,2)=.5d0*SQRT(13.d0/6.d0)
  C_KUB(7,6)=.5d0*SQRT(11.d0/6.d0)
  C_KUB(8,0)=.125d0*SQRT(33.d0)
  C_KUB(8,4)=.25d0*SQRT(7.d0/3.d0)
  C_KUB(8,8)=.125d0*SQRT(65.d0/3.d0)
  C_KUB(9,2)=.25d0*SQRT(3.d0)
  C_KUB(9,4)=.5d0*SQRT(17.d0/6.d0)
  C_KUB(9,6)=-.25d0*SQRT(13.d0)
  C_KUB(9,8)=-.5d0*SQRT(7.d0/6.d0)
  C_KUB(10,0)=.125*SQRT(65.D0/6.D0)
  C_KUB(10,2)=.125*SQRT(247.D0/6.D0)
  C_KUB(10,4)=-.25*SQRT(11.D0/2.D0)
  C_KUB(10,6)=0.0625d0*SQRT(19.D0/3.D0)
  C_KUB(10,8)=-.125*SQRT(187.D0/6.D0)
  C_KUB(10,10)=-.0625d0*SQRT(85.d0)

  CALL OPNFS(DEFFN,ERRFN,FL)
  CALL ERRFLG(ERRFN,'Error in Vorb')
  nfac=20
  call fac
!  Y=1.0D0                                                           
!  DO I=1,49
!    jJ=2*I-1
!    FCT(jJ)=Y
!    Y=Y*I
!  END DO
 ! find out whether spin is up or down or if calculation is non spin-polarized
  spin(-1)='down'
  spin(1) ='up  '
  spin(2) ='dnup'
  spin(0) ='para'
  nup=-1
  OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=999)
10 CONTINUE
  READ (1,*,END=20,ERR=999) IUNIT,FNAME,STATUS,FORM,IRECL
  if(iunit.eq.12) then
     do i=80,1,-1
        if(fname(i:i).ne.' ') then
           if(fname(i:i).eq.'p') then
              nup=1
           else if(fname(i:i).eq.'b') then
              nup=0
           endif
           goto 10
        endif
     enddo
  endif
  if(iunit.eq.11) nup=2
  goto 10
20 continue
  !
  amod(1)='LDAU'
  amod(2)='magn'
  amod(3)='Bext'
  chldau(0)=' AMF'
  chldau(1)=' FLL'
  chldau(2)=' MFH'
  chldau(3)=' MIX'
  chldau(4)=' ExE'
  !
  call init(chmix)
  allocate(vorb(natorb,3,-3:3,-3:3), vold(natorb,3,-3:3,-3:3))
  write(6,*)' natom',(mult1(iatom(ii)),ii=1,natorb)
  ! find out whether old potential is present (then iold=1)
  iold=0
  read(13,*,end=99)nmodold
  if(nmod.ne.nmodold)then
     write(6,*)' Error: attempt to mix different Vorb',nmod,nmodold
     stop
  endif
  iold=1
99 continue
  if(iold.eq.0)then
     write(6,*)' No old potential found'
  else
     write(6,104) chmix,amix
  endif
  !
  ! dir. cos. of M (nmod=2) or Bext (nmod=3) in global coordinate system
  if((nmod.eq.2).or.(nmod.eq.3)) then
     x(1)=sin(theta)*cos(phi)
     x(2)=sin(theta)*sin(phi)
     x(3)=cos(theta)
  endif
  !
  do i=1,natorb
     iat=iatom(i)
     mult=mult1(iat)
     !       write(12,*)iat,nlorb(i),' atom type, number of L'
     !!       if(iold.eq.1)read(13,*)iatold
     ! calculate orbitally dependent potential
     if((nmod.eq.2).or.(nmod.eq.3)) then
        ! dir. cos of M (nmod=2) or Bext (nmod=3) in local coordinate system
        call rotate(x,rotloc(1,1,iat),xloc)
        thetaloc=ACOS(xloc(3))
        xx=sqrt(xloc(1)**2+xloc(2)**2)
        if (xx.LT.1D-5) then
           philoc=0.D0
        else
           xab=xloc(1)/XX
           philoc=acos(xloc(1)/XX)
           if(abs(xloc(2)).gt.1D-5) &
                philoc=philoc*xloc(2)/abs(xloc(2))
        endif
        thloc=thetaloc*180./pi
        phloc=philoc*180./pi
        write(6,212)amod(nmod),(xloc(k),k=1,3),thloc,phloc
     endif
     !
     do nl=1,nlorb(i)
        !!         if(iold.eq.1)read(13,*)Lold
        L=lorb(nl,i)
        if(nmod.eq.1)then
         if(nldau.eq.4)then   ! Exact exchange
           do ll=0,6,2
            fupka(ll)=fupk(ll,i)
           enddo
           fxc=fxck(1,i)
           exc=espxc(i)
! calculate <m1|Vxc(LL,MM)|m2> for LL>0 (k>1)
           if(nsp.ne.0)then
            do m1=-l,l
             do m2=-l,l
              vxcl(m1,m2)=(0.d0,0.d0)
              do 30 k=2,nlm1(i)
               ll=abs(lm1(1,k))
               mm=lm1(2,k)

               IF (IATNR(IAT) .LT. 0) THEN    !non-cubic case

               sgMM=(-1.d0)**MM
               tclLl=(2.d0*l+1.d0)*sqrt((2.d0*ll+1.d0)/(4.d0*pi))*t3j(l,ll,l,0,0,0)
               t3p=(t3j(L,LL,L,-M1,MM,M2)+t3j(L,LL,L,-M1,-MM,M2))/sqrt(2.d0)
               t3m=(t3j(L,LL,L,-M1,MM,M2)-t3j(L,LL,L,-M1,-MM,M2))/sqrt(2.d0)
               if(mm.eq.0) then         ! m=0
                   vxcL(M1,M2)=vxcL(M1,M2)+fxck(k,i)*(-1.d0)**M1*tcLLL*t3p/sqrt(2.d0)
               else  
                if(sgMM.ge.0.d0)then    !even m
                 if(lm1(1,k).ge.0)then   !pos l
                   vxcL(M1,M2)=vxcL(M1,M2)+fxck(k,i)*(-1.d0)**M1*tcLLL*t3p  
                 else
                   vxcL(M1,M2)=vxcL(M1,M2)+fxck(k,i)*(-1.d0)**M1*tcLLL*t3m/imag
                  endif
                else                      !odd m
                 if(lm1(1,k).ge.0)then    !pos l
                   vxcL(M1,M2)=vxcL(M1,M2)-fxck(k,i)*(-1.d0)**M1*tcLLL*t3m 
                 else
                   vxcL(M1,M2)=vxcL(M1,M2)-fxck(k,i)*(-1.d0)**M1*tcLLL*t3p/imag
                 endif
                endif
               endif

               ELSEIF (IATNR(IAT) .GT. 0) THEN    !cubic case

                 TCLLL = (2D0*L+1D0)*SQRT((2D0*LL+1D0)/(4D0*PI))*T3J(L,LL,L,0,0,0)

                 IF (LL .EQ. 3) THEN

                   T3M1 = (T3J(L,LL,L,-M1,MM,M2) - T3J(L,LL,L,-M1,-MM,M2))/SQRT(2D0)

                   VXCL(M1,M2) = VXCL(M1,M2) + FXCK(K,I)*(-1D0)**M1*TCLLL*T3M1/IMAG
 
                 ELSEIF ((LL .NE. 9) .AND. (MM .EQ. 4)) THEN

                   GOTO 30

                 ELSEIF ((MM .EQ. 6) .OR. (MM .EQ. 8) .OR. (MM .EQ. 10)) THEN

                   GOTO 30

                 ELSEIF ((LL .EQ. 4) .OR. (LL .EQ. 6) .OR. (LL .EQ. 7) .OR. (LL .EQ. 9)) THEN

                   T3P1 = (T3J(L,LL,L,-M1,MM,M2) + T3J(L,LL,L,-M1,-MM,M2))/SQRT(2D0)
                   T3P2 = (T3J(L,LL,L,-M1,MM+4,M2) + T3J(L,LL,L,-M1,-MM-4,M2))/SQRT(2D0)
                   T3M1 = (T3J(L,LL,L,-M1,MM,M2) - T3J(L,LL,L,-M1,-MM,M2))/SQRT(2D0)
                   T3M2 = (T3J(L,LL,L,-M1,MM+4,M2) - T3J(L,LL,L,-M1,-MM-4,M2))/SQRT(2D0)

                   IF (MM .EQ. 0) THEN

                     VXCL(M1,M2) = VXCL(M1,M2) + FXCK(K,I)*(-1D0)**M1*TCLLL* &
                       (C_KUB(LL,MM)*T3P1/SQRT(2D0) + C_KUB(LL,MM+4)*T3P2)

                   ELSEIF (LM1(1,K) .GT. 0) THEN

                     VXCL(M1,M2) = VXCL(M1,M2) + &
                       FXCK(K,I)*(-1D0)**M1*TCLLL*(C_KUB(LL,MM)*T3P1 + C_KUB(LL,MM+4)*T3P2)

                   ELSEIF (LM1(1,K) .LT. 0) THEN

                     VXCL(M1,M2) = VXCL(M1,M2) + &
                       FXCK(K,I)*(-1D0)**M1*TCLLL*(C_KUB(LL,MM)*T3M1 + C_KUB(LL,MM+4)*T3M2)/IMAG

                   ENDIF

                 ELSEIF ((LL .EQ. 8) .OR. (LL .EQ. 10)) THEN

                   T3P1 = (T3J(L,LL,L,-M1,MM,M2) + T3J(L,LL,L,-M1,-MM,M2))/SQRT(2D0)
                   T3P2 = (T3J(L,LL,L,-M1,MM+4,M2) + T3J(L,LL,L,-M1,-MM-4,M2))/SQRT(2D0)
                   T3P3 = (T3J(L,LL,L,-M1,MM+8,M2) + T3J(L,LL,L,-M1,-MM-8,M2))/SQRT(2D0)

                   IF (MM .EQ. 0) THEN

                     VXCL(M1,M2) = VXCL(M1,M2) + FXCK(K,I)*(-1D0)**M1*TCLLL* &
                       (C_KUB(LL,MM)*T3P1/SQRT(2D0) + C_KUB(LL,MM+4)*T3P2 + C_KUB(LL,MM+8)*T3P3)

                   ELSE

                     VXCL(M1,M2) = VXCL(M1,M2) + FXCK(K,I)*(-1D0)**M1*TCLLL* &
                       (C_KUB(LL,MM)*T3P1 + C_KUB(LL,MM+4)*T3P2 + C_KUB(LL,MM+8)*T3P3)

                   ENDIF
                 ENDIF
               ENDIF

30            continue
             enddo
            enddo
           endif
         endif               ! Exact exchange end
           Ua=U(nl,i)
           Ja=J(nl,i)
           do m1=-l,l
              do m2=-l,l
                 ndm=2
                 if(nup.eq.2)ndm=1
                 do is=1,ndm   
                    dm(m1,m2,is)=dmat(i,nl,m1,m2,is)
                 enddo
              enddo
           enddo
           IF(nldau==3.AND.nup/=2) THEN
              ! Igors functional
              av_n=0.0d0
              DO ii=-l,l
                 av_n=av_n+dm(ii,ii,1)
              ENDDO
              av_n=av_n/(2*l+1)
              ddm(-l:l,-l:l)=dm(-l:l,-l:l,1)
              DO ii=-l,l
                 ddm(ii,ii)=dm(ii,ii,1)-av_n
              ENDDO
              ddm2=0.0d0
              do ii=-l,l
                 do jj=-l,l
                    do kk=-l,l
                       ddm2(ii,jj)=ddm2(ii,jj)+ddm(ii,kk)*CONJG(ddm(jj,kk))
                    ENDDO
                 ENDDO
              ENDDO
              xnumerator=0.0d0
              DO ii=-l,l
                 xnumerator=xnumerator+ddm2(ii,ii)
              ENDDO
              alpha=xnumerator/((2*l+1)*av_n*(1.d0-av_n))
              nldau=1
              call vldau(nup,nldau,l,dm,ua,ja,vu_FLL,mult,ipr,eorb_FLL)
              nldau=0
              call vldau(nup,nldau,l,dm,ua,ja,vu_amf,mult,ipr,eorb_amf)
              vu(-3:3,-3:3)=alpha*vu_amf(-3:3,-3:3)+(1.d0-alpha)*vu_FLL(-3:3,-3:3)
              WRITE(21,*) 'Mazin interpolation V=a*Vamf+(1-a)*VFLL'
              WRITE(21,107) alpha,eorb_amf,eorb_FLL
              write(21,106) alpha*eorb_amf+(1.d0-alpha)*eorb_FLL
              nldau=3
           ELSE
              call vldau(nup,nldau,l,dm,ua,ja,vu,mult,ipr,eorb_FLL)
              write(21,106) eorb_FLL
           ENDIF
106        format(':EORB:',f13.6)
107        format(':ALPH:',f13.6,':EAMF:',f13.6,':EFLL:',f13.6)
           do m1=-l,l
              do m2=-l,l
                 vorb(i,nl,m1,m2)=vu(m1,m2)
              enddo
           enddo
        else if((nmod.eq.2).or.(nmod.eq.3))then
           call orbmom(L,thetaloc,philoc,Ldzeta,ipr)
        endif
        if(nmod.eq.2)then
           if(ncalc(i).eq.1)then
              call parop(i,nl,L,OP,ipr)
           else
              OP=pop(nl,i)
           endif
           write(6,666)i,nl,ncalc(i),l,pop(nl,i),OP
666        format(' i,nl,ncalc,l,pop(l,i)',4i3,2f10.5)
           !           write(12,101)L,OP
           ! scalar product (L*M/|M|)
           if(nmodop.eq.1)then
              xL=0.
              do k=1,3
                 xL=xL+x(k)*al(k,nl,i)
              enddo
           else
              xl=alm(nl,i)
           endif
           write(6,*)' xl=Lm',xl
           do m1=-l,l
              do m2=-l,l
                 !              Vorb(m1,m2)=-OP*xL*Ldzeta(m1,m2)
                 Vorb(i,nl,m1,m2)=OP*xL*Ldzeta(m1,m2)
              enddo
           enddo
        endif
        if(nmod.eq.3)then
           !           write(12,*)L,' L'
           ! V=-muB*|Bext|*Ldzeta(operator), dzeta || Bext
           do m1=-l,l
              do m2=-l,l
                 Vorb(i,nl,m1,m2)=-Bexten*Ldzeta(m1,m2)
              enddo
           enddo
        endif
        write(6,*)
        write(6,*)' Atom',iat,' spin ',spin(nup), &
             ' potential real part (Ry)'
        write(21,*)
        write(21,*)'       Atom',iat,' spin ',spin(nup), &
             ' Potential real part (Ry)'
        do m=-l,l
           write(6,102)m,(dble(vorb(i,nl,m,m1)),m1=-l,l)
           write(21,102)m,(dble(vorb(i,nl,m,m1)),m1=-l,l)
        enddo
        write(6,*)
        write(6,*)' Potential imaginary part (Ry)'
        write(21,*)
        write(21,*)'                      Potential imaginary part (Ry)'
        do m=-l,l
           write(6,102)m,(aimag(vorb(i,nl,m,m1)),m1=-l,l)
           write(21,102)m,(aimag(vorb(i,nl,m,m1)),m1=-l,l)
        enddo
     enddo
  enddo
  ! For OP calculate correction to total energy
  if(nmod.eq.2)then
     ! parts not diagonal in spin are neglected
     if(nup.le.1)then
        do i=1,natorb
           do nl=1,nlorb(i)
              l=lorb(nl,i)
              tr=0.
              do m=-l,l
                 do m1=-l,l
                    tr=tr+vorb(i,nl,m,m1)*dmat(i,nl,m1,m,1)
                 enddo
              enddo
              eo=-tr/2.
              eorb=mult*eo
              write(6,105)i,l,eo,eorb
              write(21,106)EORB
           enddo
        enddo
     endif
  endif
  !  read the old potential
  if(iold.eq.1)then
     do i=1,natorb
        read(13,*)iatold
        do nl=1,nlorb(i)
           read(13,*)Lold
           L=lorb(nl,i)
           do m=-l,l
              do m1=-l,l
                 read(13,103) vold(i,nl,m,m1)
              enddo
           enddo
        enddo
     enddo
     ! mix old and new potential (if old is present)
     call vorbmix(vorb,vold,natorb,nlorb,lorb,chmix,amix)
  endif
  ! write potential on file 12
  do i=1,natorb
     iat=iatom(i)
     !       mult=natom(iat)
     write(12,*)iat,nlorb(i),' atom type, number of L'
     do nl=1,nlorb(i)
        L=lorb(nl,i)
        if((nmod.eq.1).and.(nldau.ne.4))then
           Ua=U(nl,i)
           Ja=J(nl,i)
           write(12,100)L,chldau(nldau),Ua,Ja
        else
           write(12,*)L,' L'
        endif
        do m=-l,l
           do m1=-l,l
              write(12,103)vorb(i,nl,m,m1)
           enddo
        enddo
        ! new vorb data written for atom iatom(i), L=lorb
     enddo
  enddo
  ! 
100 format(i4,a4,2f9.4,' L, modus, U, J (Ry)')
101 format(i4,f9.4,' L, orb.pol.par. (Ry)')
102 format(7x,' M=',i3,7f10.5)
103 format(2f18.11)
104 format(a5,' Pratt mixing of old and new potential, coef.=',f8.4)
105 format(' atom',i4,' L=',i2,' Eorb=',f10.6,' Eorb(tot)=',f10.6) 
212 format(/1x,a4,' in local orthogonal system',3f10.5,/, &
       ' angle (M,zloc)=',f8.3,' angle (M,xloc)=',f8.3,' deg'/)
999 continue
  CALL ERRCLR(ERRFN)
  stop ' ORB   END'
end PROGRAM main
