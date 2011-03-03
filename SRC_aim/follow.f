      SUBROUTINE FOLLOW(v,npmax,srch,iatinit,iposinit,iat,ipos, &
         dmax,nstep,dir,ssave)

!
!     AIM (1999-2000) by Javier D. Fuhr and Jorge O. Sofo
!
!     Instituto Balseiro and Centro Atomico Bariloche
!     S. C. de Bariloche - Rio Negro
!     Argentina
!     e-mail: fuhr@cab.cnea.gov.ar
!             sofo@cab.cnea.gov.ar
!
!     This routine follow the gradient field until it
!       finds an atom, an pseudoatom or rho < 10^(-5) 
!
!       V: initial point in orthogonal coordinates
!       NPMAX: maximum number of division in each step
!       SRCH:  (TRUE/FALSE) check if the line is outside or 
!              inside the atomic surface.
!       IATINIT,IPOSINIT: index of initial atom
!       IAT,IPOS: index of final atom
!       DMAX:  maximum length of each step
!       NSTEP: returns the number of step needed
!     Small changes by L. D. Marks, January 2006
!
      use sphe
      use crit
      use atpos
      implicit none
      INCLUDE 'param.inc'
      real*8 v,dmax,br1,br2,br3,br4,sum1,sum2,sum3
      real*8 cth,th,ph,wcth,wph,rs
      real*8 themin,themax,phimin,phimax,h0,frmin,r0,dr0
      real*8 vt,vt1,vold,rho,grho,hrho
      real*8 r,fac,h,hold,dg,vnorm,deltar,fac2,dist,xy,xyz
      real*8 vcth,vth,vph,dth,dph,rsmed,hmax,rhomin,fac0,dgmin
      real*8 hmult,snull,facmin
      real*8 t1,t2,t3,h0old,dir,len,stpos,stgrho,stngrho
      real*8 stlen
      real*8 tt1,tt2,tt3

      integer npmax,iatinit,iposinit,iat,ipos
      integer index,nth,nph,ibest,ipbest
      integer nstep,ii,jj,np,nit,kk,ith,iph,nsi
      integer nsimax,ist
      
      logical srho,sgrho,shrho,fin,deb,srch,srchold
      logical stemp,stemp2,ldeb,ldebold,ssave,debldm

      COMMON /DEBUG/ deb
      COMMON /FOLL/ stpos(3,MAXFOLL),stgrho(3,MAXFOLL), &
           stngrho(MAXFOLL),stlen(MAXFOLL),ldeb
      COMMON /BRAV/ br1(3,3),br2(3,3),br3(3,3),br4(3,3)
      COMMON /SRF/ themin,themax,phimin,phimax,h0,frmin,r0,dr0, &
         index,nsimax
      COMMON /INTEG/ cth(NGAUSS),th(NGAUSS),ph(NGAUSS),wcth(NGAUSS), &
         wph(NGAUSS),rs(NGAUSS,NGAUSS),nth,nph
!      data hmax/3.0e0/,rhomin/1.e-5/,fac0/2.1/,dgmin/1.e-5/,hmult/2.0e0/
!      data snull/1e-6/,facmin/5.e-2/
!      data dmin/1.e-3/
!      data fac0/2.1/
      parameter (hmax=3.D0, rhomin=1.D-5, dgmin=1.D-5, fac0=2.1D0 )
      parameter (hmult=1.75D0, snull=1.D-6, facmin=5.D-2)

      DIMENSION V(3),VT(3),vt1(3),vold(3)
      dimension grho(3),hrho(3,3)
      integer igmode
      logical trust
      common /ldmopts/igmode,trust
!
      debldm=.false.
      fin=.false.
      srho=.true.
      sgrho=.true.
      shrho=.false.

      srchold=srch
      ldebold=ldeb
      h0old=h0
      len=0.0D0

      if (deb) then
        ldeb=.true.
      endif

      if(deb)write (6,101) (v(jj),jj=1,3)
 101  format(':POSINIT ',3E14.6)

!     Get the initial density
      call vgh_rho(v,rho,grho,hrho,srho,sgrho,shrho,r,iat,ipos,ist,ibest,ipbest)

      if(iat.ne.0) then
!        Are we inside an RMT ?
         if(r.lt.(frmin*rmt(abs(iatnr(iat))))) then
            fin=.true.
!            write(6,*) 'r < frmin*rmt ',r,'<', &
!                 frmin*rmt(abs(iatnr(iat))),' iat=',iat,' ipos=',ipos
            goto 24
         end if
      endif
      
      if(rho.lt.rhomin) then
        fin=.true.
        write(6,*) 'CHARGE LT rhomin ',rho,' < ',rhomin
        if((rho.lt.0.d0).and.(abs(rho).gt.rhomin)) stop 'ERROR RHO < 0 !!!'
      end if
      
 24   continue

      fac=fac0
      h=hmax

      call cputim(t1)
      nstep=0
      nsi=0
!      write(6,*)'Initial step ',h
      do while(.not.fin)
        hold=h
        dg=vnorm(grho)

        call cputim(t3)
        t2=t3-t1

        if(dg.lt.dgmin) then
          write(6,*) 'gradient < dgmin ',dg,' < ',dgmin
          fin=.true.
          iat=0
          ipos=0
          if(npc.gt.0)  call newcrit(v,pos(1,index))
 204      format(':PC',2I5,3F12.8,3E12.4,I4,2E12.4)
          goto 25
        end if
        
        if(t2.gt.2400.0) then
          write(6,*) 'TIME EXCEEDED 40 min IN FOLLOW'
          write(6,*) 'h0 =',h0,'  h =',h,'  h0old =',h0old,'  dg =',dg
          write(6,*) 'fac =',fac
          stop 'TIME EXCEEDED 40 min IN FOLLOW'
        endif

        h=h0/dg*0.5D0
        if (ldeb.or.deb) write (6,*) 'h= ',h,' h0= ',h0,' dg= ',dg
        if(h.gt.hmax) h=hmax

        h=h*fac
        if(h.gt.(hold*hmult)) h=hold*hmult

        do ii=1,3
          vold(ii)=v(ii)
        enddo

        nit=0
        hold=h
        if(deb)write(6,66)'Revised step ',hold,v
66      format(a,D12.4,3F15.9)
        call cputim(tt1)
 15     if(debldm)write(6,*)':Test ',dmax
        call onestep(v,rho,grho,h,np,npmax,deltar,dir)
        nit=nit+1
!        sum1=grho(1)*grho(1)+grho(2)*grho(2)+grho(3)*grho(3)
!        sum2=v(1)*v(1)+v(2)*v(2)+v(3)*v(3)
!        sum3=v(1)*grho(1)+v(2)*grho(2)+v(3)*grho(3)
!        sum1=sum3/sqrt(sum1*sum2)
!        write(6,*)':Angle ',sum1
        if(((np.gt.npmax).or.(deltar.gt.dmax)).and. &
             (deltar.gt.dmin)) then
!          nit=nit+1
          do ii=1,3
            v(ii)=vold(ii)
          enddo
          h=h/3.D0
          if(debldm)write(6,66)'Reducing step to ',h,hold
          goto 15
        endif
        call cputim(tt3)
        tt2=tt3-tt1
        nstep=nstep+1

        if(nstep.gt.maxfoll) then
           write (6,*) ':NSTEP gt maxfoll ',nstep,maxfoll
           stop 'More than maxfoll steps, sorry'
        endif 
        len=len+deltar
        if(debldm)write(6,*)':LDM ',h,deltar,nit
        if (ldeb.or.deb) write (6,*) 'h= ',h
        
!       Step increase, with reasonable bounds
        fac2=h/hold
        if(fac2.ge.0.99) then
           if (fac.ge.0.99) then
              fac=fac*1.2D0
           else
              fac=1.1D0
           endif
        else
          if (fac2.ge.facmin) then
            fac=fac2
          else
            fac=facmin
          endif
        endif

        call vgh_rho(v,rho,grho,hrho,srho,sgrho,shrho,r,iat,ipos,ist,ibest,ipbest)

        if (deb.or.ldeb) then
          do jj=1,3
            vt1(jj)=br4(1,jj)*v(1)+ br4(2,jj)*v(2)+br4(3,jj)*v(3)
          enddo
          write (6,*) ':POS ',len,v
          write (6,*) ':TIME ',len,nit,tt2
          write (6,*) ':RHO',len,rho
          write (6,*) ':RBPOS ',vt1
          write (6,*) ':GRAD ',len,grho,vnorm(grho)
          write (6,*) ':SP ',r,iat,ipos
        endif

        if(ssave) then
           do jj=1,3
              stpos(jj,nstep)=v(jj)
              stgrho(jj,nstep)=grho(jj)
           enddo
           stngrho(nstep)=vnorm(grho)
        endif

        if(iat.ne.0) then
           if(r.lt.(frmin*rmt(abs(iatnr(iat))))) then
              fin=.true.
              goto 25
           end if
        endif
      
        if(rho.lt.rhomin) then
          fin=.true.
          if((rho.lt.0.d0).and.(abs(rho).gt.rhomin)) stop 'ERROR RHO < 0 !!!'
          iat=0
          ipos=0
          goto 25
        end if
      
        if(npcm3.gt.0) then
          do jj=1,npc
            if (icpc(jj).eq.(-3)) then
              dist=0.d0
              do kk=1,3
                dist=dist+(pc(kk,jj)-v(kk))*(pc(kk,jj)-v(kk))
              enddo
!              dist=sqrt(dist)
              if (dist.lt.(frmin*frmin*0.01)) then
                iat=0
                ipos=0
                fin=.true.
                write (6,*) 'We are inside a pseudo-atomo'
                goto 25
              endif
            endif
          enddo
        endif

        nsi=nsi+1

        if(srch.and.(nsi.ge.nsimax)) then
          nsi=0
          ith=0
          iph=0
          do ii=1,3
            vt(ii)=v(ii)-pos(ii,iatinit)
          enddo
          xy=vt(1)*vt(1)+vt(2)*vt(2)
          xyz=xy+vt(3)*vt(3)
          xyz=sqrt(xyz)
          if (xy.lt.snull) then
            vcth=1.d0
            if(vt(3).lt.0.d0) vcth=-vcth
            vph=0.d0
          else
            vcth=vt(3)/xyz
            vph=atan2(vt(2),vt(1))
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
          endif
 13       continue
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
          endif
 14       continue
          stemp=(iph.gt.0).and.(iph.lt.nph)
          stemp=stemp.and.(ith.gt.0).and.(ith.lt.nth)
          
          if(stemp) then
            stemp2=rs(ith,iph).gt.0.d0
            stemp2=stemp2.and.(rs(ith+1,iph).gt.0.d0)
            stemp2=stemp2.and.(rs(ith+1,iph+1).gt.0.d0)
            stemp2=stemp2.and.(rs(ith,iph+1).gt.0.d0)
            if(stemp2) then
              rsmed=(rs(ith,iph)*(th(ith+1)-vth)+ &
                   rs(ith+1,iph)*(vth-th(ith)))*(ph(iph+1)-vph)
              rsmed=rsmed+(rs(ith+1,iph+1)*(vth-th(ith))+ &
                   rs(ith,iph+1)*(th(ith+1)-vth))*(vph-ph(iph))
              rsmed=rsmed/((th(ith+1)-th(ith))*(ph(iph+1)-ph(iph)))
              if(rsmed.gt.xyz) then
!               write (6,*) 'We are inside the surface'
                iat=iatinit
                ipos=iposinit
              else
!               write (6,*) 'We are outside the surface'
                iat=0
                ipos=0
              endif
              fin=.true.
              goto 25
            endif
          endif
        endif
        

 25     continue
      end do

      if(deb)write (6,102) (v(jj),jj=1,3)
 102  format(':POSOUT ',5D14.6)

      srch=srchold
      ldeb=ldebold
      h0=h0old
      
      return
      end

      subroutine init_follow(v,iat,ipos,dir,dmax)
      use atpos
      implicit none
      INCLUDE 'param.inc'
      
      real*8 v,br1,br2,br3,br4
      real*8 dir
      real*8 themin,themax,phimin,phimax,h0,frmin,r0,dr0
      real*8 stpos,stgrho,stngrho,stlen,dmax
      
      integer iat,ipos
      integer index,nsimax

      logical deb,ldeb

      COMMON /debug/ deb
      COMMON /FOLL/ stpos(3,MAXFOLL),stgrho(3,MAXFOLL), &
           stngrho(MAXFOLL),stlen(MAXFOLL),ldeb
      COMMON /BRAV/ BR1(3,3),BR2(3,3),br3(3,3),br4(3,3) 
      COMMON /SRF/ themin,themax,phimin,phimax,h0,frmin,r0,dr0, &
         index,nsimax
!      COMMON /ATPOS/ ATP(3,nnpos),NPOS,NSA,NSB,NSC

      DIMENSION v(3)
      
      read(5,*) index
      write(6,*) 'ATOM ',index
      read(5,*) h0,frmin
      write (6,*) 'H ',h0,'  frmin=',frmin
      read(5,*) NSA,NSB,NSC
      write(6,120) NSA,NSB,NSC
 120  format('0NSHELL PARAMETERS IN X, Y AND Z-DIRECTION:' &
         ,3I5)         
      call gener(BR2)

      nsimax=0
      dmax=h0

      read(5,*) v(1),v(2),v(3)
      read(5,*) dir

      ldeb=.true.
      return
      end

