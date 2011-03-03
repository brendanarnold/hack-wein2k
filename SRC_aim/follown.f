      SUBROUTINE FOLLOWN(v,npmax,srch,iatinit,iposinit,iat,ipos, &
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
!     finds an atom, an pseudoatom or rho < 10^(-5) 
!
!     Rewritten, January 2006 by L. D. Marks, L-marks@northwestern.edu
!     All this routine currently does is check for termination conditions
!     and call slowstep which does the main work.
!       
!     V: initial point in orthogonal coordinates
!     IATINIT,IPOSINIT: index of initial atom
!     IAT,IPOS: index of final atom
!     h0:  (passed by common) maximum length of each step
!     NSTEP: returns the number of step needed
!
      use sphe
      use crit
      use atpos
      implicit none
      INCLUDE 'param.inc'
      real*8 v,dmax,br1,br2,br3,br4,step(3)
      real*8 cth,th,ph,wcth,wph,rs,TrustR
      real*8 themin,themax,phimin,phimax,h0,frmin,r0,dr0
      real*8 vt,vt1,vold,rho,grho,hrho
      real*8 r,fac,h,hold,dg,vnorm,deltar,fac2,dist,xy,xyz
      real*8 vcth,vth,vph,dth,dph,rsmed,hmax,rhomin,fac0,dgmin
      real*8 hmult,snull,facmin
      real*8 t1,t2,t3,h0old,dir,len,stpos,stgrho,stngrho
      real*8 stlen
      real*8 tt1,tt2,tt3

      integer npmax,iatinit,iposinit,iat,ipos
      integer index,nth,nph,ntries,ncycles
      integer nstep,ii,jj,np,kk,ith,iph
      integer nsimax,ist,ibest,ipbest
      
      logical srho,sgrho,shrho,fin,deb,srch,srchold
      logical keepgoing,ssave,debldm,ldeb

      COMMON /DEBUG/ deb
      COMMON /FOLL/ stpos(3,MAXFOLL),stgrho(3,MAXFOLL), &
           stngrho(MAXFOLL),stlen(MAXFOLL),ldeb
      COMMON /BRAV/ br1(3,3),br2(3,3),br3(3,3),br4(3,3)
      COMMON /SRF/ themin,themax,phimin,phimax,h0,frmin,r0,dr0, &
         index,nsimax
      COMMON /INTEG/ cth(NGAUSS),th(NGAUSS),ph(NGAUSS),wcth(NGAUSS), &
         wph(NGAUSS),rs(NGAUSS,NGAUSS),nth,nph
      parameter (hmax=3.D0, rhomin=1.D-5, dgmin=1.D-6, fac0=2.1D0 )
      parameter (hmult=1.75D0, snull=1.D-6, facmin=5.D-2)

      DIMENSION V(3)
      dimension grho(3),hrho(3,3)

      debldm=.false.
      keepgoing=.true.
      len=0.0D0

      if(deb)write (6,101) (v(jj),jj=1,3)
 101  format(':POSINIT ',3E14.6)

!     Get the initial density, gradient and Hessian
      srho=.true.
      sgrho=.true.
      shrho=.true.
      call vgh_rho(v,rho,grho,hrho,srho,sgrho,shrho,r,iat,ipos,ist,ibest,ipbest)

      if(iat.ne.0) then
!        Are we inside an RMT ? If so, end the search
         if(r.lt.(frmin*rmt(abs(iatnr(iat))))) then
            goto 200
         end if
      endif
!
!       How about ibest, the closest atom (frmin could be > 1)
!       There is something wrong with the logic of ibest/ipbest so
!       this is commented out
!        if (r .lt. (frmin*rmt(abs(iatnr(ibest))))) then
 !             iat=ibest
  !            ipos=ipbest
   !           goto 200
    !    endif
!!!      
      if(rho.lt.rhomin) then
!     Essentially zero density, end the search
        write(6,*) 'CHARGE LT rhomin ',rho,' < ',rhomin
        if((rho.lt.0.d0).and.(abs(rho).gt.rhomin)) stop 'ERROR RHO < 0 !!!'
        goto 200
      end if
      
 24   continue

      call cputim(t1)
      nstep=0
!
!     Infinite loop, terminated by jumping out
      do while(keepgoing)

        dg=vnorm(grho)
        call cputim(t3)
        t2=t3-t1

        if(dg.lt.dgmin) then
!         Critical point found
          write(6,*) 'gradient < dgmin ',dg,' < ',dgmin
          iat=0
          ipos=0
          if(npc.gt.0)  call newcrit(v,pos(1,index))
 204      format(':PC',2I5,3F12.8,3E12.4,I4,2E12.4)
          goto 200
        end if
        
        if(t2.gt.2400.0) then
!         Taken too long; it is very unlikely that this will ever occur
          write(6,*) 'TIME EXCEEDED 40 min IN FOLLOW'
          write(6,*) 'h0 =',h0,'  h =',h,'  h0old =',h0old,'  dg =',dg
          write(6,*) 'fac =',fac
          stop 'TIME EXCEEDED 40 min IN FOLLOW'
        endif

66      format(a,D12.4,3F15.9)
 15     continue
!       Slowly increase the trust region radiues TrustR
!       This is the maximum distance the algorithm can move before
!       getting a new gradient/hessian
!       The idea is to move slowly at first, then go faster (but never
!       too fast). There is a trap inside stepper as well.
!       Note: the control here could also be done using len, the total
!       distance moved, perhaps
!       TrustR=h0*(1.D0-0.9D0*exp(-len*10.D0) )
!       or something similar
!
        TrustR=h0*(1.D0-0.9D0*exp(-nstep*nstep*0.01D0) )
!
!       Walk along up to TrustR as a maximum distance
        call  stepper (grho,hrho,TrustR,step)
!       Update the position
        do jj=1,3
                v(jj)=v(jj)+step(jj)
        enddo
!       Update iteration counter
        nstep=nstep+1
!       Distance moved
        deltar=step(1)*step(1)+step(2)*step(2)+step(3)*step(3)
        if(deltar .gt. 1D-16)deltar=sqrt(deltar)
!
!       Too many steps, unlikely to occur (please report)
        if(nstep.gt.maxfoll) then
           write (6,*) ':NSTEP gt maxfoll ',nstep,maxfoll
           stop 'More than maxfoll steps, sorry: please report this'
        endif 
!       Total distance moved
        len=len+deltar
!
!       New point information, gradient and density
        call vgh_rho(v,rho,grho,hrho,srho,sgrho,shrho,r,iat,ipos,ist,ibest,ipbest)

        if(iat.ne.0) then
!          Within an atom
           if(r.lt.(frmin*rmt(abs(iatnr(iat))))) then
              goto 200
           end if
        endif
!
!       How about ibest, the closest atom
!       This is not quite right, commented out at present
!        if (r .lt. (frmin*rmt(abs(iatnr(ibest))))) then
 !             iat=ibest
  !            ipos=ipbest
   !           goto 200
    !    endif
!!!      
        if(rho.lt.rhomin) then
!         Essentially zero density, jump out
          if((rho.lt.0.d0).and.(abs(rho).gt.rhomin)) stop 'ERROR RHO < 0 !!!'
          iat=0
          ipos=0
          goto 200
        end if
      
!       Pseudoatom present?
        if(npcm3.gt.0) then
          do jj=1,npc
            if (icpc(jj).eq.(-3)) then
              dist=0.d0
              do kk=1,3
                dist=dist+(pc(kk,jj)-v(kk))*(pc(kk,jj)-v(kk))
              enddo
              if (dist.lt.(frmin*frmin*0.01)) then
                iat=0
                ipos=0
                write (6,*) 'We are inside a pseudo-atomo'
                goto 200
              endif
            endif
          enddo
        endif
      end do

 200  continue
      if(deb)write (6,102) (v(jj),jj=1,3)
 102  format(':POSOUT ',5D14.6)

      return
      end
!
        subroutine stepper (G,H,TrustRin,step)
!
!       L. D. Marks, L-marks@northwestern.edu, January 2006
!       This code walks along the gradient line using existing 1st & 2nd derivative
!       information. It is not "pretty", but since this is not a slow step (compared
!       to function evaluations) it does not need to be.
!       A trap is included to ensure that it does not move too far with a reasonable
!       guess at the size of the 3rd derivative error using the largest eigenvalue
!       of the Hessian
!
        implicit real*8 (a-h,o-z)
        dimension H(3,3),G(3),step(3),sold(3),GW(3),vnudge(3)
!       Number of subdivisions, and total number of steps
        parameter (ntries=30, ncycles=45)
        logical usemethod
!       Limit such that guestimated error from 3rd derivative is < thirdlim
!       Set for a 5% error
        parameter (thirdlim=0.05D0)
!
        TrustR=TrustRin
        gg=vnorm(G)
!       Calculate the largest absolute eigenvalue
!       This seems to be fast enough
        call MaxEigen(H,EigenMax,usemethod)
!
!       The third-order term is of the order of H**2/G
!       The error is going to be about 1/6 d**3 H**2/G
!       Taking the change as H/2 d**2
!       Fractional error is ~ 1/3 d*H/G < thirdlim
!       i.e. d < 3*thirdlim*G/H
!
        if(usemethod)then
                TrustLim=3.D0*thirdlim*gg/EigenMax
!               Limit the step ?
                TrustLim=max(0.01D0,TrustLim)
                TrustR=min(TrustR,TrustLim)
        endif
!
!       Estimate the increment to take
        dt=TrustR/gg/ntries
        dmax=TrustR*TrustR
!       Initialize
        do j=1,3
                step(j)=0
                GW(j)=G(j)
        enddo
!
!       Iteration, up to ncycles
        do i=1,ncycles
                total=0.D0
!               Gradient step
                do j=1,3
                        sold(j)=step(j)
                        step(j)=step(j)+dt*GW(j)
                        total=total+step(j)*step(j)
                enddo
!       Have we moved far enough ?
                if(total .gt. dmax)then
!       Yes, use previous (undershoot)
                        do j=1,3
                                step(j)=sold(j)
                        enddo
                        return
                endif
!
!       Keep going
!       Update the derivative
                do j=1,3
                   GW(j)=g(j)+step(1)*H(1,j)+step(2)*H(2,j)+step(3)*H(3,j)
                enddo
!       Derivative zero ?
                gg=vnorm(GW)
                if(gg .lt. 1d-6)return
!
!       Keep going
        enddo
        return
        end
!
        subroutine MaxEigen(H,EigenMax,usemethod)
!       Returns largest absolute eigenvalue of the Hessian using Lapack
!       L. D. Marks, L-marks@northwestern.edu, January 2006
!
!       If something goes wrong, it returns usemethod as .false.
!
        implicit real*8 (a-h,o-z)
        dimension H(3,3)
        dimension Work(20),Eigen(3),Hwork(3,3)
        logical usemethod
        usemethod=.false.
!
!       Copy H
        do j=1,3
           do k=1,3
                Hwork(k,j)=H(k,j)
           enddo
        enddo
!       Get all the eigenvalues and vectors
        call DSYEV( 'N', 'U', 3, Hwork, 3, Eigen,Work, 20 , INFO )
!
!       Did it work
        if(info .ne.0)return
!
        EigenMax=abs(Eigen(1))
        do j=2,3
                t=abs(Eigen(j))
                if(t .gt. EigenMax)EigenMax=t
        enddo
        usemethod=.true.
        return
        end

