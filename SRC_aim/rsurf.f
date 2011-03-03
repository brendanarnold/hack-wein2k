      subroutine rsurf(rr,grho,theta,phi,rr0,iatinit,npmax,srch)
!
!     AIM (1999-2000) by Javier D. Fuhr and Jorge O. Sofo
!
!     Instituto Balseiro and Centro Atomico Bariloche
!     S. C. de Bariloche - Rio Negro
!     Argentina
!     e-mail: fuhr@cab.cnea.gov.ar
!             sofo@cab.cnea.gov.ar
!
!     This routine finds the radius RR to the atomic surface 
!       corresponding to atom IATINIT, for
!       a given direction (THETA,PHI), starting from RR0.
!     
!     If SRCH = TRUE, it use the already calculated surface to
!       determine if in each step the point is inside or outside
!       the atomic surface.
!
!     Small changes by L. D. Marks, January 2006
!     Mainly to use new follow algorithm 
!
!     The algorithm is:
!       1) At some point, find the atom/pseudoatom at which at which
!       the flux lines (gradient) from this point ends (follow/follown)
!       2) Use a bisection approach to find the largest radial value
!       where the flux terminates at the target atom
!
      use crit
      use sphe
      use atpos
      implicit none
      INCLUDE 'param.inc'

      real*8 rr,grho,theta,phi,rr0
      real*8 themin,themax,phimin,phimax,h0
      real*8 frmin,r0,dr0,cth,th,ph,wcth,wph,rs,v
      real*8 rr1,rr2,drr,dr,dir
      real*8 fac,dmax,t2,t1,ct,st,cf,sf,r
      real*8 stpos,stgrho,stngrho,stlen

      integer iatinit,npmax,nsimax,ibest,ipbest
      integer index,nth,nph,iat,ipos,iposinit,i,ii,nstep
      
      logical inter,srho,sgrho,shrho,cross,deb,srchold
      logical in,in1,in2,low,srch,debold,ldeb,ssave

      COMMON /DEBUG/ deb
      COMMON /FOLL/ stpos(3,MAXFOLL),stgrho(3,MAXFOLL), &
           stngrho(3,MAXFOLL),stlen(MAXFOLL),ldeb
!      COMMON /ATPOS/ ATP(3,nnpos),NPOS,NSA,NSB,NSC
!      COMMON /SPHE/ POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO), &
!         RNOT(NATO),ndat,IATNR(NDIF),IOP(NDIF),MULT(NATO)
!      COMMON /CRIT/ pc(3,NNPOS*NDIF),evpc(3,NNPOS*NDIF), &
!         zpc(3,3,NNPOS*NDIF),pcrb(3,NNPOS*NDIF),npc,icpc(NNPOS*NDIF), &
!         npcm3
      COMMON /SRF/ themin,themax,phimin,phimax,h0,frmin,r0,dr0, &
         index,nsimax
      COMMON /INTEG/ cth(NGAUSS),th(NGAUSS),ph(NGAUSS),wcth(NGAUSS), &
         wph(NGAUSS),rs(NGAUSS,NGAUSS),nth,nph
      DIMENSION V(3),grho(3)
      integer igmode
      logical trust
      common /ldmopts/igmode,trust
      
      dir=1.0D0
      srchold=srch
      debold=deb
      ssave=.false.

 17   ct=cos(theta)
      st=sin(theta)
      cf=cos(phi)
      sf=sin(phi)

      rr=rr0
      rr1=rr
      rr2=rr
      drr=1.0
      dr=dr0

      v(1)=pos(1,index)
      v(2)=pos(2,index)
      v(3)=pos(3,index)
      call search_at(v,r,inter,iat,ipos,ibest,ipbest)
      iposinit=ipos
!      write(6,103) iatinit,iposinit
 103  format('ATOM iat=',i4,' ipos=',i4)
      i=0

      cross=.false.
      sgrho=.true.
      srho=.true.
      shrho=.false.

      in=.true.
      low=.false.

      fac=1.1D0

      dmax=h0*2.d0

      in1=.true.
      in2=in1

      do while((drr.gt.drrmin).or.(i.lt.2))
        call cputim(t1)
        i=i+1
        v(1)=pos(1,index)+rr*st*cf 
        v(2)=pos(2,index)+rr*st*sf 
        v(3)=pos(3,index)+rr*ct 
        ldeb=.false.
        if(trust)then
        call follown(v,npmax,srch,iatinit,iposinit,iat,ipos,dmax, &
             nstep,dir,ssave)
        else
        call follow(v,npmax,srch,iatinit,iposinit,iat,ipos,dmax, &
             nstep,dir,ssave)
        endif
        call cputim(t2)
        t2=t2-t1
        if(trust)then
        write(6,1011) i,rr,iat,ipos,t2,nstep
        else
        write(6,101) i,rr,iat,ipos,t2,nstep
        endif
 101    FORMAT(':STEP ',i4,' r=',f12.8,' iat=',i4,' ipos=',i4, &
           ' time(sec)=',F10.5,' nstep=',i4)
 1011   FORMAT(':STEPN ',i4,' r=',f12.8,' iat=',i4,' ipos=',i4, &
           ' time(sec)=',F10.5,' nstep=',i4)
        if(deb) write(6,102) (ATP(ii,ipos),ii=1,3)
 102    FORMAT('DE LA CELDA ',3f10.5)

        if((iat.eq.iatinit).and.(ipos.eq.iposinit)) then
          in=.true.
        else
          in=.false.
        endif

!     
!     DETERMINA EL NUEVO RADIO
!     

        if ((i.eq.1).or.((in1.eqv.in).and.(.not.cross))) then
          if (in) then
            rr2=rr1
            rr1=rr
            rr=rr+dr
          else
            rr2=rr1
            rr1=rr
            rr=rr-dr
          end if
          if((i.gt.2).and.(dr.lt.(0.6))) then
            dr=dr*1.2
            if(deb) write(6,*) ':DR ',dr
          endif
        else
          if (.not.cross) then
            cross=.true.
            rr2=rr1
          else
            if (in2) then
              if (in) then
                rr2=rr1
              else
                in1=in2
              end if
            else
              if (in) then
                in1=in2
              else
                rr2=rr1
              end if  
            end if
          end if
          rr1=rr
          rr=(rr2+rr1)*0.5D0 !/2.0
        end if

        in2=in1
        in1=in
!        drr=abs(rr2-rr1)/rr
!       Changed to an absolute position tolerance
        drr=abs(rr2-rr1)
        if(deb) write(6,*) ':DRR ',i,rr2,rr1,drr
      end do

 33   srch=srchold
      deb=debold

      return
      end

