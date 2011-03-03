      SUBROUTINE MAIN1 ()
!
!     AIM (1999-2000) by Javier D. Fuhr and Jorge O. Sofo
!
!     Instituto Balseiro and Centro Atomico Bariloche
!     S. C. de Bariloche - Rio Negro
!     Argentina
!     e-mail: fuhr@cab.cnea.gov.ar
!             sofo@cab.cnea.gov.ar
!
!     Reads case.inaim, interpret it, and execute the
!       corresponding subroutines.
!
!     Small changes by L. D. Marks, January 2006
      use rad
      use sphe
      use atpos
      IMPLICIT NONE

      CHARACTER*4 switch

      logical ortho,deb,debold,ssave
      logical stwo,sthree,sfour,srch,readcrit,readsurf

      real*8 BR1,BR2,br3,br4
      real*8 tau,otau,oiz,invoiz,tpi,a,alpha,beta,gamma,gamma1
      real*8 v,pi,dmax,dstmax,v1,v2,vt1,vt2,dir,vi(3)

      integer iz,iord
      integer lwork,index,jj,iat,ipos,i,j
      integer iatinit,iposinit,npmax,npoints,nstep

      INCLUDE 'param.inc'
      COMMON /DEBUG/ deb
      COMMON /SYM2/ tau(3,NSYM),otau(3,NSYM),oiz(3,3,NSYM), &
         invoiz(3,3,NSYM),iz(3,3,NSYM),iord
!      COMMON /RAD/ RM(NRAD,NATO),JRI(NATO)                       
!      COMMON /ATPOS/ ATP(3,nnpos),NPOS,NSA,NSB,NSC
      COMMON /BRAV/ BR1(3,3),BR2(3,3),br3(3,3),br4(3,3)
      COMMON /LATT/ A(3),alpha,beta,gamma,gamma1,ortho
!      COMMON /SPHE/ POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO), &
!         RNOT(NATO),ndat,IATNR(NDIF),IOP(NDIF),MULT(NATO)
      DIMENSION V(3),v1(3),v2(3),vt1(3),vt2(3)
      integer igmode
      logical trust
      common /ldmopts/igmode,trust

      tpi=2.0*acos(-1.d0)
      pi=acos(-1.d0)
      
      lwork=9

      deb=.false.
      readcrit=.false.
      readsurf=.false.
!
!     default options
!     Use the trust-region mode
      trust=.true.
!     Use conventional Gauss-Legendre quadrature
      igmode=1

!     READ STRUCT AND CLMSUM FILES
!
      call init()


      READ(5,112) switch

      do while(switch.ne.'END ')

!
!     Options
        if(switch .eq. 'OLDA')then
!       Old mode
          igmode=0
          trust=.false.
          write(6,122)'Switching to old mode with midpoint and'
          write(6,122)'quadrature for both theta, phi'
          goto 13
        else if(switch .eq. 'QUAD')then
!       Full theta, phi quadrature
          write(6,122)'Turning on full quadrature for theta,phi'
          igmode=0
          goto 13
        else if(switch .eq. 'LATT')then
!       Lattice mode, equal cos(theta), phi sampling
          write(6,122)'Using lattice sampling for theta, phi'
          igmode=2
          goto 13
        endif

!
!     BUSQUEDA DE PUNTOS CRITICOS
!
        if(switch.eq.'CRIT') then
          read(5,*) index
          read(5,121) switch
          read(5,*) NSA,NSB,NSC
          write(6,120) NSA,NSB,NSC
          call gener(BR2)
          dstmax=min(20.d0,min(nsa*a(1),nsb*a(2),nsc*a(3))) 

          stwo=.false.
          sthree=.false.
          sfour=.false.
          if((switch.eq.'TWO ').or.(switch.eq.'ALL ')) stwo=.true.
          if((switch.eq.'THRE').or.(switch.eq.'ALL ')) sthree=.true.
          if((switch.eq.'FOUR')) sfour=.true.
          if(switch.eq.'ONE ') then
            read(5,*) v(1),v(2),v(3)
            do jj=1,3
              vi(jj)=br4(1,jj)*pos(1,index)+ &
                   br4(2,jj)*pos(2,index)+br4(3,jj)*pos(3,index)
            enddo
            call newcrit(v,vi)
            goto 13
          endif
          if(switch.eq.'LINE') then
             read(5,*) vt1(1),vt1(2),vt1(3)
             read(5,*) vt2(1),vt2(2),vt2(3)
             do jj=1,3
                vi(jj)=br4(1,jj)*pos(1,index)+ &
                     br4(2,jj)*pos(2,index)+br4(3,jj)*pos(3,index)
             enddo
             do i=1,3
                v1(i)=vt1(1)*br1(1,i)+vt1(2)*br1(2,i)+vt1(3)*br1(3,i)
                v2(i)=vt2(1)*br1(1,i)+vt2(2)*br1(2,i)+vt2(3)*br1(3,i)
             enddo
             read(5,*) npoints
             do i=1,npoints
                write(6,*) 'Searching from point',i
                do j=1,3
                   v(j)=v1(j)+(v2(j)-v1(j))*(i-1.)/(npoints-1.)
                enddo
                call newcrit(v,vi)
             enddo
             goto 13
          endif
          if(switch.eq.'AROU') then
             call critaround(index)
          endif
          call critics(index,stwo,sthree,sfour,dstmax)

 120      format('0NSHELL PARAMETERS IN X, Y AND Z-DIRECTION:',3I5)
 121      format(a4)
 122      format(a)
          
          goto 13
        endif

!
!     SUPERFICIE DE FLUJO CERO
!
        if(switch.eq.'SURF') then
          call init_surf(iat,readcrit)
!     deb=.true.
          call surf(iat,readsurf)
          goto 13
        endif

!
!     INTEGRA RHO
!            
        if (switch.eq.'IRHO') then
          call integrho()
          goto 13
        endif

!
!     SOLO UN R DE LA SUPERFICIE DE FLUJO 0
!            
        if (switch.eq.'RSUR') then
          write(6,*) 'SOLO UN R DE LA SUPERFICIE DE FLUJO 0'
          call rsur()
          goto 13
        endif

        if (switch.eq.'FOLL') then
          call init_follow(v,iatinit,iposinit,dir,dmax)
          npmax=17
          srch=.false.
          debold=deb
          ssave=.false.
!     deb=.true.
          call follow(v,npmax,srch,iatinit,iposinit, &
             iat,ipos,dmax,nstep,dir,ssave)
          deb=debold
          goto 13
        endif

        if (switch.eq.'DOIT') then
          call doit()
          goto 13
        endif
        
        if (switch.eq.'READ') then
          call readcs(readcrit,readsurf)
          goto 13
        endif

        if (switch.eq.'DIPO') then
          call dipole()
          goto 13
        endif

!       Skip temperature factor term
        if (switch.eq.'TEMP')goto 13
        write(*,*)'Warning: Option ',switch,' not recognized.'
        write(*,*)'Continuing...'
        write(6,*)'Warning: Option ',switch,' not recognized.'
        write(6,*)'Continuing...'

!        STOP 'OPTION NOT RECOGNIZED.'

 13     read(5,112) switch
      end do

 112  FORMAT(A4)
      END                                                   

