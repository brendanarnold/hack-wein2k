!
      SUBROUTINE maxstp(x,p,n,stpmax,nat)
      use atomgrid
      implicit double precision (a-h, o-z )
      include 'param.inc'
!
!
!   SUBR. SUGEO IS A SET-UP SUBR. FOR THE LATTICE GEOMETRY
!
!   Np  No. OF PARTICLLES IN THE UNIT CELL
!
!   IEA(N)   SEQUENCE  OF  INEQUIVALENT  ATOMS  IN  RR(3,N)
!   RR(3,N)   SL PARTICLE COORD.  -  CARTESIAN
!
!
      parameter ( nmax = 1000 )
      parameter ( nsh = 1 )
      parameter ( eps = 1.d-8 )
!
      integer qparam
      dimension x(n),p(n)
!
      DIMENSION RR(3,nmax),pp(3,nmax),rmta(nmax)
      DIMENSION IAK(nmax)
      DIMENSION GG(3,3),alat(3),alat2(3),gg0(3,3)
!
!
      CHARACTER*67       ERRMSG
      CHARACTER*4        LATTIC
      COMMON  /CHAR/     LATTIC
      SAVE    /CHAR/
      logical ortho,touch
      COMMON /ORTH/   ORTHO
      COMMON /STRUK/  AA,BB,CC,ALPHA(3),PIA(3),VOL
      real*8,allocatable ::       POS2(:,:)
      COMMON /GENER/  BR1(3,3),BR2(3,3)                                 
!------------------------------------------------------------------------
!      time = mclock()/100.d0
!
!      if( .not. ortho )then
!        stpmax=1.d0
!        return
!      endif
!     write(*,*)'  .not. ortho  in maxstp   ==> stop '
!         stop
!      endif
!
      if(.not.allocated(POS2)) allocate (pos2(3,48*nat))
      alat(1) = aa
      alat(2) = bb
      alat(3) = cc
      IF(LATTIC(1:1).EQ.'R') then
         alat2(1)=DSQRT(1.0d0/3.0d0*aa*aa+1.0d0/9.0d0*cc*cc)
         alat2(2)=alat2(1)
         alat2(3)=alat2(1)
         write(6,*)'alat',alat2(1)
         write(6,*)'alat',alat2(2)
         write(6,*)'alat',alat2(3)
      else
         alat2(1)=alat(1)
         alat2(2)=alat(2)
         alat2(3)=alat(3)
      endif
!
!.....max. step < 1/2 cell
      dlambda = 1.d39
      do 10 i=1,n
         if( dabs(p(i)) .gt. eps )then
            ia = mod(i-1,3) + 1
            test = 0.5d0 * alat2(ia) / dabs( p(i) )
            dlambda = dmin1( test, dlambda )
         endif
 10   CONTINUE
!
      tpi = 8.d0 * datan(1.d0)
      call inv(br2,gg0,det)
      do 20 i1=1,3
         do 20 i2 = 1,3
            gg(i1,i2) = gg0(i2,i1) * tpi 
 20   CONTINUE
!     write(*,500)gg
!500   format('  lattice vectors are: '/ 3(3f20.10/) )
!
!     berechnen der neuen atompositionen
      qparam=0
      index=0
      do 30 ia=1,nat
        index=index+1
        index1=index
        iak(index) = ia 
        rmta(index) = rmt(ia)
        do 40 ir=1,3
          qparam=qparam+1
          if( qparam .gt. n )then
             write(*,*)'  qparam  .gt.  n in maxstp ==>  stop '
             write(*,*)'    qparam, n = ',qparam, n
             stop
          endif
          rr(ir,index)= x(qparam) 
          pp(ir,index)= p(qparam) 
          rr(ir,index) =  &
          dmod(  rr(ir,index) + alat2(ir), alat2(ir) )
 40     CONTINUE
        do 50 im=1,mult(ia)-1
          index=index+1
          iak(index) = ia 
          rmta(index) = rmt(ia)
!     rotieren und translatieren (aber umgekehrt)
          do 60 ir =1,3
            xxx=0.0
            pp(ir,index)=0.0
            do 70 ir1=1,3
              xxx=xxx+rotij(ir1,ir,Index)*rr(ir1,INDEX1)
              pp(ir,index)=pp(ir,index) &
              +rotij(ir1,ir,Index)*pp(ir1,INDEX1)
 70         CONTINUE
!ccc old            xxx=xxx-tauij(ir,index)*alat2(ir)
            xxx=xxx+tauij(ir,index)*alat2(ir)
            rr(ir,index)=dmod( xxx + alat2(ir), alat2(ir) )
 60       CONTINUE
 50     CONTINUE
 30   CONTINUE
!
      np = index
!
!------------------------------------------------------------------------
      if( np .gt. nmax )then
         write(*,*)' np > nmax in sugeo ==> stop '
         write(*,*)' np, nmax = ',np,nmax
         stop
      endif
!
      stpmax=1.01
 200  touch=.false.
      stpmax=stpmax-0.01
      index=0
      do 100 ia=1,nat
        do 100 im=1,mult(ia)
        index=index+1
           do 100 ir=1,3
           pos2(ir,index)=(rr(ir,index)+pp(ir,index)*stpmax)/alat2(ir)
 100  CONTINUE
      write(6,*)' in pairdis  ---------- '
      write(6,*)'  rr = '
      write(6,300)((rr(i,j),i=1,3),j=1,np)
      write(6,*)'  gg = '
      write(6,300)gg
      write(6,*)'  pp = '
      write(6,300)((pp(i,j),i=1,3),j=1,np)
      write(6,*)'  pos2 = '
      write(6,300)((pos2(i,j),i=1,3),j=1,index)
 110  format(30i3)
 300  format(3f16.8)
      call nn(nat,touch,POS2,Alat,ALPHA,RMT,V,PIA,VOL,IATNR,MULT,ISPLIT, &
           stpmax)
      write(6,*)'  maxstp: ',stpmax
      if (touch) goto 200
!      call pairdis(rr,np,gg,nsh,pp,rmta,dlambda)
!      stpmax = dlambda
!
!      time = mclock()/100.d0 - time
!     write(*,*)'  time in maxstp: ',time
      RETURN
!
      END
