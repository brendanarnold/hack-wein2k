      SUBROUTINE DOIT()
!
!     AIM (1999-2000) by Javier D. Fuhr and Jorge O. Sofo
!
!     Instituto Balseiro and Centro Atomico Bariloche
!     S. C. de Bariloche - Rio Negro
!     Argentina
!     e-mail: fuhr@cab.cnea.gov.ar
!             sofo@cab.cnea.gov.ar
!
!     This routine is just for testing
!
      use sphe
      use rad
      use blank
      use atpos
      implicit none
      include 'param.inc'
      
      real*8 a,alpha,beta,gamma,gamma1
      real*8 br1,br2,br3,br4
      real*8 v,r,rho,grho,hrho,theta,phi,cth,sth,cph,sph,rr
      real*8 rhor2,drhor2,rhopp,yp1,ypn

      integer jrj,ir,ilm
      integer iat,ist,i,jatom,index,ipos

      logical srho,sgrho,shrho,inter,ortho,sy,syp,sypp

!      COMMON /SPHE/ POS(3,NDIF),RMT(NATO),ROTLOC(3,3,NATO),DX(NATO), &
!           RNOT(NATO),ndat,IATNR(NDIF),IOP(NDIF),MULT(NATO)
!      COMMON /ATPOS/ ATP(3,nnpos),NPOS,NSA,NSB,NSC
      COMMON /LATT/ A(3),alpha,beta,gamma,gamma1,ortho
      COMMON /BRAV/ BR1(3,3),BR2(3,3),br3(3,3),br4(3,3)
!      COMMON /RAD/ RM(NRAD,NATO),JRI(NATO)                       
!      COMMON // CLM(NRAD,NCOM,NATO),CLM2(NRAD,NCOM,NATO), &
!         LM(2,NCOM,NATO),LMMX(NATO)      

      DIMENSION v(3),grho(3),hrho(3,3)


      jatom=1
      ilm=1

!$$$      yp1=(clm(2,1,jatom)-clm(1,1,jatom))/(rm(2,jatom)-rm(1,jatom))
!$$$      ypn=(clm(jri(jatom),1,jatom)-clm(jri(jatom)-1,1,jatom))/
!$$$     $     (rm(jri(jatom),jatom)-rm(jri(jatom)-1,jatom))
!$$$
!$$$
!         write(6,*) 'yp1=',yp1,' ypn=',ypn
         write(6,*) 'rmlast=',rm(jri(jatom),jatom)

!          call spline(rm(1,jatom),clm(1,1,jatom),jri(jatom),yp1,ypn,
!     $     clm2(1,1,jatom))

      sy=.true.
      syp=.true.
      sypp=.true.
   
!      open(30,FILE='rhoB',STATUS='OLD')
  
      do i=1,jri(jatom)
         write(6,'('':RCLM'',2E18.8)') rm(i,jatom),clm(i,ilm,jatom)
!         write(30,4999) rm(i,jatom),clm(i,ilm,jatom)   
      enddo
       
!      close(30)
!      open(31,FILE='interB',STATUS='OLD') 
      do i=1,1000
         r=rmt(jatom)*i/1000.
         ir=1+log(r/rnot(jatom))/dx(jatom)
         ir=0
         if(ir.gt.jri(jatom)) stop 'ERROR'
         call splint(rm(1,jatom),clm(1,ilm,jatom), &
              clm2(1,ilm,jatom), &
              jri(jatom),r,ir,rhor2,drhor2,rhopp,sy,syp,sypp)
         write(6,'('':RR'',4E18.8)') r,rhor2,drhor2,rhopp
!         write(31,5000) r,rhor2,drhor2,rhopp
      enddo
!      close(31)
 4999 FORMAT(2E18.8)
 5000 FORMAT(4E18.8)
      return
      end

