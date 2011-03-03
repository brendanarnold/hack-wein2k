      SUBROUTINE FOMAI2(jatom,nr,lmmax,fvdrho)                           
!                                                                       
      use param
      use defs
      use struk; USE xa
      IMPLICIT REAL*8 (A-H,O-Z)
!
!bk   93/08/26: parameters, variables, and commons 
!     for force calculation
      real*8          fvdrho(0:3,ndif),vlm1(nrad,ncom+1)
!                                                                       
      COMMON /POTNLC/ VR(NRAD)

!
!     93/09/13: calculate integral of Veff*grad(rholm)
!     + working with real spherical y_lmp, p=+/-
!     (A12) --> fvdrho(3)
!     variables:
!     fvdrho(3) --> total force from potential*grad(density) integral
!
      ja=jatom
      nr=jri(ja)
!     loop over potential
        read(19,'(3x)')
        read(19,'(15x,i3//)') lmp
        DO 5016 ilm=0,lmp
!     read potential, (lm=(0,0) already read
!     normalize potential with vlm=V(l,m)*r*r
!     lm=(0,) normalized differently originally 
!     attention: vr(ir) is in hartree !!!!!!!!!!!
          if (ilm.eq.0) then
            DO 5017 ir=1,nr
              vlm1(ir,1)=vr(ir)*sqrt(4.d0*pi)*2.d0/r(ir)
 5017       CONTINUE
          else
            read(19,'(15x,i3,5x,i2/)') llp,mp
            read(19,'(3x,4e19.12)') ( vlm1(ir,ilm+1), ir=1,nr)
            read(19,'(/)')
          endif
 5016   CONTINUE
        read(19,'(///)')
!ccc
      call vdrho(lmmax,lm,jri(jatom),dx(jatom),r,vlm1,rholm, &
        fvdrho(1,jatom))
      RETURN
      END                                                               


