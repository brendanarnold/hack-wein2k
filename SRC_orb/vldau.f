SUBROUTINE vldau(nup,nldau,l,dmat,U,J,V,mult,ipr,eorb)
  use exex
  IMPLICIT REAL*8(A-H,O-Z)
  ! calculation of additional LDA+U potential, three versions
  ! of double counting are implemented:
  ! 1/ 'self-interaction corrected' LDA+U of Anisimov et al. 1993. 
  ! 2/ 'around the mean field' (nldau=0 in the input). Czyzyk and Sawatzky 1994.
  ! 3/  Mean field Hubbard model (nldau=2 in input) Anisimov et al. 1991
  !     this option is to be used with LDA (not LSDA)
  ! *****************************************************************
  real*8 J,trace(2),n0(2),ntot
  complex*16 dmat(-3:3,-3:3,2),V(-3:3,-3:3),czero      
  real*8 Vee(-3:3,-3:3,-3:3,-3:3)
!  COMMON/FACT/FCT(0:100),nfac   
!  nfac=20
!  call fac
  pi=acos(-1.d0)
  sqfp=sqrt(4.d0*pi)
  !**********************************************************************
  !
  czero=(0.d0,0.d0)
  ipr3j=0
  if(nup.eq.2)then
     write(6,*)' DNUP block of orbital potential'
  else
     !
     ! traces of density matrices
     do isp=1,2
        trace(isp)=0.d0
        do m=-l,l
           trace(isp)=trace(isp)+dmat(m,m,isp)
        enddo
     enddo
     ! average number of electron per state
     do isp=1,2
        n0(isp)=trace(isp)/dfloat(2*l+1)
     enddo
     ntot=(n0(1)+n0(2))/2.d0
     trtot=trace(1)+trace(2)
     ! 
     if(ipr.gt.1)then
        write(6,589)(n0(isp),isp=1,2),ntot
        write(6,*)
     endif
589  format(' nspin1',f10.5,' nspin2',f10.5,' nspintot',f10.5)
     ! DFT (conventional) double counting correction 
     ! subtract from diagonal elements 
     ! of density matrix average number of electron per state. 
     ! AMF method
     if(nldau.eq.0)then
        do isp=1,2
           do m=-l,l
              dmat(m,m,isp)=dmat(m,m,isp)-n0(isp)
           enddo
           write(6,600)isp,(dble(dmat(m,m,isp)),m=-l,l)
600        format(i3,7f11.5)
        enddo
     endif
     ! MFH method
     if(nldau.eq.2)then
        do isp=1,2
           do m=-l,l
              dmat(m,m,isp)=dmat(m,m,isp)-ntot   
           enddo
           write(6,600)isp,(dble(dmat(m,m,isp)),m=-l,l)
        enddo
     endif
     !
  endif
  !
  ! calculate the matrix Vee(M,M',M'',M''')
  call Vcalc(Vee,L,U,J)
  Eldau=0.0d0
  do m=-l,l
     do m1=-l,l
        ! double sum 
        v(m,m1)=czero
        do m2=-l,l
           do m3=-l,l
             if(nldau.eq.4)then
              ! exchange contribution, present also for cross term
               v(m,m1)=v(m,m1)-Vee(m,m2,m3,m1)*dmat(m2,m3,1)
               Eldau=Eldau-Vee(m,m2,m3,m1)*(dmat(m,m1,1))*dmat(m2,m3,1)
             else
              if(nup.eq.2)then
                 v(m,m1)=v(m,m1)-Vee(m,m2,m3,m1)*dmat(m2,m3,1)
                 ! Mean field (H-F) LDA+U energy
!!!!PB                 Eldau=Eldau-Vee(m,m2,m3,m1)*dmat(m,m1,1)*dmat(m2,m3,1)
                 Eldau=Eldau-Vee(m,m2,m3,m1)*(dmat(m,m1,1))*dmat(m2,m3,1)
              else
                 ! first term with opposite spin &
                 ! second term with the same spin as potential &
                 v(m,m1)=v(m,m1)+ &
                      Vee(m,m2,m1,m3)*dmat(m2,m3,2) + &
                      (Vee(m,m2,m1,m3)-Vee(m,m2,m3,m1))*dmat(m2,m3,1)
                 ! Mean field (H-F) LDA+U energy
!!!!PB                 Eldau=Eldau+ &
!                      Vee(m,m2,m1,m3)*dmat(m,m1,1)*dmat(m2,m3,2)+ &
!                      (Vee(m,m2,m1,m3)-Vee(m,m2,m3,m1)) &
!!!!                      *dmat(m,m1,1)*dmat(m2,m3,1)
                 Eldau=Eldau+ &
                      Vee(m,m2,m1,m3)*(dmat(m,m1,1))*dmat(m2,m3,2)+ &
                      (Vee(m,m2,m1,m3)-Vee(m,m2,m3,m1)) &
                      *(dmat(m,m1,1))*dmat(m2,m3,1)
              endif
             endif
           enddo
        enddo
     enddo
     ! double counting correction to the potential 
     if((nldau.eq.1).and.(nup.ne.2))   &
          v(m,m)=v(m,m)-U*(trtot-0.5d0)+J*(trace(1)-0.5d0)
  enddo
  Eldau=Eldau/2.d0
  Edc=0.d0
  trdmv=0.d0
!  if(nup.ne.2)then       ! for updn-correction by P.Novak 27.3.2003
     ! double counting correction to the energy 
     if((nldau.eq.1).and.(nup.ne.2)) Edc=0.5d0*U*trtot*(trtot-1.d0) - &
          0.5d0*J*(trace(1)*(trace(1)-1.d0)+trace(2)*(trace(2)-1.d0))
     !
     ! trace dmat(1)*v
     do m=-l,l
        do m1=-l,l
!!!!!PB           trdmv=trdmv+dmat(m,m1,1)*v(m1,m)
           trdmv=trdmv+dmat(m,m1,1)*conjg(v(m1,m))
           if(nup.eq.2)then                         ! added
             trdmv=trdmv+v(m,m1)*dmat(m1,m,1)       ! added
           endif                                    ! added
        enddo
     enddo
!  endif
  ! overall contribution to total energy of given spin
  if(nldau.ne.4)then
    Ecorr=Eldau-Edc/2.d0-trdmv
    write(6,107)Ecorr,mult,Eldau,Edc,trdmv
107 format(' Ecorr',f11.5,' Mult',i3,' Eldau',f11.5,' Edc',f11.5, &
       ' Tr(rho.V)',f11.5)
  else
! approximate exact exchange
     Ecorr=Eldau-Exc-trdmv       !AExEx energy correction
     vspxc  =fxc/sqfp             !EECE  'double counting' potential starts
     do m=-l,l
      v(m,m)=v(m,m)-vspxc         !spherical contribution
      if(nsp.ne.0)then           !nonspherical contribution if nsp.ne.0
       do m1=-l,l
        v(m,m1)=v(m,m1)-vxcl(m,m1) !nonspherical contribution
       enddo
      endif
     enddo
     do m=-l,l                 !mixing in the hybrid and output
      do m1=-l,l                 
       v(m,m1)=hybr*v(m,m1)     
        if(m1.eq.m)then
         write(6,667)m,m1,v(m,m1),vxcl(m,m1),vspxc
        else
         write(6,668)m,m1,v(m,m1),vxcl(m,m1)
667      format(2i3,2f12.6,3x,3f12.6,' m,m1,v,vnspxc,vspxc')
668      format(2i3,2f12.6,3x,2f12.6,' m,m1,v,vnspxc')
        endif
       enddo
      enddo
      Ecorr=(Eldau-Exc-trdmv)*hybr
      write(6,108)Ecorr,mult,Eldau,Exc,trdmv
108 format(' Ecorr',f11.5,' Mult',i3,' ExHF ',f11.5,' Exc(d.c.)',f11.5, &
       ' Tr(rho.VHF)',f11.5)
  endif           !AExEx   end
  EORB=Ecorr*mult
  write(6,106)EORB
106 format(':EORB:',f13.6)
  return
end SUBROUTINE vldau
