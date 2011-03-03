SUBROUTINE FOMAI1(latom,jatom,nemin,nemax,lmx,ihmx,alm,blm, &
     clm,aalm,bblm,cclm,tuu,tdd,tud, &
     tdu,tuu12,tuu21,tuu22,tud21,tdu12,lv,lpv, &
     mv,mpv,forout,fsph,fnsp,fsph2)
  !                                                                       
  USE param
  USE defs
  USE struk
  USE lo; USE lohelp
  USE xa;  USE atspdt
  IMPLICIT REAL*8 (A-H,O-Z)
  !
  COMPLEX*16       alm((lmax2+1)*(lmax2+1),nume),               &
       blm((lmax2+1)*(lmax2+1),nume), &
       clm((lmax2+1)*(lmax2+1),nume,nloat)              
  COMPLEX*16       aalm(3,nume,(lmax2+1)*(lmax2+1)) &
       ,bblm(3,nume,(lmax2+1)*(lmax2+1)) &
       ,cclm(3,nume,(lmax2+1)*(lmax2+1),nloat) &
       ,afac,bfac,cfacf(nloat),cfac2,efac,ffac &
       ,tuu(ngau),tdd(ngau),tud(ngau),tdu(ngau) &
       ,tuu12(ngau),tuu21(ngau),tuu22(ngau) &
       ,tud21(ngau),tdu12(ngau),kinfac1,kinfac2 &
       ,kinfac3,kinfac4
  INTEGER          lv(ngau),lpv(ngau),mv(ngau),mpv(ngau)
  LOGICAL          forout
  REAL*8           fsph(0:3,ndif),fnsp(0:3,ndif),fsph2(0:3,ndif)
!                                                                       
!.... MAIN FORCE CALCULATIONS
!bk   93/08/26: 
! OBS GM: No cross terms yet for several LO's

  sqrt2=SQRT(2.0d0)

  ia=latom
  DO num=nemin,nemax
     DO l=0,lmx
        DO m=-l,l
           ly=l*(l+1)+m+1
!
!   the following write does nothing, but fixes a pgi compiler bug 
           IF(num.EQ.0) WRITE(*,*) num,l,m                
!     (B1) fsph
           afac=2.d0*alm(ly,num)*(el(l)-e(num))+blm(ly,num)
           DO jlo=1,ilo(l)
                afac=afac +  &
                   clm(ly,num,jlo)*pi12lo(jlo,l)*(elo(l,jlo)+el(l)-2.d0*e(num))
           ENDDO
           bfac=alm(ly,num) + &
                2.d0*blm(ly,num)*(el(l)-e(num))*pei(l)
           DO jlo=1,ilo(l)
              bfac=bfac +clm(ly,num,jlo) &
                   * ( pi12lo(jlo,l) + pe12lo(jlo,l)*(elo(l,jlo)+el(l)-2.d0*e(num)))
           ENDDO
           cfacf=(0.0d0,0.0d0)
           DO jlo=1,ilo(l)
              cfacf(jlo)= blm(ly,num)*pi12lo(jlo,l) &
                        + (elo(l,jlo) + el(l)-2.d0*e(num)) & 
                           * (alm(ly,num)*pi12lo(jlo,l)+pe12lo(jlo,l)*blm(ly,num))
           ENDDO

! (B4) extra surface term for new kinetic energy operator
           kinfac1=(0.0d0,0.0d0)
           kinfac2=(0.0d0,0.0d0)
           kinfac1=alm(ly,num)* p(l) + blm(ly,num)*pe(l)
           kinfac2=alm(ly,num)*dp(l) + blm(ly,num)*dpe(l)
           DO jlo=1,ilo(l)
              kinfac1=kinfac1+clm(ly,num,jlo)*plo(jlo,l)
              kinfac2=kinfac2+clm(ly,num,jlo)*dplo(jlo,l)
           ENDDO
           DO ik=1,3
              fsph(ik,ia)=fsph(ik,ia) &
                   +weight(num)*dimag( aalm(ik,num,ly)*dconjg(afac) &
                                      +bblm(ik,num,ly)*dconjg(bfac))
              DO jlo=1,ilo(l)
                 DO jlop=1,ilo(l)
                    fsph(ik,ia)=fsph(ik,ia)+weight(num)* &
                         dimag(cclm(ik,num,ly,jlo)*(dconjg(cfacf(jlop)+ &
                    2.d0*clm(ly,num,jlo)*(elo(l,jlo)-e(num))*pr12lo(jlo,jlop,l))))
                 ENDDO
              ENDDO
              kinfac3=(0.0d0,0.0d0)
              kinfac4=(0.0d0,0.0d0)
              kinfac3=aalm(ik,num,ly)*dp(l)+bblm(ik,num,ly)*dpe(l)
              kinfac4=aalm(ik,num,ly)* p(l)+bblm(ik,num,ly)*pe(l)
              DO jlo=1,ilo(l)
                 kinfac3=kinfac3+cclm(ik,num,ly,jlo)*dplo(jlo,l)
                 kinfac4=kinfac4+cclm(ik,num,ly,jlo)*plo(jlo,l)
              ENDDO
              fsph2(ik,ia)=fsph2(ik,ia)+weight(num)* &
                   dimag(rmt(jatom)**2*( dconjg(kinfac1)*kinfac3 &
                                        -dconjg(kinfac4)*kinfac2) )
           ENDDO
        ENDDO
     ENDDO
!     (A20) fnsp
     DO ih=1,ihmx
        lpy= lpv(ih)*(lpv(ih)+1)+mpv(ih)+1
        ly = lv(ih) *( lv(ih)+1)+ mv(ih)+1
        afac= dconjg(alm(lpy,num))* tuu(ih) &
            + dconjg(blm(lpy,num))* tdu(ih) 
        DO jlo=1,ilo(lpv(ih))
           afac = afac + dconjg(clm(lpy,num,jlo))* tuu21(ih)
        ENDDO
        bfac= dconjg(alm(lpy,num))* tud(ih) &
            + dconjg(blm(lpy,num))* tdd(ih)
        DO jlo=1,ilo(lpv(ih))
           bfac = bfac + dconjg(clm(lpy,num,jlo))* tud21(ih)
        ENDDO
        cfac2= dconjg(alm(lpy,num))* tuu12(ih) &
             + dconjg(blm(lpy,num))* tdu12(ih)
        DO jlo=1,ilo(lpv(ih))
            cfac2=cfac2 + dconjg(clm(lpy,num,jlo))* tuu22(ih)
        ENDDO
        DO ik=1,3
           fnsp(ik,ia)=fnsp(ik,ia) &
                +2.d0*weight(num)*dimag( afac*aalm(ik,num,ly) &
                                        +bfac*bblm(ik,num,ly))
           DO jlop=1,ilo(lv(ih))
              fnsp(ik,ia)=fnsp(ik,ia)+2.d0*weight(num) &
                                          *dimag(cfac2*cclm(ik,num,ly,jlop))
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  RETURN 
END SUBROUTINE FOMAI1

