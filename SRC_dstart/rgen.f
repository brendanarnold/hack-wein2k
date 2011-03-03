


      subroutine rgen(rhok,zatom,zatom1,akap)

      use struct 
      use waves
      IMPLICIT REAL*8 (A-H,O-Z)
!      include 'param.inc'
      
      dimension        rhok(nwave,nat+1),akap(4,nat)
      dimension        zatom(3,nat),zatom1(3,nat)
      COMPLEX*16        RHOk,ch1,ch2,TAUP
      complex*16       cima
      INTEGER          STG
!      CHARACTER*4      LATTIC                                 
!      CHARACTER*5      MODUS                                            
!      CHARACTER*10     ANAME 
!      CHARACTER*80     TITLE
!     
!      COMMON /STRUK/  POS(3,NDIF),AA,BB,CC,ALPHA(3),RMT(NATO),V(NATO),
!     *     PIA(3),VOL,IATNR(NATO),MULT(NATO),ISPLIT(NATO)       
!      COMMON /CHAR/   TITLE,LATTIC,MODUS,ANAME(NATO)                    
      COMMON /GENER/  BR1(3,3),BR2(3,3)                                 
!      COMMON /GITT/   KZZ(3,NWAV),ABSK(NWAV),nwave
!      COMMON /SYM2/   TAU(3,NSYM),IORD,IZ(3,3,NSYM)
!      COMMON /POTNLC/ R0(NATO),DX(NATO),JRI(NATO) 
      COMMON /COM/    EMIN,EF,ELECN,XWT,NSPIN, &
           NBAND,NK,MINWAV,MAXWAV           
      COMMON /FOUR/  kx,ky,kz,NST,STG(3,NSYM),TAUP(NSYM)
      

      s(xxx)=fpi/vol*(dsin(xxx)-xxx*dcos(xxx))
      

      fpi=4.d0*acos(-1.d0)
      tpi=2.d0*acos(-1.d0)
      cima=(0.d0,1.d0)
      
!..   construction of interstitial charge density
!     by summing over tails of all atoms      
      
      do ik=1,nwave
         kx=kzz(1,ik)
         ky=kzz(2,ik)
         kz=kzz(3,ik)
         akx=kx*br1(1,1)+ky*br1(1,2)+kz*br1(1,3)
         aky=kx*br1(2,1)+ky*br1(2,2)+kz*br1(2,3)
         akz=kx*br1(3,1)+ky*br1(3,2)+kz*br1(3,3)
         akk=dsqrt(akx*akx+aky*aky+akz*akz)
         do ia=1,nat
            index=0
            do ia1=1,ia-1
               index=index+mult(ia1)
            enddo
            ch2=0.
            do ip=1,mult(ia)
               index=index+1
               tk=(kx*pos(1,index)+ky*pos(2,index)+ &
                    kz*pos(3,index)) * tpi
!               ch1=exp(-cima*tk*fpi/2.)
               ch1=dcmplx(cos(tk),-sin(tk))
               ch2=ch2+ch1
            enddo
            rhok(ik,nat+1)=rhok(ik,nat+1)+ch2*rhok(ik,ia)
         enddo
      enddo

!       normierung mit anzahl der sternmitglieder

      do ik=1,nwave
         kx=kzz(1,ik)
         ky=kzz(2,ik)
         kz=kzz(3,ik)
         akx=kx*br1(1,1)+ky*br1(1,2)+kz*br1(1,3)
         aky=kx*br1(2,1)+ky*br1(2,2)+kz*br1(2,3)
         akz=kx*br1(3,1)+ky*br1(3,2)+kz*br1(3,3)
         akk=dsqrt(akx*akx+aky*aky+akz*akz)
         call stern   
         do ia=1,nat+1
            rhok(ik,ia)=rhok(ik,ia)*nst
         enddo
      enddo

         
         
!     ...  bestimmen des tail anteils in jeder MT Kugel
!     ...  wird spaeter genutzt um mit l=0 anteil im
!     ...  MT zu vergleichen
         
         
      do ia=1,nat
         index=1
         do ia1=1,ia-1
            index=index+mult(ia1)
         enddo
         zatom1(1,ia)=0.d0
         
         
         do ik=1,nwave
            kx=kzz(1,ik)
            ky=kzz(2,ik)
            kz=kzz(3,ik)
            akx=kx*br1(1,1)+ky*br1(1,2)+kz*br1(1,3)
            aky=kx*br1(2,1)+ky*br1(2,2)+kz*br1(2,3)
            akz=kx*br1(3,1)+ky*br1(3,2)+kz*br1(3,3)
            akk=dsqrt(akx*akx+aky*aky+akz*akz)
            
            call stars1(index,ch1)
            if (akk.lt.1.e-6) then
               ss=fpi*(rmt(ia)**3)/vol/3.d0
            else
               ss=s(akk*rmt(ia))/akk/akk/akk
            endif
            zatom1(1,ia)=zatom1(1,ia)+rhok(ik,nat+1)* &
                 ch1*ss/nst
         enddo
!        write(6,102) aname(ia)
!        write(6,103) zatom1(1,ia),
!    $      zatom1(1,ia)-akap(4,ia)
      enddo
 102  format(3x,'tail charge after FT induced in ', &
           'MT of atom   ',a10)
 103  format(3x,2f10.5)
      
      return
      end
