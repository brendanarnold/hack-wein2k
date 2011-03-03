      SUBROUTINE XSPLT (ALM,BLM,clm,MULT,num)    
!                                                                       
!     DECOMPOSITION OF P CHARGE IN PX,PY,PZ                             
!                                                                       
        USE charp
      use param
      use lo; use xdos
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16       ALM,BLM,clm
      DIMENSION      ALM((lmax2+1)*(lmax2+1)),BLM((lmax2+1)*(lmax2+1))     
      DIMENSION        cLM((LMAX2+1)*(LMAX2+1),NUME,nloat)
  !                                                                       
 
!      dimension        ilo(                       
!
      common /normal/  pu1u1(0:lxdos,0:lxdos),pu1ue(0:lxdos,0:lxdos) &
                      ,pu1u2(0:lxdos,0:lxdos,nloat),pueue(0:lxdos,0:lxdos) &
                      ,pueu2(0:lxdos,0:lxdos,nloat),pu2u2(0:lxdos,nloat,0:lxdos,nloat)
!---------------------------------------------------------------------  
!
!.....double loop over l,m,l',m' (combined indices) up to l=2
      do l=0,lxdos
      do m=-l,l
      ly=l*(l+1)+m+1
      do lp=0,l
      do mp=-lp,lp
      lpy=lp*(lp+1)+mp+1

!                                                                       
      xqtl(ly,lpy,num)=xqtl(ly,lpy,num)+ &
                       (alm(ly)*conjg(alm(lpy))*pu1u1(l,lp)+ &
                        blm(ly)*conjg(blm(lpy))*pueue(l,lp)+ &
                        alm(ly)*conjg(blm(lpy))*pu1ue(l,lp)+ &
                        blm(ly)*conjg(alm(lpy))*pu1ue(lp,l)) &
                        *100.D0/MULT
           do jlo=1,ilo(l)
             do jlop=1,ilo(lp)
             xqtl(ly,lpy,num)=xqtl(ly,lpy,num)+ &
                        clm(ly,num,jlo)*conjg(clm(lpy,num,jlop))*pu2u2(l,jlo,lp,jlop) &
                        *100.D0/MULT
             enddo
           xqtl(ly,lpy,num)=xqtl(ly,lpy,num)+ &
                        (clm(ly,num,jlo)*conjg(alm(lpy))*pu1u2(lp,l,jlo)+ &
                        clm(ly,num,jlo)*conjg(blm(lpy))*pueu2(lp,l,jlo)) &
                        *100.D0/MULT
           enddo
           do jlop=1,ilo(lp)
           xqtl(ly,lpy,num)=xqtl(ly,lpy,num)+ &
                        (alm(ly)*conjg(clm(lpy,num,jlop))*pu1u2(l,lp,jlop)+ &
                        blm(ly)*conjg(clm(lpy,num,jlop))*pueu2(l,lp,jlop))  &
                        *100.D0/MULT
           enddo
!
!   if(l.eq.0.and.lp.eq.0.and.num.eq.1) then
!     write(6,*) 'xqtl',xqtl(ly,lpy,num), &
!                       alm(ly),conjg(alm(lpy)),pu1u1(l,lp), &
!                        blm(ly),conjg(blm(lpy)),pueue(l,lp), &
!                        clm(ly),conjg(clm(lpy)),pu2u2(l,lp), &
!                        alm(ly),conjg(blm(lpy)),pu1ue(l,lp), &
!                        blm(ly),conjg(alm(lpy)),pu1ue(lp,l), &
!                        alm(ly),conjg(clm(lpy)),pu1u2(l,lp), &
!                        clm(ly),conjg(alm(lpy)),pu1u2(lp,l), &
!                        blm(ly),conjg(clm(lpy)),pueu2(l,lp), &
!                        clm(ly),conjg(blm(lpy)),pueu2(lp,l),'end xqtl'
!   endif
      enddo
      enddo
      enddo
      enddo
      RETURN                                                            
      END                                                               















