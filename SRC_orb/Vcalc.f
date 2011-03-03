      subroutine Vcalc(Vee,L,U,J)
      use exex
! calculation of Vee 
! according to [1] Liechtenstein et al. Phys.Rev. B 52, R5467 (1995)
      implicit real*8(a-h,o-z)
      integer q
      real*8 J,F(0:6)
      real*8 Vee(-3:3,-3:3,-3:3,-3:3)
      common/print/ipr,ipr3j,iprv
     if(nexex.eq.0)then
! calculate Slater parameters F(k) from U,J
      F(0)=U
      if(l.eq.2)then
! d-states: ratio r42=F(4)/F(2)= 0.625 ([1])
      r42=0.625
      F(2)=14.*J/(1.+r42)
      F(4)=F(2)*r42
      write(6,100)f(0),F(2),F(4)
100   format(' Slater integrals F0, F2, F4',3f8.3,' Ry')
      endif
      if(L.eq.3)then
! f-states: see Solovyev et al. Phys. Rev. B50, 16 861 (1994), ref. 33
      r42=451./675.
      r62=1001./2025.
!      F(2)=3.*J/(1.+r42+r62)
      F(2)=3.*J/(2.d0/15.d0+1.d0/11.d0*r42+50.d0/429.d0*r62) 
      F(4)=F(2)*r42
!      F(4)=F(4)*r42
      F(6)=F(2)*r62
      write(6,102)f(0),F(2),F(4),F(6)
102   format(' Slater integrals F0, F2, F4, F(6)',4f8.3,' Ry')
      endif
     else               !EECE
      do l1=0,6,2
       f(l1)=fupka(l1)
 write(6,550)l1,f(l1)
      enddo
550  format(i4,f12.6,' k, F^k')
      endif             !EECE end
! calculate <m,m2|Vee|m1,m3> eq. 6 in [1]
      do m=-l,l
!      m=-l
         do m2=-l,l
            do m1=-l,l
!            m1=-l
               do m3=-l,l
               Vee(m,m2,m1,m3)=0.
                 do k=0,2*L,2
! calculate the ak coefficients eq. 7
                 ak=0.
                   t1=t3j(l,k,l,0,0,0)
                     do q=-k,k
                     t2=t3j(l,k,l,-m,q,m1)
                     t3=t3j(l,k,l,-m2,-q,m3)
!                    sg=(-1)**(m+q+m1)
                     sg=(-1)**(m+q+m2)
                     ak=ak+sg*t2*t3
                     ddd=  sg*t2*t3
                     enddo
                 cc=(2*L+1)**2*t1**2
                 ak=(2*L+1)**2*t1**2*ak
                Vee(m,m2,m1,m3)=Vee(m,m2,m1,m3)+ak*F(k)
      if(iprv.gt.3)write(6,101)k,ak,F(k),Vee(m,m2,m1,m3)
101   format(i3,3f12.5,' k,ak,F(k),Vee a.u.')
                enddo
                if(iprv.gt.2)write(6,103)m,m2,m1,m3,Vee(m,m2,m1,m3)
103   format(4i3,f12.5,' m,m2,m1,m3,Vee a.u.')
              enddo
             enddo
           enddo
        enddo
      return
      end













