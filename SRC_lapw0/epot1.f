      SUBROUTINE EPOT0  (myid)
      use efeld
      implicit real*8 (a-h,o-z)
      character*100 LINE
      character*80  STYLE(0:30)
        iefeld=0
        read(5,101,end=105,err=105)LINE
101    format(a)
        read(LINE,*,end=102,err=102) iefeld,refeld,wefeld
        goto 103
102    wefeld=0.3
        read(LINE,*,end=105,err=105) iefeld,refeld
103    continue
        if(mod(abs(iefeld),1000).gt.2000) then
                  CALL OUTERR('LAPW0','Too many Fourier Coefficients')
                  stop 'Error: Too many Fourier Coefficients'
        endif
!
        if(abs(refeld) .lt. 1D-10)then
                iefeld=0
                return
        endif
!       Give helpful output for the user
        mode1=abs(iefeld)/1000
        mode=mode1
        if(mode1 .ge. 40)then
                mode=mode-40
                write(6,210)
                if(myid .eq. 0)write(21,210)
        else if(mode1 .ge. 20)then
                mode=mode-20
                if(mode.ne.0)then
                        write(6,211)
                        if(myid .eq. 0)write(21,211)
                endif
        endif
        do j=0,30
           style(j)='Unknown Option'
        enddo
        style(11) ='Numerical Triangular ramp from E/2 to -E/2'
        style(1) ='=0 if |z|<lambda, otherwise E*cos((z0-lambda)*pi2/(0.5-lambda))'
        style(2) ='=0 if |Z|<lambda, otherwise E*(1.D0-cos((z0-lambda)*pi2/(0.5-lambda)))'
        style(3) ='=0 if |z|<lambda, otherwise E*(1-cos((z0-lambda)*pi/(0.5-lambda)))*0.5D0'
        style(4) ='=E if |z|<lambda, otherwise E*(0.5D0-Z0)/(0.5d0-lambda)'
        style(5) ='=0 if |z|<lambda, otherwise E*(1.D0-(0.5D0-Z0)/(0.5d0-lambda))'
        style(6) ='Form from Lozovoi et al, J. Chem Phys, 115, 1661, 2001, normalized to max of 1'
        style(7) ='Form from Lozovoi et al, J. Chem Phys, 115, 1661, 2001'
        style(8) ='=Constant E'
        style(9) ='=0 if |z|<lambda, otherwise E'
        style(10)='=E if |z|<lambda, otherwise E*((1-cos((z0-lambda)*pit/cut))*0.5D0)'
        style(0)='Analytic triangular ramp'
        style(12) ='Numerical Triangular ramp from 0 to -E'
!
!       Help output
        if(iefeld .eq. -999)then
                write(*,401)'Available E-field modes'
401             format(a)
                do j=0,12
                   write(*,400)j,style(j)
400                format('Mode ',i2,' is ',a)
                enddo
                write(*,401)'Add 20 for finite Fourier series'
                write(*,401)'Add 40 for Fourier expansion within RMTs'
                stop
        endif
!        
        if(style(mode)(1:5) .eq. 'Unkno')then
             write(6,*)'Unsupported E-field option, no E-field'
             if(myid .eq. 0)write(21,*)'Unsupported E-field option, no E-field'
             iefeld=0
             return
        endif
        write(6, 200)refeld,mod(abs(iefeld),1000)
        write(6, 201)style(mode)
        if(myid .eq. 0)then
                write(21, 200)refeld,mod(abs(iefeld),1000)
                write(21, 201)style(mode)
        endif
        if((mode .eq. 0).or.(mode .eq.8).or.(mode.eq.11))goto 300
          write(6, 202)wefeld
          if(myid .eq. 0)write(21,202)wefeld
300    continue
200    format(':EFIELD       E-field of ',F10.5,' Ry with ',i4,' Fourier coefficients requested')
201    format(':EFIELD       Form is ',a)
202    format(':EFIELD       Lambda =',F10.5)
210    format(':EFIELD       Fourier expansion used within Muffin-Tins')
211    format(':EFIELD       Limited Fourier Series used for Plane Waves')
105    CONTINUE
      return
      end

      SUBROUTINE EPOTINIT(NKK,myid)
      use efeld
      use symetr
      implicit real*8 (a-h,o-z)
      jj=0
      kefeld=0
      if(iefeld .eq. 0)return
! select all (0,0,kz)
         jj=0
         kebig=0
         DO J=1,NKK
           if(KZZ(1,J).eq.0.and.KZZ(2,J).eq.0) then
            jj=jj+1
            if(jj.gt.mod(abs(iefeld),1000)) goto 43
            jefeld(jj)=j
            kebig=max(abs(kzz(3,J)),kebig)
!            write(*,*)'Using ',jj,kzz(1,j),kzz(2,j),kzz(3,j)
           endif
         enddo
 43   continue
      kefeld=min(jj,mod(abs(iefeld),1000))
! Establish values
      iefmode=abs(iefeld)/1000
      write(6,100)kefeld,kebig
      if(myid .eq. 0)write(21,100)kefeld,kebig
100   format(':EFIELD       Number of FC in E-field ',i3,' , Largest ',i4)
      call makeback(fefeldpw,wefeld,kefeld,refeld,iefmode,myid)
      return
      end
!
      SUBROUTINE EPOT3  (v,NKK1)
!     Add the PW co-efficients
      use efeld
      implicit real*8 (a-h,o-z)
      complex*16 V(nkk1)
!     Add the E-field to the PW's
          do i=1,kefeld
             jj=jefeld(i)
             v(jj)=v(jj)+fefeldpw(i)
          enddo
      return
      end
!
      SUBROUTINE EPOT4  (v,NKK1)
!     Add the PW co-efficients
      use efeld
      implicit real*8 (a-h,o-z)
      complex*16 V(nkk1,2)
!     Add the E-field to the PW's
          do i=1,kefeld
             jj=jefeld(i)
             v(jj,1)=v(jj,1)+fefeldpw(i)
             v(jj,2)=v(jj,2)+fefeldpw(i)
          enddo
      return
      end

      SUBROUTINE EPOT1c (LM,JATOM,LMMAX,Vfield)
      use efeld
      use struct
      use densit
      use work
      use vresp1
!                                                                       
!     Projects the E-field onto spherical harmonics in the RMT's
!     LDM, November 2006
!     Note, we do need to include higher harmonics both for some of
!     the additional forms and for atoms at the origin
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
!
      COMPLEX*16     YL
!
      DIMENSION         YY(ncom,nspd)                                    
      DIMENSION         bbu(NCOM)
!
      DIMENSION  YL((lmax2+1)*(lmax2+1))
      DIMENSION  LM(2,NCOM+3,NAT),LMMAX(NAT)         
      DIMENSION  A2(3),ROTA1(3)
      DIMENSION  Vfield(ncom+3,nrad)
!
      dimension spt(3,nspd),weight(nspd),points(nspd)
!
      INDEX=1
      DO J=1,JATOM-1
        INDEX=INDEX+MULT(J)
      ENDDO
!
      LLMM=LMMAX(JATOM)
!
      LUSE=0
      do L=1,LLMM
        LUSE=MAX(LUSE,abs(lm(1,L,JATOM)))
      enddo
!
!     generate Gauss-Legendre points on sphere
!                                                                       
      call gpoint(spt,weight,npt1,LUSE)

!     setup YY array
      do 1 ind=1,npt1
         call rotate (spt(1,ind),ROTLOC(1,1,JATOM),ROTA1)
         CALL YLM (ROTA1,lmax2,YL)                                          
         DO 3 LM1=1,LLMM        
          CALL SUML(LM1,YK,YL,LM(1,1,JATOM))                             
          YY(LM1,ind)=YK
 3       CONTINUE                                                          
!
  1   CONTINUE
!                                                                       
!     LOOP OVER ALL RADIAL MESH POINTS
!
      DO ir=1,JRI(JATOM)
!
       do lm1=1,llmm
        bbu(lm1)=0.0D0
       enddo
!
       RIR=R0(JATOM)*EXP(DX(JATOM)*dble(ir-1))/a(3)
       imode=abs(iefeld)/1000
       ifourier=1
       DO k=1,npt1
        z=spt(3,k)*rir+pos(3,index)
!       LDM general case
        tmp1= backget(z,imode,wefeld,ifourier)*weight(k)
!
!       Project onto harmonics
        DO LM1=1,LLMM
         bbu (lm1)=bbu (lm1)+yy(lm1,k)*tmp1
        enddo
!
       enddo
!
!     Save
!
        DO LM1=1,LLMM
                 Vfield(LM1,ir)=bbu(LM1)*refeld
        enddo
        
      enddo
      RETURN
      END
        subroutine patchcsymm(kzz,v,nkk)
        use struct, only : NSYM
        implicit real*8 (a-h,o-z)
        complex*16 v(nkk), val1,val
        complex*16 taup(nsym)
        integer kzz(3,nkk),stg(3,nsym)
        v(1)=dble(v(1))
        do j=2,nkk-1
                val1=v(j)
                call stern(j,nst,stg,taup)
                do jj=1,nst
                   i1=stg(1,jj)
                   i2=stg(2,jj)
                   i3=stg(3,jj)
                   do k=j+1,j+50
                        j1=kzz(1,k)+i1
                        j2=kzz(2,k)+i2
                        j3=kzz(3,k)+i3
                        if((j1.eq.0).and.(j2.eq.0).and.(j3.eq.0))then
                                val=(conjg(v(k)*taup(jj))+val1)*0.5d0
                                v(j)=val
                                v(k)=conjg(val)*taup(jj)
                                goto 10
                        endif
                  enddo
                enddo
10              continue
        enddo
        return
        end
      subroutine EFORCE(FHF)
      use efeld
      use struct
!     Add the Hellman-Feynman forces due to the efeld
      implicit real*8 (a-h,o-z)
      dimension FHF(0:3,1:nat),F1(3),F2(3)
      index=1
      write(6,*)
      do j=1,nat
         z=pos(3,index)
         imode=abs(iefeld)/1000
         ifourier=1
!        LDM general case
         tmp1=zz(J)*backgetder(z)
         F1(1)=0
         F1(2)=0
         F1(3)=tmp1
         call rotate (F1,ROTLOC(1,1,J),F2)
         FHF(1:3,J)=FHF(1:3,J)+ F2(1:3)
         write(6,100)j,tmp1*1d3,index,z
         index=index+mult(j)
      enddo
      write(6,*)
100   format(':EFIELD Z Force for atom ',i3,' = ',F13.6,' mRyd/au, #',i3,' z=',F13.6)
      return
      end
