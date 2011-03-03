      subroutine sumdd(r,cth,sth,rho,drrho,drdrrho,ang, &
         dtang,dtdtang,dfdtang,dfang,dfdfang,hrho,lmmax,lm,jatom)
! Calculates cubic harmonics after Kara & Kurki-Suonio
! Acta Cryst A 1981 37 201-210
! GKHM 2/5-01
!     Small changes by L. D. Marks, January 2006
       IMPLICIT NONE
       INCLUDE 'param.inc'
       INTEGER  jatom,lmmax,lm(2,ncom,*)
       REAL*8   c1,c2,c3,fac
       REAL*8   c_kub(0:10,0:10)
       INTEGER  i,j,ia,ir
       real*8 rho(ncom),ang(ncom),drrho(ncom),drdrrho(ncom), &
            dtang(ncom),dtdtang(ncom),dfdtang(ncom),dfang(ncom), &
            dfdfang(ncom),hrho(3,3)
       real*8 r,invr,invr2,cth,sth

      invr=1.D0/r
      invr2=invr*invr
!     Was a bug here
!      do i=0,10
!       do j=0,10
      do i=1,3
        do j=1,3
         hrho(i,j)=0.0d0
       enddo
      enddo
      do i=0,10
       do j=0,10
         c_kub(i,j)=0.0d0
       enddo
      enddo
      c_kub(0,0)=1.d0
      c_kub(3,2)=1.d0
      c_kub(4,0)=.5d0*SQRT(7.d0/3.d0)
      c_kub(4,4)=.5*SQRT(5.d0/3.d0)
      c_kub(6,0)=.5d0*SQRT(.5d0)
      c_kub(6,2)=.25d0*SQRT(11.d0)
      c_kub(6,4)=-.5d0*SQRT(7.d0/2.d0)
      c_kub(6,6)=-.25d0*SQRT(5.d0)
      c_kub(7,2)=.5d0*SQRT(13.d0/6.d0)
      c_kub(7,6)=.5d0*SQRT(11.d0/6.d0)
      c_kub(8,0)=.125d0*SQRT(33.d0)
      c_kub(8,4)=.25d0*SQRT(7.d0/3.d0)
      c_kub(8,8)=.125d0*SQRT(65.d0/3.d0)
      c_kub(9,2)=.25d0*SQRT(3.d0)
      c_kub(9,4)=.5d0*SQRT(17.d0/6.d0)
      c_kub(9,6)=-.25d0*SQRT(13.d0)
      c_kub(9,8)=-.5d0*SQRT(7.d0/6.d0)
      c_kub(10,0)=.125*SQRT(65.D0/6.D0)
      c_kub(10,2)=.125*SQRT(247.D0/6.D0)
      c_kub(10,4)=-.25*SQRT(11.D0/2.D0)
      c_kub(10,6)=0.0625d0*SQRT(19.D0/3.D0)
      c_kub(10,8)=-.125*SQRT(187.D0/6.D0)
      c_kub(10,10)=-.0625d0*SQRT(85.d0)


       i=1
 1     CONTINUE
          IF(i.gt.lmmax) GOTO 4
          IF(lm(1,i,jatom).EQ.0.AND.lm(2,i,jatom).EQ.0) THEN
             fac=1.d0
             call sumdd_aux(r,cth,sth,rho,drrho,drdrrho,ang,dtang, &
                  dtdtang,dfdtang,dfang,dfdfang,hrho,lmmax,lm,jatom, &
                  i,i,fac)
             i=i+1
          ELSEIF (lm(1,i,jatom).EQ.-3.AND.lm(2,i,jatom).EQ.2) THEN  
             fac=1.d0
             call sumdd_aux(r,cth,sth,rho,drrho,drdrrho,ang,dtang, &
                  dtdtang,dfdtang,dfang,dfdfang,hrho,lmmax,lm,jatom, &
                  i,i,fac)
             i=i+1
          ELSEIF (lm(1,i,jatom).EQ.4.OR.lm(1,i,jatom).EQ.6.OR. &
                  lm(1,i,jatom).EQ.-7.OR.lm(1,i,jatom).EQ.-9) THEN  
            c1=c_kub(abs(lm(1,i,jatom)),lm(2,i,jatom))
            c2=c_kub(abs(lm(1,i,jatom)),lm(2,i,jatom)+4)
            fac=c1*c1
            ia=i
            ir=i
             call sumdd_aux(r,cth,sth,rho,drrho,drdrrho,ang,dtang, &
                  dtdtang,dfdtang,dfang,dfdfang,hrho,lmmax,lm,jatom, &
                  ia,ir,fac)
            fac=c1*c2
            ia=i
            ir=i+1
             call sumdd_aux(r,cth,sth,rho,drrho,drdrrho,ang,dtang, &
                  dtdtang,dfdtang,dfang,dfdfang,hrho,lmmax,lm,jatom, &
                  ia,ir,fac)
            fac=c2*c1
            ia=i+1
            ir=i
             call sumdd_aux(r,cth,sth,rho,drrho,drdrrho,ang,dtang, &
                  dtdtang,dfdtang,dfang,dfdfang,hrho,lmmax,lm,jatom, &
                  ia,ir,fac)
            fac=c2*c2
            ia=i+1
            ir=i+1
             call sumdd_aux(r,cth,sth,rho,drrho,drdrrho,ang,dtang, &
                  dtdtang,dfdtang,dfang,dfdfang,hrho,lmmax,lm,jatom, &
                  ia,ir,fac)
             i=i+2
          ELSEIF (lm(1,i,jatom).EQ.8.OR.lm(1,i,jatom).EQ.10) THEN 
            c1=c_kub(abs(lm(1,i,jatom)),lm(2,i,jatom))
            c2=c_kub(abs(lm(1,i,jatom)),lm(2,i,jatom)+4)
            c3=c_kub(abs(lm(1,i,jatom)),lm(2,i,jatom)+8)
            fac=c1*c1
            ia=i
            ir=i
             call sumdd_aux(r,cth,sth,rho,drrho,drdrrho,ang,dtang, &
                  dtdtang,dfdtang,dfang,dfdfang,hrho,lmmax,lm,jatom, &
                  ia,ir,fac)
            fac=c1*c2
            ia=i
            ir=i+1
             call sumdd_aux(r,cth,sth,rho,drrho,drdrrho,ang,dtang, &
                  dtdtang,dfdtang,dfang,dfdfang,hrho,lmmax,lm,jatom, &
                  ia,ir,fac)
            fac=c1*c3
            ia=i
            ir=i+2
             call sumdd_aux(r,cth,sth,rho,drrho,drdrrho,ang,dtang, &
                  dtdtang,dfdtang,dfang,dfdfang,hrho,lmmax,lm,jatom, &
                  ia,ir,fac)
            fac=c2*c1
            ia=i+1
            ir=i
             call sumdd_aux(r,cth,sth,rho,drrho,drdrrho,ang,dtang, &
                  dtdtang,dfdtang,dfang,dfdfang,hrho,lmmax,lm,jatom, &
                  ia,ir,fac)
            fac=c2*c2
            ia=i+1
            ir=i+1
             call sumdd_aux(r,cth,sth,rho,drrho,drdrrho,ang,dtang, &
                  dtdtang,dfdtang,dfang,dfdfang,hrho,lmmax,lm,jatom, &
                  ia,ir,fac)
            fac=c2*c3
            ia=i+1
            ir=i+2
             call sumdd_aux(r,cth,sth,rho,drrho,drdrrho,ang,dtang, &
                  dtdtang,dfdtang,dfang,dfdfang,hrho,lmmax,lm,jatom, &
                  ia,ir,fac)
            fac=c3*c1
            ia=i+2
            ir=i
             call sumdd_aux(r,cth,sth,rho,drrho,drdrrho,ang,dtang, &
                  dtdtang,dfdtang,dfang,dfdfang,hrho,lmmax,lm,jatom, &
                  ia,ir,fac)
            fac=c3*c2
            ia=i+2
            ir=i+1
             call sumdd_aux(r,cth,sth,rho,drrho,drdrrho,ang,dtang, &
                  dtdtang,dfdtang,dfang,dfdfang,hrho,lmmax,lm,jatom, &
                  ia,ir,fac)
            fac=c3*c3
            ia=i+2
            ir=i+2
             call sumdd_aux(r,cth,sth,rho,drrho,drdrrho,ang,dtang, &
                  dtdtang,dfdtang,dfang,dfdfang,hrho,lmmax,lm,jatom, &
                  ia,ir,fac)
             i=i+3
          ELSE
             WRITE(6,*) 'UNCORRECT LM LIST FOR CUBIC STRUCTURE'
             WRITE(6,'(a2,i4,a3,i4)') 'L=',lm(1,i,jatom), &
                   ' M=',lm(2,i,jatom)
             STOP
          ENDIF
       GOTO 1
 4     CONTINUE 
       END


      subroutine sumdd_aux(r,cth,sth,rho,drrho,drdrrho,ang, &
         dtang,dtdtang,dfdtang,dfang,dfdfang,hrho,lmmax,lm, &
           jatom,ia,ir,fac)
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INTEGER  jatom,lmmax,lm(2,ncom,*)
      REAL*8   fac
      INTEGER  ia,ir
      real*8 rho(ncom),ang(ncom),drrho(ncom),drdrrho(ncom), &
           dtang(ncom),dtdtang(ncom),dfdtang(ncom),dfang(ncom), &
           dfdfang(ncom),hrho(3,3)
       real*8 r,invr,invr2,cth,sth
      
      invr=1.D0/r
      invr2=invr*invr
      hrho(1,1) = hrho(1,1)+fac*(drdrrho(ir)*ang(ia))
      hrho(1,2) = hrho(1,2)+fac*(-rho(ir)*dtang(ia)*invr2+ &
           drrho(ir)*dtang(ia)*invr)
      hrho(1,3) = hrho(1,3)+fac*((-rho(ir)*dfang(ia)*invr2+ &
           drrho(ir)*dfang(ia)*invr)/sth)
      hrho(2,2) = hrho(2,2)+fac*(rho(ir)*dtdtang(ia)*invr2+ &
           drrho(ir)*ang(ia)*invr)
      hrho(2,3) = hrho(2,3)+fac*(rho(ir)*(-cth*dfang(ia)+ &
           sth*dfdtang(ia))/(sth*sth)*invr2)
      hrho(3,3) = hrho(3,3)+fac*(rho(ir)*(dfdfang(ia)+ &
           cth*sth*dtang(ia))/(sth*sth)*invr2+ &
           drrho(ir)*ang(ia)*invr)

      return
      end

