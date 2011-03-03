      SUBROUTINE CHARGE(CHG,IR,R,JATOM,IATNR)                         
      use blank
      use rad
      IMPLICIT NONE

      INCLUDE 'param.inc'

      COMPLEX*16 YL,dtyl,dtdtyl,dfdtyl,IMAG,IMAG1

      real*8 minu,rho(ncom),ang(ncom),sqin2,r,chg
      real*8 rhop,rhopp

      integer lmmx,ilm,l,m,idm,idp,iatnr,ir,jatom,idx,jrj

      logical sy,syp,sypp
      
      COMMON /YLMS/ yl((lmax2+4)*(lmax2+4)),dtyl((lmax2+3)*(lmax2+3)), &
         dtdtyl((lmax2+2)*(lmax2+2)),dfdtyl((lmax2+2)*(lmax2+2))
!     Small changes by L. D. Marks, January 2006

      parameter (imag=(0.0D0,1.0D0), SQIN2=0.707106781186547D0)                    
      jrj=jri(jatom)
      sy=.true.
      syp=.false.
      sypp=.false.

      LMMX=LMMAX(JATOM)                                                 
      DO  ILM=1,LMMX                                                  
!        CALL RADIAL(CLM(1,ILM,JATOM),RHO(ILM),R,IR,JATOM)
        call splint(rm(1,jatom),clm(1,ilm,jatom),clm2(1,ilm,jatom), &
           jrj,r,ir,rho(ilm),rhop,rhopp,sy,syp,sypp)
        rho(ilm)=rho(ilm)/(r*r)
        L=IABS(LM(1,ILM,JATOM))
        M=LM(2,ILM,JATOM)
        MINU=1.D0
        IMAG1=(1.D0,0.D0)
        IF(LM(1,ILM,JATOM).LT.0) THEN
          IMAG1=-IMAG
          MINU=-1.D0
        END IF
        IF(MOD(M,2).EQ.1) THEN
          IMAG1=-IMAG1
          MINU=-MINU
        END IF
        IF(M.NE.0) GOTO 2
        IDX=L*(L+1)+1
        ANG(ILM)=YL(IDX)
        GOTO 1
 2      IDM=L*(L+1)+1
        IDP=IDM+M
        IDM=IDM-M
        ANG(ILM)=(YL(IDP)+MINU*YL(IDM))*SQIN2*IMAG1
!        write(6,*) 'ang(ilm) ',ilm,idp,yl(idp),idm,yl(idm),ang(ilm)
 1      CONTINUE                                                          
      end do
!     
      IF(IATNR.GT.0) CALL SUM(RHO,ANG,CHG,LMMX,lm,jatom)                    
      IF(IATNR.GT.0) RETURN                                             
!     
      CHG=0.0                                                           
      DO  ILM=1,LMMX 
        CHG=CHG+RHO(ILM)*ANG(ILM)
!     write(6,*) ilm,chg,rho(ilm),ang(ilm)
      end do
      if (chg.lt.0.d0) then
        write (6,*) 'NEGATIVE CHARGE!!! r=',r,' COMPONENTS:'
        write (6,*) 'ilm  L  M  ang(ilm) rho(ilm)'
        do ilm=1,lmmx
          write (6,*) ilm,IABS(LM(1,ILM,JATOM)),LM(2,ILM,JATOM), &
             ang(ilm),rho(ilm)
        end do
      end if
      RETURN
      END                                                               

