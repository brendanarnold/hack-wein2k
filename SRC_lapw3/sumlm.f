SUBROUTINE SUMLM(LMMAX,TCC,YL,F,FA,lm,FFX)
  ! Calculates cubic harmonics after Kara & Kurki-Suonio
  ! Acta Cryst A 1981 37 201-210
  ! GKHM 26/11-02
  IMPLICIT NONE
  INCLUDE 'param.inc'
  INTEGER      :: lmmax,lm(2,ncom)
  INTEGER      :: i,ll,mm,ll2,mm2,ll3,mm3
  INTEGER      :: ly1,ly2,ly21,ly22,ly31,ly32
  REAL*8      :: c_kub(0:10,0:10)
  REAL*8      :: tcc(ncom),ffx(ncom),f,fa,f1,minu,fa1,c1,c2,c3
  COMPLEX*16  :: YL((lmax2+1)*(lmax2+1)),ANG,IMAG,f2,fac,ang1,ang2,ang3

  IMAG=(0.0d0,1.0d0)                                            
  f=0.0d0
  FA=0.0d0                                                           
  ffx=0.0d0

  c_kub(1:10,1:10)=0.0d0
  c_kub(0,0)=1.d0
  c_kub(3,2)=1.d0
  c_kub(4,0)=.5d0*SQRT(7.d0/3.d0)
  c_kub(4,4)=.5d0*SQRT(5.d0/3.d0)
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
  c_kub(10,0)=.125d0*SQRT(65.D0/6.D0)
  c_kub(10,2)=.125d0*SQRT(247.D0/6.D0)
  c_kub(10,4)=-.25d0*SQRT(11.D0/2.D0)
  c_kub(10,6)=0.0625d0*SQRT(19.D0/3.D0)
  c_kub(10,8)=-.125d0*SQRT(187.D0/6.D0)
  c_kub(10,10)=-.0625d0*SQRT(85.d0)

  i=1
  DO WHILE (i<=lmmax)
     !1 CONTINUE
     !  IF(i.gt.lmmax) GOTO 4
     ll=ABS(lm(1,i))
     mm=lm(2,i)
     ly1=ll*(ll+1)+mm+1
     ly2=ll*(ll+1)-mm+1
     fac=1.d0/SQRT(2.d0) ! No "(-1.d0)**mm" because mm always even
     minu=1.d0
     IF(lm(1,i)<0) THEN
        fac=fac*(-imag)
        minu=-1.d0
     ENDIF
     IF(ll.EQ.0.AND.mm.EQ.0) THEN
        f1=TCC(i)*YL(1)
        f=f+f1
        ffx(i)=f1
        i=i+1
     ELSEIF (ll.EQ.3.AND.mm.EQ.2) THEN  
        ANG=fac*(YL(ly1)+minu*YL(ly2))
        fa1=TCC(i)*ANG
        FA=fa+fa1                                                     
        FFX(i)=FA1                                                         
        i=i+1
     ELSEIF (ll==4.OR.ll==6.OR.ll==7.OR.ll==9) THEN
        ll2=ABS(lm(1,i+1)) ; mm2=lm(2,i+1) 
        ly21=ll2*(ll2+1)+mm2+1 ; ly22=ll2*(ll2+1)-mm2+1
        c1=c_kub(ll,mm) ; c2=c_kub(ll2,mm2)

        IF(mm==0) THEN
           ANG1=YL(ly1)
        ELSE
           ANG1=fac*(YL(ly1)+minu*YL(ly2)) 
        ENDIF
        ANG2=fac*(YL(ly21)+minu*YL(ly22))

        F1=c1*TCC(i)+c2*TCC(i+1)
        F2=c1*ang1+c2*ANG2
        IF(lm(1,i)>0) THEN
           F=F+F1*F2                                                         
        ELSE
           FA=FA+F1*F2 
        ENDIF
        FFX(i)=F1*F2                                                      
        FFX(i+1)=0.d0
        i=i+2
     ELSEIF (ll==8.OR.ll==10) THEN 
        ll2=ABS(lm(1,i+1)) ; mm2=lm(2,i+1) 
        ll3=ABS(lm(1,i+2)) ; mm3=lm(2,i+2) 
        ly21=ll2*(ll2+1)+mm2+1 ; ly22=ll2*(ll2+1)-mm2+1
        ly31=ll3*(ll3+1)+mm3+1 ; ly32=ll3*(ll3+1)-mm3+1

        c1=c_kub(ll,mm) ; c2=c_kub(ll2,mm2) ; c3=c_kub(ll3,mm3)
        IF(mm==0) THEN
           ANG1=YL(ly1)
        ELSE
           ANG1=fac*(YL(ly1)+minu*YL(ly2)) 
        ENDIF
        ANG2=fac*(YL(ly21)+minu*YL(ly22)) 
        ANG3=fac*(YL(ly31)+minu*YL(ly32))

        f1 = c1*tcc(i) + c2*tcc(i+1) + c3*tcc(i+2)
        f2 = c1*ang1 + c2*ang2 + c3*ang3 
        IF(lm(1,i)>0) THEN
           F=F+F1*F2                                                         
        ELSE
           FA=FA+F1*F2
        ENDIF
        FFX(i)=F1*F2                                                      
        FFX(i+1)=0.d0
        FFX(i+2)=0.d0
        i=i+3
     ELSE
        WRITE(6,*) 'UNCORRECT LM LIST FOR CUBIC SITE'
        WRITE(6,'(a2,i4,a3,i4)') 'L=',lm(1,i), &
             ' M=',lm(2,i)
        STOP 'UNCORRECT LM LIST FOR CUBIC SITE'
     ENDIF
  ENDDO
  !  GOTO 1
  !4 CONTINUE 


END SUBROUTINE SUMLM
