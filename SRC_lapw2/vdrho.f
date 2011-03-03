      subroutine vdrho(lmmax,lm,jri,dx,r,vlm,clmnew,fvdrho)
      use param
      IMPLICIT REAL*8 (A-H,O-Z)

!      parameter (lmm0=6)

      REAL*8 I1,I2,MINUS,INTFELD1,INTFELD2   

      LOGICAL hl1,hl2,hl3,hl4
      dimension lm(2,0:NCOM-1),clmnew(NRAD,0:NCOM-1),r(NRAD)
      dimension vlm(NRAD,0:NCOM-1)
      dimension fvdrho(3)
!      DIMENSION intfeld1(0:lmm0,0:lmm0,-1:1,0:lmm0,0:lmm0,-1:1)
!      DIMENSION intfeld2(0:lmm0,0:lmm0,-1:1,0:lmm0,0:lmm0,-1:1) 
      DIMENSION intfeld1(0:lmax2,0:lmax2,-1:1,0:lmax2,0:lmax2,-1:1)
      DIMENSION intfeld2(0:lmax2,0:lmax2,-1:1,0:lmax2,0:lmax2,-1:1) 
      dimension value(NRAD),derv(NRAD)
      PARAMETER (ONE=1.D0)
      PARAMETER (TWO=2.D0)
      PARAMETER (THREE=3.D0) 
      PARAMETER (HALF=0.5D0)
      PARAMETER (MINUS=-1.D0)
!*******************Inline FUNKTIONEN***************************

       D2(L) = ((TWO*L+THREE)*(TWO*L+ONE))
       VFA(L,M) = SQRT((L+M+ONE)*(L+M+TWO)/D2(L))
       VFC(L,M) = SQRT((L-M+ONE)*(L-M+TWO)/D2(L))
       VFF(L,M) = SQRT((L+M+ONE)*(L-M+ONE)/D2(L))

      SQHALF=SQRT(HALF) 
!*****************integralfelder FUER LM-KOMBINATION ANLEGEN********  
!      do jatom=1,nat  
      jatom=1    
      forvrz=0.D0
      forvry=0.D0
      forvrx=0.D0

      DO 5000 indexl=0,lmax2
         DO 5001 indexm=0,lmax2
            DO 5002 indexk=-1,+1
               DO 5003 indexlh=0,lmax2
                 DO 5004 indexmh=0,lmax2
                    DO 5005 indexkh=-1,+1
      intfeld1(indexl,indexm,indexk,indexlh,indexmh,indexkh)=0.D0
      intfeld2(indexl,indexm,indexk,indexlh,indexmh,indexkh)=0.D0 
 5005              CONTINUE
 5004            CONTINUE
 5003         CONTINUE
 5002       CONTINUE
 5001     CONTINUE
 5000  CONTINUE
        la8max=iabs(lm(1,lmmax-1))
        lmmaxns=lmmax-1
        DO 5006 lm1=0,lmmaxns-1 
            indexl=abs(lm(1,lm1)) 
            indexm=lm(2,lm1)
            indexk=1
            if (lm(1,lm1).lt.0) indexk=-1 
            lplone=indexl+1
            mmione=indexm-1
            mplone=indexm+1
          DO 5007 lm11=lm1+1,lmmaxns 
            indexkh=1
            if (lm(1,lm11).lt.0) indexkh=-1
            hl1=(abs(lm(1,lm11)).eq.lplone) 
            hl2=(lm(2,lm11).eq.indexm)
            hl3=(lm(2,lm11).eq.mmione)
            hl4=(lm(2,lm11).eq.mplone)
        if (hl1.and.(hl2.or.hl3.or.hl4)) then
                
            DO 811 J=1,JRI
              VALUE(J)=CLMNEW(J,lm11)/r(j)/r(j)
 811        CONTINUE 
            call dfrad(r,value,derv,jri) 
!            CALL DERIV(R,DX,
!     &                 VALUE,1,JRI,DERV)
 
            DO 812 J=1,JRI          
            VALUE(J)=VLM(J,lm1)*CLMNEW(J,lm11) &
                     /R(J)
 812        CONTINUE  

                call charge(r,1.d0,1.d0,value,dx,jri,i2)

            DO 813 J=1,JRI
             VALUE(J)=VLM(J,lm1)*R(J)*R(J)*DERV(J)
 813        CONTINUE

                call charge(r,1.d0,1.d0,value,dx,jri,i1)

            intmindex=lm(2,lm11)
      intfeld1(indexl,indexm,indexk,lplone,intmindex,indexkh)=i1 
      intfeld2(indexl,indexm,indexk,lplone,intmindex,indexkh)=i2
!      write(*,*)lm1,lm11,indexl,indexm,indexk,lplone,intmindex,indexkh

!************************************************************** 
     

            DO 814 J=1,JRI
              VALUE(J)=CLMNEW(J,lm1)/r(j)/r(j)
 814        CONTINUE 

            call dfrad(r,value,derv,jri) 
!            CALL DERIV(R,DX,
!     &                 VALUE,1,JRI,DERV)
 
            DO 815 J=1,JRI          
              VALUE(J)=VLM(J,lm11)*CLMNEW(J,lm1) &
                       /R(J)
 815       CONTINUE  

                 call charge(r,1.d0,1.d0,value,dx,jri,i2)

            DO 816 J=1,JRI
             VALUE(J)=VLM(J,lm11) &
                      *R(J)*R(J)*DERV(J)
 816        CONTINUE

                 call charge(r,1.d0,1.d0,value,dx,jri,i1)

      intfeld1(lplone,intmindex,indexkh,indexl,indexm,indexk)=i1 
      intfeld2(lplone,intmindex,indexkh,indexl,indexm,indexk)=i2
!      write(*,*)lm1,lm11,lplone,intmindex,indexkh,indexl,indexm,indexk

             
!****************************************************************

       endif   
 5007 CONTINUE
 5006 CONTINUE 
 

      DO 5008 l=0,lA8max-1
       DO 5009 mz=-l,l 

       m=abs(mz)                               

       
       fakt=half
       if (m.eq.0) fakt=one   

       forvrz=forvrz+vff(l,mz)*fakt* &
          (  (two+l)*(intfeld2(l,m,+1,l+1,m,+1)+ &
                      intfeld2(l,m,-1,l+1,m,-1)) &
                  -l*(intfeld2(l+1,m,+1,l,m,+1)+ &
                      intfeld2(l+1,m,-1,l,m,-1)) &
                    +(intfeld1(l,m,+1,l+1,m,+1)+  &
                      intfeld1(l,m,-1,l+1,m,-1)) &
                    +(intfeld1(l+1,m,+1,l,m,+1)+ &
                      intfeld1(l+1,m,-1,l,m,-1))   )  


       faktmo=half
       faktpo=minus*half 
  
       if (mz.eq.0) then 
        faktpo=minus*sqhalf
        faktmo=sqhalf
       endif
        
       if (mz.gt.1) faktmo=minus*half
       if (mz.eq.1) faktmo=minus*sqhalf

       if (mz.lt.(-1)) faktpo=half
       if (mz.eq.(-1)) faktpo=sqhalf

       forvrx=forvrx-vfa(l,mz)*half*faktpo* &
        (  (two+l)*(intfeld2(l,abs(mz),+1,l+1,abs(mz+1),+1)+ &
                    intfeld2(l,abs(mz),-1,l+1,abs(mz+1),-1)) &
                 -l*(intfeld2(l+1,abs(mz+1),+1,l,abs(mz),+1)+ &
                     intfeld2(l+1,abs(mz+1),-1,l,abs(mz),-1)) &
                   +(intfeld1(l,abs(mz),+1,l+1,abs(mz+1),+1)+ &
                     intfeld1(l,abs(mz),-1,l+1,abs(mz+1),-1)) &
                   +(intfeld1(l+1,abs(mz+1),+1,l,abs(mz),+1)+ &
                     intfeld1(l+1,abs(mz+1),-1,l,abs(mz),-1))   )  

       forvrx=forvrx+vfc(l,mz)*half*faktmo* &
         (  (two+l)*(intfeld2(l,abs(mz),+1,l+1,abs(mz-1),+1)+ &
                    intfeld2(l,abs(mz),-1,l+1,abs(mz-1),-1)) &
                 -l*(intfeld2(l+1,abs(mz-1),+1,l,abs(mz),+1)+ &
                    intfeld2(l+1,abs(mz-1),-1,l,abs(mz),-1)) &
                   +(intfeld1(l,abs(mz),+1,l+1,abs(mz-1),+1)+ &
                     intfeld1(l,abs(mz),-1,l+1,abs(mz-1),-1)) &
                   +(intfeld1(l+1,abs(mz-1),+1,l,abs(mz),+1)+ &
                     intfeld1(l+1,abs(mz-1),-1,l,abs(mz),-1))   ) 


       faktmo=half
       if (mz.eq.1) faktmo=sqhalf
       if (mz.eq.0) faktmo=sqhalf

 
       forvry= forvry+vfc(l,mz)*half*faktmo* &
          (  (l+two)*(intfeld2(l,abs(mz),+1,l+1,abs(mz-1),-1)- &
                      intfeld2(l,abs(mz),-1,l+1,abs(mz-1),+1)) &
                    +(intfeld1(l,abs(mz),+1,l+1,abs(mz-1),-1)- &
                      intfeld1(l,abs(mz),-1,l+1,abs(mz-1),+1)) )   
                                                             

       forvry= forvry+vfc(l,mz)*half*faktmo* &
          (  l*(intfeld2(l+1,abs(mz-1),+1,l,abs(mz),-1)- &
                intfeld2(l+1,abs(mz-1),-1,l,abs(mz),+1)) &
              -(intfeld1(l+1,abs(mz-1),+1,l,abs(mz),-1)- &
                intfeld1(l+1,abs(mz-1),-1,l,abs(mz),+1)) )


       faktpo=half
       if (mz.eq.(-1)) faktpo=sqhalf
       if (mz.eq.0) faktpo=sqhalf

       forvry= forvry+vfa(l,mz)*half*faktpo* &
          (  (l+two)*(intfeld2(l,abs(mz),+1,l+1,abs(mz+1),-1)- &
                      intfeld2(l,abs(mz),-1,l+1,abs(mz+1),+1)) &
                    +(intfeld1(l,abs(mz),+1,l+1,abs(mz+1),-1)- &
                      intfeld1(l,abs(mz),-1,l+1,abs(mz+1),+1)) )   
                                                             

       forvry= forvry+vfa(l,mz)*half*faktpo* &
          (  l*(intfeld2(l+1,abs(mz+1),+1,l,abs(mz),-1)- &
                intfeld2(l+1,abs(mz+1),-1,l,abs(mz),+1)) &
              -(intfeld1(l+1,abs(mz+1),+1,l,abs(mz),-1)- &
                intfeld1(l+1,abs(mz+1),-1,l,abs(mz),+1)) )
       
!       write(*,*) l,mz,forvrx,forvry,forvrz
 5009  CONTINUE
 5008  CONTINUE  
      
       fvdrho(1)=forvrx
       fvdrho(2)=forvry
       fvdrho(3)=forvrz
!      enddo
       return
       end
