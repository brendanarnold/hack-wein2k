      subroutine ykav(j,jatom,yka,lmmtmx,lmmult,lmmax1)
!
!.... performs star average for ylm
!
      use struct
      use symetr
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
      complex*16 yka(ncom+3),taup(nsym),cphs
      complex*16 yll(nsym,(lmax2+1)*(lmax2+1))
      integer          lmmtmx(nat),lmmult(2,ncom+3,nat)
      dimension av(3),rotv2(3),rotv3(3,nsym)
!
      pi=acos(-1.d0)
      if(j.eq.1) then
         pi=acos(-1.d0)
         sqfp=sqrt(4.d0*pi)
         sqfpi=1.d0/sqfp
         do 37 lmx=1,ncom+3
 37      yka(lmx)=dcmplx(sqfpi,0.d0)
         return
      end if
!
      index=1
      do 40 latom=2,jatom
 40   index=index+mult(latom-1)
      CALL STERN(J,IND,KKK,TAUP)                                     
      DO 41 LM1=1,NCOM+3                                             
      YKA(LM1)=(0.0D0,0.0D0)                               
 41   CONTINUE                                                        
          LLMM=lmmtmx(JATOM)                                    
          if(iatnr(jatom).gt.0) llmm=lmmax1
!          write(6,*) 'ykav',llmm
!          if(iatnr(jatom).gt.0) then
!             if(lmmax(jatom).eq.6) then
!             llmm=lmmax(jatom)+3
!             else if(lmmax(jatom).eq.1) then
!             llmm=1
!             else
!             llmm=lmmax(jatom)+2
!             end if
!          end if
      DO 42 JJ=1,IND                                                  
          AV(1)=KKK(1,JJ)                                               
          AV(2)=KKK(2,JJ)                                               
          AV(3)=KKK(3,JJ)                                               
          ROTV2(1)=BR1(1,1)*av(1)+BR1(1,2)*av(2)+BR1(1,3)*av(3)   
          ROTV2(2)=BR1(2,1)*av(1)+BR1(2,2)*av(2)+BR1(2,3)*av(3)   
          ROTV2(3)=BR1(3,1)*av(1)+BR1(3,2)*av(2)+BR1(3,3)*av(3)   
          CALL ROTATE (ROTV2,ROTLOC(1,1,JATOM),ROTV3(1,jj))        
  42   continue

          CALL YLM_3(ind,ROTV3,nsym,lmax2,YLL,nsym)                              
          pi2=2.d0*pi          
          do 43 jj=1,ind

          ARG1=KKK(1,JJ)*POS(1,INDEX)                        
          ARG2=KKK(2,JJ)*POS(2,INDEX)                        
          ARG3=KKK(3,JJ)*POS(3,INDEX)    
          arg1=(arg1+arg2+arg3)*pi2                    
!          CPHS=EXP( (0.d0,1.d0)*2.d0*PI*(ARG1+ARG2+ARG3) )            
          CPHS=DCMPLX(cos(arg1),sin(arg1))*taup(jj)/ind               
          DO 942 LM1=1,llmm                                      
             LL=IABS( lmmult(1,LM1,jatom) )                         
             MM=lmmult(2,LM1,JATOM)                                 
             LTEST=LL*(LL+1)+MM+1                               
             YKA(LM1)=YKA(LM1) + CPHS*           &
             CONJG( YLL(jj,LTEST) )   
 942      continue    
  43   CONTINUE                                                       
!
      return
      end

      SUBROUTINE YLM_3(n,V,ldv,LMAX,Y,ldy)
      use struct, only : nsym
      INCLUDE 'param.inc'
      INTEGER            LMAX
      DOUBLE PRECISION   V(3,ldv)
      COMPLEX*16         Y(ldy,(lmax+1)**2),hy
!
      INTEGER            I2L, I4L2, INDEX, INDEX2, L, M, MSIGN
      DOUBLE PRECISION   A, B, C, AB, ABC, ABMAX, ABCMAX
      DOUBLE PRECISION   D4LL1C, D2L13, PI
      DOUBLE PRECISION   COSTh(nsym),SINTh(nsym),COSPh(nsym),SINPh(nsym)
      DOUBLE PRECISION   TEMP1, TEMP2, TEMP3
      DOUBLE PRECISION   YLLR, YLL1R, YL1L1R, YLMR, YLLRT
      DOUBLE PRECISION   YLLI, YLL1I, YL1L1I, YLMI
!
      PI = (4.0D+0)*ATAN(1.0D+0)
      rt3=sqrt(3.d0)
      rt15=sqrt(1.5d0)
!
!        Y(0,0)
!
      YLLR = 1.0D+0/SQRT(4.0D+0*PI)
      YLLI = 0.0D+0
      hy   = DCMPLX(YLLR,YLLI)
!
!        continue only if spherical harmonics for (L .GT. 0) are desired
!
      IF (LMAX .LE. 0)  then
         do i=1,n
            y(i,1)=hy
         enddo
         goto 999
      endif
!
!        calculate sin(Phi), cos(Phi), sin(Theta), cos(Theta)
!        Theta, Phi ... polar angles of vector V
!
      call cart_2_kugel(n,v,ldv,costh,sinth,cosph,sinph)
!
!        Y(1,0)
!
      yllrt=rt3*yllr
      do i2=1,n
      y(i2,1)=hy
!      y(i2,3) = DCMPLX(SQRT(3.0D+0)*YLLR*COSTh(i2),0.0D+0)
      y(i2,3) = DCMPLX(YLLRT*COSTh(i2),0.0D+0)
!
!        Y(1,1) ( = -DCONJG(Y(1,-1)))
!
      TEMP1 = -RT15*YLLR*SINTh(i2)
      y(i2,4) = DCMPLX(TEMP1*COSPh(i2),TEMP1*SINPh(i2))
      y(i2,2) = -DCONJG(y(i2,4))
      enddo
!
      DO 20 L = 2, LMAX
         INDEX  = L*L+1
         INDEX2 = INDEX + 2*L
         MSIGN  = 1 - 2*MOD(L,2)
!
!        YLL = Y(L,L) = f(Y(L-1,L-1)) ... Formula 1
!
         TEMP1 = -SQRT(DBLE(2*L+1)/DBLE(2*L))
      do i3=1,n
         YL1L1R = DBLE(y(i3,INDEX-1))
         YL1L1I = DIMAG(y(i3,INDEX-1))
!        TEMP1 = -SQRT(DBLE(2*L+1)/DBLE(2*L))*SINTh(i3)
         YLLR = TEMP1*SINTh(i3)*(COSPh(i3)*YL1L1R - SINPh(i3)*YL1L1I)
         YLLI = TEMP1*SINTh(i3)*(COSPh(i3)*YL1L1I + SINPh(i3)*YL1L1R)
         y(i3,INDEX2) = DCMPLX(YLLR,YLLI)
         y(i3,INDEX)  = MSIGN*DCONJG(y(i3,INDEX2))
      enddo
         INDEX2 = INDEX2 - 1
         INDEX  = INDEX  + 1
!
!        YLL1 = Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
!               (the coefficient for Y(L-2,L-1) in Formula 2 is zero)
!
         TEMP2 = SQRT(DBLE(2*L+1))
      do i4=1,n
!        TEMP2 = SQRT(DBLE(2*L+1))*COSTh(i4)
         YL1L1R = DBLE(y(i4,INDEX-2))
         YL1L1I = DIMAG(y(i4,INDEX-2))
         YLL1R = TEMP2*COSTh(i4)*YL1L1R
         YLL1I = TEMP2*COSTh(i4)*YL1L1I
         y(i4,INDEX2) = DCMPLX(YLL1R,YLL1I)
         y(i4,INDEX)  = -MSIGN*DCONJG(y(i4,INDEX2))
      enddo
         INDEX2 = INDEX2 - 1
         INDEX  = INDEX  + 1
!
         I4L2 = INDEX2 - 4*L + 2
         I2L  = INDEX2 - 2*L
         D4LL1C = SQRT(DBLE(4*L*L-1))
!        D4LL1C = COSTh(i3)*SQRT(DBLE(4*L*L-1))
         D2L13  = -SQRT(DBLE(2*L+1)/DBLE(2*L-3))
!
         DO 10 M = L - 2, 0, -1
!
!        YLM = Y(L,M) = f(Y(L-2,M),Y(L-1,M)) ... Formula 2
!
            TEMP1 = 1.0D+0/SQRT(DBLE((L+M)*(L-M)))
            TEMP2 = D4LL1C*TEMP1
            TEMP3 = D2L13*SQRT(DBLE((L+M-1)*(L-M-1)))*TEMP1
            do i5=1,n
            YLMR=TEMP2*COSTh(i5)*DBLE(y(i5,I2L))+TEMP3*DBLE(y(i5,I4L2))
!.......................................................................
           YLMI=TEMP2*COSTh(i5)*DIMAG(y(i5,I2L))+TEMP3*DIMAG(y(i5,I4L2))
            y(i5,INDEX2) = DCMPLX(YLMR,YLMI)
            y(i5,INDEX)  = MSIGN*DCONJG(y(i5,INDEX2))
            enddo
!
            MSIGN  = -MSIGN
            INDEX2 = INDEX2 - 1
            INDEX  = INDEX  + 1
            I4L2   = I4L2   - 1
            I2L    = I2L    - 1
   10    CONTINUE
   20 CONTINUE
!
  999 RETURN
!
!        End of 'YLM'
!
      END

      subroutine cart_2_kugel(n,v,ldv,costh,sinth,cosph,sinph)
      real*8 costh(n),sinth(n),cosph(n),sinph(n)
      real*8 v(3,ldv)
      real*8 abmax,ab,a,b,abcmax,c,abc
      do i1=1,n
      ABMAX  = MAX(ABS(v(1,i1)),ABS(v(2,i1)))
      IF (ABMAX .GT. 0.0D+0) THEN
         abmax=1.d0/abmax
         A = v(1,i1)*ABMAX
         B = v(2,i1)*ABMAX
         AB = 1.d0/SQRT(A*A+B*B)
         COSPh(i1) = A*AB
         SINPh(i1) = B*AB
      ELSE
         COSPh(i1) = 1.0D+0
         SINPh(i1) = 0.0D+0
      ENDIF
      ABCMAX = MAX(ABMAX,ABS(v(3,i1)))
      IF (ABCMAX .GT. 0.0D+0) THEN
         abcmax=1.d0/abcmax
         A = v(1,i1)*ABCMAX
         B = v(2,i1)*ABCMAX
         C = v(3,i1)*ABCMAX
         AB = A*A + B*B
         ABC = 1.d0/SQRT(AB + C*C)
         COSTh(i1) = C*ABC
         SINTh(i1) = SQRT(AB)*ABC
      ELSE
         COSTh(i1) = 1.0D+0
         SINTh(i1) = 0.0D+0
      ENDIF
      enddo
      return
      end
