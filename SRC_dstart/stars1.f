      subroutine stars1(index,cq)

! symmetrisiert die fourierkomponente der ladungsdichte
!  ueber den stern von G

      use struct
      implicit real*8(a-h,o-z)
!	include 'param.inc'
      
      complex*16   cq,imag,taup,chi(nsym)
      dimension    ind(nsym)
      integer      stg
!      COMMON /STRUK/  POS(3,ndif),AA,BB,CC,ALPHA(3),RMT(NATO),V(NATO),
!     *                PIA(3),VOL,IATNR(NATO),MULT(NATO),ISPLIT(NATO)       
!      COMMON /SYM2/   TAU(3,NSYM),IORD,IZ(3,3,NSYM)
      COMMON /FOUR/   iG(3),NST,STG(3,NSYM),TAUP(NSYM)


	cq=0.d0
      tpi= 2.d0 * acos(-1.d0)
      imag=(0.d0,1.d0)
      NST=0
      DO 1 I=1,IORD
         TK=0.
         tk1=0.
          DO 2 J=1,3
             tk1=tk1+tau(j,i)*ig(j)
            it=0
            DO 3 L=1,3
               it=IZ(J,L,I)*ig(l)+it
 3          CONTINUE
            STG(J,I)=it
 2       CONTINUE
         IF(NST.EQ.0) GOTO 7
         DO 4 M=1,NST
            DO 5 J=1,3
               IF(STG(J,M).NE.STG(J,I)) GOTO 4
 5          CONTINUE
!          taup(M)=taup(m)+exp(imag*tpi*tk1)
          taup(M)=taup(M)+DCMPLX(cos(tk1*tpi),sin(tk1*tpi))
          ind(M)=ind(m)+1
          goto 1
 4       CONTINUE
 7       NST=NST+1
         tk=0.
         DO 6 J=1,3
            STG(J,NST)=STG(J,I)
            tk=tk+stg(j,nst)*pos(j,index)
 6       CONTINUE
         ind(nst)=1
         taup(nst)=dcmplx(cos(tpi*tk1),sin(tpi*tk1))
         chi(nst)=dcmplx(cos(tk*tpi),sin(tk*tpi))
!         taup(nst)=exp(imag*tpi*tk1)
!         chi(nst)=exp(imag*tk*tpi)
 1    CONTINUE

      do inst=1,nst
         cq=cq +chi(inst)*taup(inst)/ind(inst)
      enddo
      return
      end
