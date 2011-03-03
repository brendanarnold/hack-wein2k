SUBROUTINE lomain(nemin,nemax,lfirst,latom,n,jatom, &
     alm,blm,clm,aalm,bblm,cclm,force,forcea)
  !                                                                       
  USE param
  USE defs
  USE struk
  USE lo
  USE xa; USE xa3; USE com
  IMPLICIT NONE
  
  INTEGER, INTENT(in) :: nemin,nemax,lfirst,latom,n,jatom
  COMPLEX*16       YL((LMAX2+1)*(LMAX2+1)), &
       ALM((LMAX2+1)*(LMAX2+1),nume),               &
       BLM((LMAX2+1)*(LMAX2+1),nume), &
       cLM((LMAX2+1)*(LMAX2+1),nume,nloat)              
  COMPLEX*16       PHSHEL                     
  
  !     parameters, variables, and commons for force calculation
  complex*16 &
       ayp,byp,cyp &
       ,aalm(3,nume,(LMAX2+1)*(LMAX2+1)) &
       ,bblm(3,nume,(LMAX2+1)*(LMAX2+1)) &
       ,cclm(3,nume,(LMAX2+1)*(LMAX2+1),nloat)
  logical     force,forcea
  INTEGER             :: i,l,m,m1,jlo,jneq,num,index
  real*8              :: k(3),krot(3),twopi,arg1,arg2,arg3,argt
  
  !                                                                       
  !.initialise a,b,c of lo                                      
  TWOPI=TWO*PI
  i=nlov
  DO L=0,LoMAX
     !        if (.not.loor(l)) goto 10
     do jlo=1,ilo(l)
        do jneq=1,mult(jatom)
           DO M1=-l,+l 
              i=i+1
              BK(1)=BKXlo(I)
              BK(2)=BKYlo(I)
              BK(3)=BKZlo(I)
              CALL ROTATE (BK,ROTIJ(1,1,LATOM),BKROT)
              BK(1)=BKROT(1)*BR1(1,1)+BKROT(2)*BR1(1,2) &
                   +BKROT(3)*BR1(1,3)   
              BK(2)=BKROT(1)*BR1(2,1)+BKROT(2)*BR1(2,2) &
                   +BKROT(3)*BR1(2,3)   
              BK(3)=BKROT(1)*BR1(3,1)+BKROT(2)*BR1(3,2) &
                   +BKROT(3)*BR1(3,3)   
              if(force.and.forcea) then
                 k(1)=kxlo(i)
                 k(2)=kylo(i)
                 k(3)=kzlo(i)
                 call rotate(k,rotij(1,1,latom),krot)
                 k(1)=krot(1)*br1(1,1) &
                      +krot(2)*br1(1,2)+krot(3)*br1(1,3)
                 k(2)=krot(1)*br1(2,1) &
                      +krot(2)*br1(2,2)+krot(3)*br1(2,3)
                 k(3)=krot(1)*br1(3,1) &
                      +krot(2)*br1(3,2)+krot(3)*br1(3,3)
                 CALL ROTATE (K,ROTLOC(1,1,JATOM),KROt)
                 k(1)=krot(1)
                 k(2)=krot(2)
                 k(3)=krot(3)
              endif
              CALL ROTATE (BK,ROTLOC(1,1,JATOM),BKRLOC)
              CALL YLM (BKRLOC,LOMAX,YL)
              ARG1=BKROT(1)*POS(1,LFIRST)*TWOPI
              ARG2=BKROT(2)*POS(2,LFIRST)*TWOPI
              ARG3=BKROT(3)*POS(3,LFIRST)*TWOPI
              ARGT=(BKXlo(I)*TAUIJ(1,LATOM)+BKYlo(I)*TAUIJ(2,LATOM)+ &
                   BKZlo(I)*TAUIJ(3,LATOM))*TWOPI
              PHSHEL=EXP( IMAG*(ARG1+ARG2+ARG3+ARGT) )
              DO NUM=NEMIN,NEMAX
                 PHS(NUM)=PHSHEL * A_lo(I,NUM)         
              ENDDO
              DO M=-l,+l
                 index=l*(l+1)+m+1 
                 DO NUM=NEMIN,NEMAX
                    ALM(index,num)=ALM(INDEX,num)+ &
                         ALo(l,jlo)*conjg(YL(INDEX))*PHS(NUM)
                    BLM(INDEX,num)=BLM(INDEX,num)+ &
                         BLo(l,jlo)*conjg(YL(INDEX))*PHS(NUM)
                    clm(index,num,jlo)=clm(index,num,jlo)+ &
                         clo(l,jlo)*CONJG(yl(index))*phs(num)
                 ENDDO
                 if (force.and.forcea) then
                    do num=nemin,nemax
                       ayp=alo(l,jlo)*conjg(yl(index))*phs(num)
                       byp=blo(l,jlo)*conjg(yl(index))*phs(num)
                       cyp=clo(l,jlo)*conjg(yl(index))*phs(num)
                       
                       aalm(1,num,index)=aalm(1,num,index)+ayp*k(1)
                       aalm(2,num,index)=aalm(2,num,index)+ayp*k(2)
                       aalm(3,num,index)=aalm(3,num,index)+ayp*k(3)
                       
                       bblm(1,num,index)=bblm(1,num,index)+byp*k(1)
                       bblm(2,num,index)=bblm(2,num,index)+byp*k(2)
                       bblm(3,num,index)=bblm(3,num,index)+byp*k(3)
                       
                       cclm(1,num,index,jlo)=cclm(1,num,index,jlo)+cyp*k(1)
                       cclm(2,num,index,jlo)=cclm(2,num,index,jlo)+cyp*k(2)
                       cclm(3,num,index,jlo)=cclm(3,num,index,jlo)+cyp*k(3)
                    ENDDO
                 endif
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  return                   
END SUBROUTINE lomain
