      SUBROUTINE LoMAIN(nemin,nemax,lfirst,latom,n,jatom, &
                        mu,LC,iso)
      USE param
      USE struct
      USE sym2
      USE abc
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16       YL((LMAX2+1)*(LMAX2+1))
      COMPLEX*16       PHSHEL,IMAG,PH_SPIN(2)
      DIMENSION   BK(3) 
      COMPLEX*16,ALLOCATABLE   :: phs(:)  
      COMMON /XA/     R(NRAD),bK1(3),BKROT(3),BKRLOC(3)
      COMMON /GENER/  BR1(3,3),BR2(3,3)                                 
      logical         loor(0:lomax),lapw(0:lomax)                                   
      common /loabc/  alo(0:lomax,2,nloat,nrf)
      common /lolog/  nlo,nlov,nlon,loor,ilo(0:lomax),lapw 

      DATA            CZERO/(0.0D0,0.0D0)/,IMAG/(0.0D0,1.0D0)/         
                                                                       
!------------------------------------------------------------------     
      PI=ACOS(-1.0D0)                                                   
      TWOPI=2.D0*PI
      ALLOCATE(phs(nume))
      i=n-(nlo+nlon)                                                      
      DO 10 L=0,LOMAX
       do 10 jlo=1,ilo(l)
        do 20 jneq=1,mult(jatom)
          DO 25 M1=-l,+l                                                    
          i=i+1  
	if (.not.(L.eq.LC)) goto 25
        BK(1)=BKX(I)
        BK(2)=BKY(I)
        BK(3)=BKZ(I)
        CALL ROTATE (BK,TMAT(1,1,mu),BKROT)
        BK(1)=BKROT(1)*BR1(1,1)+BKROT(2)*BR1(1,2)+BKROT(3)*BR1(1,3)
        BK(2)=BKROT(1)*BR1(2,1)+BKROT(2)*BR1(2,2)+BKROT(3)*BR1(2,3)
        BK(3)=BKROT(1)*BR1(3,1)+BKROT(2)*BR1(3,2)+BKROT(3)*BR1(3,3)
        CALL ROTATE (BK,ROTLOC(1,1,JATOM),BKRLOC)
        CALL YLM (BKRLOC,LMAX2,YL)
        ARG1=BKROT(1)*POS(1,LFIRST)*TWOPI
        ARG2=BKROT(2)*POS(2,LFIRST)*TWOPI
        ARG3=BKROT(3)*POS(3,LFIRST)*TWOPI
        ARGT=(BKX(I)*TAU(1,MU)+BKY(I)*TAU(2,MU)+ &
              BKZ(I)*TAU(3,MU))*TWOPI
        PHSHEL=EXP(IMAG*(ARG1+ARG2+ARG3+ARGT))
	  DO 45 is=1,iso
	  PH_SPIN(IS)=EXP((2*is-3)*IMAG*PHASE(MU)/2)
          DO 50 NUM=NEMIN,NEMAX                                           
            PHS(NUM)=PHSHEL*A(I,NUM,is)*PH_SPIN(IS)         
  50      CONTINUE                                                        
          DO 30 M=1,2*L+1    
          ind_yl=M+L*L                                                
            DO 40 NUM=NEMIN,NEMAX    
             DO 40 irf=1,nrf                                       
      ALM(m,num,irf,is)=ALM(m,num,irf,is) &
        +ALo(l,is,jlo,irf)*dconjg(YL(IND_YL))*PHS(NUM)
  40         CONTINUE
  30        CONTINUE     
  45      CONTINUE
  25    CONTINUE                                                          
  20    CONTINUE                                                          
  10  CONTINUE                                 
      return                   
      END        
