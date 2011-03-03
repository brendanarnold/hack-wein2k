      SUBROUTINE L2MAIN(cmplx)
      USE param
      USE struct
      USE case
      USE sym2
      USE com
      USE abc
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER *4     MODUS                                          
      CHARACTER *10    BNAME                                
      CHARACTER *80    VECFN,FNAME
      LOGICAL          notcalc,cmplx 
      REAL*8,ALLOCATABLE       :: a_real(:),fj(:,:),dfj(:,:)                            
      DIMENSION ilo(0:lomax)
      COMPLEX*16      YL((LMAX2+1)*(LMAX2+1))
      COMPLEX*16      PHSHEL,CFAC,IMAG,CZERO,PH_SPIN(2)
      COMMON /CHAR/   MODUS                    
      COMMON /GENER/  BR1(3,3),BR2(3,3)
      COMMON /ATSPDT/  P(0:LMAX2,2,nrf),DP(0:LMAX2,2,nrf)                            
      COMMON /RADFU/   RF1(NRAD,0:LMAX2,2,nrf),RF2(NRAD,0:LMAX2,2,nrf)
      COMMON /RINTEG/  RI_MAT(0:lmax2,nrf,nrf,2,2)
      COMMON /XA/     R(NRAD),BK(3),BKROT(3),BKRLOC(3)
      logical          loor(0:lomax),lapw(0:lmax2)   
      common /loabc/   alo(0:lomax,2,nloat,nrf)                                           
      common /lolog/   nlo,nlov,nlon,loor,ilo,lapw  
      COMMON /PROC/    VECFN(2)
      COMMON /IPROC/   IPROC
!...............................
      complex*16 h_yl(2*LMAX2+1,iblock)
      complex*16 h_alyl(2*LMAX2+1,iblock,2)
      complex*16 h_blyl(2*LMAX2+1,iblock,2)
!...............................
                                                                       
      DATA             CZERO/(0.0D0,0.0D0)/,IMAG/(0.0D0,1.0D0)/         
                                                                       
!------------------------------------------------------------------     
        assign 2022 to iform1
 2022 FORMAT(3X,4E19.12) 
      PI=ACOS(-1.0D0)                                                   
      TWOPI=2.D0*PI                                                     
      EDELT=0.0D0                                                       
      TEST=0.0D0                                                        
      ETOT=0.0D0
      DO JK=1,NKPT
	NB(JK)=0
      END DO
   do i=1,iord
        do ii=1,3
        do jj=1,3
        tmat(ii,jj,i)=iz(ii,jj,i)
        end do
        end do
   end do
   ALLOCATE(a_real(nmat),fj(0:lmax2,nmat),dfj(0:lmax2,nmat))

!     ---------------------------------                                 
!     START LOOP FOR ALL ATOMS BY DO 50                                 
!     ---------------------------------                                 
      do is=1,iso
      READ(17+is,2032) ISCF 
      end do                                               
      LFIRST=1 
      icount=1                                                  
      DO 50 JATOM=1,NAT 
        write(6,*)'JATOM=',JATOM
        write(6,808)((ROTLOC(ii,jj,jatom),jj=1,3),ii=1,3)
 808    format(3f12.6)
	notcalc=.true.
	if (jatom.eq.iatom(icount)) then
	notcalc=.false.
	icase=icount
	icount=icount+1
	end if

      LMX=3                                                         
        IF(JATOM.GT.1) LFIRST=LFIRST + MULT(JATOM-1)                      
        ITAP=30+JATOM                                                     
                                                                       
! CALCULATE RADIAL FUNCTIONS U(R), UE(R), ...                       

        do is=1,iso
	itape=8+is
        jtape=17+is
	rewind(itape)
        CALL ATPAR(JATOM,LFIRST,itape,jtape,is,iso)
        rewind(itape)
	end do
	if (notcalc) goto 50
        FAC=4.0D0*PI*RMT(JATOM)**2/SQRT(VOL)
        L=LCASE(ICASE)
        CALL ROUT(l,iso)
	CFAC=IMAG**L
!................................

        KKK=0     

!.....READ IN K POINT AND BASIS VECTORS                                 
        iloop=0
    
 788    continue
        if (IPROC.GT.0) then
         iloop=iloop+1
         write(6,*)'ILOOP ',iloop
	 do is=1,iso
	 itape=8+is
         close(itape)
         call mknam(FNAME,VECFN(is),ILOOP)
         write(6,*)'opening ',FNAME
         open(itape,FILE=FNAME,STATUS='old',FORM='unformatted',ERR=950)
	 end do
        endif

	DO I=1,NAT
        do is=1,iso
        itape=8+is
	   READ(itape) EMIST
	   READ(itape) EMIST
        enddo
        ENDDO

   4    continue
        call inispl
!......reading from vector for both spins
        DO 17 is=1,iso
        itape=8+is
        READ(itape,END=998) S,T,Z,BNAME,N,NE     
        if (is.eq.1) then
	KKK=KKK+1
        WRITE(29,205) S,T,Z,N,NE,BNAME    
	end if   
 205    FORMAT(/,1X,' K-POINT:',3F8.4,1X,I5,I4,2X,A10)                 
        IF(icase.EQ.1.and.is.eq.1) NK=NK+1                                      
        READ(itape) (KX(I),KY(I),KZ(I),I=1,N)                           
        DO 5 I=1,N                                                      
          BKX(I)=(S+KX(I))                                              
          BKY(I)=(T+KY(I))                                              
          BKZ(I)=(Z+KZ(I)) 
 5      CONTINUE                   
      CALL HARMON(N,BKX,BKY,BKZ,LMAX2,FJ,DFJ,RMT(JATOM))                       
      NEMIN=1                                                           
      NEMAX=0 
   14 READ(itape) NUM,E(NUM) 
      if (.not.cmplx) then
	READ(itape) (A_real(I),I=1,N)    
        do i=1,n
	a(i,num,is)=(1.d0,0.d0)*a_real(i)
        end do 
      else                                      
      READ(itape) (A(I,NUM,is),I=1,N) 
      end if   
      IF(E(NUM).LT.EMIN) NEMIN=NEMIN+1                                  
      IF(E(NUM).LT.EMAX) NEMAX=NEMAX+1                              
      IF(NUM.EQ.NE) GOTO 16                                             
      GOTO 14                                                           
   16 CONTINUE            
      NB(KKK)=NEMAX-NEMIN+1
 17   CONTINUE
!.......reading from vect end
      if (MODUS.eq.'TOTA'.or.MODUS.eq.'SPIN') goto 888
!.......sum over sym. operations
      DO 777 MU=1,IORD
      DO IS=1,iso
      PH_SPIN(is)=EXP((2*is-3)*IMAG*PHASE(MU)/2)
      DO 9 I=1,2*LMAX2+1                                            
       DO 9 NUM=NEMIN,NEMAX   
         DO 9 irf=1,nrf                                          
         ALM(I,num,irf,is)=CZERO                                               
    9 CONTINUE 
      END DO

      DO 120 iI=1,N-(nlo+nlon+nlov),iblock
      i3=0
      do 121 i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
      i3=i3+1
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

         DO  M=1,2*L+1                                                    
         IND_YL=M+L*L
         h_yl(M,i3)=conjg(yl(ind_yl))*phshel
         END DO
 
  121   continue

	DO 122 IS=1,iso
           i3=0
           do  i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
           i3=i3+1

	   DO M=1,2*L+1
          if (lapw(l)) then
          h_ALYL(m,i3,is)=(DFJ(L,I)*P(l,is,2)-FJ(L,I)*DP(l,is,2))* &
                      h_yl(M,i3)*ph_spin(is) 
          h_BLYL(m,i3,is)=(FJ(L,I)*DP(l,is,1)-DFJ(L,I)*P(l,is,1))* &
                       h_yl(M,i3)*ph_spin(is) 
          else
          h_ALYL(m,i3,is)=FJ(L,I)/P(l,is,1)/RMT(JATOM)**2* &
                        h_yl(M,i3)*ph_spin(is)
          h_BLYL(m,i3,is)=(0.d0,0.d0)
          end if
           END DO
           enddo

           ibb=min(iblock,N-(nlo+nlon+nlov)-ii+1)
           lda=2*LMAX2+1
           ldc=lda
           ldb=nmat
	   
           call zgemm('N','N',2*l+1,nemax-nemin+1,ibb,(1.d0,0.d0), &
         h_alyl(1,1,is),lda,a(ii,nemin,is),ldb,(1.d0,0.d0), &
         alm(1,nemin,1,is),ldc)
           call zgemm('N','N',2*l+1,nemax-nemin+1,ibb,(1.d0,0.d0), &
         h_blyl(1,1,is),lda,a(ii,nemin,is),ldb,(1.d0,0.d0), &
         alm(1,nemin,2,is),ldc)
           
  122 CONTINUE
  120 CONTINUE 
      if (nlo.ne.0) then
        call &
        lomain (nemin,nemax,lfirst,lfirst,n,jatom, &
                        mu,L,iso)
 	end if
        DO 230 M=1,2*L+1                                                    
          DO 892 NUM=NEMIN,NEMAX       
	   DO 892 is=1, iso     
            DO 892 irf=1,nrf                             
        ALM(M,NUM,irf,is)=ALM(M,NUM,irf,is)*FAC*CFAC
  892     CONTINUE                                                        
  230   CONTINUE         
         
         CALL XSPLT(nemin,nemax,lcase(icase),MU,iso)
 777  CONTINUE
         CALL QTL(nemin,nemax,icase,iso)
 888  CONTINUE
	 CALL PSPLIT(nemin,nemax,icase,iso)
      GOTO 4
	

 998  CONTINUE

      if (iloop.lt.iproc) goto 788

      CALL OUTP(ICASE)   
       CLOSE(15)
       CLOSE(29)
  50  CONTINUE                                                          
                                                                       
! ....END LOOP OVER ALL ATOMS                                           
      REWIND(itape)       
      
      RETURN                    
 950  CALL OUTERR('l2main','error reading parallel vectors')
      STOP 'L2main - Error'
!
 2032 FORMAT(50X,I2,//)                                                        
      END                                                               
