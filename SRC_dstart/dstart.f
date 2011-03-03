	program dstart
!                                                                       
!                S T A R T  D E N S I T Y   G E N E R A T I O N         
!                                ACTUAL VERSION                         
!
      use struct
      use rotmat
      use waves
      IMPLICIT REAL*8 (A-H,O-Z)
!      include 'param.inc'
      parameter (npt=nrad+ipinst)
!
      real*8, allocatable :: akap(:,:),rhoat(:,:),rholm(:,:,:)
      real*8, allocatable :: zatom(:,:),zatom1(:,:),fj(:)
      complex*16,allocatable :: rhorad(:,:,:),rhok(:,:)
      complex*16,allocatable :: cstar(:,:,:),cy(:)
      integer,allocatable :: lm(:,:,:),lmmax(:),nptat(:) 
      dimension        vv(3),vvrot(3)
!      dimension        fj(ncom)
!      dimension        lm(2,ncom,nato)
!.....max L=10 in LM-list
!      dimension        lm0(2,121)
!      dimension        lmmax(nato),nptat(nato)
!      dimension        r(npt),rhoat(npt,nato),rhok(NWAV,NATO+1)
      dimension        r(npt)
!      dimension        rholm(npt,nato,ncom),inst(nwav)
!      dimension        inst(nwav)
!      dimension        zatom(3,nato),zatom1(3,nato)
      dimension        stg(3,nsym),ig(3),ind(nsym) 
!      complex*16       ch1,ch2,cstar(ncom,nato,nsym),cima
      complex*16       ch1,ch2,cima,tmpcom
!      COMPLEX*16       cy((lmax2+1)*(lmax2+1)),taup(nsym)
      COMPLEX*16       taup(nsym)
!      CHARACTER*4      LATTIC,IREL,RELA,cform                           
!      CHARACTER*5      MODUS                                            
!      CHARACTER*10     ANAME 
      CHARACTER*11     STATUS,FORM                                      
      CHARACTER*67       ERRMSG
      CHARACTER*80       DEFFN, ERRFN
      CHARACTER*80     FNAME,FNAME1
!      LOGICAL          REL,save,ortho                                         
      LOGICAL          save,ortho,fftstop                                         
!                                                                       
!      common /xa/     rhok,rhorad,rholm
!      COMMON /STRUK/  POS(3,NDIF),AA,BB,CC,ALPHA(3),RMT(NATO),V(NATO),
!     *                PIA(3),VOL,IATNR(NATO),MULT(NATO),ISPLIT(NATO)    
      COMMON /COM/    EMIN,EF,ELECN,XWT,NSPIN,           &
                      NBAND,NK,MINWAV,MAXWAV           
!      COMMON /CHAR/   TITLE,LATTIC,MODUS,ANAME(NATO)                    
!      COMMON /CHAR/   MODUS                   
      COMMON /GENER/  BR1(3,3), &
                      BR2(3,3)                                 
!      COMMON /GITT/   KZZ(3,NWAV),ABSK(NWAV),nwave
!      COMMON /SYM2/   TAU(3,NSYM),IORD,IZ(3,3,NSYM)
      COMMON /ORTH/   ORTHO
!      COMMON /POTNLC/ R0(NATO),DX(NATO),JRI(NATO) 
!      COMMON /ROTMAT/ ROTLOC(3,3,NATO),ROTIJ(3,3,NDIF),TAUIJ(3,NDIF)    
!      COMMON /ROTMAT/ ROTIJ(3,3,NDIF),TAUIJ(3,NDIF)    
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
       fftstop=.false.
       api=acos(-1.d0)
       tpi=2.d0*api
      cima=(0.d0,1.d0)
      CALL GTFNAM(DEFFN,ERRFN)
      CALL ERRFLG(ERRFN,'Error in DSTART')
      OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
   10 CONTINUE
         READ (1,*,END=20,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL
         OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
         if(iunit.eq.4) fftstop=.true.
      GOTO 10
   20 CONTINUE
      CLOSE (1)
!     
!.....DEFINE CONSTANTS                                                  
!      DO 100 I=1,NATO                                                   
!      MULT(I)=0                                                         
!      V(I)=0.D0                                                         
!  100 RMT(I)=0.D0                                                       
!                                                                       
!C.....READ STRUCT 
                                                      
      call read_struct
!c allocate nato and ndif arrays
      allocate ( akap(4,nat) ) 
      allocate ( lmmax(nat),nptat(nat),rhoat(npt,nat) )
      allocate ( zatom(3,nat),zatom1(3,nat) )

      allocate ( ROTIJ(3,3,NDIF),TAUIJ(3,NDIF) )
      allocate ( lm(2,ncom,nat) )
!
      do i=1,nat
      zatom(1,i)=zz(i)
      enddo
        assign 2021 to iform1
        assign 2071 to iform2
        assign 2076 to iform3
 2022 FORMAT(3X,5E13.7)                                                 
 2021 FORMAT(3X,4E19.12)                                                
 2070 FORMAT(3X,3I5,2E15.7)                                             
 2071 FORMAT(3X,3I5,2E19.12)                                            
 2075 FORMAT(3X,3I5,4E15.7)                                             
 2076 FORMAT(3X,3I5,2E19.12,2e11.3)                                     
!                                                                       
      write(6,2031)
      CALL LATGEN 


!     generation des inversen gitters

!       einlesen von rmt*kmax aus eingangsdaten von lapw1

	read(17,*)
	read(17,*) rkm
	write(6,1018) rkm
 1018	format(3x,' rmt(min)*kmax = ',f10.5)

      NOTDEF=999999
!      KMAX(1)=KMAX1
!      KMAX(2)=KMAX2
!      KMAX(3)=KMAX3

!
!.....FIND MAX RECIPROCAL LATTICE VECTORS
      KXMAX=0
      KYMAX=0
      KZMAX=0


      rrmt=9999
      jatom=nat
      do 130 i1=1,jatom
         rrmt=dmin1(rrmt,rmt(i1))
 130  CONTINUE
	gmin=rkm/rrmt*2.



!       kxmax=rkm * aa/rrmt/2./api+1
!       kymax=rkm * bb/rrmt/2./api+1
!       kzmax=rkm * cc/rrmt/2./api+1
      lmax2=0
      ncom1=1
      read(15,*) 
      read(15,*) 
      read(15,*) 
	do 140 ia=1,jatom
!bk940511
      READ(15,'(121(i3,i2))') ( (LM(J,JLM,ia),J=1,2), JLM=1,NCOM)
!     NEG L MEANS NEGATIVE SPHERICAL HARMONIC COMB (SEE KURKI-SUONIO)   
      DO 224 JLM=2,NCOM                                                 
      lmax2=max(iabs(lm(1,jlm,ia)),lmax2)
  224 IF(LM(1,JLM,ia).EQ.0) GOTO 198                                    
  198 LMMAX(ia)=JLM-1   
      if(ncom1.lt.jlm-1) ncom1=jlm-1
 140  CONTINUE
      allocate ( cy((lmax2+1)*(lmax2+1)) )
      allocate ( rholm(npt,nat,ncom1+10), cstar(ncom1,nat,nsym) )
      allocate ( rhorad(nrad,ncom1,nat), fj(ncom1+6) )

!       read cutoff-parameter for density
	rewind 15
	read(15,*)
	read(15,*)
	read(15,*)
	do 150 ia=1,jatom
	  read(15,*)
 150  CONTINUE
      gmax=12.0
      read(15,*,end=26) gmax1
      if((gmax1.gt.6.d0).and.(gmax1.lt.30.1d0)) gmax=gmax1
   26 continue
!      save=.false.
      if (ortho.or.lattic(1:1).eq.'R') then
         call deter (gmax,pia,BR2,kxmax,kymax,kzmax,lattic)
      else
        kxmax=gmax*aa/2./api+1
        kymax=gmax*bb/2./api+1
        kzmax=gmax*cc/2./api+1
      endif
         if (lattic(1:1).eq.'R') then
            KXMAX=max(kxmax,int(GMAX*AA/2.d0/PI+1))
            KYMAX=max(kymax,int(GMAX*BB/2.d0/PI+1))
            KZMAX=max(kzmax,int(GMAX*CC/6.d0/PI+1))
         endif
  IFF1=(kxmax+1)*2
  IFF2=(kymax+1)*2
  IFF3=(kzmax+1)*2
  IFF1T=IFF1 
  IFF2T=IFF2
  IFF3T=IFF3
  CALL IFFLIM(IFF1T,IFF1)
  CALL IFFLIM(IFF2T,IFF2)
  CALL IFFLIM(IFF3T,IFF3)

!
!      if(kxmax.lt.kmax1.and.kymax.lt.kmax2.and.kzmax.lt.kmax3) 
!     *   save=.true.
!.....GENERATE LIST OF RECIPROCAL LATTICE VECTORS KZZ
!     THIS LIST CONTAINS H,K,L VALUES FORBIDDEN BY GENERAL COND. TOO

	write(6,*) 'gmin',gmin
	write(6,*) 'gmax',gmax
      if(gmin.gt.gmax) then
          write(6,'(":WARN  :  GMAX  .lt. GMIN; You should increase GMAX in case.in2")')
      endif
   
          call cputim(time0)
      CALL RECPR (INDMAX,KXMAX,KYMAX,KZMAX &
        ,gmin,gmax,nwave0)
          call cputim(time1)
!      nwave=nwave0
      WRITE(6,205)NWAVE0
      WRITE(6,204)NWAVE,kxmax,kymax,kzmax
  206 FORMAT(' MINIMUM FFT-mesh in LAPW0',3i5)
  205 FORMAT(I10,' FOURIER COEFFICIENTS CALCULATED UP TO GMIN')
  204 FORMAT(I10,' FOURIER COEFFICIENTS CALCULATED UP TO GMAX',/ &
       ' LARGEST COMPONENTS:',3I4)
!...   nur gittervektoren mit betrag < kmax
      write(6,206) iff1,iff2,iff3
!     write new case.in0_std
      read(14,'(6x,i2)') ixc
      write(13,'("TOT",i5,"    (5...CA-LDA, 13...PBE-GGA, 38...WC-GGA)")') ixc
      write(13,'("NR2V      IFFT      (R2V)")')
      write(13,'(3i4,f8.2,"    min IFFT-parameters, enhancement factor")') iff1,iff2,iff3,2.0
      close(13)
!
      if(fftstop) then  
        CALL ERRCLR(ERRFN)
        stop 'DSTART FFTSTOP'
      endif
      
      allocate ( rhok(nwave,nat+1) )
      RHOK=(0.0D0,0.0D0)

      index=0
      do 300 ia=1,jatom
!                                                                       
 1005 FORMAT(/,3X,'LMMAX=',I3,'  LM = ',121(i3,I2),/)                      
 1990 FORMAT(3X,'ATOMNUMBER =',I3,5X,10A4)                             

      WRITE(6,1005)  LMMAX(ia),((LM(J,JLM,ia),J=1,2), JLM=1,LMMAX(ia))
!      WRITE(21,1005) LMMAX(ia),((LM(J,JLM,ia),J=1,2), JLM=1,LMMAX(ia))
      IMAX=JRI(ia)                                                   
      DO 299 I=1,NPT                                           
  299 R(I)=R0(ia)*EXP((I-1)*DX(ia))                            

	DO  J=1,LMMAX(ia)
	  do  j1=1,jatom
            DO  I=1,IMAX
	      RHOLM(I,j1,J)=0.d0
            enddo
	  enddo
	enddo
	

      jspin=1
      IT=80
      DO ISPIN=1,JSPIN
         IT=IT+1
             read(IT,1991) nptat(ia)
 1991	     format(43x,i5)
             if(nptat(ia).gt.npt) stop 'ipinst too small'
          do id=1,5
             read(IT,*) 
          enddo 

          read(IT,iform1) (rhoat(KK,ia),KK=1,nptat(ia))
          do id=1,6
             read(IT,*) 
          enddo 
       enddo 
!.....fix for zero density
           do KK=1,nptat(ia)
           if(rhoat(KK,ia).lt.1.d-15) rhoat(KK,ia)=1.d-15
           enddo
!     bestimmung von MT und tail charge

      call  somm1(r0(ia),rhoat(1,ia),dx(ia),zamt,1,imax)

      zatom(2,ia)=zamt
!cccc      zatom(2,ia)=zamt*2
      call  somm2(r,rhoat(1,ia),dx(ia),zamt,imax,nptat(ia))
      zatom(3,ia)=zamt
!cccc      zatom(3,ia)=zamt*2

      write(6,900) aname(ia)
!       write(66,900) aname(ia)
!       write(66,1910) zatom(1,ia),zatom(2,ia),zatom(3,ia), &
!            zatom(2,ia)+zatom(3,ia)
      write(6,1910) zatom(1,ia),zatom(2,ia),zatom(3,ia), &
            zatom(2,ia)+zatom(3,ia)

 900  format(//'Charge decomposition for atom  ',a10,/ &
           ,3x,'total, MT, tail, MT+tail charge')
1910  format(4f10.5)

300   continue       
      write(6,*)
      write(6,*)
      write(6,*)

!     bestimmung der Fourierkonstanten der tail-Ladungsdichte
!     
!     die ladungsdichte des tails wird in die mt-kugel fortgesetzt,
!     i.e. es wird eine kappe aufgesetzt, um die fouriertransformation
!     besser konvergent zu machen
!     die kappe hat die form b*exp(-a*r)/(eps+exp(-a*r))
 
       do ia=1,nat
          imax=jri(ia)
          DO I=1,NPT                                            
             R(I)=R0(ia)*EXP((I-1)*DX(ia))                            
          enddo
          
!..   wert und erste ableitung der atomdichte an rmt
          
          gg1=rhoat(imax,ia)/r(imax)/r(imax)
          akap1=gg1
          gg1=rhoat(imax-1,ia)/r(imax-1)/r(imax-1)
          gg2=rhoat(imax+1,ia)/r(imax+1)/r(imax+1)
          akap2=(gg2-gg1)/2.D0/dx(ia)/r(imax)


!...  ueberhoehungsfaktor der kappe
          akap(3,ia)=1+1.5D0*rmt(ia)*rmt(ia)
          aksi=akap(3,ia)
          
!..   parameter fuer die kappe      
          afak=-rmt(ia)*akap2/akap1*aksi/(aksi-1)
          ar=afak*(1-exp(-afak)) 
          it=0
 360      it=it+1
          ar1=afak*(1-exp(-ar))
          if (abs(ar-ar1).lt.1.e-5) goto 361
          if (it.gt.1000) stop 'probleme mit log. ableitung'
          ar=ar1
          goto 360
 361      continue
          akap(1,ia)=ar/rmt(ia)
          eps=(aksi-1)/(exp(ar)-aksi)
          aa0=akap1*(eps*exp(ar)+1)
          akap(2,ia)=aa0
          akap(3,ia)=eps   

          write(6,1014) (akap(jj,ia),jj=1,3)      
 1014     format(3x,'Kappe: a,b,eps:  ',3f10.5) 
          
         
!     berechnung der fourierkomponenten fuer  kappe und tail          
!     die berechnung fuer den tail steht hier, fuer die kappe
!     in kapp

          do 301 ik=1,nwave0
             kx=kzz(1,ik)
             ky=kzz(2,ik)
             kz=kzz(3,ik)
             akx=kx*br1(1,1)+ky*br1(1,2)+kz*br1(1,3)
             aky=kx*br1(2,1)+ky*br1(2,2)+kz*br1(2,3)
             akz=kx*br1(3,1)+ky*br1(3,2)+kz*br1(3,3)
             akk=dsqrt(akx*akx+aky*aky+akz*akz)
             
             if (akk.lt.1.d-6) then
!cccc                hint=zatom(1,ia)-zatom(2,ia)
                hint=zatom(3,ia)
                call kapp(akk,akap(1,ia),rmt(ia),hint1)
                hint=hint+hint1
!cccc                hint=hint+hint1*2.
                rhok(ik,ia)=rhok(ik,ia)+hint
                goto 301
             endif
             
             rr1=r(imax)
             gg1=rhoat(imax,ia)/rr1/rr1
             g1=dlog(gg1)
             hint=0.
             
             do 303 ir=imax+1,nptat(ia)
                rr2=r(ir)
                g2=dlog(rhoat(ir,ia)/rr2/rr2)  
                ac=-(g2-g1)/(rr2-rr1)
                bc=dexp(((g1+g2)+ac*(rr1+rr2))*0.5D0)
                
                akr=akk*rr2
                akn=ac*ac+akk*akk
                ak1=-rr2*dexp(-ac*rr2)/akn
                ak2=-dexp(-ac*rr2)/(akn*akn)
                h1=ak1*(ac*dsin(akr)+akk*dcos(akr))
                h1=h1+ak2*((ac*ac-akk*akk)*dsin(akr)+ &
                     2.D0*ac*akk*dcos(akr))
                
                
                
                akr=akk*rr1
                akn=ac*ac+akk*akk
                ak1=-rr1*dexp(-ac*rr1)/akn
                ak2=-dexp(-ac*rr1)/(akn*akn)
                h1=h1-ak1*(ac*dsin(akr)+akk*dcos(akr))
                h1=h1-ak2*((ac*ac-akk*akk)*dsin(akr)+ &
                     2.D0*ac*akk*dcos(akr))
                
                
                h1=h1*bc
                hint=hint+h1
                rr1=rr2
                g1=g2
 303         continue              
             
             
!     addition einer kappe
             
             call kapp(akk,akap(1,ia),rmt(ia),hint1)
             hint=hint+hint1
             
!     
             
             
!     factor 2 wegen spinentartung                  
             rhok(ik,ia)=rhok(ik,ia)+hint/akk
!ccccc             rhok(ik,ia)=rhok(ik,ia)+2.*hint/akk
 301      continue
       enddo
           call cputim(time2)
      

!  ... generation der interstialen ladungsdichte 
       write(6,2033)
      call rgen(rhok,zatom,zatom1,akap)
           call cputim(time2a)
            print*, 'time in rgen',time2a-time2      

!      die tailladungsdichte ist jetzt im K-raum vorhanden. 
!.     jezt wird der Beitrag dieser tails in jeder  einzelnen MT 
!.     bestimmt und nach (l,m) entwickelt	





!...    generation der lm-abhaengigen sachen
      lmxxx=0
      do ia=1,jatom
               do ll=1,lmmax(ia)
                lmxxx=max(lmxxx,iabs(lm(1,ll,ia)))
!         lmxxx=max(lmxxx,iabs(lm(1,lmmax(ia),ia)))
               enddo   
      enddo   
      
        volin=1./vol
        do 305 ik=1,nwave0
           kx=kzz(1,ik)
           ky=kzz(2,ik)
           kz=kzz(3,ik)
           akx=kx*br1(1,1)+ky*br1(1,2)+kz*br1(1,3)
           aky=kx*br1(2,1)+ky*br1(2,2)+kz*br1(2,3)
           akz=kx*br1(3,1)+ky*br1(3,2)+kz*br1(3,3)
           akk=dsqrt(akx*akx+aky*aky+akz*akz)
           ig(1)=kx
           ig(2)=ky
           ig(3)=kz
           
           do ia=1,jatom
              lmxx=iabs(lm(1,lmmax(ia),ia))
              do il=1,lmmax(ia)
                 do ir=1,jri(ia)
                    rhorad(ir,il,ia)=0.
                 enddo
              enddo
           enddo
!..     radialer +  nur l-abhaengiger teil
!..     radialer anteil auch mit r*r multipliziert 
!...    atomschleife notwendig da versch atome versch mesh haben
           do ia=1,jatom
              lmxx=max(iabs(lm(1,lmmax(ia),ia)),6)
              do  ir=1,jri(ia)
                 rr1=R0(ia)*EXP((Ir-1)*DX(ia))
                 akr=akk*rr1
                 call sphbes(lmxx,akr,fj)
!             write(61,*) ik,ir,fj(1),akr                 
                 do ill=1,lmmax(ia)
                    il=iabs(lm(1,ill,ia))+1
                    rhorad(ir,ill,ia)=rhorad(ir,ill,ia)+ &
                         rhok(ik,nat+1)*fj(il)*rr1*rr1
!ccc     $                   rhok(ik,jatom+1)*fj(il)*rr1*rr1
                 enddo
              enddo
           enddo
!...    starsummation ueber ylm(g)*exp(iGpos)

         do i=1,iord
            do j=1,jatom
               do j1=1,lmmax(j)
                  cstar(j1,j,i)=0.
               enddo
            enddo
         enddo

      NST=0
      DO 501 I=1,IORD
	 tk=0.
         DO 502 J=1,3                                                   
            tk=tk + tau(j,i)*ig(j) 
            jK=0
            DO 503 L=1,3                                                
               jK=IZ(J,L,I)*iG(L)+jK
 503        CONTINUE
            STG(J,I)=jK
 502     CONTINUE
         tk=tk*tpi
         IF(NST.EQ.0) GOTO 507
         DO 504 M=1,NST
            DO 505 J=1,3
               IF(STG(J,M).NE.STG(J,I)) GOTO 504
 505        CONTINUE
           taup(M)=taup(M)+dcmplx(cos(tk),sin(tk)) 
           IND(M)=IND(M)+1                                             
            GOTO 501
 504     CONTINUE
 507     NST=NST+1
         vv(1)=0.
         vv(2)=0.
         vv(3)=0.
         DO 506 J=1,3
            STG(J,NST)=STG(J,I)                                         
             vv(1)=vv(1)+stg(j,nst)*br1(1,j)
             vv(2)=vv(2)+stg(j,nst)*br1(2,j)
             vv(3)=vv(3)+stg(j,nst)*br1(3,j)
 506     CONTINUE
         taup(nst)=dcmplx(cos(tk),sin(tk))
         ind(nst)=1
         do ia=1,jatom
              CALL ROTATE (vv,rotloc(1,1,ia),vvrot)
         call ylm(vvrot,lmxxx,cy)
            index=1
            do ia1=1,ia-1
               index=index+mult(ia1)
            enddo
            tk=0.
            do i1=1,3
               tk=tk+stg(i1,nst)*pos(i1,index)
            enddo
               tk=tk*tpi
               tmpcom=dcmplx(cos(tk),sin(tk))
            do ill=1,lmmax(ia)
               il=iabs(lm(1,ill,ia))
               ipp=isign(1,lm(1,ill,ia))
               im=lm(2,ill,ia)
               ii=il*(il+1)+1+im
               ch1=cima**il
	       if (im.eq.0) then 
		  cyy=cy(ii)
	       else
		  if (ipp.gt.0) then
		     cyy=((-1)**im)*dble(cy(ii))*sqrt(2.d0)
		  endif
		  if (ipp.lt.0) then
		     cyy=((-1)**im)*dimag(cy(ii))*sqrt(2.d0)
		  endif
	       endif
               cstar(ill,ia,nst)=cstar(ill,ia,nst)+ &
                    tmpcom*cyy*ch1
            enddo
         enddo
         
 501  CONTINUE


           tmp=4.d0*api*volin
           do ia=1,jatom
              do ill=1,lmmax(ia)
                 ch1=0.
                 DO Ist=1,NST
                    ch1=ch1+cstar(ill,ia,ist)* &
                        taup(ist)/ind(ist) 		 
                 enddo
                 ch1=ch1*tmp/nst
                 do ir=1,jri(ia)
                    ch2=ch1*rhorad(ir,ill,ia)
                   rholm(ir,ia,ill)=rholm(ir,ia,ill)+dble(ch2)
!           if(ill.eq.1) write(61,*) ir,ch1,ch2
                 enddo
              enddo
           enddo
305    continue
!...   tail-charge im mt-radius
       sqrfpii=1.d0/sqrt(4.d0*api)
       do ia=1,nat
          do ir=1,jri(ia)
             rr1=R0(ia)*EXP((Ir-1)*DX(ia))
             tmp=exp(-akap(1,ia)*rr1)
             akap1=akap(2,ia)*tmp/ &
                  (akap(3,ia)+tmp)*rr1*rr1
             rholm(ir,ia,1)=rholm(ir,ia,1)-akap1*sqrfpii
!cccc             rholm(ir,ia,1)=rholm(ir,ia,1)-2.*akap1/sqrt(4.*api)
          enddo
       enddo


       write(6,2033)
       do ia=1,jatom
          call somm1(r0(ia),rholm(1,ia,1),dx(ia),zmt,1,jri(ia))
 
         write(6,1920) aname(ia),zmt*sqrt(4.d0*api)

       enddo
1920   format(3x,'tail charge in MT for atom   ',a10,3x,f10.5)


!       die ladungsdichte im MT = atomdichte dort + tailladungsdichte 
! ...   addition der mt ladungsdichte

        do ia=1,nat
           do ir=1,jri(ia)
              rholm(ir,ia,1)=rholm(ir,ia,1)*sqrt(4.*api) &
                   +rhoat(ir,ia)
!cccc     $             +2.*rhoat(ir,ia)
           enddo
        enddo

!c..      test of charge neutrality

          z1ges=0.
          z2ges=0.
          zges=0.
          do ia=1,nat
          call somm1(r0(ia),rholm(1,ia,1),dx(ia),zmt,1,jri(ia))
          z1ges=z1ges+zmt*mult(ia)
          z2ges=z2ges+zatom1(1,ia)*mult(ia)
          zges=zges+zatom(1,ia)*mult(ia)
          enddo
          z2ges=-z2ges+rhok(1,nat+1)
          ztges=z1ges+z2ges
          dz=zges-ztges
          write(6,2033)
          write(6,1058)
          write(6,1057) zges,ztges,dz
 1057     format(3x,3f10.5)
 1058     format(3x,'test of charge neutrality', &
               ' ch. requ., act. ch., diff  ')


!...   ausschreiben

           call cputim(time3)
        write(6,*) 'time recpr',time1-time0
        write(6,*) 'time rhok ',time2-time1
        write(6,*) 'time rholm',time3-time2



        iscf=0
      WRITE(51,1970) ISCF                                               
      WRITE(51,78) TITLE,LATTIC,AA,BB,CC,JRI(1)                         
      WRITE(51,77)                                
      DO 560 JATOM=1,NAT                                                
      WRITE(51,1990) JATOM                                           
      WRITE(51,2771) LMMAX(jatom)                                    
      DO 565 LM1=1,LMMAX(jatom)                                      
      WRITE(51,2011) LM(1,LM1,jatom),LM(2,LM1,JATOM)                
      WRITE(51,iform1) ( rholm(J,JATOM,lm1), J=1,JRI(JATOM) )    
 565  WRITE(51,2031)                                              
 560  WRITE(51,2033)                                                    
!       
      WRITE(51,2051)                                                    
      WRITE(51,1980)                                                    
      WRITE(51,2061) nwave                                             
      DO 5155 J=1,nwave0                                              
 5155 WRITE(51,iform3) (kzz(JX,J),JX=1,3),rhok(J,nat+1)/vol
      aaaa=0.d0
      DO J=nwave0+1,nwave                                              
      WRITE(51,iform3) (kzz(JX,J),JX=1,3),aaaa,aaaa
      enddo
 77   FORMAT(1X,9(F9.7))                                                
 78   FORMAT(1X,A20,A4,3F10.6,I5)
 1970 FORMAT(3X,'             TOTAL CHARGE DENSITY GENERATED ',        &
       'BY',I3,'. ITERATION (NORM: CLM=CLM*R*R)')                 
 1980  FORMAT(3X)                                                       
 2771  FORMAT(3X,'NUMBER OF LM',I3//)                                  
 2011  FORMAT(3X,'CLM(R) FOR L',I3,3X,'M=',I2/)                        
 2031  FORMAT(/)                                                        
 2033  FORMAT(///)                                                      
 2051 FORMAT(3X,'MIXED TOTAL CHARGE DENSITY IN INTERSTITIAL',            &
            3X,'NEW TOTAL CHARGE DENSITY')                              
 2061 FORMAT(1X,'        ',I10,1X,'NUMBER OF PW')                                     

      CALL ERRCLR(ERRFN)
      stop 'DSTART ENDS'
!
!        error handling
!
  910 INFO = 1
!
!        'lapw2.def' couldn't be opened
!
      WRITE (ERRMSG,9000) FNAME
      CALL OUTERR('DSTART',ERRMSG)
      GOTO 999
  920 INFO = 2
!
!        file FNAME couldn't be opened
!
      WRITE (ERRMSG,9010) IUNIT
      CALL OUTERR('DSTART',ERRMSG)
      WRITE (ERRMSG,9020) FNAME
      CALL OUTERR('DSTART',ERRMSG)
      WRITE (ERRMSG,9030) STATUS, FORM
      CALL OUTERR('DSTART',ERRMSG)
      GOTO 999
  960 INFO = 7
!        Error reading file 'lapw2.def'
!
      WRITE (ERRMSG,9040) FNAME
      CALL OUTERR('DSTART',ERRMSG)
      GOTO 999
  999 STOP 'DSTART - Error'
!
 9000 FORMAT('can''t open definition file ',A40)
 9010 FORMAT('can''t open unit: ',I2)
 9020 FORMAT('       filename: ',A50)
 9030 FORMAT('         status: ',A,'  form: ',A)
 9040 FORMAT('Error reading file: ',A47)
      END
