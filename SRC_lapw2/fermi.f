      SUBROUTINE FERMI(nbmax)
!                                                                      
!     NELEC IS THE NUMBER OF ELECTRONS IN THIS COMPOUND                
!     EF IS THE FERMI-ENERGY TO BE CALCULATED                          
!                                                                      
        USE param; USE parallel
        USE bandm
        USE char
        USE kpp1
        USE xa2
        USE com,only : weigh,ef,elecn,xwt,jspin=>nspin,nat,nband,rel,nk,nb,minwav,maxwav,init_com
      use struk
      IMPLICIT REAL*8 (A-H,O-Z)
!
      CHARACTER *10    KNAME
      CHARACTER *67    ERRMSG                                          
      CHARACTER *80    FNAME,VECFN,VECFND,enefn,enefnd
      REAL*8,ALLOCATABLE  :: e_local(:,:)

      common /mkef/delef,ts2
      COMMON /PROC/    VECFN,VECFND,enefn,enefnd
      COMMON /IPROC/   IPROC 
!--------------------------------------------------------------------- 
!find nkpt, nmat and nume in energy file
      JSPIN=1                                                         
      ITAP=30
     k=0  
      iloop=0
 778  continue
      if(IPROC.GT.0) then
         iloop=iloop+1
         close(29)
         close(30)
         call mknam(FNAME,enefn,ILOOP)
         iunit=30
         open(30,FILE=FNAME,STATUS='old',ERR=950)
         call mknam(FNAME,enefnd,ILOOP)
         iunit=29
         open(29,FILE=FNAME,STATUS='unknown',ERR=950)
      endif

      DO I=1,NAT 
         READ(30,'(f9.5)') EMIST
         READ(30,'(f9.5)') EMIST
      ENDDO
      DO I=1,NAT                                                  
         READ(29,'(f9.5)',END=1004) EMIST
         READ(29,'(f9.5)',END=1004) EMIST
      ENDDO
      JSPIN=2                                                          
 1004 CONTINUE                                                         
       DO
         READ(30,'(3e19.12,a10,2i6)',IOSTAT=ios) SS,T,ZZ,KNAME,N,NEn
         IF (ios /= 0) EXIT
         k=k+1
         nmat=MAX(n,nmat)
         nume=MAX(nen,nume)
         DO ii=1,nen
            READ(30,*) NUM,E1
         ENDDO
         IF(jspin.EQ.2) THEN
            READ(29,'(3e19.12,a10,2i6)',IOSTAT=ios) SS,T,ZZ,KNAME,N,NEn
            nmat=MAX(n,nmat)
            nume=MAX(nen,nume)
            DO ii=1,nen
               READ(29,*) NUM,E1
            ENDDO
         ENDIF
      ENDDO
      IF(ILOOP.LT.IPROC) goto 778
      nkpt=k+1
      REWIND(29)
      REWIND(30)

      allocate(e_local(2*nkpt,nume))
      CALL init_bandm(nume)
      CALL init_com(nkpt,nume)
      CALL init_kpp(nkpt)
      CALL init_xa2(nume,nkpt)
!je
!je   Gaussverschmierungsmethode
!je
        if(efmod.eq.'GAUSS') then
           ef=ef-int(ef)
           if(ef.gt.0.1d0) ef=0.002d0
           IF(myid.EQ.0) write(6,'(23h   GAUSS-SMEARING WITH ,f10.5,4h Ry )') ef
           IF(myid.EQ.0) write(21,'(//,27h       GAUSS-SMEARING WITH ,f10.5,4h Ry )') &
                ef
           call fermig(weigh,elecn,nat,ef,nk,e)
           ef=ef+0.5d0
           return
        end if
!je
!je
!     Added TEMPS mode
      if(efmod(1:4).eq.'TEMP') then
!     temp. smearing
         delef=0.d0
         etemp=ef-int(ef)
         if(etemp.gt.0.1d0) etemp=0.002d0
         IF(myid.EQ.0) write(6,'(23h   TEMP.-SMEARING WITH ,f10.5,4h Ry )') etemp
         IF(myid.EQ.0) write(21,'(//,27h       TEMP.-SMEARING WITH ,f10.5,4h Ry )')  &
              etemp
!     Pass efmod to fermi5
         call fermi5(etemp,e,efmod)
         return
      else if(efmod.eq.'TETRA') then	
!     tetraeder method of bloechl (for ef.lt.100  corr is used!!)
         call fermi_tetra(nbmax)
         ts2=0.d0
         return	
      end if
      ENORM=2.0D0                                                      
      JSPIN=1                                                          
      K=0                                                              
      TEST=9.8D0                                                       
!     for babi test1 changed to lower value, otherwise no degeneracy fo
!      TEST1=2.D-10                                                    
      TEST1=2.D-8                                                      
      MINWAV=18000                                                     
      MAXWAV=0                                                         
!                                                                      
!.....MORE DIMENS.REPRESENTATION,IF DELTA E LESS THEN TEST1            
      ELECN=ELECN-1.D-10                                               
      INDEX=0                                                          
      SUMW=0.0D0                                                       

!para begin
! ensure to mimick the proper vector file
! if running on multiple processors 
      iloop=0
	iloop0=0
 777  continue
      if(IPROC.GT.0) then
         iloop=iloop+1
	 iloop0=iloop0+1
  
         close(29)
         close(30)
         call mknam(FNAME,enefn,ILOOP)
         iunit=30
         open(30,FILE=FNAME,STATUS='old',ERR=950)
         call mknam(FNAME,enefnd,ILOOP)
         iunit=29
         open(29,FILE=FNAME,STATUS='unknown',ERR=950)
      endif
!para end
      DO I=1,NAT 
         READ(30,'(f9.5)') EMIST
         READ(30,'(f9.5)') EMIST
      ENDDO
      DO I=1,NAT                                                  
         READ(29,'(f9.5)',END=1002) EMIST
         READ(29,'(f9.5)',END=1002) EMIST
      ENDDO
      JSPIN=2                                                          
 1002 CONTINUE                                                         
!                                                                      
    4 K=K+1
      if(iloop0.ne.0) KPP(ILOOP0)=K
!para begin
! testing
!      write(*,*)'reading k=',K
!para end

      IF(K.GT.2*NKPT) GOTO 900
      READ(ITAP,'(3e19.12,a10,2i6,f5.1)',END=999) SS,T,ZZ,KNAME,N,NE(K),WEIGHT(K)
      nmat=MAX(n,nmat)
      if(ne(k).gt.nume) GOTO 920
      INDEX=INDEX+NE(K)
      SUMW=SUMW+WEIGHT(K)
      IF(N.GT.MAXWAV) MAXWAV=N                                         
      IF(N.LT.MINWAV) MINWAV=N                                         
!      READ(ITAP) (KX(I),KY(I),KZ(I),I=1,N)                             
 14   READ(ITAP,*) NUM,E1                                                
      E_local(K,NUM)=E1                                                      
      if(itap.eq.30.and.e1.gt.ebmax(num)) ebmax(num)=e1
      if(itap.eq.30.and.e1.lt.ebmin(num)) ebmin(num)=e1
!      READ(ITAP,iostat=ios) (A(I),I=1,N)
      IF(NUM.EQ.NE(K)) GOTO 4                                          
      GOTO 14                                                          
!                                                                      
! ... CALCULATE WEIGHT FOR K-POINTS                                    
!                                                                      
     

!para begin
! we still have vector files to read

!  999 K=K-1
 999  CONTINUE
      K=K-1
      write(6,*)'skipping last k point'
      IF(ILOOP.LT.IPROC) goto 777


!para end
!************************ needed for weight-file **********************
      NK=K
!************************ needed for weight-file **********************

      IF(ITAP.EQ.30.AND.JSPIN.EQ.2) THEN                               
         ITAP=29
         SUMUP=SUMW
!para begin
         iloop=0
   	
         IF(IPROC.EQ.0) GOTO 4 
         GOTO 777
!para end

      ELSE IF(ITAP.EQ.29.AND.JSPIN.EQ.2) THEN                           
         DIF=2*SUMUP-SUMW                                              
         IF(ABS(DIF).GT.TEST1) GOTO 930
      END IF                                                           
!                                                                      
      DO 8 KK=1,K                                                      
      WEIGHT(KK)=WEIGHT(KK)*ENORM/SUMW                                 
      NNN=NE(KK)                                                       
      DO 8 NN=1,NNN                                                    
   8  WEIGH(KK,NN)=WEIGHT(KK)                                          
!                                                                      
      M=0                                                              
      EMIN=10.0D0                                                      
      ELN=0.0D0                                                        
!                                                                      
! ....SORTIEREN DER ENERGIEN NACH IHRER GROESSE                        
!     AUFSUMMIEREN BIS ELN=NELEC, BESTIMMT E-FERMI                     
!                                                                      
    9 KK=0                                                             
   10 KK=KK+1                                                          
      IF(KK.GT.K) GOTO 12                                              
      NN=0                                                             
   11 NN=NN+1                                                          
      IF(EMIN.GT.E_local(KK,NN)) GOTO 15                                     
   13 IF(NN.EQ.NE(KK)) GOTO 10                                         
      GOTO 11                                                          
!                                                                      
! ....SET NEW EMIN,STORE INDEX KK AND NN                               
!                                                                      
   15 EMIN=E_local(KK,NN)                                                    
      KNOT=KK                                                          
      NNOT=NN                                                          
      GOTO 13                                                          
!                                                                      
! ....SAVE EMIN, EVICT THIS ENERGY FROM LIST, SUM UNTILL EF IS REACHED 
!                                                                      
   12 M=M+1                                                            
      ELN=ELN+WEIGHT(KNOT)                                             
      IF(ELN.GT.ELECN) GOTO 16                                         
      E_local(KNOT,NNOT)=E_local(KNOT,NNOT)+10.d0
      EMIN=EMIN+10.0d0                                                   
      IF(INDEX.EQ.M) GOTO 940                                          
      GOTO 9                                                           
!                                                                      
! ....SAVE E-FERMI, REWIND TAPES AND RETURN                            
!                                                                      
   16 CONTINUE
      IF(EFMOD.eq.'ALL  ') GOTO 18                                     
      EF=EMIN+TEST1                                            
!      EF=EMIN+0.0000000001D0                                          
!                                                                      
! ....DETERMINATE THE NUMBER OF THE IRREDUZIBLE-REPRESENTATION OF THE  
! ....EIGENVALUE AT E-FERMI AND SET WEIGHT SO, THAT ELN=ELECN          
!                                                                      
      TEST2=ABS(ELN-ELECN)                                             
      SUB=WEIGHT(KNOT)-TEST2                                           
      E1=E_local(KNOT,NNOT)                                                  
      NNOT0=NNOT-1                                                     
      NOT01=NNOT-2                                                     
      NNOT1=NNOT+1                                                     
      NNOT2=NNOT+2                                                     
      IF(NNOT.EQ.1) GOTO 20                                            
      E0=E_local(KNOT,NNOT0)-10.d0                                           
      IF(ABS(E0-E1).GT.TEST1) GOTO 20                                  
      IF(NNOT.GT.2) THEN                                               
      E01=E_local(KNOT,NOT01)-10.d0                                          
!                                                                      
      IF(ABS(E0-E01).LT.TEST1) GOTO 25                                 
      END IF                                                           
      IF(NNOT1.GT.NE(KNOT)) GOTO 26                                    
      E2=E_local(KNOT,NNOT1)                                                 
      IF(ABS(E1-E2).GT.TEST1) GOTO 26                                  
      GOTO 27                                                          
!                                                                      
   20 IF(NNOT1.GT.NE(KNOT)) GOTO 22                                    
      E2=E_local(KNOT,NNOT1)                                                 
      IF(ABS(E1-E2).GT.TEST1) GOTO 22                                  
!                                                                      
      IF(NNOT2.GT.NE(KNOT)) GOTO 21                                    
      E3=E_local(KNOT,NNOT2)                                                 
      IF(ABS(E3-E1).GT.TEST1) GOTO 21                                  
!                                                                      
! ....3-DIM.REPRESENTATION, NNOT,NNOT+1,NNOT+2                         
!                                                                      
      WEIGH(KNOT,NNOT)=SUB/3.0D0                                       
      WEIGH(KNOT,NNOT1)=SUB/3.0D0                                      
      WEIGH(KNOT,NNOT2)=SUB/3.0D0                                      
      GOTO 18                                                           
!                                                                      
! ....2-DIM.REPR.                                                      
   21 WEIGH(KNOT,NNOT)=SUB/2.0D0                                       
      WEIGH(KNOT,NNOT1)=SUB/2.0D0                                      
      GOTO 18                                                           
!                                                                      
! ....ONE-DIM.REPR.                                                    
   22 WEIGH(KNOT,NNOT)=SUB                                             
      GOTO 18                                                           
!                                                                      
! ....3-DIM.REPRESENTATION, NNOT,NNOT-1,NNOT-2                         
   25 WEIGH(KNOT,NNOT)=(SUB+2*WEIGHT(KNOT))/3.0D0                      
      WEIGH(KNOT,NNOT0)=WEIGH(KNOT,NNOT)                               
      WEIGH(KNOT,NOT01)=WEIGH(KNOT,NNOT)                               
      GOTO 18                                                           
!                                                                      
! ....2-DIM.REPR., NNOT,NNOT-1                                         
   26 WEIGH(KNOT,NNOT)=(SUB+WEIGHT(KNOT))/2.0D0                        
      WEIGH(KNOT,NNOT0)=WEIGH(KNOT,NNOT)                               
      GOTO 18                                                           
!                                                                      
! ....3-DIM.REPR., NNOT-1,NNOT,NNOT+1                                  
   27 WEIGH(KNOT,NNOT)=(SUB+WEIGHT(KNOT))/3.0D0                        
      WEIGH(KNOT,NNOT0)=WEIGH(KNOT,NNOT)                               
      WEIGH(KNOT,NNOT1)=WEIGH(KNOT,NNOT)                               
18    CONTINUE
      DO j=1,nume
         DO i=1,2*nkpt
            e(i+(j-1)*2*nkpt)=e_local(i,j)
         ENDDO
      ENDDO
      RETURN
!
!        Error messages
!
 900  WRITE (ERRMSG,9000) NKPT
      CALL OUTERR('FERMI',ERRMSG)
      STOP  'FERMI - Error'
 920  WRITE (ERRMSG,9020) nume,ne(k) 
      CALL OUTERR('FERMI',ERRMSG)
      STOP  'FERMI - Error'
 930  CALL OUTERR('FERMI',' SUMUP AND SUMDN IN FERMI NOT EQUAL')
      CALL OUTERR('FERMI',' CHECK UP AND DOWN INPUTS OF LAPW1 ')
      WRITE (ERRMSG,9030) SUMUP,SUMW,DIF             
      CALL OUTERR('FERMI',ERRMSG)
      STOP  'FERMI - Error'
 940  WRITE (ERRMSG,9041) ELECN
      CALL OUTERR('FERMI',ERRMSG)
      WRITE (ERRMSG,9040) ELECN,ELN,EMIN,INDEX
      CALL OUTERR('FERMI',ERRMSG)
      CALL OUTERR('FERMI','    INCREASE ENERGY WINDOW IN CASE.IN1')
      CALL OUTERR('FERMI',' OR INCREASE PARAMETER NUME ')
      CALL OUTERR('FERMI',' OR DECREASE NE IN CASE.IN2 ')
      STOP  'FERMI - Error'
 950  WRITE (ERRMSG,9060) IUNIT
      CALL OUTERR('LAPW2',ERRMSG)
      WRITE (ERRMSG,9070) FNAME
      CALL OUTERR('LAPW2',ERRMSG)
!      WRITE (ERRMSG,9080) STATUS, FORM
!      CALL OUTERR('LAPW2',ERRMSG)
      STOP 'FERMI - Error'
!
 9000 FORMAT ('NUMBER OF K-POINTS .GT. NKPT =',i6)
 9020 FORMAT ('nume,ne(k)',2i6)
 9030 FORMAT ('SUMUP,SUMW,DIF', 3f10.7)
 9040 FORMAT ('ELECN,ELN,EMIN,INDEX', 3f9.5,i5)
 9041 FORMAT (' NOT ENOUGH EIGENVALUES FOR',f4.0,' ELECTRONS')
 9050 FORMAT (8F10.7) 
 9060 FORMAT('can''t open unit: ',I2)
 9070 FORMAT('       filename: ',A50)
 9080 FORMAT('         status: ',A,'  form: ',A)
                                              
!
  199 FORMAT(3F10.5,3X,A10,2I5,F5.2)                                   
  170 FORMAT(F5.2)                                                     
  201 FORMAT(8(3I3))                                                   
  202 FORMAT(I5,F10.5)                                                 
  203 FORMAT(8F9.6)                                                    
      END

      SUBROUTINE fermi_tetra(nbmax)
!                                                                      
!     NELEC IS THE NUMBER OF ELECTRONS IN THIS COMPOUND                
!     EF IS THE FERMI-ENERGY TO BE CALCULATED                          
!                                                                      
        USE parallel
        USE param
        USE bandm
        USE kpp1
        USE xa2,only : weight,e,ne
        USE com,only : weigh,ef,elecn,xwt,jspin=>nspin,nat,nband,rel,nk,nb,minwav,maxwav
      use struk
      IMPLICIT REAL*8 (A-H,O-Z)                                        
      PARAMETER (nw=250000)                                            
      CHARACTER *10    KNAME                                
      CHARACTER *67    ERRMSG                                          
      CHARACTER *80    FNAME,VECFN,VECFND,enefn,enefnd
      LOGICAL          EFERM
      INTEGER          :: iw(nw)
      REAL*8           :: eb(nume,nkpt,2)
      INTEGER          NEHELP(NKPT,2)
      common /correct/ icor

!para begin
      COMMON /PROC/    VECFN,VECFND,enefn,enefnd
      COMMON /IPROC/   IPROC
!para end

!--------------------------------------------------------------------- 
!   
!     icor switches non-linear correction  on/off  
      icor=1                      
      if(ef.gt.100.d0) icor=0                                          
      IF(myid.eq.0) write(6,*) '  BZ-integration with TETRA-program. icor=:',icor
      do 2 j=1,nume*2*nkpt
 2    e(j)=3.0d0
      do 3 i1=1,2
        do 3 i=1,nkpt
          do 3 j=1,nume
 3    eb(j,i,i1)=3.0d0
      MINWAV=18000                                                     
      MAXWAV=0                                                         
      nemax=0
      ITAP=30                                                          
      JSPIN=1                                                          
      K=0                                                              
      TEST=8.9D0                                                       
!     for babi test1 changed to lower value, otherwise no degeneracy fo
!      TEST1=2.D-10                                                    
      TEST1=2.D-8                                                      
!                                                                      
!.....MORE DIMENS.REPRESENTATION,IF DELTA E LESS THEN TEST1            
      ELECN=ELECN-1.D-10                                               

!para begin
! ensure to mimick the proper vector file
! if running on multiple processors
      k1=0
      ispin=1
      iloop=0
	iloop0=0	
 777  continue
      if(IPROC.GT.0) then
         iloop=iloop+1
	 iloop0=iloop0+1

         close(29)
         close(30)
         call mknam(FNAME,enefn,ILOOP)
         open(30,FILE=FNAME,STATUS='old',ERR=950)
         call mknam(FNAME,enefnd,ILOOP)
         open(29,FILE=FNAME,STATUS='unknown',ERR=950)
      endif
!para end

!.....READ FROM TAPE 10,WRITTEN BY LAPW1                               
      DO I=1,NAT 
         READ(30,'(f9.5)') EMIST
         READ(30,'(f9.5)') EMIST
      ENDDO
      DO I=1,NAT                                                  
         READ(29,'(f9.5)',END=1002) EMIST
         READ(29,'(f9.5)',END=1002) EMIST
      ENDDO
      JSPIN=2                                                          
 1002 CONTINUE
!                                                                      

    4 K=K+1                                                            
      if(iloop0.ne.0) KPP(ILOOP0)=K
!para begin
! testing
!      write(*,*)'reading k=',K,itap,ispin,iloop,iloop0
!para end

      IF(K.GT.2*NKPT) GOTO 900                                           
      READ(ITAP,'(3e19.12,a10,2i6,f5.1)',END=999) SS,T,ZZ,KNAME,N,NEHELP(K,ispin),WEI
      nmat=MAX(n,nmat)

!para begin
      NE(K+k1)=NEHELP(K,ispin)
!para end       	
      if(nehelp(k,ispin).gt.nume) GOTO 920
      if(nemax.lt.nehelp(k,ispin)) nemax=nehelp(k,ispin)
      IF(N.GT.MAXWAV) MAXWAV=N                                         
      IF(N.LT.MINWAV) MINWAV=N                                         
   14 READ(ITAP,*) NUM,E1                                                
      Eb(num,K,ispin)=E1                                               
      if(itap.eq.30.and.(e1.gt.ebmax(num))) ebmax(num)=e1
      if(itap.eq.30.and.e1.lt.ebmin(num)) ebmin(num)=e1
!      READ(ITAP) (A(I),I=1,N)                                          
      IF(NUM.EQ.NEHELP(K,ispin)) GOTO 4
      GOTO 14                                                          
!                                                                      
! ... CALCULATE WEIGHT FOR K-POINTS                                    
!                                                                      
!para begin
! we still have vector files to read

!  999 K=K-1
 999  CONTINUE
      K=K-1
!      write(6,*)'skipping last k point'
      IF(ILOOP.LT.IPROC) goto 777


!para end
        
!****************** needed for weight-file
       nk=k
!******************                                                   
      IF(ITAP.EQ.30.AND.JSPIN.EQ.2) THEN                               
         ITAP=29
         k1=k  
         ispin=2                                                       
!para begin
!      k=0
!      GOTO 4  
         iloop=0
         K=0
         IF(IPROC.EQ.0) GOTO 4
         GOTO 777
!para end
                                                          
      ELSE IF(ITAP.EQ.29.AND.JSPIN.EQ.2) THEN 
         IF(k1.ne.k) GOTO 930
      END IF                                                           
!     
!      write(*,*)  '2*nume,k,elecn/2.d0,nw',2*nume,k,elecn/2.d0,nw
!      write(*,1234) ((e(i1,i2),i1=1,2*nume),i2=1,k)
 1234 format(20(10f8.4,/))
      IF(myid.EQ.0) write(6,*)'call eord...' 
      call eord(e,eb,nemax,k,jspin) 
      IF(myid.eq.0) write(6,*)'call dos...' 
      call dos(nemax*jspin,k,e,weight,elecn/2.d0*jspin,d,ef,iw,nw) 
      IF(myid.eq.0) write(6,*)'call eweigh...'
      call eweigh(ef,weigh,weight,nemax,k,jspin,nehelp,eb,nbmax)
!                                                                      
!     REWIND TAPES AND RETURN                             
!                                                                      
      RETURN                                             
!
!        Error messages
!
 900  WRITE (ERRMSG,9000) NKPT
      CALL OUTERR('FERMI',ERRMSG)
      STOP  'FERMI - Error'
 920  WRITE (ERRMSG,9020) nume,nehelp(k,ispin) 
      CALL OUTERR('FERMI',ERRMSG)
      STOP  'FERMI - Error'
 930  CALL OUTERR('FERMI',' # of k-points in up and down not equal:')
      WRITE (ERRMSG,9030) k1,k
      CALL OUTERR('FERMI',ERRMSG)
      STOP  'FERMI - Error'
 950  WRITE (ERRMSG,9060) 
      CALL OUTERR('LAPW2',ERRMSG)
      WRITE (ERRMSG,9070) FNAME
      CALL OUTERR('LAPW2',ERRMSG)
      STOP 'FERMI - Error'
!
 9000 FORMAT ('NUMBER OF K-POINTS .GT. NKPT =',i6)
 9020 FORMAT ('nume,ne',2i6)
 9030 FORMAT ('k1, k', 2i6,' check INPUTS OF LAPW1 ')
 9060 FORMAT('can''t open unit: ')
 9070 FORMAT('       filename: ',A50)
!                                                                      
      END                                                              

      subroutine eord(e,eb,nemax,k,jspin) 
        USE param
      IMPLICIT REAL*8 (A-H,O-Z)                                        
      dimension eb(nume,nkpt,2),e(nemax*jspin,k)
      do 10 k1=1,k
      ne2=0
      do 10 ispin=1,jspin
      do 10 ne1=1,nemax
      ne2=ne2+1
   10 e(ne2,k1)=eb(ne1,k1,ispin)
      return
      end

      subroutine eweigh(ef,weigh,weight,nemax,k,jspin,ne,eb,nbmax)       
        USE param
        USE parallel
        USE bandm
      IMPLICIT REAL*8 (A-H,O-Z)                                        
      dimension eb(nume,nkpt,2)
      dimension weigh(2*nkpt,nume),weight(nemax*jspin,k),ne(k)
      TEST1=2.D-8                                                      
      nbmax=0  
      emax=-10.d0                                                      
      DO 8 KK=1,K                                                      
      NNN=NE(KK)             
      NNN=NEmax             
!      write(6,50) ((kk,nn,eb(nn,kk,jspin1),weight(nn+(jspin1-1)*nnn,kk),nn=1,nnn),jspin1=1,jspin)                      
  50  format((2i5,2f10.6))
      DO 8 NN=1,NNN   
       do jspin1=1,jspin
      if(abs(weight(nn+(jspin1-1)*nnn,kk)).gt.test1) then
         nbmax=max(nbmax,nn)
         emax=max(emax,eb(nn,kk,jspin1))                      
      end if
       enddo
   8  WEIGH(KK,NN)=WEIGHT(nn,KK) *2.0d0 / jspin                        
      if(myid.eq.0) write(6,*) '  number of occupied bands:',nbmax
      if(myid.eq.0) write(6,*) '  highest energy:',emax
!!!for nonmagnetic case check insulator case
!!!      if(jspin.eq.1) then
!      if(abs(ebmax(nbmax)-emax).lt.1.d-5.and.  &
!      if(abs(ef.gt.emax).and.  &
      if(ef.gt.emax.and.  &
         abs(emax-ebmin(nbmax+1)).gt.1.d-3) then
        if(myid.eq.0) write(6,*) 'insulator !'
        if((ef-emax).lt.1.d-4.or.(ef-emax).ge.1.d-4) then
          if(myid.eq.0) write(6,*) 'EF-inconsistency corrected'
          if(myid.eq.0) write(21,888) 
 888      format('       Insulator, EF-inconsistency corrected')
          ef=emax
        endif
      endif
!!!      endif
      ef=ef+0.5d0            
      return
      end
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
!     ----                                                         ----
!     ----  BLOCK KINTR                                            ----
!     ----  CALCULATION OF SAMPLING WEIGHTS                        ----
!     ----                                                         ----
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
!     .....................................................DOS.........
      SUBROUTINE DOS(NB,NKP,EB,WGHT,RNTOT,D,EF,W,NWX)                  
        USE parallel
!     **                                                              *
!     **  CALCULATES THE SAMPLING WEIGHTS FROM TETRAHEDRON INTEGRATION*
!     **                                                              *
!     **  INPUT :                                                     *
!     **    NB          NUMBER OF BANDS                               *
!     **    NKP         NUMBER OF IRREDUCIBLE K-POINTS                *
!     **    EB          ENERGIES ( DETERMINE FERMI SURFACE )          *
!     **    RNTOT       NUMBER OF OCCUPIED STATES                     *
!     **    W           (INTEGER) WORK ARRAY                          *
!     **    NWX         LENGTH OF WROK ARRAY W                        *
!     **  OUTPUT :                                                    *
!     **    WGHT        SAMPLING WEIGHTS                              *
!     **    EF          FERMI LEVEL                                   *
!     **                                                              *
!     **  AUTHOR : PETER E. BLOECHL                                   *
!     **                                                              *
!     **  SUBROUTINES USED:                                           *
!     **  EFERMI,SAMFAC,DEF0,TETR0,INITDR,TETR1,TOTNOS,EFI,WEIGHT
!     **  DRVAL
!     **                                                           **KI
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)                              
      IMPLICIT INTEGER (O)                                             
      DIMENSION EB(NB*NKP),WGHT(NB*NKP)                                
      CHARACTER *67    ERRMSG
      INTEGER W(NWX)                                                   
      DATA MWRIT/100/ICHECK/1/TOLMAX/1.D-4/                            
!     -----------------------------------------------------------------
!     --  CALCULATE FERMI LEVEL                                       -
!     -----------------------------------------------------------------
!      write(*,*) NB,NKP,RNTOT,NWX
!      write(*,*) eb(1)
!           write(6,*) 'before efermi in dos'
      CALL EFERMI(RNTOT,EF,TOLMAX,NB,NKP,EB,W,NWX)                     
                                                                       
      IF(myid.eq.0) write(6,*) '  FERMI ENERGY AT ',EF                                   
!     -----------------------------------------------------------------
!     --  CALCULATE WEIGHTS                                           -
!     -----------------------------------------------------------------
      CALL SAMFAC(NB,NKP,EB,EF,WGHT,W,NWX)                             
                                                                       
      IF(ICHECK.EQ.0) RETURN                                           
!     -----------------------------------------------------------------
!     --  CHECK WHETHER SUMRULE IS FULLFILLED                         -
!     -----------------------------------------------------------------
      SUM=0.D0                                                         
      DO 100 I=1,NB*NKP                                                
      SUM=SUM+WGHT(I)                                                  
100   CONTINUE                                                         
      IF(DABS(SUM-RNTOT).GT.TOLMAX) THEN                               
      IF(myid.eq.0) WRITE (6,9000) SUM,RNTOT
      IF(myid.eq.0) write(21,9001) SUM,RNTOT
!        GOTO 900
      ELSE                                                             
        IF(myid.eq.0) write(6,*) '  SUM RULE OF TETRAHEDRON INTEGRATION CHECKED '
      END IF                                                           
      RETURN                                                           
  900 CALL OUTERR('FERMI',' INTEGRATION FAILED.....STOP IN DOS')
      WRITE (ERRMSG,9000) SUM,RNTOT
      CALL OUTERR('FERMI',ERRMSG)
      STOP  'FERMI - Error'
 9000 FORMAT(' RESULT OF INTEGRATION: ',f10.5, '; SHOULD BE: ',f10.5)
 9001 FORMAT(':WARN : RESULT OF INTEGRATION: ',f10.5, '; SHOULD BE: ' &
       ,f10.5)
      END                                                              
!     .....................................................EFERMI......
      SUBROUTINE EFERMI(RNTOT,EF,TOLMAX,NB,NKP,EB,W,NWX)               
!     **                                                              *
!     **  CALCUALTES THE FERMILEVEL BY INTEGRATION OF THE TOTAL       *
!     **  DENSITY OF STATES                                           *
!     **                                                              *
!     **  INPUT :                                                     *
!     **    RNTOT       NUMBER OF OCCUPIED STATES                     *
!     **    NB          NUMBER OF BANDS                               *
!     **    NKP         NUMBER OF IRREDUCIBLE K-POINTS                *
!     **    EB          ENERGIES ( DETERMINE FERMI SURFACE )          *
!     **    W           (INTEGER) WORK ARRAY                          *
!     **    NWX         LENGTH OF WROK ARRAY W                        *
!     **    TOLMAX      TOLERANCE IN THE NUMBER OF STATES AT EF       *
!     **  OUTPUT :                                                    *
!     **    EF          FERMI LEVEL                                   *
!     **                                                              *
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)                              
      IMPLICIT INTEGER (O)                                             
      DIMENSION EB(NB,NKP),E(4),IKP(4)                                 
      INTEGER W(NWX)                                                   
      CHARACTER*67 ERRMSG
      DATA NP/1000/                                                    
      CALL DEF0(NWX)                                                   
!     -----------------------------------------------------------------
!     --  FIND EMIN EMAX     (ENERGYBANDS ARE ASSUMED                 -
!     --                      TO BE ORDERED WITH RESPECT TO SYMMETRY  -
!     -----------------------------------------------------------------
      EMIN=EB(1,1)                                                     
      EMAX=EB(NB,1)                                                    
      DO 100 IK=1,NKP                                                  
      DO 100 IB=1,NB                                                   
      EMIN=DMIN1(EMIN,EB(IB,IK))
      EMAX=DMAX1(EMAX,EB(IB,IK))
100   CONTINUE      
      DE=(EMAX-EMIN)/DBLE(NP-1)
      EMAX=EMAX+DE  
      EMIN=EMIN-DE 
      CALL DEFDR(ONOS,NP)
      CALL DEFDR(OSOS,NP)
      CALL TETR0(VOL1,NKP,NTET,MWRIT) 
      CALL DEFI(OWORK,5*MWRIT)                                         
      ILOOP=0                                                          
1000  CONTINUE                                                         
      ILOOP=ILOOP+1                                                    
!      write(*,*) 'in efermi',iloop
      CALL INITDR(W(ONOS),NP)                                          
      CALL INITDR(W(OSOS),NP)                                          
      INIT=1                                                           
      SUM=0.D0                                                         
      DO 200 ITET=1,NTET                                               
      CALL TETR1(INIT,ITET,IWGHT,IKP,MWRIT,W(OWORK))                   
      VOL=VOL1*DBLE(IWGHT)                                             
      SUM=SUM+VOL                                                      
      DO 220 IB=1,NB                                                   
      DO 210 I=1,4                                                     
      E(I)=EB(IB,IKP(I))                                               
210   CONTINUE                                                         
!c      if(iloop.eq.4)write(*,*)'totnos',itet,ntet,ib,nb,(e(i),i=1,4)
      CALL TOTNOS(VOL,E,EMIN,EMAX,NP,W(ONOS),W(OSOS))                  
220   CONTINUE                                                         
200   CONTINUE                                                         
      IF(DABS(SUM-1.D0).GT.1.D-5) GOTO 900
!     -----------------------------------------------------------------
!     --  GET FERMI LEVEL                                             -
!     -----------------------------------------------------------------
      tol=tolmax
!         write(*,*) 'in efi'
      CALL EFI(RNTOT,TOL,EF,EMIN,EMAX,NP,W(ONOS),W(OSOS))              
!         write(*,*) 'after efi',tol,tolmax,emax,emin,np,rntot
                                                                       
!     -----------------------------------------------------------------
!     --  CHECK ACCURACY AND RESTART IF NECCESARY                     -
!     -----------------------------------------------------------------
      IF(TOL.GT.TOLMAX) THEN                                           
        ESTEP=(EMAX-EMIN)/DBLE(NP-1)     
        IP=1+(EF-EMIN)/ESTEP                                           
!          write(*,*) ip,DRVAL(W(ONOS),IP),DRVAL(W(ONOS),IP+1)
        EMIN=EMIN+ESTEP*DBLE(IP-1)                                     
        EMAX=EMIN+ESTEP                                                
         if(estep.lt.1.d-8) then
         write(*,*) 'WARNING: EF not accurate, new emin,emax,NE-min,', &
         'NE-max',emin,emax,DRVAL(W(ONOS),IP),DRVAL(W(ONOS),IP+1)
          ef=(emin+emax)/2.d0
         goto 2000
         endif
        IF(RNTOT-DRVAL(W(ONOS),IP).LE.TOLMAX) THEN                     
          EF=EMIN                                                      
          GOTO 2000                                                    
        ELSE IF(DRVAL(W(ONOS),IP+1)-RNTOT.LE.TOLMAX) THEN              
          EF=EMAX                                                      
          GOTO 2000                                                    
        END IF                                                         
        IF(ILOOP.GT.5) GOTO 910
        write(6,*) 'ILOOP ',ILOOP                                          
        GOTO 1000                                                      
      END IF                                                           
2000  CONTINUE                                                         
      CALL RLSE(ONOS)                                                  
      RETURN                                                           
  900 CALL OUTERR('FERMI',' TETRAHEDRA DO NOT FILL VOLUME')
      CALL OUTERR('FERMI',' STOP IN EFERMI')
      WRITE (ERRMSG,9000) SUM
      CALL OUTERR('FERMI',ERRMSG)
      STOP  'FERMI - Error'
  910 CALL OUTERR('FERMI',' CANNOT FIND FERMI LEVEL')
      CALL OUTERR('FERMI',' STOP IN EFERMI')
      WRITE (ERRMSG,9010) TOL
      CALL OUTERR('FERMI',ERRMSG)
      WRITE (ERRMSG,9020) EMIN,EMAX
      CALL OUTERR('FERMI',ERRMSG)
      WRITE (ERRMSG,9030) DRVAL(W(ONOS),IP),DRVAL(W(ONOS),IP+1)
      CALL OUTERR('FERMI',ERRMSG)
      STOP  'FERMI - Error'
 9000 FORMAT(' SUM ',f10.5, ' SHOULD BE 1.')
 9010 FORMAT(' TOL ',f10.5)
 9020 FORMAT(' EMIN ',f10.5, ' EMAX ',f10.5)
 9030 FORMAT(' NOS(EMIN) ',f10.5, ' NOS(EMAX) ',f10.5)
      END                                                              
!     .....................................................EFI.........
      SUBROUTINE EFI(RNTOT,TOL,EFERMI,EMIN,EMAX,NP,NOS,SOS)            
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                              
      DOUBLE PRECISION NOS(NP),NOSUP,NOSLOW,NOSIP                      
      CHARACTER*67 ERRMSG
      DIMENSION SOS(NP)                                                
      ADD=0.D0                                                         
      DO 100 I=1,NP                                                    
      ADD=ADD+SOS(I)                                                   
      NOS(I)=NOS(I)+ADD                                                
100   CONTINUE                                                         
      IF(NOS(1).GT.RNTOT+.5D0*TOL.OR.NOS(NP).LT.RNTOT-.5D0*TOL) GOTO 900
      IPUP=NP                                                          
      IPLOW=1                                                          
      nosup=nos(ipup)
      noslow=nos(iplow)
      DO 200 IFIND=1,NP                                                
      IP=IPLOW+0.5*(IPUP-IPLOW)                                        
      NOSIP=NOS(IP)                                                    
      IF(RNTOT-NOSIP.GT.0.d0) THEN                                        
        IPLOW=IP                                                       
        NOSLOW=NOSIP                                                   
      ELSE                                                             
        IPUP=IP                                                        
        NOSUP=NOSIP                                                    
      END IF                                                           
      IF(IPUP.EQ.IPLOW+1) GOTO 300                                     
200   CONTINUE                                                         
      GOTO 910
300   CONTINUE                                                         
      TOL=NOSUP-NOSLOW                                                 
      ESTEP=(EMAX-EMIN)/DBLE(NP-1)                                     
      ELOW=EMIN+DBLE(IPLOW-1)*ESTEP                                    
      DNOS=NOSUP-NOSLOW                                                
      IF(DNOS.NE.0.D0) THEN                                            
        EFERMI=ELOW+(RNTOT-NOSLOW)/(NOSUP-NOSLOW)*ESTEP                
      ELSE                                                             
        EFERMI=ELOW                                                    
      END IF                                                           
      IF(EFERMI-ELOW.LT.0) PRINT*,'ERROR IN EFI '                      
      RETURN                                                           
  900 CALL OUTERR('FERMI','EFERMI OUT OF ENERGY RANGE')
      CALL OUTERR('FERMI','STOP IN EFI')
      WRITE (ERRMSG,9000) EMIN
      CALL OUTERR('FERMI',ERRMSG)
      WRITE (ERRMSG,9010) NOS(1)
      CALL OUTERR('FERMI',ERRMSG)
      WRITE (ERRMSG,9020) EMAX
      CALL OUTERR('FERMI',ERRMSG)
      WRITE (ERRMSG,9030) NOS(NP)
      CALL OUTERR('FERMI',ERRMSG)
      WRITE (ERRMSG,9040) ADD
      CALL OUTERR('FERMI',ERRMSG)
      WRITE (ERRMSG,9050) (SOS(I),I=100,1000,100)
      CALL OUTERR('FERMI',ERRMSG)
      WRITE (ERRMSG,9060) (NOS(I),I=100,1000,100)
      CALL OUTERR('FERMI',ERRMSG)
      STOP  'FERMI - Error'
  910 CALL OUTERR('FERMI','EFERMI NOT FOUND')
      CALL OUTERR('FERMI','STOP IN EFI')
      STOP  'FERMI - Error'
 9000 FORMAT('ENERGY OF LOWER BOUND                 :',f10.5)
 9010 FORMAT('NUMBER OF STATES AT THE LOWER BOUND   :',f10.5)
 9020 FORMAT('ENERGY OF UPPER BOUND                 :',f10.5)
 9030 FORMAT('NUMBER OF STATES AT THE UPPER BOUND   :',f10.5)
 9040 FORMAT('ADD ',f10.5)
 9050 FORMAT('SOS ',10f5.3)
 9060 FORMAT('NOS ',10f5.5)
      END                                                              
!     .....................................................TOTNOS......
      SUBROUTINE TOTNOS(VOL,E,EMIN,EMAX,NP,NOS,SOS)                    
!     **                                                              *
!     **  CALCULATES THE INTEGRATED DOS                               *
!     **  FOR ONE TETRAHEDRON ON THE ENERGYMESH                       *
!     **  INPUT :                                                     *
!     **    VOL         WEIGHT OF THIS TETRAHEDRON                    *
!     **    E           ENERGIES AT THE EDGEPOINTS                    *
!     **    EMIN        MINIMUM VALUE OF ENERGY MESH FOR NOS AND SOS  *
!     **    EMAX        MAXIMUM VALUE OF ENERGY MESH FOR NOS AND SOS  *
!     **    NP          NUMBER OF POINTS ON THE ENERGY MESH           *
!     **  OUTPUT:                                                     *
!     **    SOS         > NOS(E)+ SUM OVER E: SOS(E) =                *
!     **    NOS         > NUMBER OF STATES BELOW E                    *
!     **                                                              *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                              
      DOUBLE PRECISION NOS                                             
      DIMENSION E(4)                                                   
      DIMENSION NOS(NP),SOS(NP)                                        
!     -----------------------------------------------------------------
!     --  INTEGRATION WITHOUT FERMISURFACE                            -
!     -----------------------------------------------------------------
      X=DMIN1(E(1),E(2),E(3),E(4))                                     
      IF(X.GE.EMAX) THEN                                               
        RETURN                                                         
      END IF                                                           
      X=DMAX1(E(1),E(2),E(3),E(4))                                     
      IF(X.LE.EMIN) THEN                                               
        SOS(1)=SOS(1)+VOL                                              
        RETURN                                                         
      END IF                                                           
!     -----------------------------------------------------------------
!     --  ORDER ENERGIES                                              -
!     -----------------------------------------------------------------
      DO 100 I=1,3                                                     
      DO 100 J=I+1,4                                                   
      SVAR1=DMIN1(E(I),E(J))                                           
      SVAR2=DMAX1(E(I),E(J))                                           
      E(I)=SVAR1                                                       
      E(J)=SVAR2                                                       
100   CONTINUE                                                         
!     -----------------------------------------------------------------
!     --  CALCULATE UNCORRECTED INTEGRAL AS MEANVALUE                 -
!     -----------------------------------------------------------------
      E21=E(2)-E(1)                                                    
      if(e21.lt.1.d-10) e21=1.d-10
      E31=E(3)-E(1)                                                    
      E41=E(4)-E(1)                                                    
      E32=E(3)-E(2)                                                    
      E42=E(4)-E(2)                                                    
      E43=E(4)-E(3)                                                    
      ESTEP=(EMAX-EMIN)/DBLE(NP-1)                                     
      IMIN=IDINT(2.D0+(E(1)-EMIN)/ESTEP)                               
      IMIN=MAX0(1,IMIN)                                                
      IMAX=IDINT(1.D0+(E(2)-EMIN)/ESTEP)                               
      IMAX=MIN0(NP,IMAX)                                               
      EN=EMIN+ESTEP*(IMIN-1)                                           
      IF(IMAX.GE.IMIN) THEN                                            
        A=VOL/(E21*E31*E41)                                            
        DO 200 I=IMIN,IMAX                                             
        NOS(I)=NOS(I)+A*(EN-E(1))**3                                   
        EN=EN+ESTEP                                                    
200     CONTINUE                                                       
      END IF                                                           
      IMIN=MAX0(1,IMAX+1)                                              
      IMAX=INT(1.D0+(E(3)-EMIN)/ESTEP)                                 
      IMAX=MIN0(NP,IMAX)                                               
      IF(IMAX.GE.IMIN) THEN                                            
        A=VOL*E21**2/(E31*E41)                                         
        B=3.D0*VOL*E21/(E31*E41)                                       
        C=3.D0*VOL/(E31*E41)                                           
        D=-VOL/(E32*E41*E31*E42)*(E31+E42)                             
        DO 300 I=IMIN,IMAX                                             
        DE=EN-E(2)                                                     
!       NOS(I)=NOS(I)+A+B*DE+C*DE**2+D*DE**3                           
        NOS(I)=NOS(I)+A+DE*(B+DE*(C+D*DE))                             
        EN=EN+ESTEP                                                    
300     CONTINUE                                                       
      END IF                                                           
      IMIN=MAX0(1,IMAX+1)                                              
      IMAX=INT(1.D0+(E(4)-EMIN)/ESTEP)                                 
      IMAX=MIN0(NP,IMAX)                                               
      IF(E43.GT.0.D0) THEN                                             
        A=VOL                                                          
        D=VOL/(E41*E42*E43)                                            
        DO 400 I=IMIN,IMAX                                             
        NOS(I)=NOS(I)+A+D*(EN-E(4))**3                                 
        EN=EN+ESTEP                                                    
400     CONTINUE                                                       
      END IF                                                           
      IMIN=MAX0(1,IMAX+1)                                              
      IF(IMIN.GT.NP) RETURN                                            
      SOS(IMIN)=SOS(IMIN)+VOL                                          
      RETURN                                                           
      END                                                              
!     .....................................................SAMFAC......
      SUBROUTINE SAMFAC(NB,NKP,EB,EF,WGHT,W,NWX)                       
!     **                                                              *
!     **  CALCULATES SAMPLING WEIGHTS                                 *
!     **  INPUT :                                                     *
!     **    NB          NUMBER OF BANDS                               *
!     **    NKP         NUMBER OF K-POINTS                            *
!     **    EF          FERMI LEVEL                                   *
!     **    W           INTEGER WORK ARRAY                            *
!     **    NWX         LENTH OF WORK ARRAY                           *
!     **  OUTPUT :                                                    *
!     **    WGHT        SAMPLING WEIGHTS                              *
!     **                                                              *
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)                              
      IMPLICIT INTEGER (O)                                             
      DIMENSION EB(NB,NKP),WGHT(NB,NKP)                                
      DIMENSION E(4),WGHT0(4),IKP(4)                                   
      INTEGER W(NWX)                                                   
      CALL DEF0(NWX)                                                   
      CALL INITDR(WGHT,NB*NKP)                                         
      CALL TETR0(VOL0,NKP,NTET,MWRIT) 
      CALL DEFI(OWORK,5*MWRIT)                                         
      INIT=1                                                           
      DO 100 ITET=1,NTET                                               
      CALL TETR1(INIT,ITET,IWGHT,IKP,MWRIT,W(OWORK))                   
      VOL=VOL0*DBLE(IWGHT)                                             
      DO 100 IB=1,NB                                                   
      DO 110 I=1,4                                                     
      E(I)=EB(IB,IKP(I))                                               
110   CONTINUE                                                         
      CALL INITDR(WGHT0,4)                                             
      CALL WEIGHT(VOL,E,EF,WGHT0)                                      
      DO 120 I=1,4                                                     
      WGHT(IB,IKP(I))=WGHT(IB,IKP(I))+WGHT0(I)                         
120   CONTINUE                                                         
100   CONTINUE                                                         
      CALL RLSE(OWORK)                                                 
      RETURN                                                           
      END                                                              
!     .....................................................WHEIGT......
      SUBROUTINE WEIGHT(VOL,E,EF,WGHT)                                 
!     **                                                              *
!     **  CALCULATES THE WEIGHTS FOR TETRAHEDRON-SAMPLING             *
!     **  CORRESPONDING TO INTEGRATION OVER ONE TETRAHEDRON           *
!     **                                                              *
!     **  CORRECTION FOR THE NONLINEAR SHAPE INCLUDED IF ICOR=1       *
!     **                                                              *
!     **  AUTHOR : P.BLOECHL                                          *
!     **                                                              *
!     **    VOL.........VOLUME OF THIS TETRAHEDRON                    *
!     **    EF..........FERMI ENERGY                                  *
!     **    D...........KT (NOT USED)                                 *
!     **    E...........ENERGIES AT THE EDGEPOINTS                    *
!     **                                                              *
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)                              
      DIMENSION E(4),WGHT(4)                                           
      DIMENSION FA(4),FB(4),INDEX(4)                                   
      common /correct/ icor
!      DATA ICOR/1/                                                    
!     -----------------------------------------------------------------
!     --  INTEGRATION WITHOUT FERMISURFACE                            -
!     -----------------------------------------------------------------
      X=DMIN1(E(1),E(2),E(3),E(4))                                     
      IF(X.GE.EF) THEN                                                 
        CALL INITDR(WGHT,4)                                            
        RETURN                                                         
      END IF                                                           
      X=DMAX1(E(1),E(2),E(3),E(4))                                     
      IF(X.LE.EF) THEN                                                 
        VPRIME=.25D0*VOL                                               
        DO 10 I=1,4                                                    
        WGHT(I)=VPRIME                                                 
10      CONTINUE                                                       
        RETURN                                                         
      END IF                                                           
!     -----------------------------------------------------------------
!     --  ORDER ENERGIES                                              -
!     -----------------------------------------------------------------
!     -- INDEX HOLDS THE ORIGINAL POSITION OF THE ENERGIES AND WEIGHTS 
      DO 120 I=1,4                                                     
      INDEX(I)=I                                                       
120   CONTINUE                                                         
      DO 100 I=1,3                                                     
      IP=I                                                             
      DO 110 J=I+1,4                                                   
      IF(E(IP).GT.E(J)) IP=J                                           
110   CONTINUE                                                         
      IF(IP.GT.I) THEN                                                 
        X=E(IP)                                                        
        E(IP)=E(I)                                                     
        E(I)=X                                                         
        K=INDEX(IP)                                                    
        INDEX(IP)=INDEX(I)                                             
        INDEX(I)=K                                                     
      END IF                                                           
100   CONTINUE                                                         
!     -----------------------------------------------------------------
!     --  CALCULATE UNCORRECTED INTEGRAL AS MEANVALUE                 -
!     -----------------------------------------------------------------
      E21=E(2)-E(1)                                                    
      E31=E(3)-E(1)                                                    
      E41=E(4)-E(1)                                                    
      E32=E(3)-E(2)                                                    
      E42=E(4)-E(2)                                                    
      E43=E(4)-E(3)                                                    
      CALL INITDR(WGHT,4)                                              
      IF(EF.GT.E(1).AND.EF.LE.E(2)) THEN                               
        DE=EF-E(1)                                                     
        VPRIME=.25D0*VOL*DE**3/(E21*E31*E41)                           
        WGHT(1)=VPRIME*(4.D0-DE/E21-DE/E31-DE/E41)                     
        WGHT(2)=VPRIME*DE/E21                                          
        WGHT(3)=VPRIME*DE/E31                                          
        WGHT(4)=VPRIME*DE/E41                                          
!       ------  PARAMETERS FOR CORRECION                               
        DOS=3.D0*VPRIME*4.D0/(EF-E(1))                                 
      ELSE IF(EF.GT.E(2).AND.EF.LT.E(3)) THEN                          
        DE1=EF-E(1)                                                    
        DE2=EF-E(2)                                                    
        DE3=E(3)-EF                                                    
        DE4=E(4)-EF                                                    
!       ------  TETRAHEDRON X1,X2,X13',X14'                            
        VPRIME=VOL*DE1**2/(E41*E31)*.25D0                              
        WGHT(2)=VPRIME                                                 
        WGHT(3)=VPRIME*(DE1/E31)                                       
        WGHT(4)=VPRIME*(DE1/E41)                                       
        WGHT(1)=VPRIME*(3.D0-DE1/E41-DE1/E31)                          
!       ------  TETRAHEDRON X2,X13',X23',X14'                          
        VPRIME=.25D0*VOL*DE2*DE3*DE1/(E32*E31*E41)                     
        WGHT(1)=WGHT(1)+VPRIME*(2.D0-DE1/E31-DE1/E41)                  
        WGHT(2)=WGHT(2)+VPRIME*(2.D0-DE2/E32)                          
        WGHT(3)=WGHT(3)+VPRIME*(DE2/E32+DE1/E31)                       
        WGHT(4)=WGHT(4)+VPRIME*(DE1/E41)                               
!       ------  TETRAHEDRON X2,X23',X24',X14'                          
        VPRIME=.25D0*VOL*DE2**2*DE4/(E42*E32*E41)                      
        WGHT(1)=WGHT(1)+VPRIME*(1.D0-DE1/E41)                          
        WGHT(2)=WGHT(2)+VPRIME*(3.D0-DE2/E32-DE2/E42)                  
        WGHT(3)=WGHT(3)+VPRIME*(DE2/E32)                               
        WGHT(4)=WGHT(4)+VPRIME*(DE2/E42+DE1/E41)                       
!       ------  DOS=A+B*(EF-E2)+C*(EF-E2)**2                           
        DA=3.D0*VOL*E21/(E31*E41)                                      
        DB=6.D0*VOL/(E31*E41)                                          
        DC=-3.D0*VOL/(E32*E41*E31*E42)*(E31+E42)                       
        DOS=DA+DB*DE2+DC*DE2**2                                        
      ELSE IF(EF.GE.E(3).AND.EF.LT.E(4)) THEN                          
        DE=E(4)-EF                                                     
        VPRIME=.25D0*VOL*DE**3/(E41*E42*E43)                           
        VOL14=.25D0*VOL                                                
        WGHT(1)=VOL14-VPRIME*DE/E41                                    
        WGHT(2)=VOL14-VPRIME*DE/E42                                    
        WGHT(3)=VOL14-VPRIME*DE/E43                                    
        WGHT(4)=VOL14-VPRIME*(4.D0-DE/E41-DE/E42-DE/E43)               
!       ------  PARAMETERS FOR CORRECION                               
        DOS=3.D0*VPRIME*4.D0/(E(4)-EF)                                 
      ELSE                                                             
        GOTO 900
      END IF                                                           
!     -----------------------------------------------------------------
!     --  ADD CORRECTION FOR QUADRATIC DEVIATION                      -
!     -----------------------------------------------------------------
      IF(ICOR.EQ.1) THEN                                               
        DO 500 M=1,4                                                   
        DO 510 N=1,4                                                   
        WGHT(M)=WGHT(M)+.25D0*(E(N)-E(M))*DOS*.1D0                     
510     CONTINUE                                                       
500     CONTINUE                                                       
      END IF                                                           
!     -----------------------------------------------------------------
!     --  REORDER WEIGHTS                                             -
!     -----------------------------------------------------------------
      DO 600 I=1,4                                                     
      FA(INDEX(I))=WGHT(I)                                             
      FB(INDEX(I))=E(I)                                                
600   CONTINUE                                                         
      DO 610 I=1,4                                                     
      WGHT(I)=FA(I)                                                    
      E(I)=FB(I)                                                       
610   CONTINUE                                                         
      RETURN                                                           
  900 CALL OUTERR('FERMI','ERROR IN TETINT')
      STOP  'FERMI - Error'
      END                                                              
!     .....................................................TETR0.......
      SUBROUTINE TETR0(VOL,NKP,NTET,MWRIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      COMMON /IPROC/ IPROC                       
!     READ(15,REC=1)NKP,NTET,VOL,MWRIT,NREC                            
      REWIND 14                                                        
      READ(14,1234)NKP1,NTET,VOL,MWRIT,NREC
      if(nkp1.ne.nkp) goto 900
 1234 format(2i10,e20.12,2i10) 
      RETURN                                                           
 900  CALL OUTERR('FERMI', &
                  'number of k-points inconsistent when reading kgen')
      CALL OUTERR('FERMI','check IN1 and KGEN files!')
      STOP  'FERMI - Error'
      END                                                              
!     .....................................................TETR1.......
      SUBROUTINE TETR1(INIT,ITET,IWGHT,IKP,MWRIT,IWORK)                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                              
      DIMENSION IKP(4),IWORK(5*MWRIT)                                  
      SAVE IPOS,ntet                                                   
      IF(INIT.EQ.1) THEN                                               
!       READ(15,REC=1)NKP,NTET,V,MWRIT,NREC                            
        REWIND 14                                                      
        READ(14,1234)NKP,NTET,V,MWRIT,NREC                             
 1234 format(2i10,e20.12,2i10) 
        INIT=0                                                         
        IPOS=0                                                         
      END IF                                                           
      IREC=(ITET-1)/MWRIT+1                                            
      IF(ITET.GT.NTET) GOTO 900
      IF(IREC.NE.IPOS) THEN                                            
!       READ(15,REC=IREC+1)IWORK                                       
        READ(14,1235)IWORK                                             
 1235 format(6i10) 
!       ---------------------------------------------------------------
        NLEFT=NTET-(IREC-1)*MWRIT                                      
        NLEFT=5*MIN0(MWRIT,NLEFT)                                      
!       ---------------------------------------------------------------
        IPOS=IREC                                                      
      END IF                                                           
      IP=5*(ITET-1-(IPOS-1)*MWRIT)                                     
      IWGHT=IWORK(IP+1)                                                
      DO 100 I=1,4                                                     
      IKP(I)=IWORK(IP+1+I)                                             
100   CONTINUE                                                         
      RETURN                                                           
  900 CALL OUTERR('FERMI','ASK FOR NONEXISTING TETRAHEDRON')
      CALL OUTERR('FERMI','STOP IN TETR1')
      STOP  'FERMI - Error'
      END                                                              
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
!     ----  BLOCK MIXPIC                                            ---
!     ----  MIXED PICKLES                                           ---
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
!     .................................................................
      SUBROUTINE DEF0(NMAX)                                            
!     **                                                              *
!     **  ALLOCATES SPACE ON A (INTEGER WORK ARRAY)                   *
!     **  USE DEF0 BEFORE USE TO GIVE AVAILABLE SPACE ON WORK ARRAY   *
!     **  USE DEFI TO ALLOCATE SPACE FOR INTEGER ARRAYS               *
!     **  USE DEFDR TO ALLOCATE SPACE FOR INTEGER ARRAYS              *
!     **  USE RLSE TO RELEASE SPACE                                   *
!     **                                                              *
!     **  INPUT:                                                      *
!     **    NMAX        LENGTH OF INTEGER WORK ARRAY                  *
!     **    LENG        NUMBER OF ELEMENTS IN THE ARRAY TO BE         *
!     **                MAPPED ONTO THE WORK ARRAY (ENTRY: DEFI,DEFDR)*
!     **    ONAME       ALL ARRAYS FROM POINTER ONAME ARE DROPPED     *
!     **                (ENTRY: RLSE)                                 *
!     **  OUTPUT :                                                    *
!     **    ONAME       POINTER OF ARRAY TO BE ALLOCATED              *
!     **                                          (ENTRY: DEFI,DEFDR) *
!     **  REMARKS :                                                   *
!     **    AN INTEGER NUMBER IS ASSUMED TO HAVE 4 BYTES              *
!     **    A DOUBLE PRECISION NUMBER IS ASSUMED TO HAVE 8 BYTES      *
!     **                                                              *
      IMPLICIT INTEGER (O)                                             
      SAVE OMAX,OMAXX                                                  
      OMAX=1                                                           
      OMAXX=NMAX                                                       
      RETURN                                                           
!     =================================================================
      ENTRY DEFI(ONAME,LENG)                                           
      ONAME=OMAX                                                       
      OMAX=OMAX+LENG                                                   
      IF(OMAX.GT.OMAXX) GOTO 9999                                      
      RETURN                                                           
!     =================================================================
      ENTRY DEFDR(ONAME,LENG)                                          
      ONAME=OMAX                                                       
      OMAX=OMAX+LENG*2                                                 
      IF(OMAX.GT.OMAXX) GOTO 9999                                      
      RETURN                                                           
!     =================================================================
      ENTRY RLSE(ONAME)                                                
      IF(ONAME.LE.0.OR.ONAME.GT.OMAX) GOTO 9999                        
      OMAX=ONAME                                                       
      RETURN                                                           
9999  CONTINUE                                                         
      CALL OUTERR('FERMI','ERROR IN DEF0')
      STOP  'FERMI - Error'
      END                                                              
!     .................................................................
      FUNCTION DRVAL(NAME,INDEX)                                       
      DOUBLE PRECISION NAME(INDEX),drval                               
      DRVAL=NAME(INDEX)                                                
      RETURN                                                           
      END                                                              
!     .................................................................
      SUBROUTINE INITI(NAME,LENG)                                      
      DIMENSION NAME(LENG)                                             
      DO 100 I=1,LENG                                                  
      NAME(I)=0                                                        
100   CONTINUE                                                         
      RETURN                                                           
      END                                                              
!     .................................................................
      SUBROUTINE INITDR(NAME,LENG)                                     
      DOUBLE PRECISION NAME(LENG)                                      
      DO 100 I=1,LENG                                                  
      NAME(I)=0.D0                                                     
100   CONTINUE                                                         
      RETURN                                                           
      END                                                              
!
!je
!je
      SUBROUTINE FERMIG(weigh,elecn,nat,ef,nk,e)
!
!     J.E. 26.07.93
!
      USE param; USE parallel
      USE bandm
      use defs
      USE kpp1
      use struk
      USE xa2, only: weight,ne
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER *10    KNAME                         
      CHARACTER *80    FNAME,VECFN,VECFND,enefn,enefnd
      CHARACTER *67    ERRMSG
     
!                                                                      
      REAL*8 ztot,delta,efermi, emin, emax,nel, nel0, dnnt, ef0, ef1, de, wt, fac
      integer i, j, nnt, nc, ifl  
      REAL*8          :: e(2*nkpt,nume)
!                                                                      
      DIMENSION      WEIGH(2*NKPT,NUME),NB(2*NKPT)          
      common /mkef/delef,ts2      
      COMMON /PROC/    VECFN,VECFND,enefn,enefnd
      COMMON /IPROC/   IPROC 
!--------------------------------------------------------------------- 
      ztot=ELECN          
      DEL1 = EF
      nel = 0
      ITAP=30                                                          
      JSPIN=1                                                          
      TEST1=2.D-10 
      K=0                                                              
      MAXWAV=0   
      MINWAV=18000   
      SUMW=0.0D0              

!para begin
! ensure to mimick the proper vector file
! if running on multiple processors
      iloop=0
	iloop0=0
 777  continue
      if(IPROC.GT.0) then
         iloop=iloop+1
         iloop0=iloop0+1

         close(29)
         close(30)
         call mknam(FNAME,enefn,ILOOP)
         open(30,FILE=FNAME,STATUS='old',ERR=950)
         call mknam(FNAME,enefnd,ILOOP)
         open(29,FILE=FNAME,STATUS='unknown',ERR=950)
      endif
!para end
                                          
!                                                                      
!.....READ FROM TAPE 10,WRITTEN BY LAPW1                               
      DO I=1,NAT 
         READ(30,'(f9.5)') EMIST
         READ(30,'(f9.5)') EMIST
      ENDDO
      DO I=1,NAT                                                  
         READ(29,'(f9.5)',END=1002) EMIST
         READ(29,'(f9.5)',END=1002) EMIST
      ENDDO
      JSPIN=2                                                          
 1002 CONTINUE                                                         
!                                                                      
    4 K=K+1
      if(iloop0.ne.0)    KPP(ILOOP0)=K
!para begin
! testing
!      write(*,*)'reading k=',K
!para end
      IF(K.GT.2*NKPT) THEN
      WRITE(6,*) ' NUMBER OF K-POINTS .GT. NKPT =',NKPT
      WRITE(21,*) ' NUMBER OF K-POINTS .GT. NKPT =',NKPT
      STOP  'NKPT TOO SMALL'
      END IF                                                           
      READ(ITAP,'(3e19.12,a10,2i6,f5.1)',END=999) SS,T,ZZ,KNAME,N,NE(K),WEIGHT(K)               
      nmat=MAX(n,nmat)
      IF(NE(K).GT.NUME) THEN
         WRITE(6,*) 'NUME LT NE',NUME,NE(K)
         STOP 'NUME'
      END IF 
      SUMW=SUMW+WEIGHT(K)                                              
      IF(N.GT.MAXWAV) MAXWAV=N                                         
      IF(N.LT.MINWAV) MINWAV=N                                         
!      READ(ITAP) (KX(I),KY(I),KZ(I),I=1,N)                             
   14 READ(ITAP,*) NUM,E1                                                
      E(K,NUM)=E1                                                      
      if(itap.eq.30.and.e1.gt.ebmax(num)) ebmax(num)=e1
      if(itap.eq.30.and.e1.lt.ebmin(num)) ebmin(num)=e1
!      READ(ITAP) (A(I),I=1,N)                                          
      IF(NUM.EQ.NE(K)) GOTO 4                                          
      GOTO 14                                                          
!                                                                      
! ... CALCULATE WEIGHT FOR K-POINTS                                    
!                                                                      
!para begin
! we still have vector files to read

!  999 K=K-1
 999  CONTINUE
      K=K-1
!      write(6,*)'skipping last k point'
      nk=k                                                    
      IF(ILOOP.LT.IPROC) goto 777
!      write(6,*)'KPP(',ILOOP0,')=',KPP(ILOOP0)
!para end
                                                            
      IF(ITAP.EQ.30.AND.JSPIN.EQ.2) THEN                               
      ITAP=29
      SUMUP=SUMW                                                       
!para begin
         iloop=0
         IF(IPROC.EQ.0) GOTO 4
         GOTO 777
!para end                                                           
      ELSE IF(ITAP.EQ.29.AND.JSPIN.EQ.2) THEN                           
      nk=k/2                                                    
      DIF=2*SUMUP-SUMW                                                 
      IF(ABS(DIF).GT.TEST1) THEN                                       
         WRITE(6,*) ' SUMUP AND SUMDN IN FERMI NOT EQUAL, CHECK UP AND'
         WRITE(6,*) ' DOWN INPUTS OF LAPW1 ',SUMUP,SUMW,DIF            
      STOP                                                             
      END IF                                                           
      END IF                                                           
!                                                                      
!
!     Symmetriegewichte normieren 
!
      DO 8 KK=1,K                                                      
      WEIGHT(KK)=WEIGHT(KK)/SUMW                                  
   8  CONTINUE  
!
      DO 121 I=1,K
      DO 121 J=1,NE(I)
         WEIGH(I,J)=0.D0
 121  CONTINUE
!
!
!     guess Fermi-Energy and intitialize Gaussian Smearing
!
      delta = 1.E-6
      nnt = ztot/2.
      dnnt = ztot/2. - nnt
      if ((ztot-1.d0).lt.0.00001) then
         ef0 = E(1,1)
      else
         ef0 = E(1,nnt) + dnnt*(E(1,nnt+1) - E(1,nnt))
      endif
      efermi = ef0
      ifl = 0 
      nc = 0
!     Start loop 
 35   nel = 0.
      nc = nc +1
!     Test if maximun number of iterations is reached
      if (nc .ge. 100) then
        WRITE(6,*) 'FERMILEVEL NOT CONVERGED'
        STOP 'FERMIG'
      endif
!
      emax = ef0 - 1000.
      emin = 1000.
!     Loop over k-Points
      do 40 i = 1, K
!     loop over Bands
         do 40 j = 1, NE(i)
            if (E(i,j) .lt. emin) emin = E(i,j)
            de = (E(i,j) - efermi)/del1
            if (de .lt. -3.) then
              wt = 2.
            else if (de .lt. 0.) then
!e            wt = 1 + erfc(-de)
              wt = 2 - erfc(-de)
            else if (de.le.3.) then
!e            wt = 1 -erfc(de)
              wt = erfc(de)
	    else
              wt = 0.
            endif
            nel = nel + wt * WEIGHT(i)
            WEIGH(i,j) = wt*WEIGHT(i)
            if (wt .ge. 1.d-5) then
             if (emax .lt. E(i,j)) emax = E(i,j)
	    endif
 40   continue
!     special treatment of the first iteration step
      if (ifl .eq. 0) then
         ifl = 1
         efermi = efermi + 0.01
         nel0 = nel
         go to 35
      endif
!     all other steps
      if (abs(ztot-nel) .gt. delta) then
       ef1 = efermi
       fac = (nel0 - nel)/(ztot-nel)
       if (abs(fac) .le. 1.E-1) go to 41 
       efermi = efermi + (ef0 -efermi)/fac
       ef0 = ef1
       nel0 = nel
       go to 35
 41    ef0 = ef1
       nel0 = nel
       if ((ztot - nel) .le. 0.) then
          efermi = efermi - 0.01
       else
          efermi = efermi + 0.01
       endif
       go to 35
      endif
      emax = emax + 0.0001
      emin = emin - 0.0001
!
      EF=efermi
!
!     Caluculate energycorrection caused by gaussian smearing
!     Start loop over k-points (for spinpol.systems divide by 2
      eta=0.d0
      do 60 i = 1, K
!     Start loops over bands
         do 60 j = 1,NE(i)
           de = (E(i,j) - efermi)/del1
           de = de * de
           if (de.lt.15.) then
              eta = eta + 0.5*del1*EXP(-de)*WEIGHT(i)
           endif
!     end loops
 60   continue
      eta = -eta*2./sqrt(pi) / jspin
      if (myid.eq.0) then
      write(6,154) eta/2.d0 
      write(21,154) eta /2.d0
      endif
154   format(/,9x,'BANDENERGY CORRECTION:',f10.6)
      ts2=eta/2.d0
 70   continue
!
      return
!
 950  CONTINUE
      CALL OUTERR('LAPW2','error opening unit')
      WRITE (ERRMSG,9070) FNAME
      CALL OUTERR('LAPW2',ERRMSG)
      STOP 'FERMI - Error'
 9070 FORMAT('       filename: ',A50)
      end
!
!
      double precision function erfc(xx)
!
!
!     sandia mathematical program library
!     applied mathematics division 2613
!     sandia laboratories
!     albuquerque, new mexico  87185
!     control data 6600/7600  version 7.2  may 1978
!
!  Last modification: MOD_DATE
!
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!  * the primary document for the library of which this routine is     
!  * part is sand77-1441.                                              
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!
!     written by j.e. vogel from approximations derived by w.j. cody .
!
!     abstract
!
!         erfc(x) computes 2.0/sqrt(pi) times the integral from x to
!         infinity of exp(-x**2). this is done using rational approx-
!         imations.  eleven correct significant figures are provided.
!
!     description of parameters
!
!         x may be any real value
!
!     erfc is documented completely in sc-m-70-275.
!
!     .. Scalar Arguments ..
      double precision xx
!     ..
!     .. Local Scalars ..
      double precision a, r, sqpi, x, x2, xi2
      integer i
!     ..
!     .. Local Arrays ..
      double precision p1(4), p2(6), p3(4), q1(4), q2(6), q3(4)
!     ..
!     .. Intrinsic Functions ..
      intrinsic abs, exp
!     ..
!     .. Data statements ..
      data (p1(i),i=1,4)/242.6679552305318D0, 21.97926161829415D0, &
           6.996383488619136D0, -3.560984370181539D-02/, &
           (q1(i),i=1,4)/215.0588758698612D0, 91.16490540451490D0, &
           15.08279763040779D0, 1.0D0/, (p2(i),i=1,6)/22.898992851659D0, &
           26.094746956075D0, 14.571898596926D0, 4.2677201070898D0, &
           .56437160686381D0, -6.0858151959688D-06/, &
           (q2(i),i=1,6)/22.898985749891D0, 51.933570687552D0, &
           50.273202863803D0, 26.288795758761D0, 7.5688482293618D0, &
           1.0D0/, (p3(i),i=1,4)/-1.21308276389978D-2, &
           -.1199039552681460D0, -.243911029488626D0, &
           -3.24319519277746D-2/, (q3(i),i=1,4)/4.30026643452770D-02, &
           .489552441961437D0, 1.43771227937118D0, 1.0D0/
      data sqpi/.564189583547756D0/
!     ..
!
      x = abs(xx)
      x2 = x*x
      if (xx .lt. -6.0D0) then
         erfc = 2.0D0
      else if (x .gt. 25.8D0) then
         erfc = 0.0D0
      else if (x .gt. 4.0D0) then
         xi2 = 1.D0/x2
         r = xi2*(p3(1)+xi2*(p3(2)+xi2*(p3(3)+xi2*p3(4))))/ &
             (q3(1)+xi2*(q3(2)+xi2*(q3(3)+xi2*q3(4))))
         a = exp(-x2)*(sqpi+r)/x
         if (xx .lt. 0.0D0) then
            erfc = 2.0D0 - a
         else
            erfc = a
         end if
      else if (x .gt. .46875D0) then
         a = exp(-x2)*(p2(1)+x*(p2(2)+x*(p2(3)+x*(p2(4)+x*(p2(5)+ &
             x*p2(6))))))
         a = a/(q2(1)+x*(q2(2)+x*(q2(3)+x*(q2(4)+x*(q2(5)+x*q2(6))))))
         if (xx .le. 0.0D0) then
            erfc = 2.0D0 - a
         else
            erfc = a
         end if
      else
         a = x*(p1(1)+x2*(p1(2)+x2*(p1(3)+x2*p1(4))))
         a = a/(q1(1)+x2*(q1(2)+x2*(q1(3)+x2*q1(4))))
         if (xx .lt. 0.D0) a = -a
         erfc = 1.0D0 - a
      end if
!
      return
!
      end
!
