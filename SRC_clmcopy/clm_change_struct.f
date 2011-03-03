!$hp9000_800 intrinsics on
      PROGRAM clmcopy                                                     
! writes new clmsum file when struct and struct_new files have different jri
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
      CHARACTER*4    IREL,LATTIC,cform                                 
      CHARACTER*5    vech(20,3)                                                
      CHARACTER*10   ANAME                                              
      CHARACTER*11   STATUS,FORM                                        
      CHARACTER*20   TITLE                                              
      CHARACTER*77   FNAME
      CHARACTER*79   MARGN
      CHARACTER*91   MARGN91
      LOGICAL        jatch2,forces
!
      integer,allocatable :: jatch(:,:),lmch(:,:,:,:),lmchm(:)             
      complex*16 tk,imag                                
      COMPLEX*16,allocatable ::  ROKVL(:,:)                
      INTEGer,allocatable :: KVECVL(:,:)
      real*8,allocatable :: CLMNEW(:,:,:,:),fhelp(:,:),fsuhelp(:,:),clmfac(:,:)
      INTEGer,allocatable :: LM(:,:,:),LMMAX(:),JRI(:),JRInew(:),MULT(:)
!
      integer msym(3,3),krot(3),krotst(3,nsym)
      real*8  taum(3)
      DIMENSION TAU(3,NSYM),fakc(20,2)
      COMMON /SYMM/  LSYM(3,3,NSYM),IORD,INUM(NSYM)                     
      DATA      tiny/1.0d-15/
!---------------------------------------------------------------------  
!                                                                       
      call getarg(2,fname)
      if(fname.eq.'      ') call getarg(1,fname)
      OPEN(1,FILE=fname,STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM, &
      ERR=8002)
      GOTO 8003
 8000 WRITE(*,*) ' ERROR IN OPENING MIXER.DEF !!!!'
      STOP 'MIXER.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, &
      '  FORM:',FORM
      STOP 'OPEN FAILED'
 8001 CONTINUE
      CLOSE (1)

      tpi=2.d0*acos(-1.d0)
      imag=(0.d0,1.d0)

!.....READ TAPE20=POTE                                                  
!                                                                       
      READ(20,1000) TITLE                                               
      READ(20,1010) LATTIC,NAT,cform,IREL              
      READ(21,1000) TITLE                                               
      READ(21,1010) LATTIC,NAT,cform,IREL              
      nato=nat                       
      allocate ( CLMNEW(NRAD,NCOM,NATO,2),fhelp(nato,4),fsuhelp(nato,4))
      allocate ( LM(2,NCOM,NATO),LMMAX(NATO),JRI(NATO),JRInew(NATO),MULT(NATO))
      allocate ( jatch(nato,2),lmch(nato,ncom,2,2),lmchm(nato),clmfac(nato,ncom))                                              
      READ(20,1030) AA,BB,CC                                            
      READ(21,1030) AA,BB,CC                                            
!      if(cform.eq.'NEW ') then
        assign 2021 to iform1
        assign 2071 to iform2
        assign 2076 to iform3
!      else
!        assign 2020 to iform1
!        assign 2070 to iform2
!        assign 2075 to iform3
!      end if                                    
 2020 FORMAT(3X,5E13.7)                                                 
 2021 FORMAT(3X,4E19.12)                                                 
 2070 FORMAT(3X,3I5,2E15.7)                                             
 2071 FORMAT(3X,3I5,2E19.12)                                             
 2075 FORMAT(3X,3I5,4E15.7)                                             
 2076 FORMAT(3X,3I5,2E19.12,2e11.3)                                             
      INDEX=0                                                           
!     INDEX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE    
!     NOTEQUIVALENT ATOMS                                               
!                                                                       
      DO 5 JATOM = 1,NAT                                                
         INDEX=INDEX+1                                                  
         READ(20,1040) IATNR,POS1,POS2,             &
                       POS3,MULT(JATOM),ISPLIT            
         READ(21,1040) IATNR,POS1,POS2,             &
                       POS3,MULT(JATOM),ISPLIT            
            MULTND=MULT(JATOM)-1                                        
            DO 10 M=1,MULTND                                            
            INDEX=INDEX+1                                               
            READ(20,1041) IATNR,POS1,POS2,          &
                          POS3                                   
            READ(21,1041) IATNR,POS1,POS2,          &
                          POS3                                   
 10         CONTINUE                                                    
         READ(20,1050) ANAME,JRI(JATOM),RNOT,RMTT,         &
                       ZZ                                        
         READ(20,*)                                                     
         READ(20,*)                                                     
         READ(20,*)                                                     
         READ(21,1050) ANAME,JRInew(JATOM),RNOT,RMTT,         &
                       ZZ                                        
         READ(21,*)                                                     
         READ(21,*)                                                     
         READ(21,*)                                                     
 5    CONTINUE                                                          
!                                                                       
!                                                                       
!.....LOOK, WHICH CHARGE-DENSITY TAPES EXIST                            
!                                                                       
      READ(17,2032,IOSTAT=ITAPE)                                        
      write(18,2032,IOSTAT=ITAPE)                                        
!                                                                       
      DO 15 JATOM=1,NAT                                                 
            READ(17,1980)                                               
            READ(17,2000) LMMAX(JATOM)                                  
               if(lmmax(jatom).gt.ncom) stop 'ncom too small'
            DO 20 LM1=1,LMMAX(JATOM)                                    
            READ(17,2010) LM(1,LM1,JATOM),LM(2,LM1,JATOM)               
            READ(17,iform1) (CLMNEW(J,LM1,JATOM,1),J=1,JRI(JATOM))
 20         READ(17,2031)                                               
            READ(17,2033)                                               
 15      continue
!
! calculate clmfiles
!
      DO 16 JATOM=1,NAT    
        jato=jatom
        write(18,1980)                                               
        write(18,2001) LMMAX(JATOm)                                  
        DO 21 LM1=1,LMMAX(JATOm)                                    
          write(18,2011) LM(1,LM1,JATOm),LM(2,LM1,JATOm)               
            WRITE(18,IFORM1) (CLMNEW(J,LM1,JATOm,1),J=1,JRInew(JATOm))
            write(18,2031)                                               
 21         continue
            write(18,2033)                                               
 16   continue                      
!
         READ(17,1980)                                                  
         READ(17,2060) NKKVL                                            
         write(18,1980)                                                  
         write(18,2061) NKKVL   
         nwav=nkkvl
      allocate (  ROKVL(NWAV,2))                
      allocate ( KVECVL(3,NWAV))

! reading fourier coeff.
         DO  J=1,NKKVL                                               
           READ(17,iform2)  (KVECVL(JX,J),JX=1,3),ROKVL(J,1)              
           write(18,iform2)  (KVECVL(JX,J),JX=1,3),ROKVL(J,1)              
         enddo       
      STOP 'clmcopy END'                                                  
!                                                                       
!                                                                       
 77   FORMAT(1X,9(F9.7))                                                
 78   FORMAT(1X,A20,A4,3F10.6,I5)                                       
 700  FORMAT(A79)                                                    
 791  FORMAT(A91)                                                    
 710  FORMAT(4x,i2)                                                    
 720  FORMAT(a4,i2.2,a8,i3,a62)
 730  FORMAT(a4,i2.2,a29,i3,a41) 
 740  FORMAT(a4,i2.2,a25,i3,a45) 
 750  FORMAT(a4,i2.2,a72) 
 751  FORMAT(a4,i2.2,a84) 
 760  FORMAT(a4,i2.2,a20,i3,a48) 
 770  FORMAT(15X,4f15.3)
 780  FORMAT (7x,'VALENCE-FORCE IN mRy = |F|',8x,'Fx',13x, &
              'Fy',13x,'Fz')
 790  FORMAT (':FVA',i2.2,':',1x,i2,'.ATOM',4f15.3)
 1790  FORMAT (':FSU',i2.2,':',1x,i2,'.ATOM',4f15.3)
 1000 FORMAT(A20)                                                       
 1010 FORMAT(A4,23X,I3,1x,a4,/,13X,A4)                                        
 1030 FORMAT(3F10.7)                                                    
 1040 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
 1041 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
 1055 FORMAT(I4)                                                        
 1056 FORMAT(3(3I2,F5.2,/),I8)                                          
 1980 FORMAT(3X)                                                        
 1990 FORMAT(3X,'ATOM NUMBER ',I3,5X,10A4)                             
 2000 FORMAT(15X,I3//)                                                  
 2001 FORMAT(3X,'NUMBER OF LM',I3//)                                   
 2010 FORMAT(15X,I3,5X,I2/)                                             
 2011 FORMAT(3X,'CLM(R) FOR L',I3,3X,'M=',I2/)                         
 2031 FORMAT(/)                                                         
 2032 FORMAT(//)                                                        
 2033 FORMAT(///)                                                       
 2060 FORMAT(/,13X,I6)                                                  
 2061 FORMAT(/,1X,'NUMBER OF PW',I6)                                     
 3010 format(3a4)
 3011 format(3a4,2f4.1)
      END                                                               
