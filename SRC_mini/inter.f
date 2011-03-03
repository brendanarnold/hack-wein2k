      subroutine inter(alph,iterho,nat,itecon)
!                                                                       
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
!     (ifft**3+(ifft+1)**3)*2 must be gt.(2*NWAV+NRAD*NCOM*NATO*2)*8
      PARAMETER (IFFT1=90)
      PARAMETER (IFFT2=90)
      PARAMETER (IFFT3=90)
!
      DIMENSION      alph(nato+1,4,2)
      CHARACTER*19     nkktext
      CHARACTER*79   MARGN                                             
      CHARACTER*67   ERRMSG
      COMPLEX*16     ROKMIX,ROKNEW,ROKOLD,ROKOLD1,ROKOLD2
      complex*16 cfft,ust
!                                                                       
      COMMON /SYMETR/  TAU(3,NSYM),KVCNEW(3,NWAV),IZ(3,3,NSYM),             &
                       INUM(NSYM),IORD                                  
      COMMON /GENER/   BR1(3,3),BR2(3,3),avec(3,3)                    
      COMMON /XA/ CLMOLD(NRAD,NCOM,NATO,2),CLMOLD1(NRAD,NCOM,NATO,2), &
                CLMOLD2(NRAD,NCOM,NATO,2),CLMNEW(NRAD,NCOM,NATO,2), &
                ROKNEW(NWAV,2),ROKOLD(NWAV,2), &
                ROKOLD1(NWAV,2),ROKOLD2(NWAV,2),rokmix(NWAV,2),             &
                R(NRAD),RNOT(NATO),              &
                DX(NATO),ZZ(NATO),JRI(NATO),                       &
                LM(2,NCOM,NATO),LMMAX(NATO)   
      COMMON /XA2/ cfft(IFFT1*IFFT2*IFFT3), &
            UST((IFFT1+1)*(IFFT2+1)*(IFFT3+1))
!
!        Common blocks
!
!
!        LATTIC  - lattice type
!        NAME(i) - name of atom i in the unit-cell (e.g. 'Titanium')
!
      CHARACTER*4        LATTIC
      CHARACTER*10       NAME(NATO)
      COMMON  /CHAR/     LATTIC, NAME
      SAVE    /CHAR/
!
!        ALAT(1:3)   - lattice constants (a,b,c)
!        ALPHA(1:3)  - angles of the unit-cell's unit-vectors
!        IATNR(i)    - atom index of (inequivalent) atom i
!                      also indicates cubic and non-cubic symmetries
!                      IATNR(i) .GT. 0 ... cubic symmetry
!                      IATNR(i) .LT. 0 ... non-cubic symmetry
!        MULT(i)     - number of equivalent atoms for inequiv. atom i
!        PIA(1:3)    - reciprocal lattice constants (2*pi/a, 2*pi/b,
!                      2*pi/c)
!        POS(1:3,i)  - position vector of atom i (including equivalent
!                      atoms) in the unit-cell
!        RMT(i)      - muffin tin radius of atom i
!        V(i)        - relative muffin tin spherevolume for atom i
!        VI          - volume of the direct lattice unit-cell
!
      INTEGER            IATNR(NATO), MULT(NATO)
      DOUBLE PRECISION   VI
      DOUBLE PRECISION   ALAT(3), ALPHA(3), PIA(3), POS(3,NDIF)
      DOUBLE PRECISION   RMT(NATO), V(NATO)
      COMMON  /STRUK/    POS, ALAT, ALPHA, RMT, V, PIA, VI, IATNR, MULT
      SAVE    /STRUK/
      dimension qpw(2),qel(nato,2)
!---------------------------------------------------------------------  
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!.....LOOK, WHICH CHARGE-DENSITY TAPES EXIST                            
!.....READ PRATT-FACTORS; INPUT RHO=RHONEW*PRATT + (1-PRATT)*RHOOLD     
!                                                                       
      JSPIN=1                                                           
      READ(10,2032,IOSTAT=ITAPE)                                        
      READ(60,2032,IOSTAT=ITAPE)                                        
      IF(ITAPE.EQ.0) JSPIN=2
!                                                                       
!.....MIX DENSITY IN SPHERES                                            
!                                                                       
      DO 15 JATOM=1,NAT   
         do 16 jrj=1,jri(jatom)
 16      r(jrj)=rnot(jatom)*exp(dx(jatom)*(jrj-1))
!
         READ(10,1980)                                                    
         READ(10,2000) LMMAX(jatom)
         DO 40 LM1=1,LMMAX(jatom)                                             
            READ(10,2010) LM(1,LM1,JATOM),LM(2,LM1,JATOM) 
            READ(10,2021) (CLMNEW(J,LM1,JATOM,1),J=1,JRI(JATOM))       
 40      READ(10,2031)                                                  
         READ(10,2033)                                                  
!                                                                       
         IF(JSPIN.EQ.2) THEN                                              
            READ(60,1980)                                                    
            READ(60,2000) LMMAX1                                             
            DO 540 LM1=1,LMMAX1                                            
               READ(60,2031)                                               
               READ(60,2021) (CLMNEW(J,LM1,JATOM,2),J=1,JRI(JATOM))     
540         READ(60,2031)                                                  
            READ(60,2033)                                                  
         END IF                                                           
  15  CONTINUE 
!                                                                  
!.....READ DENSITY IN INTERSTITIAL                                       
!                                                                       
      READ(10,1980)                                                     
!      READ(10,2060) NKKNEW 
  read(10,'(/,a19)') nkktext
  read(nkktext,'(9x,i10)',err=6767) nkknew
  goto 6768
 6767 read(nkktext,'(13x,i6)') nkknew
 6768 continue
      if(nkknew.gt.nwav) goto 900                                             
      DO 150 J=1,NKKNEW                                              
 150  READ(10,2071) (KVCNEW(JX,J),JX=1,3),ROKNEW(J,1)                
      IF(JSPIN.EQ.2) THEN                                               
         READ(60,1980)                                                     
!         READ(60,2060) NKKNEW                                              
  read(60,'(/,a19)') nkktext
  read(nkktext,'(9x,i10)',err=6769) nkknew
  goto 6770
 6769 read(nkktext,'(13x,i6)') nkknew
 6770 continue
         DO 5150 J=1,NKKNEW                                             
 5150    READ(60,2071) (KVCNEW(JX,J),JX=1,3),ROKNEW(J,2)                
      ENDIF                                                            
!                                                                       
!                                                                       
!..................................................................
!.....integrate charges, check normalization of charge density
!.....calc iff parameters
      call iffpar(ifft1,ifft2,ifft3,iff1,iff2,iff3,nkknew)
!.....calc stepfunction
      CALL REAN0(NKKnew,KvcNEW,NAT,IFF1,IFF2,IFF3,UST)
!
      write(6,*) '      CHARGES OF NEW CHARGE DENSITY'
      margn='N'
      call integr(nkknew,kvcNEW,roknew,nat,jspin,jri,rnot, &
        dx,r,clmnew,iff1,iff2,iff3,qpw,qel,cfft,ust,margn)
      margn='01'
      call normal(roknew,nkknew,jspin,nat,mult,zz,qpw,qel,qmx,margn)
!
!     MIX DENSITIES
!
      CALL READ3(nkknew,iterho,nat,jspin)
      if (iterho.lt.itecon) then
         write(6,*) 'NO INTERPOLATION, TO LESS ITERATIONS', &
                     itecon,iterho
         goto 8000
      endif
      call integr(nkknew,kvcNEW,rokold,nat,jspin,jri,rnot, &
        dx,r,clmold,iff1,iff2,iff3,qpw,qel,cfft,ust,margn)
      margn='01'
      call normal(rokold,nkknew,jspin,nat,mult,zz,qpw,qel,qmx,margn)
      if (itecon.lt.2) goto 8100
      call integr(nkknew,kvcNEW,rokold1,nat,jspin,jri,rnot, &
        dx,r,clmold1,iff1,iff2,iff3,qpw,qel,cfft,ust,margn)
      margn='02'
      call normal(rokold1,nkknew,jspin,nat,mult,zz,qpw,qel,qmx,margn)
      if (itecon.lt.3) goto 8100
      call integr(nkknew,kvcNEW,rokold2,nat,jspin,jri,rnot, &
        dx,r,clmold2,iff1,iff2,iff3,qpw,qel,cfft,ust,margn)
      margn='02'
      call normal(rokold2,nkknew,jspin,nat,mult,zz,qpw,qel,qmx,margn)
 8100 do 44 jatom=1,nat
      do 44 lm1=1,lmmax(jatom)
      DO 44 J=1,JRI(JATOM)                                        
 44   CLMNEW(J,LM1,JATOM,1)=CLMNEW(J,LM1,JATOM,1)*alph(jatom,1,1)+   &
                            CLMOLD(J,LM1,JATOM,1)*alph(jatom,2,1)+ &
                            CLMOLD1(J,LM1,JATOM,1)*alph(jatom,3,1)+ &
                            CLMOLD2(J,LM1,JATOM,1)*alph(jatom,4,1) 
      IF(JSPIN.EQ.2) THEN      
         do 544 jatom=1,nat
         do 544 lm1=1,lmmax(jatom)
         DO 544 J=1,JRI(JATOM)                                       
544      CLMNEW(J,LM1,JATOM,2)=CLMNEW(J,LM1,JATOM,2)*alph(jatom,1,2)+   &
                               CLMOLD(J,LM1,JATOM,2)*alph(jatom,2,2)+ &
                               CLMOLD1(J,LM1,JATOM,2)*alph(jatom,3,2)+ &
                               CLMOLD2(J,LM1,JATOM,2)*alph(jatom,4,2) 
      ENDIF                                                       
      DO 135 J=1,NKKNEW                                                 
      DO 135 ISPIN=1,JSPIN                                              
         ROKMIX(J,ISPIN)=ROKNEW(J,ISPIN)*alph(nat+1,1,ispin)+ &
                         ROKOLD(J,ISPIN)*alph(nat+1,2,ispin)+ &
                         ROKOLD1(J,ISPIN)*alph(nat+1,3,ispin)+ &
                         ROKOLD2(J,ISPIN)*alph(nat+1,4,ispin) 
 135  CONTINUE                                                          
!
      write(6,*) '      CHARGES OF MIXED CHARGE DENSITY'
         margn='C'
      call integr(nkknew,kvcNEW,rokmix,nat,jspin,jri,rnot, &
        dx,r,clmnew,iff1,iff2,iff3,qpw,qel,cfft,ust,margn)
      margn='02'
      call normal(rokmix,nkknew,jspin,nat,mult,zz,qpw,qel,qmx,margn)
!.....................................................................
!.....WRITE NEW INPUT DENSITY TO TAPE51,TAPE 52, TAPE11                 
!                                                                       
!     REWIND 11                                                         
      WRITE(51,1970) 1                                               
      WRITE(51,78) LATTIC,ALAT(1),ALAT(2),ALAT(3),JRI(1)
      WRITE(51,77)                                
      DO 560 JATOM=1,NAT                                                
         WRITE(51,1990) JATOM                                           
         WRITE(51,2001) LMMAX(JATOM)                                    
         DO 565 LM1=1,LMMAX(JATOM)                                      
            WRITE(51,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)              
            WRITE(51,2021) ( CLMNEW(J,LM1,JATOM,1), J=1,JRI(JATOM) )    
565         WRITE(51,2031)                                              
560   WRITE(51,2033)                                                    
!                                                                       
      WRITE(51,2051)                                                    
      WRITE(51,1980)                                                    
      WRITE(51,2061) NKKNEW                                             
      DO 5155 J=1,NKKNEW                                                
5155  WRITE(51,2076) (KVCNEW(JX,J),JX=1,3),ROKMIX(J,1),ROKNEW(J,1)      
!....................................................................
              IF(JSPIN.EQ.2) THEN                                           
!
      WRITE(52,1970) 1                                               
      WRITE(52,78) LATTIC,ALAT(1),ALAT(2),ALAT(3),JRI(1)
      WRITE(52,77)                               
      DO 660 JATOM=1,NAT                                                
         WRITE(52,1990) JATOM                                           
         WRITE(52,2001) LMMAX(JATOM)                                    
         DO 665 LM1=1,LMMAX(JATOM)                                      
            WRITE(52,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)              
            WRITE(52,2021) ( CLMNEW(J,LM1,JATOM,2), J=1,JRI(JATOM) )    
665         WRITE(52,2031)                                              
660   WRITE(52,2033)                                                    
!                                                                       
      WRITE(52,2051)                                                    
      WRITE(52,1980)                                                    
      WRITE(52,2061) NKKNEW                                             
      DO 6155 J=1,NKKNEW                                                
6155  WRITE(52,2076) (KVCNEW(JX,J),JX=1,3),ROKMIX(J,2),ROKNEW(J,2)      
!                                                                       
      WRITE(11,1970) 1                                               
      WRITE(11,78) LATTIC,ALAT(1),ALAT(2),ALAT(3),JRI(1)
      WRITE(11,77)                               
      DO 60 JATOM=1,NAT                                                 
         WRITE(11,1990) JATOM                                           
         WRITE(11,2001) LMMAX(JATOM)                                    
         DO 65 LM1=1,LMMAX(JATOM)                                       
            WRITE(11,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)              
      DO 63 J=1,JRI(JATOM)                                              
 63   CLMNEW(J,LM1,JATOM,1)=CLMNEW(J,LM1,JATOM,1)+CLMNEW(J,LM1,JATOM,2) 
            WRITE(11,2021) ( CLMNEW(J,LM1,JATOM,1), J=1,JRI(JATOM) )    
 65         WRITE(11,2031)                                              
 60   WRITE(11,2033)                                                    
!                                                                       
      WRITE(11,2051)                                                    
      WRITE(11,1980)                                                    
      WRITE(11,2061) NKKNEW                                             
      DO 155 J=1,NKKNEW                                                 
      ROKMIX(J,1)=ROKMIX(J,1)+ROKMIX(J,2)                               
      ROKNEW(J,1)=ROKNEW(J,1)+ROKNEW(J,2)                               
 155  WRITE(11,2076) (KVCNEW(JX,J),JX=1,3),ROKMIX(J,1),ROKNEW(J,1)      
      END IF                                                            
! ..................................................................
!                                                                       
 8000 RETURN
!
!        error handling
!
  900 INFO = 0
!
!        NWAV too small
!
      CALL OUTERR('MIXER','nwav too small')
      GOTO 999
!
  940 INFO = 4
      CALL OUTERR('MIXER','Maximum number of K-points exceeded.')
      GOTO 999
!                                                                       
  999 STOP 'INTER - Error'
!                                                                       
 77   FORMAT(1X,9(F9.7))                                                
 78   FORMAT(1X,20X,A4,3F10.6,I5)                                       
 700  FORMAT(I3,A76)                                                    
 701  FORMAT(I3,A4)                                                     
 702  FORMAT(A5,13X,A60)                                                
 703  FORMAT(A8,F8.5,A60)                                               
 704  FORMAT('CONTINUE',F8.5,A60)                                       
 705  FORMAT(F10.3,6X,A60)                                              
 743  FORMAT(3X,A76)                                                    
 1000 FORMAT(A20)                                                       
 1040 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
 1041 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
 1055 FORMAT(I4)                                                        
 1061 FORMAT(7X,A5,' MIXING SCHEME WITH',F6.3)                          
 1070 FORMAT(I3,'. ITERATION',//)                                       
 1970 FORMAT(3X,'             TOTAL CHARGE DENSITY GENERATED ',        &
             'BY',I3,'. ITERATION ')                 
 1980 FORMAT(3X)                                                        
 1990 FORMAT(3X,'ATOM NUMBER=',I3,5X,10A4)                             
 2000 FORMAT(15X,I3//)                                                  
 2001 FORMAT(3X,'NUMBER OF LM',I3//)                                   
 2010 FORMAT(15X,I3,5X,I2/)                                             
 2011 FORMAT(3X,'CLM(R) FOR L',I3,3X,'M=',I2/)                         
 2021 FORMAT(3X,4E19.12)
 2031 FORMAT(/)                                                         
 2032 FORMAT(//)                                                        
 2033 FORMAT(///)                                                       
 2051 FORMAT(3X,'MIXED TOTAL CHARGE DENSITY IN INTERSTITIAL',            &
            3X,'NEW TOTAL CHARGE DENSITY')                              
 2060 FORMAT(/,13X,I6)                                                  
 2061 FORMAT(1X,'        ',I10,1X,'NUMBER OF PW')                                     
 2071 FORMAT(3X,3I5,2E19.12)
 2076 FORMAT(3X,3I5,2E19.12,2e11.3)
 5020 FORMAT(6F10.7)
 9050 FORMAT('MULT(',I3,'=',I3,', INDEX=',I3)
      END                                                               
