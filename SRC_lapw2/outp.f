      SUBROUTINE OUTP(eferm,so,nspin1,fnamehelp,doos)
!                                                                       
!     READS PARTIAL CHARGES FROM TAPES WRITTEN BY L2MAIN, WRITES        
!     PARTIAL CHARGES  IN COMPACT FORM ON OUTPUT AND ON TAPE QTL        
!
!     Modified by LDM to include DOOS option (density of occupied states)                                                                       
      USE param
      use defs
      USE char
      use struk; USE com
      USE xdos
      IMPLICIT REAL*8 (A-H,O-Z)

      CHARACTER*5      AQTL                                       
      CHARACTER*7      LL,LDUMMY                                        
      CHARACTER*45,allocatable  ::     TEXT(:)
      CHARACTER*67     ERRMSG                                            
      CHARACTER *80    fnamehelp,fnamehelp2
      LOGICAL          so,doos                                              
!                                                                       
      DIMENSION       QTL1(3,22),LL(22)      
      INTEGER,allocatable   ::  ILINP(:),ILOUT(:),IEND(:),ISTART(:)
      INTEGER,allocatable   ::  IEND1(:),ISTAR1(:)
      real*8,allocatable    ::  sonorm(:,:)
      lxdos2=(lxdos+1)*(lxdos+1)
!--------------------------------------------------------------------   
! 
      allocate (TEXT(NAT),ILINP(NAT),ILOUT(NAT),IEND(NAT),ISTART(NAT))
      allocate (IEND1(NAT),ISTAR1(NAT))
      kso=0
      if(so) kso=2                                                             
      AQTL='QTL  '                                                      
      DO 19 JATOM=1,NAT                                                 
         JSPLIT=ISPLIT(JATOM)                                           
         IF(JSPLIT.EQ.0) THEN                                           
            ILINP(JATOM) =7                                             
            IEND(JATOM)  =3                                             
            ISTART(JATOM)=4                                             
            IEND1(JATOM) = ILINP(JATOM)                                 
            ISTAR1(JATOM)=ILINP(JATOM)+1                                
            ILOUT(JATOM) =4                                             
            text(jatom)='tot,0,1,2,3' 
         ELSE IF(JSPLIT.EQ.99) THEN                                      
            ILINP(JATOM) =7                                             
            IEND(JATOM)  =3                                             
            ISTART(JATOM)=4                                             
            IEND1(JATOM) = ILINP(JATOM)                                 
            ISTAR1(JATOM)=ILINP(JATOM)+1                                
            ILOUT(JATOM) =4                                             
            text(jatom)='tot,0,1,2,3,xdos(i,j),j=1,i),i=1,lxdos2)' 
         ELSE IF(JSPLIT.EQ.88) THEN                                      
            ILINP(JATOM) =7                                             
            IEND(JATOM)  =3                                             
            ISTART(JATOM)=4                                             
            IEND1(JATOM) = ILINP(JATOM)                                 
            ISTAR1(JATOM)=ILINP(JATOM)+1                                
            ILOUT(JATOM) =4                                             
            text(jatom)='tot,0,1,2,3,xdos(i,i),i=1,lxdos2)' 
         ELSE IF(JSPLIT.EQ.1) THEN                                      
            ILINP(JATOM) =9                                             
            IEND(JATOM)  =2                                             
            IEND1(JATOM) = ILINP(JATOM)                                 
            ISTAR1(JATOM)=ILINP(JATOM)+1                                
            ISTART(JATOM)=5                                             
            ILOUT(JATOM) =6                                             
            text(jatom)='tot,0,1,PZ,PX+PY,2,3' 
         ELSE IF(JSPLIT.EQ.2) THEN                                      
               ILINP(JATOM) =9                                          
               IEND(JATOM)  =3                                          
               ISTART(JATOM)=6                                          
               IEND1(JATOM) = ILINP(JATOM)                              
               ISTAR1(JATOM)=ILINP(JATOM)+1                             
               ILOUT(JATOM) =6                                          
               text(jatom)='tot,0,1,2,D-eg,D-t2g,3' 
            ELSE IF(JSPLIT.EQ.3) THEN                                   
               ILINP(JATOM) =10                                         
               IEND(JATOM)  =3                                          
               ISTART(JATOM)=7                                          
               IEND1(JATOM) = ILINP(JATOM)                              
               ISTAR1(JATOM)=ILINP(JATOM)+1                             
               ILOUT(JATOM) =7                                          
               text(jatom)='tot,0,1,2,DZ2,DX2Y2+DXY,DXZ+DYZ,3' 
            ELSE IF(JSPLIT.EQ.-2) THEN                                  
               ILINP(JATOM) =13                                         
               IEND(JATOM)  =2                                          
               ISTART(JATOM)=5                                          
               IEND1(JATOM) = 5                              
               ISTAR1(JATOM)=10                             
               ILOUT(JATOM) =10                                         
               text(jatom)='tot,0,1,PZ,PX+PY,2,DZ2,DXY,DX2Y2,DXZ+DYZ,3' 
            ELSE IF(JSPLIT.EQ.6) THEN                                   
               ILINP(JATOM) =10                                         
               IEND(JATOM)  =2                                          
               ISTART(JATOM)=6                                          
               IEND1(JATOM) = ILINP(JATOM)                              
               ISTAR1(JATOM)=ILINP(JATOM)+1                             
               ILOUT(JATOM) =7                                          
               text(jatom)='tot,0,1,PX,PY,PZ,2,3' 
            ELSE IF(JSPLIT.EQ.5) THEN                                   
               ILINP(JATOM) =12                                         
               IEND(JATOM)  =3                                          
               ISTART(JATOM)=9                                          
               IEND1(JATOM) = ILINP(JATOM)                              
               ISTAR1(JATOM)=ILINP(JATOM)+1                             
               ILOUT(JATOM) =9                                          
               text(jatom)='tot,0,1,2,DZ2,DX2Y2,DXY,DXZ,DYZ,3' 
            ELSE IF(JSPLIT.EQ.4) THEN                                   
               ILINP(JATOM) =12                                         
               IEND(JATOM)  =2                                          
               ISTART(JATOM)=5                                          
               IEND1(JATOM) =5                                          
               ISTAR1(JATOM)=9                                          
               ILOUT(JATOM) =9                                          
               text(jatom)='tot,0,1,PZ,PX+PY,2,DZ2,DX2Y2+DXY,DXZ+DYZ,3' 
            ELSE IF(JSPLIT.EQ.8) THEN                                   
               ILINP(JATOM) =15                                         
               IEND(JATOM)  =2                                          
               ISTART(JATOM)=6                                          
               IEND1(JATOM) =6                                          
               ISTAR1(JATOM)=12                                         
               ILOUT(JATOM) =12                                         
               text(jatom)='tot,0,1,PX,PY,PZ,2,DZ2,DX2Y2,DXY,DXZ,DYZ,3' 

            else IF(JSPLIT.EQ.15) THEN                                        
!   not installed yet !!
            ILINP(JATOM) =15                                             
            IEND(JATOM)  =2                                             
            ISTART(JATOM)=6                                             
            IEND1(JATOM) = 6                                 
            ISTAR1(JATOM)=14                                
            ILOUT(JATOM) =12                                             
         ENDIF                                                          
 19   CONTINUE                                                          
!                                                                       
!     ILINP :  DEFINES THE NUMBER OF QTL'S AND SPLIT TYPES WRITTEN      
!              IN L2MAIN TO TAPES 31,32,...,38                          
! =  7 :  L=0,1,2,3,4,5,6                                 (ISPLIT=0)    
! =  9 :  L=0,1,PZ,(PX,PY),2,3,4,5,6                      (ISPLIT=1,HCP)
! =  9    L=0,1,2,EG,T2G,3,4,5,6                          (ISPLIT=2,CUB)
! = 10 :  L=0,1,2,DZ2,(DXY,DX2Y2),(DXZ,DYZ),3,4,5,6       (ISPLIT=3,HCP)
! = 12 :  L=0,1,PZ,(PX,PY),2,DZ2,(DXY,DX2Y2),(DXZ,DYZ),3,4,5,6(ISPLIT=4)
! = 10 :  L=0,1,PX,PY,PZ,2,3,4,5,6                           (ISPLIT=6) 
! = 11 :  L=0,1,PZ,(PX,PY),2,DZ2,DXY,DX2Y2,(DXZ,DYZ),3,4,5,6 (ISPLIT=-2)
! = 12 :  L=0,1,2,DZ2,DX2Y2,DXY,DXZ,DYZ,3,4,5,6              (ISPLIT=5) 
! = 15 :  L=0,1,PX,PY,PZ,2,DZ2,DX2Y2,DXY,DXZ,DYZ,3,4,5,6     (ISPLIT=8) 
! = 22 :  L=0,1,PX,PY,PZ,2,DZ2,DX2Y2,DXY,DXZ,DYZ,3,4,5,6     (ISPLIT=15) 
!                                                                       
!     IEND  :  DEFINES THE LAST  UNSPLITTED QTL VALUE BEFORE SPLIT      
!              EG: 0,1,2,EG,T2G,3,4,5,6  ===> IEND=3   (L=2)            
!     ISTART:  DEFINES THE FIRST UNSPLITTED QTL VALUE AFTER  SPLIT      
!              EG: 0,1,2,EG,T2G,3,4,5,6  ===> ISTART=6 (L=3)            
!                                                                       
!     ILOUT :  DEFINES THE NUMBER OF QTL'S AND SPLIT TYPES WRITTEN      
!              TO TAPE QTL                                              
!       = 4 :  OUTPUT UNTIL F (ISPLIT=0)                                
!       = 6 :  OUTPUT UNTIL F (ISPLIT=1,2)                              
!       = 7 :  OUTPUT UNTIL F (ISPLIT=3,6)                              
!       =10 :  OUTPUT UNTIL F (ISPLIT=-2)                               
!       = 9 :  OUTPUT UNTIL F (ISPLIT=5)                                
!                                                                       
!    ===>  MAX 22 L,M COMPONENTS ON TAPE 31
!    ===>  MAX 12 L,M COMPONENTS ON OUTPUT                              
!                                                                       
      CLOSE(31)
      DO jatom=1,nat
         CALL open_helpfile(fnamehelp,jatom,fnamehelp2)
         OPEN(1000+jatom,file=fnamehelp2,form='formatted')
      ENDDO
      nband=10000
      nbandmax=0
      DO 10 JK=1,NK    
         if(nband.gt.nb(jk)) nband=nb(jk)                                     
         if(nbandmax.lt.nb(jk)) nbandmax=nb(jk)              
         DO 20 JB=1,NB(JK)                                              
            TC=0.                                                       
            DO 30 JATOM=1,NAT                                           
               TCATOM=0.
               IF(JB.EQ.1) READ(1000+jatom,1000) R,S,T,NV,NE,KNAME           
               READ(1000+jatom,1010) NUM,E,W                                 
!cccc               IF(JB.NE.NUM) GOTO 900
               DO 40 JL=1,ILINP(JATOM)                                  
                  READ(1000+jatom,1020) LL(JL),( QTL1(J,JL),J=1,3 )          
 40            CONTINUE                                                 
               DO 50 JL=1,IEND(JATOM)                                   
                  TCATOM=TCATOM + QTL1(1,JL)*MULT(JATOM)                
                  TC    =TC     + QTL1(1,JL)*MULT(JATOM)                
 50            CONTINUE                                                 
               DO 55 JL=ISTART(JATOM),IEND1(JATOM)                      
                  TCATOM=TCATOM + QTL1(1,JL)*MULT(JATOM)                
                  TC    =TC     + QTL1(1,JL)*MULT(JATOM)                
 55            CONTINUE                                                 
               DO 60 JL=ISTAR1(JATOM),ILINP(JATOM)                      
                  TCATOM=TCATOM + QTL1(1,JL)*MULT(JATOM)                
                  TC    =TC     + QTL1(1,JL)*MULT(JATOM)                
 60            CONTINUE                                                 
!                                                                       
               IF(JATOM.EQ.1 .AND. JB.EQ.1) WRITE(6,1001) R,S,T,         &
                                                          NV,NE,KNAME   
               IF(JATOM.EQ.1) WRITE(6,1011) NUM,E,W                     
               WRITE(6,1030) ( LL(JL),            JL=1,ILOUT(JATOM) )   
               WRITE(6,1035) (        QTL1(1,JL), JL=1,ILOUT(JATOM) )   
               WRITE(6,1040) (        QTL1(2,JL), JL=1,ILOUT(JATOM) )   
               WRITE(6,1041) (        QTL1(3,JL), JL=1,ILOUT(JATOM) )   
               IF(MODUS.EQ.AQTL) THEN                                   
                 IF(JATOM.EQ.1.AND.JB.EQ.1) WRITE(15,1001) R,S,T,        &
                                                           NV,NE,KNAME  
                 IF(JATOM.EQ.1)             WRITE(15,1011) NUM,E,W      
                 WRITE(15,1050) TCATOM,( QTL1(1,JL), JL=1,ILOUT(JATOM) )
!
                 if(isplit(jatom).eq.99) then
                   read(1000+jatom,207) ((xqtl(i,j,1),j=1,i),i=1,lxdos2)
                   WRITE(15,208) ((xqtl(i,j,1),j=1,i),i=1,lxdos2)
  207              FORMAT(7X,(10F10.5))                             
  208              FORMAT((10F10.5))                             
                 ENDIF                                                    
                 if(isplit(jatom).eq.88) then
                   read(1000+jatom,207) (xqtl(i,i,1),i=1,lxdos2)
                   WRITE(15,208) (xqtl(i,i,1),i=1,lxdos2)
                 ENDIF                                                    
               ENDIF                                                    
 30         CONTINUE                                                    
          TCOUT=100.D0-TC                                               
          WRITE(6,1060) TCOUT                                           
          IF(MODUS.EQ.AQTL) WRITE(15,1061) TCOUT                        
 20      CONTINUE                                                       
 10   CONTINUE                                                          
      write(6,*) '     NBAND in QTL-file:',nband
      write(21,*) '     NBAND in QTL-file:',nband


!!$!                                                                       
!!$!.....CHECK IF EOF IS REALLY REACHED                                    
      DO 100 JATOM=1,NAT
         itape=1000+jatom                                                
         READ(1000+jatom,1020,END=100) LDUMMY,QDUMY,QDUMY1,QDUMY2            
         WRITE (6,9000) ITAPE,LDUMMY,QDUMY,QDUMY1,QDUMY2                 
!         STOP
         GOTO 910
 100  CONTINUE                                                          
      IF(MODUS.NE.AQTL) RETURN                                          

! Read norm file in spinpol. so case
      if(kso.eq.2.and.nspin1.eq.2) then
      rewind 12
      allocate (sonorm(nbandmax,nk))
      DO JK=1,NK    
         read(12,'(4e20.11)')  (sonorm(jb,jk),  JB=1,NB(JK))
         write(6,*) 'jk',jk
         write(6,'(4e20.11)')  (sonorm(jb,jk),  JB=1,NB(JK))
      enddo                                              
      endif
!                                                                       
!.....WRITE QTL TO FILE 16 = QTL                                        
      NAT1=NAT + 1                                                      
      WRITE(16,2000) TITLE                                              
      WRITE(16,2010) AA,BB,CC,EFerm                                        
      WRITE(16,2020) MINWAV,MAXWAV,NSPIN1,NAT,kso                            
      DO 105 JATOM=1,NAT                                                
      WRITE(16,2030) JATOM,MULT(JATOM),isplit(JATOM),text(jatom)        
 105  CONTINUE                                                          
!                                                                       
      DO 110 JBAND=1,NBAND                                              
         WRITE(16,2040) JBAND                                           
         REWIND 15                                                      
         DO 120 JK=1,NK 
            IF(NB(NK).LT.NBAND) WRITE(*,*) ' FOR K-POINT ',NK,' ONLY ', &
            NB(NK),' INSTEAD OF ',NBAND,' BANDS ARE FOUND !!!!'              
            READ(15,1000) R,S,T,NV,NE,KNAME
            DO 125 JB=1,NB(JK)                                          
               READ(15,1010) NUM,E,W
!              Added scaling term for doos
               IF(JB.EQ.1)WZero=W   
               DO 130 JATOM=1,NAT                                       
                  READ(15,1050) TCATOM,( QTL1(1,JL), JL=1,ILOUT(JATOM) )
                  RScale = 1.D0/100.D0
                  if(DOOS)RScale=RScale*W/WZero    
!                  TCATOM=TCATOM/100.D0                                  
                  TCATOM=TCATOM*RScale
                  DO 135 JL=1,ILOUT(JATOM)                              
!                    QTL1(1,JL)=QTL1(1,JL)/100.D0                       
                     QTL1(1,JL)=QTL1(1,JL)*RScale   
 135              CONTINUE                                              
                  IF(jb.EQ.JBAND) WRITE(16,2050) E,JATOM,TCATOM,        &
                                   ( QTL1(1,JL), JL=1,ILOUT(JATOM) )    
!
                 if(isplit(jatom).eq.99) then
                   read(15,208) ((xqtl(i,j,1),j=1,i),i=1,lxdos2)
                   IF(jb.EQ.JBAND) WRITE(16,208)  &
                      ((xqtl(i,j,1)/100.d0,j=1,i),i=1,lxdos2)
                 ENDIF                                                    
                 if(isplit(jatom).eq.88) then
                   read(15,208) (xqtl(i,i,1),i=1,lxdos2)
                   IF(jb.EQ.JBAND) WRITE(16,208)  &
                      (dble(xqtl(i,i,1))/100.d0,i=1,lxdos2)
                 ENDIF                                                    
!
 130           CONTINUE                                                 
               READ(15,1061) TCOUT                                      
               TCOUT=TCOUT/100.D0                                       
! for SO set tcout to norm of this spin
      if(kso.eq.2.and.nspin1.eq.2) tcout=sonorm(jb,jk)
               IF(jb.EQ.JBAND) WRITE(16,2050) E,NAT1,TCOUT             
 125        CONTINUE                                                    
 120     CONTINUE                                                       
 110  CONTINUE                                                          
      deallocate (TEXT,ILINP,ILOUT,IEND,ISTART)
      deallocate (IEND1,ISTAR1)
      RETURN                                                            
!                                                                       
!        Error messages
!
  900 CALL OUTERR('OUTP','JB NE NUM.')
      STOP 'OUTP - Error'
  910 CALL OUTERR('OUTP','EOF NOT REACHED.')
      WRITE (ERRMSG,9000) ITAPE,LDUMMY,QDUMY,QDUMY1,QDUMY2                 
      CALL OUTERR('OUTP',ERRMSG)
      STOP 'OUTP - Error'
!
 9000 FORMAT('FOR TAPE ',I4,12X,A6,'QINSID=',F10.5,'Q(U)=',F10.5,     &
             'Q(UE)=',F10.5)                                            
 1000 FORMAT(/,10X,3F8.4,1X,I5,I4,2X,A10)                                 
 1001 FORMAT(/,1X,' K-POINT:',3F8.4,1X,I5,I4,2X,A10)                      
 1010 FORMAT(7X,I4,4X,F9.5,9X,F10.7)                                    
 1011 FORMAT(1X,' BAND#',I4,'  E=',F9.5,'  WEIGHT=',F10.7)             
 1020 FORMAT(1X,A6,F10.5,2X,2F10.5)                                     
 1030 FORMAT(1X,'       ',12(2X,A7) )                                   
 1035 FORMAT(1X,'QINSID:',12(2X,F7.4) )                                 
 1040 FORMAT(1X,'Q(U)  :',12(2X,F7.4) )                                 
 1041 FORMAT(1X,'Q(UE) :',12(2X,F7.4))                                  
 1050 FORMAT(1X,13F10.5)                                                
 1060 FORMAT(1X,'QOUT  :',F14.4)                                         
 1061 FORMAT(8X,F15.5)                                                  
 2000 FORMAT(A80,/)                                                       
 2010 FORMAT(1X,'LATTICE CONST.=',3F8.4,3X,'FERMI ENERGY=',F10.5)       
 2020 FORMAT(I5,' < NMAT <',I5,3X,'SPIN=',I1,3X,'NAT=',I3, &
             6x,'SO',i2)         
 2030 FORMAT(1X,'JATOM',I3,2X,'MULT=',I2,2x,'ISPLIT=',i2,2x,a45)              
 2040 FORMAT(1X,'BAND:',I4)                                          
 2050 FORMAT(F10.5,I3,F8.5,3X,12F8.5)                                
      END                                                               
