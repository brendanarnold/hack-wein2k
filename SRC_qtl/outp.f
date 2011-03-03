      SUBROUTINE OUTP(ICASE)
      USE param         
      USE struct                                    
      USE case
      USE com
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*4      MODUS                                           
      CHARACTER*7      LL,LDUMMY                                        
      CHARACTER*10     KNAME                                            
      CHARACTER*45     TEXT                                            
      CHARACTER*67     ERRMSG                                            
!      
      COMMON /CHAR/   MODUS
      DIMENSION       QTL1(NDIM2),LL(NDIM2)     
!--------------------------------------------------------------------
      mw=0
	write(6,*)'NUMBER OF K-POINTS:',nkpt
	do i=1,nkpt
	write(6,*)i,nb(i),'bands'
	end do
	text='projected density of states'   
	 ITAPE=29
	 IF (MODUS.EQ.'SUMA') THEN
	 IL=ISUM(0,ICASE)
	 ELSEIF (MODUS.EQ.'FULL') THEN
	 IL= (2*LCASE(ICASE)+1)*iso    
         ELSEIF (MODUS.EQ.'TOTA') THEN
	 IL=1
	 ELSEIF (MODUS.EQ.'SPIN') THEN
         IL=2
         ELSE
         STOP 'MODUS UNKNOWN'                                      
	 END IF
      dum=0.000
      rewind(ITAPE) 
      rewind(15)
      nband=1000
      DO 11 MK=1,nkpt  
         DO 20 JB=1,NB(MK)
         if(nband.gt.nb(mk)) nband=nb(mk)
               IF(JB.EQ.1) THEN
               READ(ITAPE,1000) R,S,T,iNV,iNE,KNAME   
               END IF       
               READ(ITAPE,1010) NUM,E,W
               DO JL=1,IL 
               READ(ITAPE,1020)LL(JL),QTL1(JL)     
               END DO  
            IF(JB.EQ.1) WRITE(15,1001) R,S,T,iNV,iNE,KNAME
               WRITE(15,1011) NUM,E,W      
               WRITE(15,1050)(QTL1(JL), JL=1,IL)
 20      CONTINUE    
 11     CONTINUE       
!                                                                       
!.....CHECK IF EOF IS REALLY REACHED                                    
         READ(ITAPE,1020,END=100) LDUMMY,QDUMY,QDUMY1,QDUMY2            
         GOTO 910
 100  CONTINUE                                                          
!.....WRITE QTL TO FILE 16 = QTL                                        
      WRITE(30+ICASE,2000) TITLE                                              
      WRITE(30+ICASE,2010) AA,BB,CC,ef                                        
      WRITE(30+ICASE,2020) MW,MW,2,1                           
      WRITE(30+ICASE,2030)1,MULT(1),isplit(1),text        
      DO 110 JBAND=1,NBAND                                              
         WRITE(30+ICASE,2040) JBAND                                           
         REWIND 15                                                      
         DO 120 MK=1,nkpt 
            IF(NB(nkpt).LT.NBAND) WRITE(*,*) ' FOR K-POINT ',nkpt,' ONLY ', &
            NB(nkpt),' INSTEAD OF ',NBAND,' BANDS ARE FOUND !!!!'                                                 
            READ(15,1000) R,S,T,NV,NE,KNAME
            DO 125 JB=1,NB(MK)                                          
               READ(15,1010) NUM,E,W   
                  READ(15,1050)(QTL1(JL),JL=1,IL)
                  IF(jb.EQ.JBAND) WRITE(30+ICASE,2050) E,1,        &
                               1.0,(QTL1(JL),JL=1,IL)    
               IF(jb.EQ.JBAND) WRITE(30+ICASE,2050) E             
 125        CONTINUE                                                    
 120     CONTINUE                                                       
 110  CONTINUE                                                          
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
 9000 FORMAT('FOR TAPE ',I2,12X,A6,'QINSID=',F10.5,'Q(U)=',F10.5,     &
             'Q(UE)=',F10.5)                                            
 1000 FORMAT(/,10X,3F8.4,1X,I5,I4,2X,A10)                                 
 1001 FORMAT(/,1X,' K-POINT:',3F8.4,1X,I5,I4,2X,A10)                      
 1010 FORMAT(8X,I3,4X,F9.5,9X,F10.7)                                    
 1011 FORMAT(1X,' BAND #',I3,'  E=',F9.5,'  WEIGHT=',F10.7)             
 1020 FORMAT(1X,A6,e17.8)                                     
 1030 FORMAT(1X,'       ',12(2X,A7) )                                   
 1035 FORMAT(1X,'QINSID:',12(2X,F7.4) )                                 
 1040 FORMAT(1X,'Q(U)  :',12(2X,F7.4) )                                 
 1041 FORMAT(1X,'Q(UE) :',12(2X,F7.4))                                  
 1050 FORMAT(1X,12e17.8)                                                
 1060 FORMAT(1X,'QOUT  :',F7.4)                                         
 1061 FORMAT(8X,F10.5)                                                  
 2000 FORMAT(A80,/)                                                       
 2010 FORMAT(1X,'LATTICE CONST.=',3F8.5,3X,'FERMI ENERGY=',F10.5)       
 2020 FORMAT(1X,I4,' < NMAT < ',I4,3X,'SPIN=',I1,3X,'NATO',I3)         
 2030 FORMAT(1X,'JATOM',I3,2X,'MULT=',I2,2x,'ISPLIT=',i2,2x,a45)              
 2040 FORMAT(1X,'BAND:',1X,I3)                                          
 2050 FORMAT(F10.5,1X,I2,F8.5,3X,14f8.5)                                
      END                                                               
