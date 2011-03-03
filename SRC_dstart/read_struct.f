      subroutine read_struct
      use struct
	use reallocate
      RELA='RELA'                                     
!
      READ(20,1000) TITLE                                               
      READ(20,1010) LATTIC,NAT,CFORM,IREL                                     
!      IF (NAT .GT. NATO) stop 'nato too small'
!     READ IN LATTICE CONSTANTS                                         
      READ(20,1020) AA,BB,CC,ALPHA(1),ALPHA(2),ALPHA(3)
      IF(ABS(ALPHA(1)).LT.1.D-5) ALPHA(1)=90.0D0
      IF(ABS(ALPHA(2)).LT.1.D-5) ALPHA(2)=90.0D0
      IF(ABS(ALPHA(3)).LT.1.D-5) ALPHA(3)=90.0D0
      IF(IREL.EQ.'RELA') REL=.TRUE.                                     
      IF(IREL.EQ.'NREL') REL=.FALSE.                                    
      allocate ( rmt(nat), v(nat), iatnr(nat), mult(nat), isplit(nat) )
      allocate ( R0(nat),DX(nat),JRI(nat),ANAME(nat),zz(nat) )
      allocate ( ROTLOC(3,3,NAT) )
! allocate pos to max value (48*nat), reduce later
      allocate ( POS(3,48*nat) )
! 
!     INDEX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE    
!     NOTEQUIVALENT ATOMS                                               
      INDEX=0                                                           
      DO 50 JATOM = 1,NAT                                               
         INDEX=INDEX+1                                                  
!         IF (INDEX .GT. NDIF) stop 'ndif too small'
         READ(20,1030) IATNR(JATOM),( POS(J,INDEX),J=1,3 ),MULT(JATOM),  &
                       ISPLIT(JATOM)                                    
         IF (MULT(JATOM).EQ.0) THEN                                     
            WRITE(6,1041) JATOM,INDEX,MULT(JATOM)                       
            STOP ' DSTART: MULT EQ 0'                                    
         ENDIF                                                          
            DO 55 M=1,MULT(JATOM)-1                                     
               INDEX=INDEX+1                                            
!         IF (INDEX .GT. NDIF) stop 'ndif too small'
               READ(20,1031) IATNR(JATOM),( POS(J,INDEX),J=1,3)         
 55         CONTINUE                                                    
            READ(20,1050) ANAME(JATOM),JRI(JATOM),R0(JATOM),RMT(JATOM) &
      	               ,zz(jatom) 	 
         DX(JATOM)=LOG(RMT(JATOM)/R0(JATOM)) / (JRI(JATOM)-1)           
         RMT(JATOM)=R0(JATOM)*EXP( DX(JATOM)*(JRI(JATOM)-1) )           
         READ(20,1051) ((ROTLOC(I,J,JATOM),I=1,3),J=1,3)                
 50   CONTINUE                                                          
         ndif=index
! reduce pos to ndif
	call doreallocate(pos, 3, ndif)

      WRITE(6,800)                                                      
      WRITE(6,805)  TITLE                                               
      WRITE(6,810)  LATTIC                                              
      WRITE(6,820)  AA,BB,CC                                            
      WRITE(6,840)  NAT                                                 
      WRITE(6,850)  IREL                                                

      READ(20,10) IORD                                                  
      DO 5 J=1,IORD                                                     
 5    READ(20,11) ( (Iz(J1,J2,J),J1=1,3),TAU(J2,J), J2=1,3 )          
!.....READ IN SYMMETRY OPERATIONS AND NONPRIMITIVE TRANSLATIONS         
 10   FORMAT(I4)                                                        
 11   FORMAT(3(3I2,F11.8/))                                              
!---------------------------------------------------------------------  


 1000 FORMAT(A80)                                                       
 1010 FORMAT(A4,24X,I2,1X,A4,/,13X,A4,18X,A4)                                 
 1020 FORMAT(6F10.7,10X,F10.7)                                          
 1030 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
 1031 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1041 FORMAT(///,3X,'ERROR IN DSTART : MULT(JATOM)=0 ...',/, 20X &
        ,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)           
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
 1051 FORMAT(20X,3F10.8)                                             
 800  FORMAT(////,30X,50(1H-),/,33X,'S T R U C T U R A L   ',            &
             'I N F O R M A T I O N',/,30X,50(1H-),//)                  
 805  FORMAT(3X,'SUBSTANCE',20X,'= ',A80,/)                             
 810  FORMAT(3X,'LATTICE',22X,'= ',A4)                                  
 820  FORMAT(3X,'LATTICE CONSTANTS ARE',8X,'= ',3F12.7)                 
 840  FORMAT(3X,'NUMBER OF ATOMS IN UNITCELL  = ',I3)                   
 850  FORMAT(3X,'MODE OF CALCULATION IS',7X,'= ',A4)                    

      end
