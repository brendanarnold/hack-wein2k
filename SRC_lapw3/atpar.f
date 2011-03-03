      SUBROUTINE ATPAR(NAT)                                             
!-----X----X----X----X----X----X----X----X----X----X----X----X----X----X
!-----X LINEAR EXPANDED ORBITAL APW ROUTINE--------               X----X
!-----X SET UP THE RADIAL INTEGRAL PARAMETERS.                    X----X
!-----X SPIN-ORBITLESS RELATIVISTIC EQUATIONS USED.               X----X
!-----X L IS 7  (NT)                                              X----X
!-----X REQUIRED SUBROUTINES ARE OUTWIN AND RINT13.               X----X
!-----X SET UP FOR 4 ATOMS, AND CALL 4 TIMES                      X----X
!-----X                                      P.BLAHA              X----X
!-----X----X----X----X----X----X----X----X----X----X----X----X----X----X
      use atomgrid
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      CHARACTER*10 ANAME                                                
      COMMON /STRUK/ AA,BB,CC,alpha(3),PIA(3),VOL,ndif
      COMMON /POTNLC/ VR(NRAD)
!      COMMON /CUBWK1/ A(885),B(885),AP(885),BP(885),AE(885),BE(885)     
    1 FORMAT(A10,I5,5X,2E10.5,I5,F5.2)                                  
    2 FORMAT(5E14.7)                                                    
   12 FORMAT(8I2)                                                       
   13 FORMAT(/,' AT.NR.:',I2,5X,A10,5X,'POSITION:.',3F8.3,10X,         &
      'MULTIPLICITY: ',I5)                                              
!                                                                       
!     INDEX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE    
!     NOTEQUIVALENT ATOMS                                               
      INDEX=0                                                           
      DO 50 JATOM = 1,NAT                                               
         INDEX=INDEX+1                                                  
         READ(20,1030) IATNR(JATOM),( POS(INDEX,j),J=1,3 ),MULT(JATOM),  &
                       ISPLIT(JATOM)                                    
         IF (MULT(JATOM).EQ.0) THEN                                     
            WRITE(66,1040) JATOM,INDEX,MULT(JATOM)                      
            STOP ' LAPW2: MULT EQ 0'                                    
         ENDIF                                                          
            DO 55 M=1,MULT(JATOM)-1                                     
               INDEX=INDEX+1                                            
               READ(20,1031) IATNR(JATOM),( POS(INDEX,j),J=1,3)         
 55         CONTINUE                                                    
         READ(20,1050) ANAME,JRI(JATOM),Rnot(JATOM),RMT(JATOM)     
         DX(JATOM)=LOG(RMT(JATOM)/Rnot(JATOM)) / (JRI(JATOM)-1)           
         READ(20,1051) ((ROTLOC(I,J,JATOM),I=1,3),J=1,3)                
 1030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
 1031 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1040 FORMAT(///,3X,'ERROR IN LAPW2 : MULT(JATOM)=0 ...',                &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)           
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
 1051    FORMAT(20X,3F10.8)                                             
 50   CONTINUE                                                          
      ndif=index
!                                                                       
   51 CONTINUE                                                          
!                                                                       
      RETURN                                                            
      END                                                               
