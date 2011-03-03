      SUBROUTINE SCFANA(NAT,JSPIN,FORSUM,ESUM,force,iscf,cut)
!
!     read scf file and sums up total energy and forces
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
!
      CHARACTER*79   MARGN                                             
!
!        Common blocks
!
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
!
!
!
      LOGICAL   FORCE(NATO),cut
      DIMENSION ESX(25,NATO),FORSUM(NATO,4),FHELP(4)       
!
!---------------------------------------------------------------------  
!                                                                       
!
! INITIALIZE FORCE AND TOTAL ENERGY PARAMETER
!
      ESUM=0.d0                                                           
      EDENS=0.d0                                                          
      DO 100 JATOM=1,NAT
      do 101 j=1,4
         FORSUM(JATOM,j)=0.d0 
         FORCE(JATOM)=.false.
  101 CONTINUE
      DO 110 I=1,25
            ESX(I,JATOM)=0.d0
  110    CONTINUE
  100 CONTINUE

!                                                                       
!.....READ TAPE22=SCFDATA UNTIL BOTTOM                                  
!                                                                       
      ISCF=0                                                            
  10  READ(22,700,END=20) MARGN
      IF(MARGN(2:4).EQ.'ITE')GO TO 30                                          
      WRITE(21,700)MARGN                                               
      GO TO 10                                                          
  30  READ(MARGN,710) ISCF 
      WRITE(21,700)MARGN
      GOTO 10                                                           
  20  CONTINUE                                                          
!
!..... READ RIGHT SCF SECTION
!
      REWIND 22                                                         
 60   READ(22,700,END=40) MARGN
      IF (MARGN(2:4).NE.'ITE') GO TO 60                                        
      READ(MARGN,710) ISCFT 
      IF (ISCFT.NE.ISCF) GO TO 60                                        
 50   READ(22,700,END=40) MARGN
!
!..... READ TOTAL ENERGY PARTS
!
      IF (MARGN(2:4).EQ.'WAR') cut=.true.
      IF (MARGN(2:4).EQ.'DEN') READ(MARGN,720)EDENS
      IF (MARGN(2:4).EQ.'SUM') THEN                                         
         READ(MARGN,720) ES1 
         ESUM=ESUM+ES1                                                     
      ENDIF                                                             
      IF (MARGN(2:4).EQ.'1S ') THEN                                         
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(1,JATOM)=ESX(1,JATOM)+ESXX                                 
      ELSE IF(MARGN(2:4).EQ.'2S ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(2,JATOM) =ESX(2,JATOM)+ESXX                                
      ELSE IF(MARGN(2:4).EQ.'2P ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(3,JATOM) =ESX(3,JATOM)+ESXX                                
      ELSE IF(MARGN(2:4).EQ.'2PP') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(4,JATOM) =ESX(4,JATOM)+ESXX                                
      ELSE IF(MARGN(2:4).EQ.'3S ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(5,JATOM) =ESX(5,JATOM)+ESXX                                
      ELSE IF(MARGN(2:4).EQ.'3P ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(6,JATOM) =ESX(6,JATOM)+ESXX                                
      ELSE IF(MARGN(2:4).EQ.'3PP') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(7,JATOM) =ESX(7,JATOM)+ESXX                                
      ELSE IF(MARGN(2:4).EQ.'3D ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(8,JATOM) =ESX(8,JATOM)+ESXX                                
      ELSE IF(MARGN(2:4).EQ.'3DD') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(9,JATOM) =ESX(9,JATOM)+ESXX                                
      ELSE IF(MARGN(2:4).EQ.'4S ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(10,JATOM) =ESX(10,JATOM)+ESXX                              
      ELSE IF(MARGN(2:4).EQ.'4P ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(11,JATOM) =ESX(11,JATOM)+ESXX                              
      ELSE IF(MARGN(2:4).EQ.'4PP') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(12,JATOM) =ESX(12,JATOM)+ESXX                              
      ELSE IF(MARGN(2:4).EQ.'4D ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(13,JATOM) =ESX(13,JATOM)+ESXX                              
      ELSE IF(MARGN(2:4).EQ.'4DD') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(14,JATOM) =ESX(14,JATOM)+ESXX                              
      ELSE IF(MARGN(2:4).EQ.'4F ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(15,JATOM) =ESX(15,JATOM)+ESXX                              
      ELSE IF(MARGN(2:4).EQ.'4FF') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(16,JATOM) =ESX(16,JATOM)+ESXX                              
      ELSE IF(MARGN(2:4).EQ.'5S ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(17,JATOM) =ESX(17,JATOM)+ESXX                              
      ELSE IF(MARGN(2:4).EQ.'5P ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(18,JATOM) =ESX(18,JATOM)+ESXX                              
      ELSE IF(MARGN(2:4).EQ.'5PP') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(19,JATOM) =ESX(19,JATOM)+ESXX                              
      ELSE IF(MARGN(2:4).EQ.'5D ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(20,JATOM) =ESX(20,JATOM)+ESXX                              
      ELSE IF(MARGN(2:4).EQ.'5DD') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(21,JATOM) =ESX(21,JATOM)+ESXX                              
      ELSE IF(MARGN(2:4).EQ.'6S ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(22,JATOM) =ESX(22,JATOM)+ESXX                              
      ELSE IF(MARGN(2:4).EQ.'6P ') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(23,JATOM) =ESX(23,JATOM)+ESXX                              
      ELSE IF(MARGN(2:4).EQ.'6PP') THEN                                    
         READ(MARGN,710) JATOM
         READ(MARGN,730) ESXX                                                
         ESX(24,JATOM) =ESX(24,JATOM)+ESXX                              
      END IF                                                            
!
!..... READ FORCE PARTS
!
      IF (MARGN(2:4).EQ.'FHF') THEN                                         
         READ(MARGN,710) JATOM
         READ(MARGN,740) (FHELP(j),j=1,4)
         FORCE(JATOM)=.true.
      do 200 j=2,4
         FORSUM(JATOM,j)=FORSUM(JATOM,j)+FHELP(j)
 200  CONTINUE
      END IF                                                            
      IF (MARGN(2:4).EQ.'FVA') THEN                                         
         READ(MARGN,710) JATOM
         READ(MARGN,740) (FHELP(j),j=1,4)
      do 210 j=2,4
         FORSUM(JATOM,j)=FORSUM(JATOM,j)+FHELP(j)                               
 210  CONTINUE
      END IF                                                            
      IF (MARGN(2:4).EQ.'FCO') THEN                                         
         READ(MARGN,710) JATOM
         READ(MARGN,740) (FHELP(j),j=1,4)
      do 220 j=2,4
         FORSUM(JATOM,j)=FORSUM(JATOM,j)+FHELP(j)/jspin
 220  CONTINUE
      END IF                                                            
      GO TO 50                                                          
 40   CONTINUE                                                          
! ..................................................................
!
! CALCULATE TOTAL ENERGY
!
      ESUM=ESUM+EDENS                                                   
      DO 11 JATOM=1,NAT                                                 
      FORSUM(JATOM,1)=sqrt &
            (forsum(jatom,2)**2 &
            +forsum(jatom,3)**2 &
            +forsum(jatom,4)**2)
      DO 9 J=1,25                                                       
 9    ESUM=ESUM+ESX(J,JATOM)*2.*MULT(JATOM)/JSPIN                       
      ESUM=ESUM+ESX(3,JATOM)*2.*MULT(JATOM)/JSPIN                       
      ESUM=ESUM+ESX(6,JATOM)*2.*MULT(JATOM)/JSPIN                       
      ESUM=ESUM+ESX(9,JATOM)*2.*MULT(JATOM)/JSPIN                       
      ESUM=ESUM+ESX(11,JATOM)*2.*MULT(JATOM)/JSPIN                      
      ESUM=ESUM+ESX(14,JATOM)*2.*MULT(JATOM)/JSPIN                      
      ESUM=ESUM+ESX(18,JATOM)*2.*MULT(JATOM)/JSPIN                      
      ESUM=ESUM+ESX(21,JATOM)*2.*MULT(JATOM)/JSPIN                      
      ESUM=ESUM+ESX(23,JATOM)*2.*MULT(JATOM)/JSPIN                      
      ESUM=ESUM+ESX(8,JATOM)*4.*MULT(JATOM)/JSPIN                       
      ESUM=ESUM+ESX(13,JATOM)*4.*MULT(JATOM)/JSPIN                      
      ESUM=ESUM+ESX(16,JATOM)*4.*MULT(JATOM)/JSPIN                      
      ESUM=ESUM+ESX(20,JATOM)*4.*MULT(JATOM)/JSPIN                      
 11   ESUM=ESUM+ESX(15,JATOM)*6.*MULT(JATOM)/JSPIN                      
!                                                                       
      RETURN
!
!
!
 700  FORMAT(A79)                                                    
 710  FORMAT(4X,i2)                                              
 720  FORMAT(30X,f20.6)                                              
 730  FORMAT(20X,f20.6)                                             
 740  FORMAT(15X,4f15.3)                                             
!
      END 

