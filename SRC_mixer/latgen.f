      SUBROUTINE LATGEN(NAT)                                            
!                                                                       
!.....SET UP THE LATTICE GENERATION                                     
!     BR1(3,3) IS A BRAVAIS MATRIX, WHICH TRANSFORMS                    
!     RECIPROCAL LATTICE VECTORS (FROM LAPW1) IN  KARTESIAN             
!     COORDINATES.                                                      
!                                                                       
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
!
      COMMON /GENER/   BR1(3,3),BR2(3,3),avec(3,3)                            
!        Common blocks
!
!
!        LATTIC  - lattice type
!                  'P' ...... primitive latice (cubic, tetragonal,
!                             orthorhombic, monoclinic, triclin)
!                  'FC' ..... face centered
!                  'BC' ..... body centered
!                  'HEX' .... hexagonal
!                  'CXY' .... c-base centered (only orthorombic)
!                  'CYZ' .... a-base centered (only orthorombic)
!                  'CXZ' .... b-base centered (only orthorombic)
!                  'R' ...... rhomboedral
!
      CHARACTER*4        LATTIC
      COMMON  /CHAR/     LATTIC
      SAVE    /CHAR/
!
!        ORTH - .TRUE. for orthogonal lattice
!               .FALSE. otherwise
!
      LOGICAL            ORTHO
      COMMON  /ORTH/     ORTHO
      SAVE    /ORTH/
!
!        ALAT(1:3)   - lattice constants (a,b,c)
!        ALPHA(1:3)  - angles of the unit-cell's unit-vectors
!                      also indicates cubic and non-cubic symmetries
!        PIA(1:3)    - reciprocal lattice constants (2*pi/a, 2*pi/b,
!                      2*pi/c)
!        VI          - volume of the direct lattice unit-cell
!
      DOUBLE PRECISION   VI 
      DOUBLE PRECISION   ALAT(3), ALPHA(3), PIA(3)
      COMMON  /STRUK/    ALAT, ALPHA, PIA, VI
      SAVE    /STRUK/
!
!---------------------------------------------------------------------  
!                      
      PI=ACOS(-1.0D0)                                                   
      SQRT3=SQRT(3.D0)
      ALPHA(1)=ALPHA(1)*PI/180.0D0                                             
      ALPHA(2)=ALPHA(2)*PI/180.0D0                                           
      ALPHA(3)=ALPHA(3)*PI/180.0D0                                            
      PIA(1)=2.D0*PI/ALAT(1)                                                 
      PIA(2)=2.D0*PI/ALAT(2)                                                 
      PIA(3)=2.D0*PI/ALAT(3)                                                 
      IF(LATTIC(1:1).EQ.'H') GOTO 10                                    
      IF(LATTIC(1:1).EQ.'S') GOTO 20                                    
      IF(LATTIC(1:1).EQ.'P') GOTO 20                                    
      IF(LATTIC(1:1).EQ.'F') GOTO 30                                    
      IF(LATTIC(1:1).EQ.'B') GOTO 40                                    
      IF(LATTIC(1:1).EQ.'C') GOTO 50                                    
      IF(LATTIC(1:1).EQ.'R') GOTO 60                                    
      GOTO 900
!                                                                       
!.....HEXAGONAL LATTICE                                                 
 10   CONTINUE                                                          
      BR1(1,1)=2.D0/SQRT3*PIA(1)                                        
      BR1(1,2)=1.D0/SQRT3*PIA(2)                                        
      BR1(1,3)=0.0D0                                                    
      BR1(2,1)=0.0D0                                                    
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=0.0D0                                                    
      BR1(3,1)=0.0D0                                                    
      BR1(3,2)=0.0D0                                                    
      BR1(3,3)=PIA(3)                                                   
!                                                                       
      BR2(1,1)=2.D0/SQRT3*PIA(1)                                        
      BR2(1,2)=1.D0/SQRT3*PIA(2)                                        
      BR2(1,3)=0.0D0                                                    
      BR2(2,1)=0.0D0                                                    
      BR2(2,2)=PIA(2)                                                   
      BR2(2,3)=0.0D0                                                    
      BR2(3,1)=0.0D0                                                    
      BR2(3,2)=0.0D0                                                    
      BR2(3,3)=PIA(3)                                                   
!                                                                       
      RVFAC=2.D0/SQRT(3.D0)                                             
      ORTHO=.FALSE.                                             
      GOTO 100                                                          
!                                                                       
!.....RHOMBOHEDRAL CASE                                                
 60   BR1(1,1)=1.D0/SQRT(3.D0)*PIA(1)                                          
      BR1(1,2)=1.D0/SQRT(3.D0)*PIA(1)                                          
      BR1(1,3)=-2.d0/sqrt(3.d0)*PIA(1)                                         
      BR1(2,1)=-1.0*PIA(2)                                                    
      BR1(2,2)=1.0*PIA(2)                                                      
      BR1(2,3)=0.0*PIA(2)                                                      
      BR1(3,1)=1.0*PIA(3)                                                      
      BR1(3,2)=1.0*PIA(3)                                                      
      BR1(3,3)=1.0*PIA(3)                                                      
!
      BR2(1,1)=1.D0/SQRT(3.D0)*PIA(1)                                          
      BR2(1,2)=1.D0/SQRT(3.D0)*PIA(1)                                          
      BR2(1,3)=-2.d0/sqrt(3.d0)*PIA(1)                                         
      BR2(2,1)=-1.0*PIA(2)                                                     
      BR2(2,2)=1.0*PIA(2)                                                      
      BR2(2,3)=0.0*PIA(2)                                                      
      BR2(3,1)=1.0*PIA(3)                                                      
      BR2(3,2)=1.0*PIA(3)                                                      
      BR2(3,3)=1.0*PIA(3)                                                      
      RVFAC=6.D0/SQRT(3.D0)
      ORTHO=.FALSE.                                             
      GOTO 100                                                          
!                                                                       
!.....FC LATTICE                                                        
 30   CONTINUE                                                          
      BR1(1,1)=PIA(1)                                                   
      BR1(1,2)=0.0D0                                                    
      BR1(1,3)=0.0D0                                                    
      BR1(2,1)=0.0D0                                                    
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=0.0D0                                                    
      BR1(3,2)=0.0D0                                                    
      BR1(3,1)=0.0D0                                                    
      BR1(3,3)=PIA(3)                                                   
!                    
!
      BR2(1,1)=-PIA(1)                                                  
      BR2(1,2)= PIA(1)                                                  
      BR2(1,3)= PIA(1)                                                  
      BR2(2,1)= PIA(2)                                                  
      BR2(2,2)=-PIA(2)                                                  
      BR2(2,3)= PIA(2)                                                  
      BR2(3,1)= PIA(3)                                                  
      BR2(3,2)= PIA(3)                                                  
      BR2(3,3)=-PIA(3)                                                  
!                                                                       
      RVFAC=4.D0                                                        
      ORTHO=.TRUE.                                             
      GOTO 100                                                          
!                                                                       
!.....BC LATTICE                                                        
 40   CONTINUE                                                          
      BR1(1,1)=PIA(1)                                                   
      BR1(1,2)=0.0D0                                                    
      BR1(1,3)=0.0D0                                                    
      BR1(2,1)=0.0D0                                                    
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=0.0D0                                                    
      BR1(3,1)=0.0D0                                                    
      BR1(3,2)=0.0D0                                                    
      BR1(3,3)=PIA(3)                                                   
!                                                                       
      BR2(1,1)= 0.0D0                                                    
      BR2(1,2)= PIA(1)                                                  
      BR2(1,3)= PIA(1)                                                  
      BR2(2,1)= PIA(2)                                                  
      BR2(2,2)= 0.0D0                                                    
      BR2(2,3)= PIA(2)                                                  
      BR2(3,1)= PIA(3)                                                  
      BR2(3,2)= PIA(3)                                                  
      BR2(3,3)= 0.0D0
!
      RVFAC=2.D0
      ORTHO=.TRUE.                                             
      GOTO 100                                                    
!             
 50   CONTINUE                                                          
      IF(LATTIC(2:3).EQ.'XZ') GOTO 51                                    
      IF(LATTIC(2:3).EQ.'YZ') GOTO 52                                    
!.....CXY LATTICE                                                          
      BR1(1,1)=PIA(1)                                                   
      BR1(1,2)=0.0D0                                                    
      BR1(1,3)=0.0D0                                                    
      BR1(2,1)=0.0D0                                                    
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=0.0D0                                                    
      BR1(3,1)=0.0D0                                                    
      BR1(3,2)=0.0D0                                                    
      BR1(3,3)=PIA(3)                                                   
!                                                                       
      BR2(1,1)= PIA(1)                                                    
      BR2(1,2)= PIA(1)                                                  
      BR2(1,3)= 0.0D0                                                  
      BR2(2,1)=-PIA(2)                                                  
      BR2(2,2)= PIA(2)                                                    
      BR2(2,3)= 0.0D0                                                 
      BR2(3,1)= 0.0D0                                                  
      BR2(3,2)= 0.0D0                                                 
      BR2(3,3)= PIA(3)                                                    
!                                                                       
      RVFAC=2.D0                                                        
      ORTHO=.TRUE.                                             
      GOTO 100                                                          
!                                                                       
!.....CXZ CASE (CXZ LATTICE BUILD UP)                                     
 51   CONTINUE                                     
!.....CXZ ORTHOROMBIC CASE
      IF(ABS(ALPHA(3)-PI/2.0D0).LT.0.0001) then
         BR1(1,1)=PIA(1)                                                   
         BR1(1,2)=0.0D0                                                    
         BR1(1,3)=0.0D0                                                    
         BR1(2,1)=0.0D0                                                    
         BR1(2,2)=PIA(2)                                                   
         BR1(2,3)=0.0D0                                                    
         BR1(3,1)=0.0D0                                                    
         BR1(3,2)=0.0D0                                                    
         BR1(3,3)=PIA(3)                                                   
!                                                                       
         BR2(1,1)= PIA(1)                                                   
         BR2(1,2)= 0.0                                                   
         BR2(1,3)= PIA(1)                                                      
         BR2(2,1)= 0.0                                                      
         BR2(2,2)= PIA(2)                                                     
         BR2(2,3)= 0.0                                                     
         BR2(3,1)=-PIA(3)                                                     
         BR2(3,2)= 0.0                                                     
         BR2(3,3)= PIA(3)                                                     
!                                                                       
         RVFAC=2.0                                                         
         ORTHO=.TRUE.                                             
         GOTO 100                                                          
      ELSE
!.....CXZ MONOCLINIC CASE
         write(*,*) 'gamma not equal 90'
         SINAB=SIN(ALPHA(3))
         COSAB=COS(ALPHA(3))
!
         BR1(1,1)= PIA(1)/SINAB
         BR1(1,2)= -PIA(2)*COSAB/SINAB
         BR1(1,3)= 0.0
         BR1(2,1)= 0.0
         BR1(2,2)= PIA(2)
         BR1(2,3)= 0.0
         BR1(3,1)= 0.0
         BR1(3,2)= 0.0
         BR1(3,3)= PIA(3)
!
         BR2(1,1)= PIA(1)/SINAB
         BR2(1,2)= -PIA(2)*COSAB/SINAB
         BR2(1,3)= PIA(1)/SINAB
         BR2(2,1)= 0.0
         BR2(2,2)= PIA(2)
         BR2(2,3)= 0.0
         BR2(3,1)=-PIA(3)
         BR2(3,2)= 0.0
         BR2(3,3)= PIA(3)
!
         RVFAC=2.0/SINAB
         ORTHO=.FALSE.
         GOTO 100
      ENDIF
!                                                                       
!.....CYZ CASE (CYZ LATTICE BUILD UP)                                     
 52   CONTINUE                                     
      BR1(1,1)=PIA(1)                                                   
      BR1(1,2)=0.0D0                                                    
      BR1(1,3)=0.0D0                                                    
      BR1(2,1)=0.0D0                                                    
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=0.0D0                                                    
      BR1(3,1)=0.0D0                                                    
      BR1(3,2)=0.0D0                                                    
      BR1(3,3)=PIA(3)                                                   
!                                                                       
      BR2(1,1)= PIA(1)                                                      
      BR2(1,2)= 0.0                                                   
      BR2(1,3)= 0.0                                                      
      BR2(2,1)= 0.0                                                      
      BR2(2,2)= PIA(2)                                                     
      BR2(2,3)= PIA(2)                                                     
      BR2(3,1)= 0.0                                                     
      BR2(3,2)=-PIA(3)                                                     
      BR2(3,3)= PIA(3)                                                     
!                                                                       
      RVFAC=2.0                                                         
      ORTHO=.TRUE.                                             
      GOTO 100                                                          
!                                                                       
!
!     PRIMITIV CASE
!
  20  CONTINUE
      SINBC=SIN(ALPHA(1))
      COSAB=COS(ALPHA(3))
      COSAC=COS(ALPHA(2))
      COSBC=COS(ALPHA(1))
      WURZEL=SQRT(SINBC**2-COSAC**2-COSAB**2+2*COSBC*COSAC*COSAB)
      BR2(1,1)= SINBC/WURZEL*PIA(1)
      BR2(1,2)= (-COSAB+COSBC*COSAC)/(SINBC*WURZEL)*PIA(2)
      BR2(1,3)= (COSBC*COSAB-COSAC)/(SINBC*WURZEL)*PIA(3)
      BR2(2,1)= 0.0
      BR2(2,2)= PIA(2)/SINBC
      BR2(2,3)= -PIA(3)*COSBC/SINBC
      BR2(3,1)= 0.0
      BR2(3,2)= 0.0
      BR2(3,3)= PIA(3)
!
      BR1(1,1)= SINBC/WURZEL*PIA(1)
      BR1(1,2)= (-COSAB+COSBC*COSAC)/(SINBC*WURZEL)*PIA(2)
      BR1(1,3)= (COSBC*COSAB-COSAC)/(SINBC*WURZEL)*PIA(3)
      BR1(2,1)= 0.0
      BR1(2,2)= PIA(2)/SINBC
      BR1(2,3)= -PIA(3)*COSBC/SINBC
      BR1(3,1)= 0.0
      BR1(3,2)= 0.0
      BR1(3,3)= PIA(3)
!
      RVFAC= 1.d0/WURZEL
      ORTHO=.TRUE.
      if(abs(alpha(1)-pi/2.d0).gt.0.0001) ortho=.false.
      if(abs(alpha(2)-pi/2.d0).gt.0.0001) ortho=.false.
      if(abs(alpha(3)-pi/2.d0).gt.0.0001) ortho=.false.
!
      GOTO 100
!                                                                       
!                                                                       
!.....DEFINE VOLUME OF UNIT CELL                                        
 100  CONTINUE                                                          
!c     call gbass(br2,avec)
      VI=ALAT(1)*ALAT(2)*ALAT(3)/RVFAC  
!                                                                       
!.....DEFINE ROTATION MATRICES IN NONSYMMORPHIC CASE                    
!c      CALL ROTDEF (NAT)                                                 
!
      RETURN
!
!        Error messages
!
  900 CALL OUTERR('LATGEN','wrong lattice.')
      STOP 'LATGEN - Error'
!
!        End of 'LATGEN'
!
      END    

