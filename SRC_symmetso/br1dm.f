      SUBROUTINE BR1DM(BR1,lattic,alat,alpha)   
!                                                                       
!     BR1(3,3)  : TRANSFORMS INTEGER RECIPROCAL LATTICE VECTORS AS      
!                 GIVEN IN THE VECTORLIST OF LAPW1  ( GENERATED IN      
!                 COORS, TRANSFORMED IN BASISO, AND WRITTEN OUT IN      
!                 WFTAPE) INTO CARTESIAN SYSTEM                         
!                 CIAL COORDINATE SYSTEM ( IN UNITS OF 2 PI / A )       
!                 TO CARTESIAN SYSTEM                                   
!     Taken from LAPWDM
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      character*4 lattic
      dimension alat(3),alpha(3)
      dimension br1(3,3),pia(3)
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
!
!        Error: wrong lattice, stop execution
!
         GOTO 900
!                                                                       
!.....HEXAGONAL LATTICE                                                 
 10   CONTINUE                                                          
      BR1(1,1)=2.D0/SQRT3*PIA(1)                                        
      BR1(1,2)=1.D0/SQRT3*PIA(1)                                        
      BR1(1,3)=0.0D0                                                    
      BR1(2,1)=0.0D0                                                    
      BR1(2,2)=PIA(2)                                                   
      BR1(2,3)=0.0D0                                                    
      BR1(3,1)=0.0D0                                                    
      BR1(3,2)=0.0D0                                                    
      BR1(3,3)=PIA(3)                                                   
!                                                                       
      GOTO 100                                                          
!.....RHOMBOHEDRAL CASE                                                    
 60   BR1(1,1)=1.D0/SQRT(3.D0)*PIA(1)                                          
      BR1(1,2)=1.D0/SQRT(3.D0)*PIA(1)                                          
      BR1(1,3)=-2.d0/sqrt(3.d0)*PIA(1)                                         
      BR1(2,1)=-1.0d0*PIA(2)                                                  
      BR1(2,2)=1.0d0*PIA(2)                                                    
      BR1(2,3)=0.0d0*PIA(2)                                                    
      BR1(3,1)=1.0d0*PIA(3)                                                    
      BR1(3,2)=1.0d0*PIA(3)                                                    
      BR1(3,3)=1.0d0*PIA(3)                                                    
!
      GOTO 100                                                          
!                                                                       
!.....PRIMITIVE LATTICE                                                 
!                                                                       
  20  CONTINUE
      SINBC=SIN(ALPHA(1))
      COSAB=COS(ALPHA(3))
      COSAC=COS(ALPHA(2))
      COSBC=COS(ALPHA(1))
      WURZEL=SQRT(SINBC**2-COSAC**2-COSAB**2+2*COSBC*COSAC*COSAB)
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
      GOTO 100                                                          
!                                                                       
!.....CXZ CASE (CXZ LATTICE BUILD UP)                                     
 51   CONTINUE                                     
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
         GOTO 100                                                          
      ELSE
!.....CXZ MONOCLINIC CASE 
         write(*,*) '  gamma not equal 90'
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
      GOTO 100                                                          
 100  CONTINUE                                                          
      do i=1,3
!     write(6,654)(br1(i,j),j=1,3)
654   format(3f10.5,3x,3f10.5)
      enddo
!
      RETURN
!
!        Error messages
!
  900 write(6,*)'wrong lattice.'
!
      END                                                               
