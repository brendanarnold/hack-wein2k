!                                                                       
      SUBROUTINE DIRLAT (nat,ortho,alpha,beta,gamma)                          
!                                                                       
!     LATGEN GENERATES THE BRAVAIS MATRIX BR2(3,3), WHICH TRANSFORMS    
!     A RECIPROCAL LATTICE VECTOR OF A SPECIAL COORDINATE SYSTEM (IN    
!     UNITS OF 2 PI/A(I)) INTO CARTESIAN SYSTEM                         
!     Convention:  R_i = (i,*)
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
                                                                        
      LOGICAL           ORTHO
      CHARACTER*4       LATTIC                                          
!                                                                       
      COMMON /CHAR   /  LATTIC                                          
      COMMON /GENER  /  BR2(3,3)                                        
!-----------------------------------------------------------------------
!
      pi=4.d0*atan(1.d0)
      IF(LATTIC(1:1).EQ.'H') GOTO 10                                    
      IF(LATTIC(1:1).EQ.'S') GOTO 20                                    
      IF(LATTIC(1:1).EQ.'P') GOTO 20                                    
      IF(LATTIC(1:1).EQ.'B') GOTO 30                                    
      IF(LATTIC(1:1).EQ.'F') GOTO 40                                    
      IF(LATTIC(1:1).EQ.'C') GOTO 50                                    
      IF(LATTIC(1:1).EQ.'R') GOTO 80                                    
      STOP 'LATTIC WRONG'                                               
!                                                                       
!.....HEXAGONAL CASE                                                    
 10   BR2(1,1)=SQRT(3.d0)/2.d0                                              
      BR2(1,2)=-.5d0                                                      
      BR2(1,3)=0.0d0                                                      
      BR2(2,1)=0.0d0                                                     
      BR2(2,2)=1.0d0                                                      
      BR2(2,3)=0.0d0                                                      
      BR2(3,1)=0.0d0                                                      
      BR2(3,2)=0.0d0                                                      
      BR2(3,3)=1.d0                                                       
      ORTHO=.FALSE.                                             
      GOTO 100                                                          
!                                                                       
!.....PRIMITIVE LATTICE CASE 
 20   if(beta.ne.90.d0) then
! 	write(*,*) 'beta ne 90'
!	goto 21
      endif    
      if(alpha.ne.90.d0) then
!       write(*,*) 'alpha ne 90'
!        goto 22
      endif    
      BR2(1,1)=1.0d0*sin(gamma/180.d0*pi)*sin(beta/180.d0*pi)                
      BR2(1,2)=1.0d0*cos(gamma/180.d0*pi)*sin(beta/180.d0*pi)                 
      BR2(1,3)=1.0d0*cos(beta/ 180.d0*pi)                                   
      BR2(2,1)=0.0d0                                                      
      BR2(2,2)=1.0d0*sin(alpha/180.d0*pi)
      BR2(2,3)=1.0d0*cos(alpha/180.d0*pi)
      BR2(3,1)=0.0d0                                                      
      BR2(3,2)=0.0d0                                                      
      BR2(3,3)=1.0d0                                                      
      ORTHO=.TRUE. 
      if(gamma.ne.90.d0) ortho=.false.                                      
      if(beta.ne.90.d0) ortho=.false.                                      
      if(alpha.ne.90.d0) ortho=.false.                                      
      GOTO 100     
 21   BR2(1,1)=1.0d0*sin(beta/180.d0*pi)
      BR2(1,2)=0.0d0
      BR2(1,3)=1.0d0*cos(beta/180.d0*pi)
      BR2(2,1)=0.0d0         
      BR2(2,2)=1.0d0          
      BR2(2,3)=0.0d0          
      BR2(3,1)=0.0d0          
      BR2(3,2)=0.0d0          
      BR2(3,3)=1.0d0          
      ortho=.false.
      GOTO 100     
!
 22   BR2(1,1)=1.0d0
      BR2(1,2)=0.0d0
      BR2(1,3)=0.0d0
      BR2(2,1)=0.0d0
      BR2(2,2)=1.0d0*sin(alpha/180.d0*pi)
      BR2(2,3)=1.0d0*cos(alpha/180.d0*pi)
      BR2(3,1)=0.0d0
      BR2(3,2)=0.0d0
      BR2(3,3)=1.0d0
      ortho=.false.
      GOTO 100
!                                                                       
!.....BC CASE (DIRECT LATTICE)                                          
 30   CONTINUE                                                          
      BR2(1,1)=-0.5d0                                                     
      BR2(1,2)=0.5d0                                                      
      BR2(1,3)=0.5d0                                                      
      BR2(2,1)=0.5d0                                                     
      BR2(2,2)=-0.5d0                                                     
      BR2(2,3)=0.5d0                                                      
      BR2(3,1)=0.5d0                                                      
      BR2(3,2)=0.5d0                                                      
      BR2(3,3)=-0.5d0                                                     
      ORTHO=.TRUE.                                             
      GOTO 100                                                          
!                                                                       
!.....FC CASE (DIRECT LATTICE)                                          
 40   CONTINUE                                                          
      BR2(1,1)=0.0d0                                                      
      BR2(1,2)=0.5d0                                                      
      BR2(1,3)=0.5d0                                                      
      BR2(2,1)=0.5d0                                                      
      BR2(2,2)=0.0d0                                                      
      BR2(2,3)=0.5d0                                                      
      BR2(3,1)=0.5d0                                                      
      BR2(3,2)=0.5d0                                                      
      BR2(3,3)=0.0d0                                                      
      ORTHO=.TRUE.                                             
      GOTO 100                                                          
!                                                                       
!.....CXY  CASE (DIRECT LATTICE)                                          
 50   CONTINUE                                                          
      IF(LATTIC(2:3).EQ.'XZ') GOTO 60                                    
      IF(LATTIC(2:3).EQ.'YZ') GOTO 70                                    
      BR2(1,1)=0.5d0                                                      
      BR2(1,2)=-0.5d0                                                     
      BR2(1,3)=0.0d0                                                      
      BR2(2,1)=0.5d0                                                      
      BR2(2,2)=0.5d0                                                      
      BR2(2,3)=0.0d0                                                      
      BR2(3,1)=0.0d0                                                      
      BR2(3,2)=0.0d0                                                      
      BR2(3,3)=1.0d0                                                      
      ORTHO=.TRUE.                                             
      GOTO 100                                                          
!                                                                       
!.....CXZ  CASE (DIRECT LATTICE)                                          
 60   CONTINUE 
!.....CXZ ORTHOROMBIC CASE
      if(gamma.eq.90.d0) then
         BR2(1,1)=0.5d0                                                      
         BR2(1,2)=0.0d0                                                     
         BR2(1,3)=-0.5d0                                                      
         BR2(2,1)=0.0d0                                                      
         BR2(2,2)=1.0d0                                                      
         BR2(2,3)=0.0d0                                                      
         BR2(3,1)=0.5d0                                                     
         BR2(3,2)=0.0d0                                                      
         BR2(3,3)=0.5d0                                                      
         ORTHO=.TRUE.                                             
         GOTO 100                                                          
      ELSE
!.....CXZ MONOCLINIC CASE
         write(*,*) 'gamma not equal 90'
         SINAB=SIN(gamma/180.d0*pi)
         COSAB=COS(gamma/180.d0*pi)
!
         BR2(1,1)=0.5d0*sinab                                                
         BR2(1,2)=0.5d0*cosab                                               
         BR2(1,3)=-0.5d0                                                      
         BR2(2,1)=0.0d0                                                      
         BR2(2,2)=1.0d0                                                      
         BR2(2,3)=0.0d0                                                      
         BR2(3,1)=0.5d0*sinab                                               
         BR2(3,2)=0.5d0*cosab                                                
         BR2(3,3)=0.5d0                                                      
         ORTHO=.FALSE.
         GOTO 100
      ENDIF

!                                                                       
!.....CYZ  CASE (DIRECT LATTICE)                                          
 70   CONTINUE                                                          
      BR2(1,1)=1.0d0                                                      
      BR2(1,2)=0.0d0                                                     
      BR2(1,3)=0.0d0                                                      
      BR2(2,1)=0.0d0                                                      
      BR2(2,2)=0.5d0                                                      
      BR2(2,3)=0.5d0                                                      
      BR2(3,1)=0.0d0                                                      
      BR2(3,2)=-0.5d0                                                     
      BR2(3,3)=0.5d0                                                      
      ORTHO=.TRUE.                                             
      GOTO 100                                                          
!.....RHOMBOHEDRAL CASE
 80   BR2(1,1)=1/2.d0/sqrt(3.d0)
      BR2(1,2)=-1/2.d0                                                     
      BR2(1,3)=1/3.d0                                                      
      BR2(2,1)=1/2.d0/SQRT(3.d0)                                          
      BR2(2,2)=1*0.5d0                                                
      BR2(2,3)=1/3.d0                                                      
      BR2(3,1)=-1/SQRT(3.d0)                                         
      BR2(3,2)=0.d0                                                
      BR2(3,3)=1/3.d0                                                      
      ORTHO=.FALSE.                                             
      GOTO 100                                                         
!                                                                       
 100  CONTINUE                                                          
      write(66,*) 'Bravais Matrix:'
      write(66,999) br2
 999  format(3f15.5)
!                                                                       
      RETURN                                                            
      END                                                               

