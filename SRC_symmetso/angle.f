      SUBROUTINE ANGLE(XMS,alat,alpha,lattic,theta,phi)
      IMPLICIT REAL*8 (A-H,O-Z)
      character*4 lattic
      logical     ortho
      dimension alat(3),alpha(3),br2(3,3)
!*******************************************************************
!
      DIMENSION XMS(3),dif(3),help(3)

      do i=1,3
      dif(i)=xms(i)
      enddo
      pi=4.d0*atan(1.d0)
      alpha0=alpha(1)
      beta0=alpha(2)
      gamma=alpha(3)
      gamma1=alpha(3)*pi/180.d0
      beta1=alpha(2)*pi/180.d0
      alpha1=alpha(1)*pi/180.d0
      cosg1=0.
      IF(LATTIC(1:1).EQ.'H') GOTO 10                                    
      IF(LATTIC(1:1).EQ.'S') GOTO 20                                    
      IF(LATTIC(1:1).EQ.'P') GOTO 20                                    
      IF(LATTIC(1:1).EQ.'B') GOTO 30                                    
      IF(LATTIC(1:1).EQ.'F') GOTO 40                                    
      IF(LATTIC(1:1).EQ.'C') GOTO 50                                    
      IF(LATTIC(1:1).EQ.'R') GOTO 80                                    
      write(6,*)lattic
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
 20   continue
!
      cosg1=(cos(gamma1)-cos(alpha1)*cos(beta1))/sin(alpha1)/sin(beta1)
      gamma0=acos(cosg1)
      BR2(1,1)=1.0d0*sin(gamma0)*sin(beta1)                
      BR2(1,2)=1.0d0*cos(gamma0)*sin(beta1)                 
      BR2(1,3)=1.0d0*cos(beta1)                                   
      BR2(2,1)=0.0d0                                                      
      BR2(2,2)=1.0d0*sin(alpha1)
      BR2(2,3)=1.0d0*cos(alpha1)
      BR2(3,1)=0.0d0                                                      
      BR2(3,2)=0.0d0                                                      
      BR2(3,3)=1.0d0                                                      
      ORTHO=.TRUE. 
      if(alpha(3).ne.90.d0) ortho=.false.                                      
      if(alpha(2).ne.90.d0) ortho=.false.                                      
      if(alpha(1).ne.90.d0) ortho=.false.                                      
        write(*,*) alpha0,beta0,gamma0,ortho,br2
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
         SINAB=SIN(gamma1)
         COSAB=COS(gamma1)
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
 999  format(3f15.5)
!                                                                       
      cosgam=cos(alpha(3)/180.d0*pi)
      singam=sin(alpha(3)/180.d0*pi)           
      IF (.not.ortho) THEN                                       
        help(1)=dif(1)  
        help(2)=dif(2)  
        help(3)=dif(3)  
      if(lattic(1:1).eq.'R') then
        dif(1)=help(1)*BR2(1,1)+help(2)*BR2(2,1)+help(3)*BR2(3,1)             
        dif(2)=help(1)*BR2(1,2)+help(2)*BR2(2,2)+help(3)*BR2(3,2)             
        dif(3)=help(1)*BR2(1,3)+help(2)*BR2(2,3)+help(3)*BR2(3,3)           
      elseif(lattic(1:3).eq.'CXZ') then
        dif(1)=help(1)*singam            
        dif(2)=(help(1)*cosgam*alat(1)+help(2)*alat(2))/alat(2)             
        dif(3)=help(3)           
      else
        dif(1)=(help(1)*BR2(1,1)*alat(1)+help(2)*BR2(2,1)*alat(2)+   &
                help(3)*BR2(3,1)*alat(3))/alat(1)             
        dif(2)=(help(1)*BR2(1,2)*alat(1)+help(2)*BR2(2,2)*alat(2)+   &
                help(3)*BR2(3,2)*alat(3))/alat(2)             
        dif(3)=(help(1)*BR2(1,3)*alat(1)+help(2)*BR2(2,3)*alat(2)+   &
                help(3)*BR2(3,3)*alat(3))/alat(3)           
      endif
      endif
      dist=0.d0                                                             
      DO  L=1,3                                                      
      DIST=DIST+dif(L)*dif(L)                                 
      enddo
     
      DIST=SQRT(DIST)                                                   
        THETA=ACOS(dif(3)/dist)
	XX=SQRT(dif(1)**2+dif(2)**2)
        IF (XX.LT.1D-5) THEN
        PHI=0.D0
        ELSE
	PHI=ACOS(dif(1)/XX)
	IF (ABS(dif(2)).GT.1D-5) PHI=PHI*dif(2)/ABS(dif(2))
	END IF
      write(6,*) 'THETA, PHI', theta,phi
      RETURN                                                            
      END                                                               
