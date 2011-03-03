      SUBROUTINE GEN_BRAV(LATTIC, ang,BR1, BR2)
      implicit real*8 (a-h,o-z)
      include 'params.inc'
      dimension br1(3,3),br2(3,3)
!     
!.... DEFINE DIRECT space BRAVAIS MATRIX  BR1(i,*)
!     BR2 is not used
!
      CHARACTER*4 LATTIC

      logical ortho

      dimension ang(3)
      alpha=ang(1)
      beta =ang(2)
      gamma=ang(3)
      pi=acos(-1.D0)
       do i=1,3
       do j=1,3
       br1(i,j)=0.d0
       br2(i,j)=0.d0
       enddo
       enddo

      if(gamma.eq.0.d0) gamma=90.d0
      gamma1=gamma
      alpha1=alpha
      beta1=beta
      alpha=alpha/180.d0*pi
      beta=beta/180.d0*pi                                          
      gamma=gamma/180.d0*pi                                          

!                                                                       
!....DEFINE DIRECT space BRAVAIS MATRIX  BR2(i,*)
      IF(LATTIC(1:1).EQ.'S'.OR.LATTIC(1:1).EQ.'P') THEN
         cosg1=(cos(gamma)-cos(alpha)*cos(beta))/sin(alpha)/sin(beta)
         gamma0=acos(cosg1)
         BR2(1,1)=ALAT(1)*sin(gamma0)*sin(beta)
         BR2(1,2)=ALAT(1)*cos(gamma0)*sin(beta)
         BR2(1,3)=ALAT(1)*cos(beta)    
         BR2(2,1)=0.0d0              
         BR2(2,2)=ALAT(2)*sin(alpha)
         BR2(2,3)=ALAT(2)*cos(alpha)
         BR2(3,1)=0.0d0              
         BR2(3,2)=0.0d0              
         BR2(3,3)=ALAT(3)
         BR1(1,1)=ALAT(1)*sin(gamma0)*sin(beta)
         BR1(1,2)=ALAT(1)*cos(gamma0)*sin(beta)
         BR1(1,3)=ALAT(1)*cos(beta)    
         BR1(2,1)=0.0d0              
         BR1(2,2)=ALAT(2)*sin(alpha)
         BR1(2,3)=ALAT(2)*cos(alpha)
         BR1(3,1)=0.0d0              
         BR1(3,2)=0.0d0              
         BR1(3,3)=ALAT(3)              
         ortho=.true.
         if(gamma1.ne.90.d0) ortho=.false.
         if(alpha1.ne.90.d0) ortho=.false.
         if(beta1.ne.90.d0) ortho=.false.
      ELSE IF(LATTIC(1:1).EQ.'F') THEN                                  
         BR2(1,1)=0.5D0*ALAT(1)
         BR2(2,1)=0.5D0*ALAT(1)
         BR2(2,2)=0.5D0*ALAT(2)
         BR2(3,2)=0.5D0*ALAT(2)
         BR2(1,3)=0.5D0*ALAT(3)
         BR2(3,3)=0.5D0*ALAT(3)
         BR1(1,1)=ALAT(1)
         BR1(2,2)=ALAT(2)
         BR1(3,3)=ALAT(3)
         ortho=.true.                                                 
      ELSE IF(LATTIC(1:1).EQ.'B') THEN                                  
         BR2(1,1)=-0.5D0*ALAT(1)
         BR2(2,1)=0.5D0*ALAT(1)
         BR2(3,1)=0.5D0*ALAT(1)
         BR2(1,2)=0.5D0*ALAT(2)
         BR2(2,2)=-0.5D0*ALAT(2)
         BR2(3,2)=0.5D0*ALAT(2)
         BR2(1,3)=0.5D0*ALAT(3)
         BR2(2,3)=0.5D0*ALAT(3)
         BR2(3,3)=-0.5D0*ALAT(3)
         BR1(1,1)=ALAT(1)
         BR1(2,2)=ALAT(2)
         BR1(3,3)=ALAT(3)
         ortho=.true.                                                 
      ELSE IF(LATTIC(1:1).EQ.'H') THEN                                  
         BR1(1,1)=SQRT(3.D0)/2.D0*ALAT(1)
         BR1(1,2)=-0.5D0*ALAT(2)
         BR1(2,2)=ALAT(2)
         BR1(3,3)=ALAT(3)
         BR2(1,1)=SQRT(3.D0)/2.D0*ALAT(1)
         BR2(1,2)=-0.5D0*ALAT(2)
         BR2(2,2)=ALAT(2)
         BR2(3,3)=ALAT(3)
         ortho=.false.                                                 
      ELSE IF(LATTIC(1:1).EQ.'R') THEN                                  
         BR1(1,1)=1.D0/SQRT(3.D0)/2.D0*ALAT(1)
         BR1(1,2)=-0.5D0*ALAT(2)
         BR1(1,3)=1.D0/3.D0*ALAT(3)
         BR1(2,1)=1.D0/SQRT(3.D0)/2.D0*ALAT(1)
         BR1(2,2)=0.5D0*ALAT(2)
         BR1(2,3)=1.D0/3.D0*ALAT(3)
         BR1(3,1)=-1.D0/SQRT(3.D0)*ALAT(1)
         BR1(3,2)=0.D0                                         
         BR1(3,3)=1.D0/3.D0*ALAT(3)
         BR2(1,1)=1.D0/SQRT(3.D0)/2.D0*ALAT(1)
         BR2(1,2)=-0.5D0*ALAT(2)
         BR2(1,3)=1.D0/3.D0*ALAT(3)
         BR2(2,1)=1.D0/SQRT(3.D0)/2.D0*ALAT(1)
         BR2(2,2)=0.5D0*ALAT(2)
         BR2(2,3)=1.D0/3.D0*ALAT(3)
         BR2(3,1)=-1.D0/SQRT(3.D0)*ALAT(1)
         BR2(3,2)=0.D0                                         
         BR2(3,3)=1.D0/3.D0*ALAT(3)
         ortho=.false.                                                 
      ELSE IF(LATTIC(1:3).EQ.'CXY') THEN 
         BR2(1,1)=0.5D0*ALAT(1)
         BR2(2,1)=0.5D0*ALAT(1)
         BR2(3,1)=0.D0                                                 
         BR2(1,2)=0.5D0*ALAT(2)
         BR2(2,2)=-0.5D0*ALAT(2)
         BR2(3,2)=0.D0                                                 
         BR2(1,3)=0.D0                                                 
         BR2(2,3)=0.D0                                                 
         BR2(3,3)=ALAT(3)                                                
         BR1(1,1)=ALAT(1)
         BR1(2,2)=ALAT(2)
         BR1(3,3)=ALAT(3)
         ortho=.true.                                                 
         do i=1,iord
                do i1=1,3
                        do i2=1,3
                                iz(i2,i1,i+iord)=iz(i2,i1,i)
                        enddo
                enddo
                tau(1,i+iord)=tau(1,i)+0.5
                tau(2,i+iord)=tau(2,i)+0.5
                tau(3,i+iord)=tau(3,i)
        enddo
        iord=iord*2
      ELSE IF(LATTIC(1:3).EQ.'CYZ') THEN 
         BR2(1,1)=ALAT(1)                                                
         BR2(2,1)=0.D0                                                 
         BR2(3,1)=0.D0                                                 
         BR2(1,2)=0.D0                                                 
         BR2(2,2)=-0.5D0*ALAT(2)
         BR2(3,2)=0.5D0*ALAT(2)
         BR2(1,3)=0.D0                                                 
         BR2(2,3)=0.5D0*ALAT(3)
         BR2(3,3)=0.5D0*ALAT(3)
         BR1(1,1)=ALAT(1)
         BR1(2,2)=ALAT(2)
         BR1(3,3)=ALAT(3)
         ortho=.true.                                                 
         do i=1,iord
                do i1=1,3
                        do i2=1,3
                                iz(i2,i1,i+iord)=iz(i2,i1,i)
                        enddo
                enddo
                tau(1,i+iord)=tau(1,i)
                tau(2,i+iord)=tau(2,i)+0.5
                tau(3,i+iord)=tau(3,i)+0.5
        enddo
        iord=iord*2
      ELSE IF(LATTIC(1:3).EQ.'CXZ') THEN                  
         BR2(1,1)=0.5d0*ALAT(1)*sin(gamma)
         BR2(1,2)=0.5d0*ALAT(1)*cos(gamma)
         BR2(1,3)=-0.5d0*ALAT(3)
         BR2(2,2)=ALAT(2)                                                
         BR2(3,1)=0.5d0*ALAT(1)*sin(gamma)
         BR2(3,2)=0.5d0*ALAT(1)*cos(gamma)
         BR2(3,3)=0.5d0*ALAT(3)
         BR1(1,1)=ALAT(1)*sin(gamma)
         br1(1,2)=ALAT(1)*cos(gamma)
         BR1(2,2)=ALAT(2)
         BR1(3,3)=ALAT(3) 
         ortho=.false.
         if(abs(gamma-asin(1.D0)).lt.1D-5)ortho=.true.                                                          
         do i=1,iord
                do i1=1,3
                        do i2=1,3
                                iz(i2,i1,i+iord)=iz(i2,i1,i)
                        enddo
                enddo
                tau(1,i+iord)=tau(1,i)+0.5
                tau(2,i+iord)=tau(2,i)
                tau(3,i+iord)=tau(3,i)+0.5
        enddo
        iord=iord*2

      ELSE                                                              
         STOP 'LATTIC NOT DEFINED'                           
      END IF                                                            
!
!     Patch truncation errors
      do i=1,3
        do j=1,3
                if(abs(br1(j,i)) .lt. 1D-12)br1(j,i)=0.D0
        enddo
      enddo
      RETURN
      END
