      SUBROUTINE GEN_BRAV(LATTIC)
      implicit none
!     
!.... DEFINE DIRECT space BRAVAIS MATRIX  BR2(i,*))
!.... BR1  IS A MATRIX WHICH TRANSFORMS REAL HEX INTO ORTHOG COORD
!.... BR4 is the inverse if br1     
!
      real*8 br1,br2,br3,br4,a,alpha,beta,gamma,gamma1,y
      real*8 beta1,alpha1,cosg1,gamma0

      integer lwork,ipiv,i,j,id,info
      
      CHARACTER*4 LATTIC

      logical ortho

      COMMON /BRAV/ BR1(3,3),BR2(3,3),br3(3,3),br4(3,3)
      COMMON /LATT/ A(3),alpha,beta,gamma,gamma1,ortho
      dimension y(3,3),ipiv(3)
!      DATA BR1/9*0.0/,BR2/9*0.0/                                        
       do i=1,3
       do j=1,3
       br1(i,j)=0.d0
       br2(i,j)=0.d0
       enddo
       enddo
      lwork=30

      if(gamma.eq.0.d0) gamma=90.d0
      gamma1=gamma
      alpha1=alpha
      beta1=beta
      alpha=alpha/180.d0*acos(-1.d0)                                          
      beta=beta/180.d0*acos(-1.d0)                                          
      gamma=gamma/180.d0*acos(-1.d0)                                          

!                                                                       
!....DEFINE DIRECT space BRAVAIS MATRIX  BR2(i,*))                            
!....BR1  IS A MATRIX WHICH TRANSFORMS REAL HEX INTO ORTHOG COORD       
      IF(LATTIC(1:1).EQ.'S'.OR.LATTIC(1:1).EQ.'P') THEN
         cosg1=(cos(gamma)-cos(alpha)*cos(beta))/sin(alpha)/sin(beta)
         gamma0=acos(cosg1)
         BR2(1,1)=A(1)*sin(gamma0)*sin(beta)
         BR2(1,2)=A(1)*cos(gamma0)*sin(beta)
         BR2(1,3)=A(1)*cos(beta)    
         BR2(2,1)=0.0d0              
         BR2(2,2)=A(2)*sin(alpha)
         BR2(2,3)=A(2)*cos(alpha)
         BR2(3,1)=0.0d0              
         BR2(3,2)=0.0d0              
         BR2(3,3)=A(3)
         BR1(1,1)=A(1)*sin(gamma0)*sin(beta)
         BR1(1,2)=A(1)*cos(gamma0)*sin(beta)
         BR1(1,3)=A(1)*cos(beta)    
         BR1(2,1)=0.0d0              
         BR1(2,2)=A(2)*sin(alpha)
         BR1(2,3)=A(2)*cos(alpha)
         BR1(3,1)=0.0d0              
         BR1(3,2)=0.0d0              
         BR1(3,3)=A(3)              
!  old version with monoclinic lattic only
!      BR2(1,1)=A(1)*sin(gamma)
!      br2(1,2)=a(1)*cos(gamma)                                               
!      BR2(2,2)=A(2)                                                     
!      BR2(3,3)=A(3)                                                     
!      BR1(1,1)=A(1)*sin(gamma)
!      br1(1,2)=a(1)*cos(gamma)                                             
!      BR1(2,2)=A(2)                                                  
!      BR1(3,3)=A(3) 
         ortho=.true.
         if(gamma1.ne.90.d0) ortho=.false.
         if(alpha1.ne.90.d0) ortho=.false.
         if(beta1.ne.90.d0) ortho=.false.
      ELSE IF(LATTIC(1:1).EQ.'F') THEN                                  
         BR2(1,1)=0.5D0*A(1)
         BR2(2,1)=0.5D0*A(1)
         BR2(2,2)=0.5D0*A(2)
         BR2(3,2)=0.5D0*A(2)
         BR2(1,3)=0.5D0*A(3)
         BR2(3,3)=0.5D0*A(3)
         BR1(1,1)=A(1)
         BR1(2,2)=A(2)
         BR1(3,3)=A(3)
         ortho=.true.                                                 
      ELSE IF(LATTIC(1:1).EQ.'B') THEN                                  
         BR2(1,1)=-0.5D0*A(1)
         BR2(2,1)=0.5D0*A(1)
         BR2(3,1)=0.5D0*A(1)
         BR2(1,2)=0.5D0*A(2)
         BR2(2,2)=-0.5D0*A(2)
         BR2(3,2)=0.5D0*A(2)
         BR2(1,3)=0.5D0*A(3)
         BR2(2,3)=0.5D0*A(3)
         BR2(3,3)=-0.5D0*A(3)
         BR1(1,1)=A(1)
         BR1(2,2)=A(2)
         BR1(3,3)=A(3)
         ortho=.true.                                                 
      ELSE IF(LATTIC(1:1).EQ.'H') THEN                                  
         BR1(1,1)=SQRT(3.D0)/2.D0*A(1)
         BR1(1,2)=-0.5D0*A(2)
         BR1(2,2)=A(2)
         BR1(3,3)=A(3)
         BR2(1,1)=SQRT(3.D0)/2.D0*A(1)
         BR2(1,2)=-0.5D0*A(2)
         BR2(2,2)=A(2)
         BR2(3,3)=A(3)
         ortho=.false.                                                 
      ELSE IF(LATTIC(1:1).EQ.'R') THEN                                  
         BR1(1,1)=1.D0/SQRT(3.D0)/2.D0*A(1)
         BR1(1,2)=-0.5D0*A(2)
         BR1(1,3)=1.D0/3.D0*A(3)
         BR1(2,1)=1.D0/SQRT(3.D0)/2.D0*A(1)
         BR1(2,2)=0.5D0*A(2)
         BR1(2,3)=1.D0/3.D0*A(3)
         BR1(3,1)=-1.D0/SQRT(3.D0)*A(1)
         BR1(3,2)=0.D0                                         
         BR1(3,3)=1.D0/3.D0*A(3)
         BR2(1,1)=1.D0/SQRT(3.D0)/2.D0*A(1)
         BR2(1,2)=-0.5D0*A(2)
         BR2(1,3)=1.D0/3.D0*A(3)
         BR2(2,1)=1.D0/SQRT(3.D0)/2.D0*A(1)
         BR2(2,2)=0.5D0*A(2)
         BR2(2,3)=1.D0/3.D0*A(3)
         BR2(3,1)=-1.D0/SQRT(3.D0)*A(1)
         BR2(3,2)=0.D0                                         
         BR2(3,3)=1.D0/3.D0*A(3)
         ortho=.false.                                                 
      ELSE IF(LATTIC(1:3).EQ.'CXY') THEN 
         BR2(1,1)=0.5D0*A(1)
         BR2(2,1)=0.5D0*A(1)
         BR2(3,1)=0.D0                                                 
         BR2(1,2)=0.5D0*A(2)
         BR2(2,2)=-0.5D0*A(2)
         BR2(3,2)=0.D0                                                 
         BR2(1,3)=0.D0                                                 
         BR2(2,3)=0.D0                                                 
         BR2(3,3)=A(3)                                                
         BR1(1,1)=A(1)
         BR1(2,2)=A(2)
         BR1(3,3)=A(3)
         ortho=.true.                                                 
      ELSE IF(LATTIC(1:3).EQ.'CYZ') THEN 
         BR2(1,1)=A(1)                                                
         BR2(2,1)=0.D0                                                 
         BR2(3,1)=0.D0                                                 
         BR2(1,2)=0.D0                                                 
         BR2(2,2)=-0.5D0*A(2)
         BR2(3,2)=0.5D0*A(2)
         BR2(1,3)=0.D0                                                 
         BR2(2,3)=0.5D0*A(3)
         BR2(3,3)=0.5D0*A(3)
         BR1(1,1)=A(1)
         BR1(2,2)=A(2)
         BR1(3,3)=A(3)
         ortho=.true.                                                 
      ELSE IF(LATTIC(1:3).EQ.'CXZ') THEN                  
         BR2(1,1)=0.5d0*A(1)*sin(gamma)
         BR2(1,2)=0.5d0*A(1)*cos(gamma)
         BR2(1,3)=-0.5d0*A(3)
         BR2(2,2)=A(2)                                                
         BR2(3,1)=0.5d0*A(1)*sin(gamma)
         BR2(3,2)=0.5d0*A(1)*cos(gamma)
         BR2(3,3)=0.5d0*A(3)
         BR1(1,1)=A(1)*sin(gamma)
         br1(1,2)=A(1)*cos(gamma)
         BR1(2,2)=A(2)
         BR1(3,3)=A(3) 
         ortho=.false.
         if(abs(gamma-asin(1.D0)).lt.1D-5)ortho=.true.                                                          
      ELSE                                                              
         STOP 'LATTIC NOT DEFINED'                           
      END IF                                                            
      call gbass(br1,br3)
      do i=1,3
        do j=1,3
          br4(i,j)=br1(i,j)
        enddo
      enddo

!$$$      call dgetrf (3,3,br4,3,ipiv,info)
!$$$      if ( info .ne. 0 ) then
!$$$        write(6,*) '   INFOtrf :', info
!$$$        STOP 'ERROR BRAVAIS'
!$$$      endif
!$$$      call dgetri (3,br4,3,ipiv,work,lwork,info)
!$$$      if ( info .ne. 0 ) then
!$$$        write(6,*) '   INFOtri :', info
!$$$        STOP 'ERROR BRAVAIS'
!$$$      endif

      call ludcmp(br4,3,3,ipiv,id,info)
      if(info.ne.0) then
        write(6,*) 'Error inverting br1:'
        do i=1,3
          write(6,*) (br4(i,j),j=1,3)
        enddo
        stop 'ERROR INVERTING BR1'
      endif

      do  i=1,3
        do j=1,3
          y(i,j)=0.D0
        enddo 
        y(i,i)=1.D0
      enddo 
      do j=1,3
        call lubksb(br4,3,3,ipiv,y(1,j))
      enddo         
      do i=1,3
        do j=1,3
          br4(i,j)=y(i,j)
        enddo
      enddo

      RETURN
      END
