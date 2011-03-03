      SUBROUTINE HSLGRS(Z,NAT)
      use atomgrid                                     
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      character*11 gitter
      DIMENSION SIG(npt),Z(NAT)                                    
      DATA MM/4/,MRM/350/                                    
      DATA PI/3.1415926538/                                             
      allocate ( RN(npt,nat),Rnot(nat),dx(nat),NRMAX(NAT))
      allocate ( RHO1(npt,NAT))                  
!      H=0.024                                                           
!      R1=0.00005                                                        
      PI=ACOS(-1.D0)                                                      
!.....file12   sigma from lstart
!.....file13   sigma from hsatom
!     
      read(12,*,IOSTAT=ITEST)
        if(ITEST.ne.0) then
        call hslgr1(Z,NAT)
        return
        end if
      rewind(12)                                                              
      DO 500 N=1,NAT                                                    
      IF(Z(N).LT.0.5) GOTO 500                                          
      READ(12,1) CASE,KK,Rnot(n),dx(n)                                               
      WRITE(6,6) CASE,KK,Z(N) 
      if (kk.gt.npt) then 
         write(6,*) 'kk is too large', kk
         stop ' hsl'
      endif
      KK1=KK-1                                                          
      MESH=KK
!                                                           
! construct log-mesh with H,R1 from lstart
!                                                                      
      Rn(1,n)=Rnot(n)                                                         
      DO 200 K=2,kk                                               
  200 Rn(K,n)=Rnot(n)*EXP(dx(n)*(K-1))                                             
!
      READ(12,2) (SIG(K),K=1,KK)                 
! ....SIG  RADIAL CHARGE DENSITY AS OBTAINED BY 'lstart'
!
      DO 300 K=1,kk                                                     
  300 RHO1(K,N)=SIG(K)/(Rn(k,n)*Rn(k,n)*4.*PI)                                     
  350 CONTINUE                                                          
      KM=K                                                             
      NRMAX(N)=KM                                                    
  500 CONTINUE                                                          
    1 FORMAT(A80/I4,2E12.6)                                                  
    2 FORMAT(5d16.9)                                                    
    6 FORMAT(1X,A80/I5,F10.0/)                                           
    9 FORMAT('0RM,RX,MM,MRM',2F12.6,2I5/'0RATIO,H,R1',3F15.8)           
      END                                                               
