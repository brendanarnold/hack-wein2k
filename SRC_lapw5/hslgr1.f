      SUBROUTINE HSLGR1(Z,NAT)       
      use atomgrid                                   
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      DIMENSION R(NPT),SIG(NPT),Z(NAT)                                    
      DATA KNMAX/500/,MM/4/,MRM/350/                                    
      DATA PI/3.1415926538/                                             
      PI=ACOS(-1.D0)                                                      
!                                                                       
      DO 500 N=1,NAT                                                    
      dx(N)=0.024                                                           
      Rnot(N)=0.0005                                                         
      IF(Z(N).LT.0.5) GOTO 500                                          
      RN(1,n)=Rnot(N)                                                          
      DO 200 K=2,KNMAX                                                  
  200 RN(K,n)=Rnot(N)*EXP(dx(N)*(K-1))                                            
      READ(13,1) CASE,KK                                                
      WRITE(6,6) CASE,KK,Z(N)                                           
      KK1=KK-1                                                          
      MESH=KK                                                           
!                                                                       
!  ...CONSTRCT R-MESH OF HERMAN-SKILLMAN-TYPE                           
!                                                                       
      ONTH=1.0/3.                                                       
      C=0.88534138/Z(N)**ONTH                                           
      DELTAR=0.0025*C                                                   
      MBLOCK=MESH/40                                                    
      I=1                                                               
      R(1)=0.                                                           
      DO 220 J=1,MBLOCK                                                 
      IMK=(J-1)*40+1                                                    
      DO 210 K=1,40                                                     
      I=I+1                                                             
      AK=K                                                              
  210 R(I)=R(IMK)+AK*DELTAR                                             
  220 DELTAR=DELTAR+DELTAR                                              
!                                                                       
      READ(13,2) (SIG(K),K=1,KK)                                        
! ....SIG  RADIAL CHARGE DENSITY AS OBTAINED BY 'ATOM'                  
      DO 300 K=1,KNMAX                                                  
      RR=RN(K,n)                                                          
      DO 250 I=2,KK                                                     
      IF(R(I).GT.RR) GOTO 260                                           
  250 CONTINUE                                                          
      GOTO 350                                                          
  260 IF(I.GE.KK1) GOTO 350                                             
      IF(I.EQ.2) I=3                                                    
      KM=K                                                              
      CALL INTERP(R(I-2),SIG(I-2),4,RR,SI,DSI,.FALSE.)                  
  300 RHO1(K,N)=SI/(RR*RR*4.*PI)                                        
  350 CONTINUE                                                          
      NRMAX(N)=KM                                                       
  500 CONTINUE                                                          
    1 FORMAT(A4/I4)                                                     
    2 FORMAT(5E16.8)                                                    
    6 FORMAT(1X,A4,I5,F10.0/)                                           
    9 FORMAT('0RM,RX,MM,MRM',2F12.6,2I5/'0RATIO,H,R1',3F15.8)           
      END                                                               
