      SUBROUTINE vorbmix(vorbNEW,vorbOLD,natorb,nlorb,lorb,chmix,qmx)
!    orbital potential MIXING USING BROYDENS METHOD                        
!    modified QMIX5 program (SRC_mixer) P.Novak, July 2001    
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
!!      PARAMETER (LABC=3)                        
!!      PARAMETER (NBSIZE=NATOrb*LABC*(2*LABC+1)*(2*LABC+1))          
!                                                                       
      COMPLEX*16 vorbnew(natorb,3,-3:3,-3:3)
      COMPLEX*16 vorbold(natorb,3,-3:3,-3:3)
!                                                                       
      COMPLEX*16,allocatable :: YB(:),PS(:,:),SB(:),FN1(:),T1(:),VT(:),UI(:)
      complex*16 zero,one
      character*5 chmix
      dimension nlorb(natorb),lorb(labc,natorb)
      NBSIZE=NATOrb*LABC*(2*LABC+1)*(2*LABC+1)
      allocate ( YB(NBSIZE),PS(NBSIZE,2),  &           
         SB(NBSIZE),FN1(NBSIZE),T1(NBSIZE),VT(NBSIZE),UI(NBSIZE))   
      if(chmix.eq.'PRATT')then
      WRITE(6,1061) chmix, QMX                                            
      WRITE(21,1061) chmix, QMX                                           
 1061 FORMAT(7X,A5,' mixing scheme with',F6.3)                          
      DO 118 N=1,NATorb
      DO 118 LH=1,Nlorb(n)
      L=lorb(lh,n)                             
      DO 118 m=-l,l  
      do 118 m1=-l,l                                
      vorbNEW(N,LH,m,m1)=vorbNEW(N,LH,m,m1)*QMX+  &                 
     vorbOLD(N,LH,m,m1)*(1.-QMX)                                    
  118 CONTINUE                                                          
      return
      endif
!
      zero=(0.,0.)
      one= (1.,0.)
      IS=0                                                 
      DO 250 N=1,NATORB
      DO 250 LH=1,NLORB(N)
      L=lorb(lh,n)                             
      DO 250 m=-l,l  
      do 250 m1=-l,l                                
      IS=IS+1                                                           
      PS(IS,1)=vorbold(n,lh,m,m1)                                       
      PS(IS,2)=vorbnew(n,lh,m,m1)                                       
  250 CONTINUE                                                          
!                                                                       
      MAXMIX=IS                                                         
!                                                                       
      READ(31,END=119)DMIX,LASTIT                                       
      WRITE(6,1061) chmix, QMX                                            
      WRITE(21,1061) chmix, QMX                                           
      DMIX=QMX                                                          
      READ(31)(FN1(K),K=1,MAXMIX)                                       
      READ(31)(SB(K),K=1,MAXMIX)                                        
      YNORM=ZERO                                                        
      DO 104 K=1,MAXMIX                                                 
      SB(K)=PS(K,1)-SB(K)                                               
      TEMP=PS(K,2)-PS(K,1)                                              
      YB(K)=TEMP-FN1(K)                                                 
      YNORM=YNORM+YB(K)**2                                              
  104 FN1(K)=TEMP                                                       
      LASTIT=LASTIT+1                                                   
      REWIND 31                                                         
      REWIND 32                                                         
      WRITE(31)DMIX,LASTIT                                              
      WRITE(31)(FN1(K),K=1,MAXMIX)                                      
      WRITE(31)(PS(K,1),K=1,MAXMIX)                                     
      DO 105 K=1,MAXMIX                                                 
      UI(K)=DMIX*YB(K)+SB(K)                                            
      VT(K)=YB(K)/YNORM                                                 
  105 CONTINUE                                                          
      IF(LASTIT.NE.2)THEN                                               
       JMAX=LASTIT-1                                                    
       DO 500 J=2,JMAX                                                   
       READ(32)(SB(K),K=1,MAXMIX)                                       
       READ(32)(T1(K),K=1,MAXMIX)                                       
       AKJ=ZERO                                                         
       DO 501 K=1,MAXMIX                                                 
  501  AKJ=AKJ+T1(K)*YB(K)                                              
       DO 500 K=1,MAXMIX                                                 
       UI(K)=UI(K)-AKJ*SB(K)                                            
  500  CONTINUE                                                         
      ENDIF                                                             
      WRITE(32)(UI(K),K=1,MAXMIX)                                       
      WRITE(32)(VT(K),K=1,MAXMIX)                                       
      REWIND 32                                                         
      XX=(ONE-DMIX)                                                     
      YY=DMIX                                                           
      DO 502 K=1,MAXMIX                                                  
  502 SB(K)=XX*PS(K,1)+YY*PS(K,2)                                       
      DO 504 J=2,LASTIT                                                  
      READ(32)(UI(K),K=1,MAXMIX)                                        
      READ(32)(VT(K),K=1,MAXMIX)                                        
      CKI=ZERO                                                          
      DO 503 K=1,MAXMIX                                                  
  503 CKI=CKI+VT(K)*FN1(K)                                              
      DO 504 K=1,MAXMIX                                                  
  504 SB(K)=SB(K)-CKI*UI(K)                                             
!                                                                       
      IS=0  
      DO 280 N=1,natorb                                                  
      DO 280 LH=1,NLORB(N)                                               
      L=lorb(lh,n)                             
      DO 280 m=-l,l  
      do 280 m1=-l,l                                
      IS=IS+1                                                           
      vorbNEW(N,LH,m,m1)=SB(IS)                                       
  280 CONTINUE                                                          
      REWIND 31                                                         
      REWIND 32                                                         
      GOTO 120                                                          
! ... PRATT MIXING FOR FIRST ITERATION                                  
  119 REWIND 31                                                         
      chmix='PRATT'
      WRITE(6,1061) chmix, QMX                                            
      WRITE(21,1061) chmix, QMX                                           
      LASTIT=1                                                          
      DMIX=QMX                                                          
      WRITE(31)DMIX,LASTIT                                              
      DO 101 K=1,MAXMIX                                                 
  101 FN1(K)=PS(K,2)-PS(K,1)                                            
      WRITE(31)(FN1(K),K=1,MAXMIX)                                      
      WRITE(31)(PS(K,1),K=1,MAXMIX)                                     
      DO 117 N=1,NATorb
      DO 117 LH=1,Nlorb(n)
      L=lorb(lh,n)                             
      DO 117 m=-l,l  
      do 117 m1=-l,l                                
      vorbNEW(N,LH,m,m1)=vorbNEW(N,LH,m,m1)*DMIX+ &                  
     vorbOLD(N,LH,m,m1)*(1.-DMIX)                                    
  117 CONTINUE                                                          
!                                                                       
  120 CONTINUE                                                          
      RETURN                                                            
      END                                                               
