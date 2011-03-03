!............................................................NOCULC.... 
!                                                                       
      SUBROUTINE NOCULC (NS, RNUMB,DENSTY,SNUMB,nlast)                                            
      INCLUDE 'param.inc'
      real*4 RNUMB(nlast,MF), DENSTY(nlast,MF)
      real*4 SNUMB(nlast,MF)                                        
      COMMON /NCOUNT/ NNOC                                              
!  THE PARTIAL DS AND NS FROM ONE BAND AND ONE MICROZONE IS CALC.       
!      REAL *4 RNUMB, DENSTY, NOS                                         
      COMMON /EME/ EMIN, EMAX, EFA, DET, NO, NFUN, NU                   
      DIMENSION X(4), FAC(MF,4), FAC1(MF), FAC2(MF), FAC3(MF)           
      DIMENSION FAC4(MF), FAC5(MF)                                      
      COMMON /EMICRO/ D(4), F(MF,4), V                                  
      EQUIVALENCE (E0,X(1)), (E1,X(2)), (E2,X(3)), (E3,X(4))            
      DATA TP25, T2, T3, TSMALL /0.25D0,2.E0,3.E0,0.0000001/            
!     **  X..EIGENVALUES, FAC..EIGENVECTORS OR ANGULAR MOMENTUM WEIGHT  
      DO 20 LL=1,4                                                      
      DO 10 NF=1,NFUN                                                   
   10 FAC(NF,LL)=F(NF,LL)                                               
   20 X(LL)=D(LL)                                                       
!     **  ORDERING ACCORDING TO INCREASING EIGENVALUES                  
      DO 50 I=1,3                                                       
      J=4-I                                                             
      MARK=0                                                            
      DO 40 K=1,J                                                       
      KP1=K+1                                                           
      IF (X(K).LE.X(KP1)) GO TO 40                                      
      DO 30 NF=1,NFUN                                                   
      B=FAC(NF,KP1)                                                     
      FAC(NF,KP1)=FAC(NF,K)                                             
   30 FAC(NF,K)=B                                                       
      B=X(KP1)                                                          
      X(KP1)=X(K)                                                       
      X(K)=B                                                            
      MARK=1                                                            
   40 CONTINUE                                                          
      IF (MARK.EQ.0) GO TO 60                                           
   50 CONTINUE                                                          
   60 CONTINUE                                                          
!     CHECK=(E1-E0)*(E2-E1)*(E3-E2)                                     
!     IF(CHECK.NE.0.0)GO TO 50                                          
!     WRITE(6,6000)E0,E1,E2,E3                                          
!6000 FORMAT(1H ,' E0,E1,E2,E3 =',4D16.8)                               
!     IF(E0.EQ.E1)E0=E0-3.0D-04                                         
!     IF(E1.EQ.E2)E1=E1-2.0D-04                                         
!     IF(E2.EQ.E3)E2=E2-1.0D-04                                         
!  50 CONTINUE                                                          
!  X(1),X(2),X(3),X(4) ARE IN INCREASING ORDER                          
      IF (X(1).GT.EMAX) RETURN                                          
!     **  INDICES OF EDGEPOINTS ON ENERGY MESH                          
      I0=(E0-EMIN+TSMALL)*EFA+NU+1                                      
      I1=(E1-EMIN+TSMALL)*EFA+NU+1                                      
      I2=(E2-EMIN+TSMALL)*EFA+NU+1                                      
      I3=(E3-EMIN+TSMALL)*EFA+NU+1                                      
      IF (NO.LT.I0) RETURN                                              
      IF (NU.GE.I3) GO TO 190                                           
      E21=E1-E0                                                         
      E31=E2-E0                                                         
      E41=E3-E0                                                         
      E32=E2-E1                                                         
      E42=E3-E1                                                         
      E43=E3-E2                                                         
      IU=MAX0(I0,NU)                                                    
      E=EMIN+(IU-NU-1)*DET                                              
      IF (IU.LT.I1) GO TO 70                                            
      IF (IU.LT.I2) GO TO 100                                           
      GO TO 140                                                         
!     **  E0<E<E1,E2,E3                                                 
   70 IF (I1.EQ.I0) GO TO 100                                           
      F1=V/(E21*E31*E41)                                                
      J=MIN0(I1-1,NO)                                                   
      IU=MAX0(I0,IU)                                                    
!     **  LOOP OVER E                                                   
      DO 90 I=IU,J                                                      
      E=E+DET                                                           
      DOS1=(E-E0)*(E-E0)*F1                                             
      DOS2=(E-E0)*DOS1                                                  
      DOS1=DOS1*T3                                                      
      DOS3=DOS2*(E-E0)*TP25                                             
!     **  LOOP OVERKOMPONENTS                                           
      DO 80 NF=1,NFUN                                                   
      F2=(FAC(NF,2)-FAC(NF,1))/E21+(FAC(NF,3)-FAC(NF,1))/E31+(FAC(NF,4)  &
       -FAC(NF,1))/E41                                                  
!     **  INTEGRATED DENSITY OF STATES                                  
      RNUMB(I,NF+NS)=RNUMB(I,NF+NS)+DOS2*FAC(NF,1)+DOS3*F2              
!     **  DENSITY OF STATES                                             
   80 DENSTY(I,NF+NS)=DENSTY(I,NF+NS)+DOS1*FAC(NF,1)+DOS2*F2            
   90 CONTINUE                                                          
  100 IF (NO.LT.I1) RETURN                                              
!     **  E0,E1<E<E2,E3                                                 
      IF (I2.EQ.I1) GO TO 140                                           
      F1=V/E42                                                          
      F2=F1/E42                                                         
      B=F1/E31                                                          
      F3=B/E32                                                          
      F4=B/E41                                                          
      J=MIN0(I2-1,NO)                                                   
      IU=MAX0(I1,IU)                                                    
      DO 130 I=IU,J                                                     
      E=E+DET                                                           
      DOS1=F1*(E-E1)*TP25                                               
      DOS2=F2*(E-E1)**2*TP25                                            
      DOS3=F3*(E2-E)*(E-E1)                                             
      DOS4=DOS3*(E-E1)                                                  
      DOS5=-DOS3*(E2-E)*TP25                                            
      DOS6=DOS5*(E-E1)                                                  
      DOS7=F4*(E3-E)*(E-E0)                                             
      DOS8=DOS7*(E-E0)                                                  
      DOS9=DOS8*TP25                                                    
      DOS10=-DOS9*(E3-E)                                                
      DO 110 NF=1,NFUN                                                  
      FAC1(NF)=(FAC(NF,3)-FAC(NF,1))/E31+(FAC(NF,4)-FAC(NF,2))/E42+(FAC  &
       (NF,3)-FAC(NF,2))/E32                                            
      FAC2(NF)=(FAC(NF,3)-FAC(NF,1))/E31+(FAC(NF,4)-FAC(NF,2))/E42+(FAC  &
       (NF,4)-FAC(NF,1))/E41                                            
      FAC3(NF)=FAC(NF,1)+T2*FAC(NF,2)+(FAC(NF,3)-FAC(NF,1))*E21/E31     
      FAC4(NF)=T2*FAC(NF,1)+FAC(NF,2)-(FAC(NF,4)-FAC(NF,2))*E21/E42     
  110 FAC5(NF)=T2*(FAC(NF,1)+FAC(NF,4))+(FAC(NF,3)-FAC(NF,1))*E41/E31   
      DO 120 NF=1,NFUN                                                  
      RNUMB(I,NF+NS)=RNUMB(I,NF+NS)+DOS1*(FAC(NF,1)+T2*FAC(NF,2)+FAC(NF, &
       3))+DOS2*(FAC(NF,4)-FAC(NF,2))+DOS5*(FAC3(NF)+FAC(NF,3))+DOS6     &
       *FAC1(NF)+DOS9*FAC5(NF)+DOS10*FAC2(NF)                           
  120 DENSTY(I,NF+NS)=DENSTY(I,NF+NS)+DOS3*FAC3(NF)+DOS4*FAC1(NF)+DOS7   &
       *FAC4(NF)+DOS8*FAC2(NF)                                          
  130 CONTINUE                                                          
  140 IF (NO.LT.I2) RETURN                                              
      IF (I3.EQ.I2) GO TO 170                                           
      F1=V/(E43*E42*E41)                                                
      J=MIN0(I3-1,NO)                                                   
      IU=MAX0(I2,IU)                                                    
      DO 160 I=IU,J                                                     
      E=E+DET                                                           
      DOS1=(E-E3)*(E-E3)*F1                                             
      DOS2=(E-E3)*DOS1                                                  
      DOS1=DOS1*T3                                                      
      DOS3=DOS2*(E-E3)*TP25                                             
      DO 150 NF=1,NFUN                                                  
      F2=(FAC(NF,4)-FAC(NF,1))/E41+(FAC(NF,4)-FAC(NF,2))/E42+(FAC(NF,4)  &
       -FAC(NF,3))/E43                                                  
      RNUMB(I,NF+NS)=RNUMB(I,NF+NS)+DOS2*FAC(NF,4)+DOS3*F2+V*(FAC(NF,1)  &
       +FAC(NF,2)+FAC(NF,3)+FAC(NF,4))*TP25                             
  150 DENSTY(I,NF+NS)=DENSTY(I,NF+NS)+DOS1*FAC(NF,4)+DOS2*F2            
  160 CONTINUE                                                          
  170 IF (NO.LT.I3) RETURN                                              
      DO 180 NF=1,NFUN                                                  
  180 SNUMB(I3,NF+NS)=SNUMB(I3,NF+NS)+V*(FAC(NF,1)+FAC(NF,2)+FAC(NF,3)   &
       +FAC(NF,4))*TP25                                                 
      RETURN                                                            
  190 NNOC=NNOC+1                                                       
      RETURN                                                            
      END                                                               
