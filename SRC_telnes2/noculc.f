!BOP
! !ROUTINE: noculc
! !INTERFACE:
      SUBROUTINE NOCULC (NS, RNUMB,DENSTY,SNUMB,nlast,d,f,v)
! !USES:
      use tetra_params
! !INPUT/OUTPUT PARAMETERS:
!   ns     :
!   rnumb  :
!   densty :
!   snumb  :
!   nlast  :
!   d      :
!   f      :
!   v      :
! !DESCRIPTION:
!  A routine from the tetra package.  Not cleaned for compatibility reasons.
!  THE PARTIAL DS AND NS FROM ONE BAND AND ONE MICROZONE IS CALC.       
! !REVISION HISTORY:
!   Taken from SRC\_tetra/tetra.f of wien2k\_04.10
!   Updated November 2004 (Kevin Jorissen)
!   Cleanup February 2005 (Kevin Jorissen)
!EOP


!KJ old statement      SUBROUTINE NOCULC (NS, RNUMB,DENSTY,SNUMB,nlast)      
! this line changed KJ - I added 3 arguments.                                                                                                                                                            
!KJ IN/OUTPUT
      implicit none !KJ
      integer ns,nlast !KJ
      real*8 RNUMB(nlast,mg), DENSTY(nlast,mg) !KJ
      real*8 SNUMB(nlast,mg)                   !KJ

!KJ LOCALS :      

! KJ I destroyed the common because it objected to my modules :
!KJ      COMMON /EMICRO/ D(4), F(MF,4), V                                  
      real*8 d(4),f(mg,4),v  ! this line added KJ
	  real*8 e0,e1,e2,e3,tp25,t2,t3,tsmall,fac,x,b,emax,emin,efa,e21,e31,e41,e32,e42,e43,e,det !KJ
	  real*8 f1,dos1,dos2,dos3,f2,f3,f4,dos4,dos5,dos6,dos7,dos8,dos9,dos10,fac1,fac2,fac3,fac4,fac5 !KJ
	  integer ll,nf,nfun,i,j,mark,k,kp1,i0,nu,i1,i2,i3,no,iu,nnoc !KJ
      COMMON /NCOUNT/ NNOC                                              
      COMMON /EME/ EMIN, EMAX, EFA, DET, NO, NFUN, NU                   
      DIMENSION X(4), FAC(mg,4), FAC1(mg), FAC2(mg), FAC3(mg)           
      DIMENSION FAC4(mg), FAC5(mg)                                      
! KJ end of my changes

      EQUIVALENCE (E0,X(1)), (E1,X(2)), (E2,X(3)), (E3,X(4))            
      DATA TP25, T2, T3, TSMALL /0.25D0,2.E0,3.E0,0.0000001/            

!     **  X..EIGENVALUES, FAC..EIGENVECTORS OR ANGULAR MOMENTUM WEIGHT  
      do LL=1,4                                                      
        do NF=1,NFUN                                                   
          FAC(NF,LL)=F(NF,LL)                                               
        enddo
        X(LL)=D(LL)                                                       
      enddo
!     **  ORDERING ACCORDING TO INCREASING EIGENVALUES                  
      do I=1,3                                                       
        J=4-I                                                             
        MARK=0                                                            
        do K=1,J                                                       
          KP1=K+1                                                           
          IF (X(K).LE.X(KP1)) cycle
          do NF=1,NFUN                                                   
            B=FAC(NF,KP1)                                                     
            FAC(NF,KP1)=FAC(NF,K)                                             
            FAC(NF,K)=B
          enddo
          B=X(KP1)                                                          
          X(KP1)=X(K)                                                       
          X(K)=B                                                            
          MARK=1                                                            
        enddo
      IF (MARK.EQ.0) exit                                           
	  enddo
                                                  
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
      do I=IU,J                                                      
        E=E+DET                                                           
        DOS1=(E-E0)*(E-E0)*F1                                             
        DOS2=(E-E0)*DOS1                                                  
        DOS1=DOS1*T3                                                      
        DOS3=DOS2*(E-E0)*TP25                                             
!     **  LOOP OVER COMPONENTS                                           
        do NF=1,NFUN                                                   
          F2=(FAC(NF,2)-FAC(NF,1))/E21+(FAC(NF,3)-FAC(NF,1))/E31+(FAC(NF,4)  &
           -FAC(NF,1))/E41                                                  
!     **  INTEGRATED DENSITY OF STATES                                  
          RNUMB(I,NF+NS)=RNUMB(I,NF+NS)+DOS2*FAC(NF,1)+DOS3*F2              
!     **  DENSITY OF STATES                                             
          DENSTY(I,NF+NS)=DENSTY(I,NF+NS)+DOS1*FAC(NF,1)+DOS2*F2
        enddo            
      enddo                                                          
  100 IF (NO.LT.I1) RETURN                                              
!     **  E0,E1<E<E2,E3                                                 
      if (I2.ne.I1) then
        F1=V/E42                                                          
        F2=F1/E42                                                         
        B=F1/E31                                                          
        F3=B/E32                                                          
        F4=B/E41                                                          
        J=MIN0(I2-1,NO)                                                   
        IU=MAX0(I1,IU)                                                    
        do I=IU,J                                                     
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
          do NF=1,NFUN                                                  
            FAC1(NF)=(FAC(NF,3)-FAC(NF,1))/E31+(FAC(NF,4)-FAC(NF,2))/E42+(FAC  &
             (NF,3)-FAC(NF,2))/E32                                            
            FAC2(NF)=(FAC(NF,3)-FAC(NF,1))/E31+(FAC(NF,4)-FAC(NF,2))/E42+(FAC  &
             (NF,4)-FAC(NF,1))/E41                                            
            FAC3(NF)=FAC(NF,1)+T2*FAC(NF,2)+(FAC(NF,3)-FAC(NF,1))*E21/E31     
            FAC4(NF)=T2*FAC(NF,1)+FAC(NF,2)-(FAC(NF,4)-FAC(NF,2))*E21/E42     
            FAC5(NF)=T2*(FAC(NF,1)+FAC(NF,4))+(FAC(NF,3)-FAC(NF,1))*E41/E31   
          enddo
          do NF=1,NFUN                                                  
            RNUMB(I,NF+NS)=RNUMB(I,NF+NS)+DOS1*(FAC(NF,1)+T2*FAC(NF,2)+FAC(NF, &
             3))+DOS2*(FAC(NF,4)-FAC(NF,2))+DOS5*(FAC3(NF)+FAC(NF,3))+DOS6     &
             *FAC1(NF)+DOS9*FAC5(NF)+DOS10*FAC2(NF)                           
            DENSTY(I,NF+NS)=DENSTY(I,NF+NS)+DOS3*FAC3(NF)+DOS4*FAC1(NF)+DOS7   &
             *FAC4(NF)+DOS8*FAC2(NF)                                          
	      enddo
        enddo
	  endif
 140     IF (NO.LT.I2) RETURN                                              
      IF (I3.EQ.I2) GO TO 170                                           
      F1=V/(E43*E42*E41)                                                
      J=MIN0(I3-1,NO)                                                   
      IU=MAX0(I2,IU)                                                    
      do I=IU,J                                                     
        E=E+DET                                                           
        DOS1=(E-E3)*(E-E3)*F1                                             
        DOS2=(E-E3)*DOS1                                                  
        DOS1=DOS1*T3                                                      
        DOS3=DOS2*(E-E3)*TP25                                             
        do NF=1,NFUN                                                  
          F2=(FAC(NF,4)-FAC(NF,1))/E41+(FAC(NF,4)-FAC(NF,2))/E42+(FAC(NF,4)  &
           -FAC(NF,3))/E43                                                  
          RNUMB(I,NF+NS)=RNUMB(I,NF+NS)+DOS2*FAC(NF,4)+DOS3*F2+V*(FAC(NF,1)  &
           +FAC(NF,2)+FAC(NF,3)+FAC(NF,4))*TP25                             
          DENSTY(I,NF+NS)=DENSTY(I,NF+NS)+DOS1*FAC(NF,4)+DOS2*F2            
        enddo                                                         
      enddo
  170 IF (NO.LT.I3) RETURN                                              
      do NF=1,NFUN                                                  
        SNUMB(I3,NF+NS)=SNUMB(I3,NF+NS)+V*(FAC(NF,1)+FAC(NF,2)+FAC(NF,3)+FAC(NF,4))*TP25                                                 
	  enddo
      RETURN                                                            
  190 NNOC=NNOC+1                
      RETURN                                                            
      END                                                               

