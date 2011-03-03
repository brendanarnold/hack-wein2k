      SUBROUTINE STERN (JJA,NST,STG,TAUP,kzz,nkk)                               
!                                                                       
!     STERN GENERATES THE STAR OF REC LATTICE VECTOR KZZ(1,JJA)         
!     THE STAR VECTORS ARE STORED IN STG, THE STAR-SIZE IN  NST         
!     IZ CONTAINS THE SYMMETRY-MATRICES                                 
!                                                                       
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
!
      INCLUDE 'param.inc'
!
      COMPLEX*16          TAUP,IMAG                                     
      INTEGER          G,STG,INDEX(NSYM)                                
      dimension KZZ(3,Nkk)
      COMMON /SYMETR/  TAU(3,NSYM),IZ(3,3,NSYM),             &
                       INUM(NSYM),IORD                                  
!                                                                       
      DIMENSION        STG(3,NSYM),TAUP(NSYM),G(3)                      
      DATA IMAG        /(0.D0,1.D0)/                                    
!-------------------------------------------------------------------    
!                                                                       
      TPI=2.D0*ACOS(-1.D0)                                              
      G(1)=KZZ(1,JJA)                                                   
      G(2)=KZZ(2,JJA)                                                   
      G(3)=KZZ(3,JJA)                                                   
      NST=0                                                             
!                                                                       
!.....START LOOP OVER ALL SYMMETRY OPERATIONS                           
      DO 1 I=1,IORD                                                     
         TK=0.                                                          
         DO 2 J=1,3                                                     
         TK=TK+TAU(J,I)*G(J)*TPI                                     
           K=0                                                         
           DO 3 L=1,3                                                  
   3       K=IZ(J,L,I)*G(L)+K                                          
           STG(J,I)=K                                                     
 2      continue

         IF(NST.EQ.0) GOTO 7                                            
!                                                                       
!........PROOF, IF THE VECTOR STG(J,I) IS A NEW STARMEMBER OR NOT       
         DO 4 M=1,NST                                                   
            DO 5 J=1,3                                                  
               IF(STG(J,M).NE.STG(J,I)) GOTO 4                          
   5        CONTINUE                                                    
!           STG(J,I) IS NOT A NEW STARMEMBER, IT IS EQUIV TO STG(J,M).  
!           BUT THE TAUPHASE OF STG(J,I) CAN BE NEW.  THEREFORE THE     
!           ALREADY DEFINED PHASE TAUP(M) IS AVERAGED WITH THE PHASE    
!           OF STG(J,M)                                                 
!            TAUP(M)=TAUP(M)+EXP(TK*IMAG)                                
            TAUP(M)=TAUP(M)+DCMPLX(cos(TK),sin(TK))
            INDEX(M)=INDEX(M)+1                                         
            GOTO 1                                                      
   4     CONTINUE                                                       
!                                                                       
!........NEW VECTOR FOUND ]]]                                           
   7     NST=NST+1                                                      
         DO 6 J=1,3                                                     
   6     STG(J,NST)=STG(J,I)                                            
!         TAUP(NST)=EXP(TK*IMAG)                                         
         TAUP(NST)=DCMPLX(cos(TK),sin(TK))
         INDEX(NST)=1                                                   
   1  CONTINUE                                                          
!                                                                       
      DO 10 I=1,NST           
      TAUP(I)=TAUP(I)/DBLE(INDEX(I))                                          
  10  continue
      RETURN                                                            
      END                                                               
