      subroutine lsqfit(x,q,nstep,jspin,nat)                            
!      LEAST SQUARE FIT  AND       
!     STORED IN ARRAY VEX.                                              
!                                                                       
!     THE  LINEAR EQUATION SYSTEM  OF THE FIT IS  SOLVED  BY  THE       
!     LINPACK-ROUTINES DSICO (FACTORISATION) AND DSISL (SOLUTION)       
!     IN D O U B L E  P R E C I S I O N  MODE. FOR SINGLE PRECISION     
!     USE THE ROUTINES SSICO AND SSISL (CDC) !                          
!                                                                       
!
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
!
      character*1 num(0:9)
      character*10 pos,cha,cdn,cup
      dimension AA(NATO*3,NATO*3),X(mxpm,maxit),weight(maxit)
      dimension B(NATO*3),Q(NATO,2,maxit)
      dimension tt0(NATO),TT1(NATO),TT2(NATO)
      DIMENSION WSPACE(NATO*3),KPVT(NATO*3)                         
!---------------------------------------------------------------------  
      data num /'0','1','2','3','4','5','6','7','8','9'/
!.....read q,x for all steps
!
!      write(6,*) 'fitpoints:'
!      read(*,*) inter0
      inter0=30
!      nat=6
!      jspin=1
      icount=0
!
!      open(9,file='co_hcp.tmpM',status='old')
!      ind=1
! 887  read(9,*,end=888)
!      do i=1,nat*3
!        read(9,*) x(i,ind)
!      enddo
!      ind=ind+1
!      goto 887
! 888   write(6,*) ind-1,' timesteps read'
!
!      do jatom=1,nat
!      iunit=9+jatom
!      pos='pos0'//num(jatom)
!      cup='cup0'//num(jatom)
!      open(iunit,file=pos,status='old')
!      open(iunit+40,file=cup,status='old')
!
!      ind=1
! 3    continue
!     read(iunit,*,end=2) (x(ind,lm1),lm1=jatom*3-2,jatom*3)
!      weight(ind)=1.d0
!      read(iunit+40,*,end=2) q(jatom,1,ind)
!      write(6,*) (x(ind,lm1),lm1=jatom*3-2,jatom*3),q(jatom,1,ind)
!      ind=ind+1
!      goto 3
! 2    nstep=ind-1
!      enddo
!
!      write(6,*) nstep,' timesteps read'
!
      LLMM=nat*3                                                 
!************************************************
      DO 4 I=1,NATO*3                                                   
      DO 4 J=1,NATO*3                                                   
 4    AA(I,J)=0.0D0                                                     
!************************************************
      write(6,*) nstep,inter0,nstep-inter0+1,' timesteps read'
      do 1 ind=nstep-inter0+1,nstep
      weight(ind)=1.d0
!                                                                       
!.....CALC MATRIX AA OF LEAST SQUARE FIT
!
         DO 5 LM1=1,LLMM                                                
         DO 5 LM2=1,LLMM                                                
          AA(LM1,LM2)=AA(LM1,LM2)+x(LM2,IND)*x(LM1,IND)*weight(ind)     
 5       CONTINUE                                                       
  1   CONTINUE                                                          
!                                                                       
!.....FACTOR MATRIX AA                                                  
!                                                                       
      CALL SSICO(AA,NATO*3,LLMM,KPVT,RCOND,WSPACE)                      
      WRITE(6,*) 'Condition of least squares matrix:',RCOND
!                            
!.....loop over atoms
!.....loop over spin
      do 99 jatom=1,nat
      do 99 ispin=1,jspin
!
      do 10 lm1=1,llmm
        b(lm1)=0.0D0
 10   continue
!
      do 77 k=nstep-inter0+1,nstep
         DO 9 LM1=1,LLMM                                                
         b(lm1)=b(lm1)+q(jatom,ispin,k)*x(lm1,k)*weight(k)
 9       continue
 77   CONTINUE
!
!******************************************************
!.....SOLVE  A.X = B                                                    
!                                                                       
      CALL SSISL(AA,NATO*3,LLMM,KPVT,B)                                 
!          WRITE(6,1540)                                                
          SIGMA=0.D0                                                    
      do 12 k=nstep-inter0+1,nstep
          TEST=0.D0                                                     
          DO 13 LM1=1,LLMM                                              
  13      TEST=TEST+b(LM1)*x(lm1,k)                                   
          IF(MOD(k,5).EQ.0) THEN                                       
            DIFVXC=q(jatom,ispin,k)-TEST                                
            WRITE(6,1550)k,q(jatom,ispin,k),TEST,DIFVXC         
            WRITE(6,*)x((jatom-1)*3+1,k),x((jatom-1)*3+2,k), &
                  x((jatom-1)*3+3,k)
 1550       FORMAT(I8,2F10.5,F15.5)                                     
          ENDIF                                                         
  12      SIGMA=SIGMA+(q(jatom,ispin,k)-TEST)**2        
          SIGMA=SQRT(SIGMA/inter0)                                      
          WRITE(6,*) 'ATOM',JATOM,'  SIGMA OF FIT FOR SPIN',ispin        &
         ,SIGMA                 
          TEST=0.D0                                                     
!.....prediction
          DO 113 LM1=1,LLMM                                             
  113     TEST=TEST+b(LM1)*x(lm1,nstep+1)     
          q(jatom,ispin,nstep+1)=test
!
          write(6,1549) jatom,ispin,test,(b(lm1),lm1=1,llmm)
 1549     format(/,'New charge for atom',i3,', spin:',i3,f10.5,/, &
          'Linear combination coefficients:',/,(3f20.10))
           WRITE(6,*)x((jatom-1)*3+1,nstep+1), &
                  x((jatom-1)*3+2,nstep+1), &
                  x((jatom-1)*3+3,nstep+1)
!
 99       continue

 999      continue
!      stop
      END
!
!                                                               
      SUBROUTINE SSICO(A,LDA,N,KPVT,RCOND,Z)                            
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER LDA,N,KPVT(1)                                             
      REAL*8 A(LDA,1),Z(1)                                              
      REAL*8 RCOND                                                      
!                                                                       
!     SSICO FACTORS A REAL SYMMETRIC MATRIX BY ELIMINATION WITH         
!     SYMMETRIC PIVOTING AND ESTIMATES THE CONDITION OF THE MATRIX.     
!                                                                       
!     IF  RCOND  IS NOT NEEDED, SSIFA IS SLIGHTLY FASTER.               
!     TO SOLVE  A*X = B , FOLLOW SSICO BY SSISL.                        
!     TO COMPUTE  INVERSE(A)*C , FOLLOW SSICO BY SSISL.                 
!     TO COMPUTE  INVERSE(A) , FOLLOW SSICO BY SSIDI.                   
!     TO COMPUTE  DETERMINANT(A) , FOLLOW SSICO BY SSIDI.               
!     TO COMPUTE  INERTIA(A), FOLLOW SSICO BY SSIDI.                    
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       dble(LDA, N)                                           
!                THE SYMMETRIC MATRIX TO BE FACTORED.                   
!                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED.         
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  A .                
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!     OUTPUT                                                            
!                                                                       
!        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH      
!                WERE USED TO OBTAIN IT.                                
!                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)     
!                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT         
!                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE            
!                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL            
!                WITH 1 BY 1 AND 2 BY 2 BLOCKS.                         
!                                                                       
!        KPVT    INTEGER(N)                                             
!                AN INTEGER VECTOR OF PIVOT INDICES.                    
!                                                                       
!        RCOND   REAL                                                   
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .        
!                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS       
!                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE             
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . 
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION     
!                           1.0 + RCOND .EQ. 1.0                        
!                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING           
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF         
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE          
!                UNDERFLOWS.                                            
!                                                                       
!        Z       dble(N)                                                
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.  
!                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS      
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT           
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .                    
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     LINPACK SSIFA                                                     
!     BLAS SAXPY,SDOT,SSCAL,SASUM                                       
!     FORTRAN ABS,MAX  ,IABS,SIGN                                       
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      REAL*8 AK,AKM1,BK,BKM1,SDOT,DENOM,EK,T                            
      REAL*8 ANORM,S,SASUM,YNORM                                        
      INTEGER I,INFO,J,JM1,K,KP,KPS,KS                                  
!                                                                       
!                                                                       
!     FIND NORM OF A USING ONLY UPPER HALF                              
!                                                                       
      DO 30 J = 1, N                                                    
         Z(J) = SASUM(J,A(1,J),1)                                       
         JM1 = J - 1                                                    
         IF (JM1 .LT. 1) GO TO 20                                       
         DO 10 I = 1, JM1                                               
            Z(I) = Z(I) + ABS(A(I,J))                                   
   10    CONTINUE                                                       
   20    CONTINUE                                                       
   30 CONTINUE                                                          
      ANORM = 0.0D0                                                     
      DO 40 J = 1, N                                                    
         ANORM = MAX  (ANORM,Z(J))                                      
   40 CONTINUE                                                          
!                                                                       
!     FACTOR                                                            
!                                                                       
      CALL SSIFA(A,LDA,N,KPVT,INFO)                                     
!                                                                       
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .              
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .         
!     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL           
!     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .                   
!     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.            
!                                                                       
!     SOLVE U*D*W = E                                                   
!                                                                       
      EK = 1.0D0                                                        
      DO 50 J = 1, N                                                    
         Z(J) = 0.0D0                                                   
   50 CONTINUE                                                          
      K = N                                                             
   60 IF (K .EQ. 0) GO TO 120                                           
         KS = 1                                                         
         IF (KPVT(K) .LT. 0) KS = 2                                     
         KP = IABS(KPVT(K))                                             
         KPS = K + 1 - KS                                               
         IF (KP .EQ. KPS) GO TO 70                                      
            T = Z(KPS)                                                  
            Z(KPS) = Z(KP)                                              
            Z(KP) = T                                                   
   70    CONTINUE                                                       
         IF (Z(K) .NE. 0.0D0) EK = SIGN(EK,Z(K))                        
         Z(K) = Z(K) + EK                                               
         CALL SAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)                          
         IF (KS .EQ. 1) GO TO 80                                        
            IF (Z(K-1) .NE. 0.0D0) EK = SIGN(EK,Z(K-1))                 
            Z(K-1) = Z(K-1) + EK                                        
            CALL SAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)                   
   80    CONTINUE                                                       
         IF (KS .EQ. 2) GO TO 100                                       
            IF (ABS(Z(K)) .LE. ABS(A(K,K))) GO TO 90                    
               S = ABS(A(K,K))/ABS(Z(K))                                
               CALL SSCAL(N,S,Z,1)                                      
               EK = S*EK                                                
   90       CONTINUE                                                    
            IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)                   
            IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0                         
         GO TO 110                                                      
  100    CONTINUE                                                       
            AK = A(K,K)/A(K-1,K)                                        
            AKM1 = A(K-1,K-1)/A(K-1,K)                                  
            BK = Z(K)/A(K-1,K)                                          
            BKM1 = Z(K-1)/A(K-1,K)                                      
            DENOM = AK*AKM1 - 1.0D0                                     
            Z(K) = (AKM1*BK - BKM1)/DENOM                               
            Z(K-1) = (AK*BKM1 - BK)/DENOM                               
  110    CONTINUE                                                       
         K = K - KS                                                     
      GO TO 60                                                          
  120 CONTINUE                                                          
      S = 1.0D0/SASUM(N,Z,1)                                            
      CALL SSCAL(N,S,Z,1)                                               
!                                                                       
!     SOLVE TRANS(U)*Y = W                                              
!                                                                       
      K = 1                                                             
  130 IF (K .GT. N) GO TO 160                                           
         KS = 1                                                         
         IF (KPVT(K) .LT. 0) KS = 2                                     
         IF (K .EQ. 1) GO TO 150                                        
            Z(K) = Z(K) + SDOT(K-1,A(1,K),1,Z(1),1)                     
            IF (KS .EQ. 2)                                               &
               Z(K+1) = Z(K+1) + SDOT(K-1,A(1,K+1),1,Z(1),1)            
            KP = IABS(KPVT(K))                                          
            IF (KP .EQ. K) GO TO 140                                    
               T = Z(K)                                                 
               Z(K) = Z(KP)                                             
               Z(KP) = T                                                
  140       CONTINUE                                                    
  150    CONTINUE                                                       
         K = K + KS                                                     
      GO TO 130                                                         
  160 CONTINUE                                                          
      S = 1.0D0/SASUM(N,Z,1)                                            
      CALL SSCAL(N,S,Z,1)                                               
!                                                                       
      YNORM = 1.0D0                                                     
!                                                                       
!     SOLVE U*D*V = Y                                                   
!                                                                       
      K = N                                                             
  170 IF (K .EQ. 0) GO TO 230                                           
         KS = 1                                                         
         IF (KPVT(K) .LT. 0) KS = 2                                     
         IF (K .EQ. KS) GO TO 190                                       
            KP = IABS(KPVT(K))                                          
            KPS = K + 1 - KS                                            
            IF (KP .EQ. KPS) GO TO 180                                  
               T = Z(KPS)                                               
               Z(KPS) = Z(KP)                                           
               Z(KP) = T                                                
  180       CONTINUE                                                    
            CALL SAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)                       
            IF (KS .EQ. 2) CALL SAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)    
  190    CONTINUE                                                       
         IF (KS .EQ. 2) GO TO 210                                       
            IF (ABS(Z(K)) .LE. ABS(A(K,K))) GO TO 200                   
               S = ABS(A(K,K))/ABS(Z(K))                                
               CALL SSCAL(N,S,Z,1)                                      
               YNORM = S*YNORM                                          
  200       CONTINUE                                                    
            IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)                   
            IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0                         
         GO TO 220                                                      
  210    CONTINUE                                                       
            AK = A(K,K)/A(K-1,K)                                        
            AKM1 = A(K-1,K-1)/A(K-1,K)                                  
            BK = Z(K)/A(K-1,K)                                          
            BKM1 = Z(K-1)/A(K-1,K)                                      
            DENOM = AK*AKM1 - 1.0D0                                     
            Z(K) = (AKM1*BK - BKM1)/DENOM                               
            Z(K-1) = (AK*BKM1 - BK)/DENOM                               
  220    CONTINUE                                                       
         K = K - KS                                                     
      GO TO 170                                                         
  230 CONTINUE                                                          
      S = 1.0D0/SASUM(N,Z,1)                                            
      CALL SSCAL(N,S,Z,1)                                               
      YNORM = S*YNORM                                                   
!                                                                       
!     SOLVE TRANS(U)*Z = V                                              
!                                                                       
      K = 1                                                             
  240 IF (K .GT. N) GO TO 270                                           
         KS = 1                                                         
         IF (KPVT(K) .LT. 0) KS = 2                                     
         IF (K .EQ. 1) GO TO 260                                        
            Z(K) = Z(K) + SDOT(K-1,A(1,K),1,Z(1),1)                     
            IF (KS .EQ. 2)                                               &
               Z(K+1) = Z(K+1) + SDOT(K-1,A(1,K+1),1,Z(1),1)            
            KP = IABS(KPVT(K))                                          
            IF (KP .EQ. K) GO TO 250                                    
               T = Z(K)                                                 
               Z(K) = Z(KP)                                             
               Z(KP) = T                                                
  250       CONTINUE                                                    
  260    CONTINUE                                                       
         K = K + KS                                                     
      GO TO 240                                                         
  270 CONTINUE                                                          
!     MAKE ZNORM = 1.0                                                  
      S = 1.0D0/SASUM(N,Z,1)                                            
      CALL SSCAL(N,S,Z,1)                                               
      YNORM = S*YNORM                                                   
!                                                                       
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM                         
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0                               
      RETURN                                                            
      END                                                               
      SUBROUTINE SSISL(A,LDA,N,KPVT,B)                                  
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER LDA,N,KPVT(1)                                             
      REAL*8 A(LDA,1),B(1)                                              
!                                                                       
!     SSISL SOLVES THE REAL SYMMETRIC SYSTEM                            
!     A * X = B                                                         
!     USING THE FACTORS COMPUTED BY SSIFA.                              
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       dble(LDA,N)                                            
!                THE OUTPUT FROM SSIFA.                                 
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  A .                
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!        KPVT    INTEGER(N)                                             
!                THE PIVOT VECTOR FROM SSIFA.                           
!                                                                       
!        B       dble(N)                                                
!                THE RIGHT HAND SIDE VECTOR.                            
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        B       THE SOLUTION VECTOR  X .                               
!                                                                       
!     ERROR CONDITION                                                   
!                                                                       
!        A DIVISION BY ZERO MAY OCCUR IF  SSICO  HAS SET RCOND .EQ. 0.0 
!        OR  SSIFA  HAS SET INFO .NE. 0  .                              
!                                                                       
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 
!     WITH  P  COLUMNS                                                  
!           CALL SSIFA(A,LDA,N,KPVT,INFO)                               
!           IF (INFO .NE. 0) GO TO ...                                  
!           DO 10 J = 1, P                                              
!              CALL SSISL(A,LDA,N,KPVT,C(1,J))                          
!        10 CONTINUE                                                    
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.            
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS SAXPY,SDOT                                                   
!     FORTRAN IABS                                                      
!                                                                       
!     INTERNAL VARIABLES.                                               
!                                                                       
      REAL*8 AK,AKM1,BK,BKM1,SDOT,DENOM,TEMP                            
      INTEGER K,KP                                                      
!                                                                       
!     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND                    
!     D INVERSE TO B.                                                   
!                                                                       
      K = N                                                             
   10 IF (K .EQ. 0) GO TO 80                                            
         IF (KPVT(K) .LT. 0) GO TO 40                                   
!                                                                       
!           1 X 1 PIVOT BLOCK.                                          
!                                                                       
            IF (K .EQ. 1) GO TO 30                                      
               KP = KPVT(K)                                             
               IF (KP .EQ. K) GO TO 20                                  
!                                                                       
!                 INTERCHANGE.                                          
!                                                                       
                  TEMP = B(K)                                           
                  B(K) = B(KP)                                          
                  B(KP) = TEMP                                          
   20          CONTINUE                                                 
!                                                                       
!              APPLY THE TRANSFORMATION.                                
!                                                                       
               CALL SAXPY(K-1,B(K),A(1,K),1,B(1),1)                     
   30       CONTINUE                                                    
!                                                                       
!           APPLY D INVERSE.                                            
!                                                                       
            B(K) = B(K)/A(K,K)                                          
            K = K - 1                                                   
         GO TO 70                                                       
   40    CONTINUE                                                       
!                                                                       
!           2 X 2 PIVOT BLOCK.                                          
!                                                                       
            IF (K .EQ. 2) GO TO 60                                      
               KP = IABS(KPVT(K))                                       
               IF (KP .EQ. K - 1) GO TO 50                              
!                                                                       
!                 INTERCHANGE.                                          
!                                                                       
                  TEMP = B(K-1)                                         
                  B(K-1) = B(KP)                                        
                  B(KP) = TEMP                                          
   50          CONTINUE                                                 
!                                                                       
!              APPLY THE TRANSFORMATION.                                
!                                                                       
               CALL SAXPY(K-2,B(K),A(1,K),1,B(1),1)                     
               CALL SAXPY(K-2,B(K-1),A(1,K-1),1,B(1),1)                 
   60       CONTINUE                                                    
!                                                                       
!           APPLY D INVERSE.                                            
!                                                                       
            AK = A(K,K)/A(K-1,K)                                        
            AKM1 = A(K-1,K-1)/A(K-1,K)                                  
            BK = B(K)/A(K-1,K)                                          
            BKM1 = B(K-1)/A(K-1,K)                                      
            DENOM = AK*AKM1 - 1.0D0                                     
            B(K) = (AKM1*BK - BKM1)/DENOM                               
            B(K-1) = (AK*BKM1 - BK)/DENOM                               
            K = K - 2                                                   
   70    CONTINUE                                                       
      GO TO 10                                                          
   80 CONTINUE                                                          
!                                                                       
!     LOOP FORWARD APPLYING THE TRANSFORMATIONS.                        
!                                                                       
      K = 1                                                             
   90 IF (K .GT. N) GO TO 160                                           
         IF (KPVT(K) .LT. 0) GO TO 120                                  
!                                                                       
!           1 X 1 PIVOT BLOCK.                                          
!                                                                       
            IF (K .EQ. 1) GO TO 110                                     
!                                                                       
!              APPLY THE TRANSFORMATION.                                
!                                                                       
               B(K) = B(K) + SDOT(K-1,A(1,K),1,B(1),1)                  
               KP = KPVT(K)                                             
               IF (KP .EQ. K) GO TO 100                                 
!                                                                       
!                 INTERCHANGE.                                          
!                                                                       
                  TEMP = B(K)                                           
                  B(K) = B(KP)                                          
                  B(KP) = TEMP                                          
  100          CONTINUE                                                 
  110       CONTINUE                                                    
            K = K + 1                                                   
         GO TO 150                                                      
  120    CONTINUE                                                       
!                                                                       
!           2 X 2 PIVOT BLOCK.                                          
!                                                                       
            IF (K .EQ. 1) GO TO 140                                     
!                                                                       
!              APPLY THE TRANSFORMATION.                                
!                                                                       
               B(K) = B(K) + SDOT(K-1,A(1,K),1,B(1),1)                  
               B(K+1) = B(K+1) + SDOT(K-1,A(1,K+1),1,B(1),1)            
               KP = IABS(KPVT(K))                                       
               IF (KP .EQ. K) GO TO 130                                 
!                                                                       
!                 INTERCHANGE.                                          
!                                                                       
                  TEMP = B(K)                                           
                  B(K) = B(KP)                                          
                  B(KP) = TEMP                                          
  130          CONTINUE                                                 
  140       CONTINUE                                                    
            K = K + 2                                                   
  150    CONTINUE                                                       
      GO TO 90                                                          
  160 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE CSPFA(AP,N,KPVT,INFO)                                  
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N,KPVT(1),INFO                                            
      COMPLEX*16 AP(1)                                                  
!                                                                       
!     CSPFA FACTORS A COMPLEX SYMMETRIC MATRIX STORED IN                
!     PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING.               
!                                                                       
!     TO SOLVE  A*X = B , FOLLOW CSPFA BY CSPSL.                        
!     TO COMPUTE  INVERSE(A)*C , FOLLOW CSPFA BY CSPSL.                 
!     TO COMPUTE  DETERMINANT(A) , FOLLOW CSPFA BY CSPDI.               
!     TO COMPUTE  INVERSE(A) , FOLLOW CSPFA BY CSPDI.                   
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        AP      COMPLEX (N*(N+1)/2)                                    
!                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE        
!                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY  
!                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .      
!                SEE COMMENTS BELOW FOR DETAILS.                        
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH      
!                WERE USED TO OBTAIN IT STORED IN PACKED FORM.          
!                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)     
!                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT         
!                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE            
!                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL            
!                WITH 1 BY 1 AND 2 BY 2 BLOCKS.                         
!                                                                       
!        KPVT    INTEGER(N)                                             
!                AN INTEGER VECTOR OF PIVOT INDICES.                    
!                                                                       
!        INFO    INTEGER                                                
!                = 0  NORMAL VALUE.                                     
!                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS      
!                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,       
!                     BUT IT DOES INDICATE THAT CSPSL OR CSPDI MAY      
!                     DIVIDE BY ZERO IF CALLED.                         
!                                                                       
!     PACKED STORAGE                                                    
!                                                                       
!          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER            
!          TRIANGLE OF A SYMMETRIC MATRIX.                              
!                                                                       
!                K = 0                                                  
!                DO 20 J = 1, N                                         
!                   DO 10 I = 1, J                                      
!                      K = K + 1                                        
!                      AP(K)  = A(I,J)                                  
!             10    CONTINUE                                            
!             20 CONTINUE                                               
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.            
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS CAXPY,CSWAP,ICAMAX                                           
!     FORTRAN ABS,DIMAG,MAX  ,REAL,SQRT                                 
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      COMPLEX*16 AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T                    
      REAL*8 ABSAKK,ALPHA,COLMAX,ROWMAX                                 
      INTEGER ICAMAX,IJ,IJJ,IK,IKM1,IM,IMAX,IMAXP1,IMIM,IMJ,IMK         
      INTEGER J,JJ,JK,JKM1,JMAX,JMIM,K,KK,KM1,KM1K,KM1KM1,KM2,KSTEP     
      LOGICAL SWAP                                                      
!                                                                       
      COMPLEX*16 ZDUM                                                   
      REAL*8 CABS1                                                      
      CABS1(ZDUM) = ABS(dble(ZDUM)) + ABS(DIMAG(ZDUM))                  
!                                                                       
!     INITIALIZE                                                        
!                                                                       
!     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.                       
      ALPHA = (1.0D0 + SQRT(17.0D0))/8.0D0                              
!                                                                       
      INFO = 0                                                          
!                                                                       
!     MAIN LOOP ON K, WHICH GOES FROM N TO 1.                           
!                                                                       
      K = N                                                             
      IK = (N*(N - 1))/2                                                
   10 CONTINUE                                                          
!                                                                       
!        LEAVE THE LOOP IF K=0 OR K=1.                                  
!                                                                       
!     ...EXIT                                                           
         IF (K .EQ. 0) GO TO 200                                        
         IF (K .GT. 1) GO TO 20                                         
            KPVT(1) = 1                                                 
            IF (CABS1(AP(1)) .EQ. 0.0D0) INFO = 1                       
!     ......EXIT                                                        
            GO TO 200                                                   
   20    CONTINUE                                                       
!                                                                       
!        THIS SECTION OF CODE DETERMINES THE KIND OF                    
!        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,            
!        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND          
!        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS                
!        REQUIRED.                                                      
!                                                                       
         KM1 = K - 1                                                    
         KK = IK + K                                                    
         ABSAKK = CABS1(AP(KK))                                         
!                                                                       
!        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN                  
!        COLUMN K.                                                      
!                                                                       
         IMAX = ICAMAX(K-1,AP(IK+1),1)                                  
         IMK = IK + IMAX                                                
         COLMAX = CABS1(AP(IMK))                                        
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30                         
            KSTEP = 1                                                   
            SWAP = .FALSE.                                              
         GO TO 90                                                       
   30    CONTINUE                                                       
!                                                                       
!           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN               
!           ROW IMAX.                                                   
!                                                                       
            ROWMAX = 0.0D0                                              
            IMAXP1 = IMAX + 1                                           
            IM = IMAX*(IMAX - 1)/2                                      
            IMJ = IM + 2*IMAX                                           
            DO 40 J = IMAXP1, K                                         
               ROWMAX = MAX  (ROWMAX,CABS1(AP(IMJ)))                    
               IMJ = IMJ + J                                            
   40       CONTINUE                                                    
            IF (IMAX .EQ. 1) GO TO 50                                   
               JMAX = ICAMAX(IMAX-1,AP(IM+1),1)                         
               JMIM = JMAX + IM                                         
               ROWMAX = MAX  (ROWMAX,CABS1(AP(JMIM)))                   
   50       CONTINUE                                                    
            IMIM = IMAX + IM                                            
            IF (CABS1(AP(IMIM)) .LT. ALPHA*ROWMAX) GO TO 60             
               KSTEP = 1                                                
               SWAP = .TRUE.                                            
            GO TO 80                                                    
   60       CONTINUE                                                    
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70      
               KSTEP = 1                                                
               SWAP = .FALSE.                                           
            GO TO 80                                                    
   70       CONTINUE                                                    
               KSTEP = 2                                                
               SWAP = IMAX .NE. KM1                                     
   80       CONTINUE                                                    
   90    CONTINUE                                                       
         IF (MAX  (ABSAKK,COLMAX) .NE. 0.0D0) GO TO 100                 
!                                                                       
!           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.           
!                                                                       
            KPVT(K) = K                                                 
            INFO = K                                                    
         GO TO 190                                                      
  100    CONTINUE                                                       
         IF (KSTEP .EQ. 2) GO TO 140                                    
!                                                                       
!           1 X 1 PIVOT BLOCK.                                          
!                                                                       
            IF (.NOT.SWAP) GO TO 120                                    
!                                                                       
!              PERFORM AN INTERCHANGE.                                  
!                                                                       
               CALL CSWAP(IMAX,AP(IM+1),1,AP(IK+1),1)                   
               IMJ = IK + IMAX                                          
               DO 110 JJ = IMAX, K                                      
                  J = K + IMAX - JJ                                     
                  JK = IK + J                                           
                  T = AP(JK)                                            
                  AP(JK) = AP(IMJ)                                      
                  AP(IMJ) = T                                           
                  IMJ = IMJ - (J - 1)                                   
  110          CONTINUE                                                 
  120       CONTINUE                                                    
!                                                                       
!           PERFORM THE ELIMINATION.                                    
!                                                                       
            IJ = IK - (K - 1)                                           
            DO 130 JJ = 1, KM1                                          
               J = K - JJ                                               
               JK = IK + J                                              
               MULK = -AP(JK)/AP(KK)                                    
               T = MULK                                                 
               CALL CAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)                    
               IJJ = IJ + J                                             
               AP(JK) = MULK                                            
               IJ = IJ - (J - 1)                                        
  130       CONTINUE                                                    
!                                                                       
!           SET THE PIVOT ARRAY.                                        
!                                                                       
            KPVT(K) = K                                                 
            IF (SWAP) KPVT(K) = IMAX                                    
         GO TO 190                                                      
  140    CONTINUE                                                       
!                                                                       
!           2 X 2 PIVOT BLOCK.                                          
!                                                                       
            KM1K = IK + K - 1                                           
            IKM1 = IK - (K - 1)                                         
            IF (.NOT.SWAP) GO TO 160                                    
!                                                                       
!              PERFORM AN INTERCHANGE.                                  
!                                                                       
               CALL CSWAP(IMAX,AP(IM+1),1,AP(IKM1+1),1)                 
               IMJ = IKM1 + IMAX                                        
               DO 150 JJ = IMAX, KM1                                    
                  J = KM1 + IMAX - JJ                                   
                  JKM1 = IKM1 + J                                       
                  T = AP(JKM1)                                          
                  AP(JKM1) = AP(IMJ)                                    
                  AP(IMJ) = T                                           
                  IMJ = IMJ - (J - 1)                                   
  150          CONTINUE                                                 
               T = AP(KM1K)                                             
               AP(KM1K) = AP(IMK)                                       
               AP(IMK) = T                                              
  160       CONTINUE                                                    
!                                                                       
!           PERFORM THE ELIMINATION.                                    
!                                                                       
            KM2 = K - 2                                                 
            IF (KM2 .EQ. 0) GO TO 180                                   
               AK = AP(KK)/AP(KM1K)                                     
               KM1KM1 = IKM1 + K - 1                                    
               AKM1 = AP(KM1KM1)/AP(KM1K)                               
               DENOM = 1.0D0 - AK*AKM1                                  
               IJ = IK - (K - 1) - (K - 2)                              
               DO 170 JJ = 1, KM2                                       
                  J = KM1 - JJ                                          
                  JK = IK + J                                           
                  BK = AP(JK)/AP(KM1K)                                  
                  JKM1 = IKM1 + J                                       
                  BKM1 = AP(JKM1)/AP(KM1K)                              
                  MULK = (AKM1*BK - BKM1)/DENOM                         
                  MULKM1 = (AK*BKM1 - BK)/DENOM                         
                  T = MULK                                              
                  CALL CAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)                 
                  T = MULKM1                                            
                  CALL CAXPY(J,T,AP(IKM1+1),1,AP(IJ+1),1)               
                  AP(JK) = MULK                                         
                  AP(JKM1) = MULKM1                                     
                  IJJ = IJ + J                                          
                  IJ = IJ - (J - 1)                                     
  170          CONTINUE                                                 
  180       CONTINUE                                                    
!                                                                       
!           SET THE PIVOT ARRAY.                                        
!                                                                       
            KPVT(K) = 1 - K                                             
            IF (SWAP) KPVT(K) = -IMAX                                   
            KPVT(K-1) = KPVT(K)                                         
  190    CONTINUE                                                       
         IK = IK - (K - 1)                                              
         IF (KSTEP .EQ. 2) IK = IK - (K - 2)                            
         K = K - KSTEP                                                  
      GO TO 10                                                          
  200 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE CSPSL(AP,N,KPVT,B)                                     
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N,KPVT(1)                                                 
      COMPLEX*16 AP(1),B(1)                                             
!                                                                       
!     CSISL SOLVES THE COMPLEX SYMMETRIC SYSTEM                         
!     A * X = B                                                         
!     USING THE FACTORS COMPUTED BY CSPFA.                              
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        AP      COMPLEX(N*(N+1)/2)                                     
!                THE OUTPUT FROM CSPFA.                                 
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!        KPVT    INTEGER(N)                                             
!                THE PIVOT VECTOR FROM CSPFA.                           
!                                                                       
!        B       COMPLEX(N)                                             
!                THE RIGHT HAND SIDE VECTOR.                            
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        B       THE SOLUTION VECTOR  X .                               
!                                                                       
!     ERROR CONDITION                                                   
!                                                                       
!        A DIVISION BY ZERO MAY OCCUR IF  CSPCO  HAS SET RCOND .EQ. 0.0 
!        OR  CSPFA  HAS SET INFO .NE. 0  .                              
!                                                                       
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 
!     WITH  P  COLUMNS                                                  
!           CALL CSPFA(AP,N,KPVT,INFO)                                  
!           IF (INFO .NE. 0) GO TO ...                                  
!           DO 10 J = 1, P                                              
!              CALL CSPSL(AP,N,KPVT,C(1,J))                             
!        10 CONTINUE                                                    
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.            
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS CAXPY,CDOTU                                                  
!     FORTRAN IABS                                                      
!                                                                       
!     INTERNAL VARIABLES.                                               
!                                                                       
      COMPLEX*16 AK,AKM1,BK,BKM1,CDOTU,DENOM,TEMP                       
      INTEGER IK,IKM1,IKP1,K,KK,KM1K,KM1KM1,KP                          
!                                                                       
!     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND                    
!     D INVERSE TO B.                                                   
!                                                                       
      K = N                                                             
      IK = (N*(N - 1))/2                                                
   10 IF (K .EQ. 0) GO TO 80                                            
         KK = IK + K                                                    
         IF (KPVT(K) .LT. 0) GO TO 40                                   
!                                                                       
!           1 X 1 PIVOT BLOCK.                                          
!                                                                       
            IF (K .EQ. 1) GO TO 30                                      
               KP = KPVT(K)                                             
               IF (KP .EQ. K) GO TO 20                                  
!                                                                       
!                 INTERCHANGE.                                          
!                                                                       
                  TEMP = B(K)                                           
                  B(K) = B(KP)                                          
                  B(KP) = TEMP                                          
   20          CONTINUE                                                 
!                                                                       
!              APPLY THE TRANSFORMATION.                                
!                                                                       
               CALL CAXPY(K-1,B(K),AP(IK+1),1,B(1),1)                   
   30       CONTINUE                                                    
!                                                                       
!           APPLY D INVERSE.                                            
!                                                                       
            B(K) = B(K)/AP(KK)                                          
            K = K - 1                                                   
            IK = IK - K                                                 
         GO TO 70                                                       
   40    CONTINUE                                                       
!                                                                       
!           2 X 2 PIVOT BLOCK.                                          
!                                                                       
            IKM1 = IK - (K - 1)                                         
            IF (K .EQ. 2) GO TO 60                                      
               KP = IABS(KPVT(K))                                       
               IF (KP .EQ. K - 1) GO TO 50                              
!                                                                       
!                 INTERCHANGE.                                          
!                                                                       
                  TEMP = B(K-1)                                         
                  B(K-1) = B(KP)                                        
                  B(KP) = TEMP                                          
   50          CONTINUE                                                 
!                                                                       
!              APPLY THE TRANSFORMATION.                                
!                                                                       
               CALL CAXPY(K-2,B(K),AP(IK+1),1,B(1),1)                   
               CALL CAXPY(K-2,B(K-1),AP(IKM1+1),1,B(1),1)               
   60       CONTINUE                                                    
!                                                                       
!           APPLY D INVERSE.                                            
!                                                                       
            KM1K = IK + K - 1                                           
            KK = IK + K                                                 
            AK = AP(KK)/AP(KM1K)                                        
            KM1KM1 = IKM1 + K - 1                                       
            AKM1 = AP(KM1KM1)/AP(KM1K)                                  
            BK = B(K)/AP(KM1K)                                          
            BKM1 = B(K-1)/AP(KM1K)                                      
            DENOM = AK*AKM1 - 1.0D0                                     
            B(K) = (AKM1*BK - BKM1)/DENOM                               
            B(K-1) = (AK*BKM1 - BK)/DENOM                               
            K = K - 2                                                   
            IK = IK - (K + 1) - K                                       
   70    CONTINUE                                                       
      GO TO 10                                                          
   80 CONTINUE                                                          
!                                                                       
!     LOOP FORWARD APPLYING THE TRANSFORMATIONS.                        
!                                                                       
      K = 1                                                             
      IK = 0                                                            
   90 IF (K .GT. N) GO TO 160                                           
         IF (KPVT(K) .LT. 0) GO TO 120                                  
!                                                                       
!           1 X 1 PIVOT BLOCK.                                          
!                                                                       
            IF (K .EQ. 1) GO TO 110                                     
!                                                                       
!              APPLY THE TRANSFORMATION.                                
!                                                                       
               B(K) = B(K) + CDOTU(K-1,AP(IK+1),1,B(1),1)               
               KP = KPVT(K)                                             
               IF (KP .EQ. K) GO TO 100                                 
!                                                                       
!                 INTERCHANGE.                                          
!                                                                       
                  TEMP = B(K)                                           
                  B(K) = B(KP)                                          
                  B(KP) = TEMP                                          
  100          CONTINUE                                                 
  110       CONTINUE                                                    
            IK = IK + K                                                 
            K = K + 1                                                   
         GO TO 150                                                      
  120    CONTINUE                                                       
!                                                                       
!           2 X 2 PIVOT BLOCK.                                          
!                                                                       
            IF (K .EQ. 1) GO TO 140                                     
!                                                                       
!              APPLY THE TRANSFORMATION.                                
!                                                                       
               B(K) = B(K) + CDOTU(K-1,AP(IK+1),1,B(1),1)               
               IKP1 = IK + K                                            
               B(K+1) = B(K+1) + CDOTU(K-1,AP(IKP1+1),1,B(1),1)         
               KP = IABS(KPVT(K))                                       
               IF (KP .EQ. K) GO TO 130                                 
!                                                                       
!                 INTERCHANGE.                                          
!                                                                       
                  TEMP = B(K)                                           
                  B(K) = B(KP)                                          
                  B(KP) = TEMP                                          
  130          CONTINUE                                                 
  140       CONTINUE                                                    
            IK = IK + K + K + 1                                         
            K = K + 2                                                   
  150    CONTINUE                                                       
      GO TO 90                                                          
  160 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)                            
!                                                                       
!     CONSTANT TIMES A VECTOR PLUS A VECTOR.                            
!     USES UNROLLED LOOP FOR INCREMENTS EQUAL TO ONE.                   
!     JACK DONGARRA, LINPACK, 3/11/78.                                  
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 SX(1),SY(1),SA                                             
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N                                 
!                                                                       
      IF(N.LE.0)RETURN                                                  
      IF (SA .EQ. 0.0) RETURN                                           
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               
!                                                                       
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                
!          NOT EQUAL TO 1                                               
!                                                                       
      IX = 1                                                            
      IY = 1                                                            
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 
      DO 10 I = 1,N                                                     
        SY(IY) = SY(IY) + SA*SX(IX)                                     
        IX = IX + INCX                                                  
        IY = IY + INCY                                                  
   10 CONTINUE                                                          
      RETURN                                                            
!                                                                       
!        CODE FOR BOTH INCREMENTS EQUAL TO 1                            
!                                                                       
!                                                                       
!        CLEAN-UP LOOP                                                  
!                                                                       
   20 M = MOD(N,4)                                                      
      IF( M .EQ. 0 ) GO TO 40                                           
      DO 30 I = 1,M                                                     
        SY(I) = SY(I) + SA*SX(I)                                        
   30 CONTINUE                                                          
      IF( N .LT. 4 ) RETURN                                             
   40 MP1 = M + 1                                                       
      DO 50 I = MP1,N,4                                                 
        SY(I) = SY(I) + SA*SX(I)                                        
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)                            
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)                            
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)                            
   50 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      FUNCTION SDOT(N,SX,INCX,SY,INCY)  
!                                                                       
!     FORMS THE DOT PRODUCT OF TWO VECTORS.                             
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.                  
!     JACK DONGARRA, LINPACK, 3/11/78.                                  
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 SX(1),SY(1),STEMP                                          
      REAL*8 SDOT
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N                                 
!                                                                       
      STEMP = 0.0D0                                                     
      SDOT = 0.0D0                                                      
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               
!                                                                       
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                
!          NOT EQUAL TO 1                                               
!                                                                       
      IX = 1                                                            
      IY = 1                                                            
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 
      DO 10 I = 1,N                                                     
        STEMP = STEMP + SX(IX)*SY(IY)                                   
        IX = IX + INCX                                                  
        IY = IY + INCY                                                  
   10 CONTINUE                                                          
      SDOT = STEMP                                                      
      RETURN                                                            
!                                                                       
!        CODE FOR BOTH INCREMENTS EQUAL TO 1                            
!                                                                       
!                                                                       
!        CLEAN-UP LOOP                                                  
!                                                                       
   20 M = MOD(N,5)                                                      
      IF( M .EQ. 0 ) GO TO 40                                           
      DO 30 I = 1,M                                                     
        STEMP = STEMP + SX(I)*SY(I)                                     
   30 CONTINUE                                                          
      IF( N .LT. 5 ) GO TO 60                                           
   40 MP1 = M + 1                                                       
      DO 50 I = MP1,N,5                                                 
        STEMP = STEMP + SX(I)*SY(I) + SX(I + 1)*SY(I + 1) +              &
         SX(I + 2)*SY(I + 2) + SX(I + 3)*SY(I + 3) + SX(I + 4)*SY(I + 4)
   50 CONTINUE                                                          
   60 SDOT = STEMP                                                      
      RETURN                                                            
      END                                                               
      SUBROUTINE  SSCAL(N,SA,SX,INCX)                                   
!                                                                       
!     SCALES A VECTOR BY A CONSTANT.                                    
!     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO 1.                     
!     JACK DONGARRA, LINPACK, 3/11/78.                                  
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 SA,SX(1)                                                   
      INTEGER I,INCX,M,MP1,N,NINCX                                      
!                                                                       
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1)GO TO 20                                             
!                                                                       
!        CODE FOR INCREMENT NOT EQUAL TO 1                              
!                                                                       
      NINCX = N*INCX                                                    
      DO 10 I = 1,NINCX,INCX                                            
        SX(I) = SA*SX(I)                                                
   10 CONTINUE                                                          
      RETURN                                                            
!                                                                       
!        CODE FOR INCREMENT EQUAL TO 1                                  
!                                                                       
!                                                                       
!        CLEAN-UP LOOP                                                  
!                                                                       
   20 M = MOD(N,5)                                                      
      IF( M .EQ. 0 ) GO TO 40                                           
      DO 30 I = 1,M                                                     
        SX(I) = SA*SX(I)                                                
   30 CONTINUE                                                          
      IF( N .LT. 5 ) RETURN                                             
   40 MP1 = M + 1                                                       
      DO 50 I = MP1,N,5                                                 
        SX(I) = SA*SX(I)                                                
        SX(I + 1) = SA*SX(I + 1)                                        
        SX(I + 2) = SA*SX(I + 2)                                        
        SX(I + 3) = SA*SX(I + 3)                                        
        SX(I + 4) = SA*SX(I + 4)                                        
   50 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      FUNCTION SASUM(N,SX,INCX) 
!                                                                       
!     TAKES THE SUM OF THE ABSOLUTE VALUES.                             
!     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.                   
!     JACK DONGARRA, LINPACK, 3/11/78.                                  
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 SX(1),STEMP                                                
      REAL*8 SASUM
      INTEGER I,INCX,M,MP1,N,NINCX                                      
!                                                                       
      SASUM = 0.0D0                                                     
      STEMP = 0.0D0                                                     
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1)GO TO 20                                             
!                                                                       
!        CODE FOR INCREMENT NOT EQUAL TO 1                              
!                                                                       
      NINCX = N*INCX                                                    
      DO 10 I = 1,NINCX,INCX                                            
        STEMP = STEMP + ABS(SX(I))                                      
   10 CONTINUE                                                          
      SASUM = STEMP                                                     
      RETURN                                                            
!                                                                       
!        CODE FOR INCREMENT EQUAL TO 1                                  
!                                                                       
!                                                                       
!        CLEAN-UP LOOP                                                  
!                                                                       
   20 M = MOD(N,6)                                                      
      IF( M .EQ. 0 ) GO TO 40                                           
      DO 30 I = 1,M                                                     
        STEMP = STEMP + ABS(SX(I))                                      
   30 CONTINUE                                                          
      IF( N .LT. 6 ) GO TO 60                                           
   40 MP1 = M + 1                                                       
      DO 50 I = MP1,N,6                                                 
        STEMP = STEMP + ABS(SX(I)) + ABS(SX(I + 1)) + ABS(SX(I + 2))     &
        + ABS(SX(I + 3)) + ABS(SX(I + 4)) + ABS(SX(I + 5))              
   50 CONTINUE                                                          
   60 SASUM = STEMP                                                     
      RETURN                                                            
      END                                                               
      SUBROUTINE CAXPY(N,CA,CX,INCX,CY,INCY)                            
!                                                                       
!     CONSTANT TIMES A VECTOR PLUS A VECTOR.                            
!     JACK DONGARRA, LINPACK, 3/11/78.                                  
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 CX(1),CY(1),CA                                         
      INTEGER I,INCX,INCY,IX,IY,N                                       
!                                                                       
      IF(N.LE.0)RETURN                                                  
      IF (ABS(dble(CA)) + ABS(DIMAG(CA)) .EQ. 0.0 ) RETURN              
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               
!                                                                       
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                
!          NOT EQUAL TO 1                                               
!                                                                       
      IX = 1                                                            
      IY = 1                                                            
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 
      DO 10 I = 1,N                                                     
        CY(IY) = CY(IY) + CA*CX(IX)                                     
        IX = IX + INCX                                                  
        IY = IY + INCY                                                  
   10 CONTINUE                                                          
      RETURN                                                            
!                                                                       
!        CODE FOR BOTH INCREMENTS EQUAL TO 1                            
!                                                                       
   20 DO 30 I = 1,N                                                     
        CY(I) = CY(I) + CA*CX(I)                                        
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE  CSWAP (N,CX,INCX,CY,INCY)                             
!                                                                       
!     INTERCHANGES TWO VECTORS.                                         
!     JACK DONGARRA, LINPACK, 3/11/78.                                  
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 CX(1),CY(1),CTEMP                                      
      INTEGER I,INCX,INCY,IX,IY,N                                       
!                                                                       
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               
!                                                                       
!       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL       
!         TO 1                                                          
!                                                                       
      IX = 1                                                            
      IY = 1                                                            
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 
      DO 10 I = 1,N                                                     
        CTEMP = CX(IX)                                                  
        CX(IX) = CY(IY)                                                 
        CY(IY) = CTEMP                                                  
        IX = IX + INCX                                                  
        IY = IY + INCY                                                  
   10 CONTINUE                                                          
      RETURN                                                            
!                                                                       
!       CODE FOR BOTH INCREMENTS EQUAL TO 1                             
   20 DO 30 I = 1,N                                                     
        CTEMP = CX(I)                                                   
        CX(I) = CY(I)                                                   
        CY(I) = CTEMP                                                   
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      INTEGER FUNCTION ICAMAX(N,CX,INCX)                                
!                                                                       
!     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.            
!     JACK DONGARRA, LINPACK, 3/11/78.                                  
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 CX(1)                                                  
      REAL*8 SMAX                                                       
      INTEGER I,INCX,IX,N                                               
      COMPLEX*16 ZDUM                                                   
      REAL*8 CABS1                                                      
      CABS1(ZDUM) = ABS(dble(ZDUM)) + ABS(DIMAG(ZDUM))                  
!                                                                       
      ICAMAX = 0                                                        
      IF( N .LT. 1 ) RETURN                                             
      ICAMAX = 1                                                        
      IF(N.EQ.1)RETURN                                                  
      IF(INCX.EQ.1)GO TO 20                                             
!                                                                       
!        CODE FOR INCREMENT NOT EQUAL TO 1                              
!                                                                       
      IX = 1                                                            
      SMAX = CABS1(CX(1))                                               
      IX = IX + INCX                                                    
      DO 10 I = 2,N                                                     
         IF(CABS1(CX(IX)).LE.SMAX) GO TO 5                              
         ICAMAX = I                                                     
         SMAX = CABS1(CX(IX))                                           
    5    IX = IX + INCX                                                 
   10 CONTINUE                                                          
      RETURN                                                            
!                                                                       
!        CODE FOR INCREMENT EQUAL TO 1                                  
!                                                                       
   20 SMAX = CABS1(CX(1))                                               
      DO 30 I = 2,N                                                     
         IF(CABS1(CX(I)).LE.SMAX) GO TO 30                              
         ICAMAX = I                                                     
         SMAX = CABS1(CX(I))                                            
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      FUNCTION CDOTU(N,CX,INCX,CY,INCY)  
!                                                                       
!     FORMS THE DOT PRODUCT OF TWO VECTORS.                             
!     JACK DONGARRA, LINPACK, 3/11/78.                                  
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 CX(1),CY(1),CTEMP                                      
      COMPLEX*16 CDOTU
      INTEGER I,INCX,INCY,IX,IY,N                                       
!                                                                       
      CTEMP = (0.0D0,0.0D0)
      CDOTU = (0.0D0,0.0d0)
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               
!                                                                       
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                
!          NOT EQUAL TO 1                                               
!                                                                       
      IX = 1                                                            
      IY = 1                                                            
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 
      DO 10 I = 1,N                                                     
        CTEMP = CTEMP + CX(IX)*CY(IY)                                   
        IX = IX + INCX                                                  
        IY = IY + INCY                                                  
   10 CONTINUE                                                          
      CDOTU = CTEMP                                                     
      RETURN                                                            
!                                                                       
!        CODE FOR BOTH INCREMENTS EQUAL TO 1                            
!                                                                       
   20 DO 30 I = 1,N                                                     
        CTEMP = CTEMP + CX(I)*CY(I)                                     
   30 CONTINUE                                                          
      CDOTU = CTEMP                                                     
      RETURN                                                            
      END                                                               
      SUBROUTINE SSIFA(A,LDA,N,KPVT,INFO)                               
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER LDA,N,KPVT(1),INFO                                        
      REAL*8 A(LDA,1)                                                   
!                                                                       
!     SSIFA FACTORS A REAL SYMMETRIC MATRIX BY ELIMINATION              
!     WITH SYMMETRIC PIVOTING.                                          
!                                                                       
!     TO SOLVE  A*X = B , FOLLOW SSIFA BY SSISL.                        
!     TO COMPUTE  INVERSE(A)*C , FOLLOW SSIFA BY SSISL.                 
!     TO COMPUTE  DETERMINANT(A) , FOLLOW SSIFA BY SSIDI.               
!     TO COMPUTE  INERTIA(A) , FOLLOW SSIFA BY SSIDI.                   
!     TO COMPUTE  INVERSE(A) , FOLLOW SSIFA BY SSIDI.                   
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       dble(LDA,N)                                            
!                THE SYMMETRIC MATRIX TO BE FACTORED.                   
!                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED.         
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  A .                
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH      
!                WERE USED TO OBTAIN IT.                                
!                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)     
!                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT         
!                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE            
!                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL            
!                WITH 1 BY 1 AND 2 BY 2 BLOCKS.                         
!                                                                       
!        KPVT    INTEGER(N)                                             
!                AN INTEGER VECTOR OF PIVOT INDICES.                    
!                                                                       
!        INFO    INTEGER                                                
!                = 0  NORMAL VALUE.                                     
!                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS      
!                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,       
!                     BUT IT DOES INDICATE THAT SSISL OR SSIDI MAY      
!                     DIVIDE BY ZERO IF CALLED.                         
!                                                                       
!     LINPACK. THIS VERSION DATED 08/14/78 .                            
!     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.            
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS SAXPY,SSWAP,ISAMAX                                           
!     FORTRAN ABS,MAX  ,SQRT                                            
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      REAL*8 AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T                        
      REAL*8 ABSAKK,ALPHA,COLMAX,ROWMAX                                 
      INTEGER IMAX,IMAXP1,J,JJ,JMAX,K,KM1,KM2,KSTEP,ISAMAX              
      LOGICAL SWAP                                                      
!                                                                       
!                                                                       
!     INITIALIZE                                                        
!                                                                       
!     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.                       
      ALPHA = (1.0D0 + SQRT(17.0D0))/8.0D0                              
!                                                                       
      INFO = 0                                                          
!                                                                       
!     MAIN LOOP ON K, WHICH GOES FROM N TO 1.                           
!                                                                       
      K = N                                                             
   10 CONTINUE                                                          
!                                                                       
!        LEAVE THE LOOP IF K=0 OR K=1.                                  
!                                                                       
!     ...EXIT                                                           
         IF (K .EQ. 0) GO TO 200                                        
         IF (K .GT. 1) GO TO 20                                         
            KPVT(1) = 1                                                 
            IF (A(1,1) .EQ. 0.0D0) INFO = 1                             
!     ......EXIT                                                        
            GO TO 200                                                   
   20    CONTINUE                                                       
!                                                                       
!        THIS SECTION OF CODE DETERMINES THE KIND OF                    
!        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,            
!        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND          
!        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS                
!        REQUIRED.                                                      
!                                                                       
         KM1 = K - 1                                                    
         ABSAKK = ABS(A(K,K))                                           
!                                                                       
!        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN                  
!        COLUMN K.                                                      
!                                                                       
         IMAX = ISAMAX(K-1,A(1,K),1)                                    
         COLMAX = ABS(A(IMAX,K))                                        
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30                         
            KSTEP = 1                                                   
            SWAP = .FALSE.                                              
         GO TO 90                                                       
   30    CONTINUE                                                       
!                                                                       
!           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN               
!           ROW IMAX.                                                   
!                                                                       
            ROWMAX = 0.0D0                                              
            IMAXP1 = IMAX + 1                                           
            DO 40 J = IMAXP1, K                                         
               ROWMAX = MAX  (ROWMAX,ABS(A(IMAX,J)))                    
   40       CONTINUE                                                    
            IF (IMAX .EQ. 1) GO TO 50                                   
               JMAX = ISAMAX(IMAX-1,A(1,IMAX),1)                        
               ROWMAX = MAX  (ROWMAX,ABS(A(JMAX,IMAX)))                 
   50       CONTINUE                                                    
            IF (ABS(A(IMAX,IMAX)) .LT. ALPHA*ROWMAX) GO TO 60           
               KSTEP = 1                                                
               SWAP = .TRUE.                                            
            GO TO 80                                                    
   60       CONTINUE                                                    
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70      
               KSTEP = 1                                                
               SWAP = .FALSE.                                           
            GO TO 80                                                    
   70       CONTINUE                                                    
               KSTEP = 2                                                
               SWAP = IMAX .NE. KM1                                     
   80       CONTINUE                                                    
   90    CONTINUE                                                       
         IF (MAX  (ABSAKK,COLMAX) .NE. 0.0D0) GO TO 100                 
!                                                                       
!           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.           
!                                                                       
            KPVT(K) = K                                                 
            INFO = K                                                    
         GO TO 190                                                      
  100    CONTINUE                                                       
         IF (KSTEP .EQ. 2) GO TO 140                                    
!                                                                       
!           1 X 1 PIVOT BLOCK.                                          
!                                                                       
            IF (.NOT.SWAP) GO TO 120                                    
!                                                                       
!              PERFORM AN INTERCHANGE.                                  
!                                                                       
               CALL SSWAP(IMAX,A(1,IMAX),1,A(1,K),1)                    
               DO 110 JJ = IMAX, K                                      
                  J = K + IMAX - JJ                                     
                  T = A(J,K)                                            
                  A(J,K) = A(IMAX,J)                                    
                  A(IMAX,J) = T                                         
  110          CONTINUE                                                 
  120       CONTINUE                                                    
!                                                                       
!           PERFORM THE ELIMINATION.                                    
!                                                                       
            DO 130 JJ = 1, KM1                                          
               J = K - JJ                                               
               MULK = -A(J,K)/A(K,K)                                    
               T = MULK                                                 
               CALL SAXPY(J,T,A(1,K),1,A(1,J),1)                        
               A(J,K) = MULK                                            
  130       CONTINUE                                                    
!                                                                       
!           SET THE PIVOT ARRAY.                                        
!                                                                       
            KPVT(K) = K                                                 
            IF (SWAP) KPVT(K) = IMAX                                    
         GO TO 190                                                      
  140    CONTINUE                                                       
!                                                                       
!           2 X 2 PIVOT BLOCK.                                          
!                                                                       
            IF (.NOT.SWAP) GO TO 160                                    
!                                                                       
!              PERFORM AN INTERCHANGE.                                  
!                                                                       
               CALL SSWAP(IMAX,A(1,IMAX),1,A(1,K-1),1)                  
               DO 150 JJ = IMAX, KM1                                    
                  J = KM1 + IMAX - JJ                                   
                  T = A(J,K-1)                                          
                  A(J,K-1) = A(IMAX,J)                                  
                  A(IMAX,J) = T                                         
  150          CONTINUE                                                 
               T = A(K-1,K)                                             
               A(K-1,K) = A(IMAX,K)                                     
               A(IMAX,K) = T                                            
  160       CONTINUE                                                    
!                                                                       
!           PERFORM THE ELIMINATION.                                    
!                                                                       
            KM2 = K - 2                                                 
            IF (KM2 .EQ. 0) GO TO 180                                   
               AK = A(K,K)/A(K-1,K)                                     
               AKM1 = A(K-1,K-1)/A(K-1,K)                               
               DENOM = 1.0D0 - AK*AKM1                                  
               DO 170 JJ = 1, KM2                                       
                  J = KM1 - JJ                                          
                  BK = A(J,K)/A(K-1,K)                                  
                  BKM1 = A(J,K-1)/A(K-1,K)                              
                  MULK = (AKM1*BK - BKM1)/DENOM                         
                  MULKM1 = (AK*BKM1 - BK)/DENOM                         
                  T = MULK                                              
                  CALL SAXPY(J,T,A(1,K),1,A(1,J),1)                     
                  T = MULKM1                                            
                  CALL SAXPY(J,T,A(1,K-1),1,A(1,J),1)                   
                  A(J,K) = MULK                                         
                  A(J,K-1) = MULKM1                                     
  170          CONTINUE                                                 
  180       CONTINUE                                                    
!                                                                       
!           SET THE PIVOT ARRAY.                                        
!                                                                       
            KPVT(K) = 1 - K                                             
            IF (SWAP) KPVT(K) = -IMAX                                   
            KPVT(K-1) = KPVT(K)                                         
  190    CONTINUE                                                       
         K = K - KSTEP                                                  
      GO TO 10                                                          
  200 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE  SSWAP (N,SX,INCX,SY,INCY)                             
!                                                                       
!     INTERCHANGES TWO VECTORS.                                         
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO 1.                    
!     JACK DONGARRA, LINPACK, 3/11/78.                                  
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 SX(1),SY(1),STEMP                                          
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N                                 
!                                                                       
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               
!                                                                       
!       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL       
!         TO 1                                                          
!                                                                       
      IX = 1                                                            
      IY = 1                                                            
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 
      DO 10 I = 1,N                                                     
        STEMP = SX(IX)                                                  
        SX(IX) = SY(IY)                                                 
        SY(IY) = STEMP                                                  
        IX = IX + INCX                                                  
        IY = IY + INCY                                                  
   10 CONTINUE                                                          
      RETURN                                                            
!                                                                       
!       CODE FOR BOTH INCREMENTS EQUAL TO 1                             
!                                                                       
!                                                                       
!       CLEAN-UP LOOP                                                   
!                                                                       
   20 M = MOD(N,3)                                                      
      IF( M .EQ. 0 ) GO TO 40                                           
      DO 30 I = 1,M                                                     
        STEMP = SX(I)                                                   
        SX(I) = SY(I)                                                   
        SY(I) = STEMP                                                   
   30 CONTINUE                                                          
      IF( N .LT. 3 ) RETURN                                             
   40 MP1 = M + 1                                                       
      DO 50 I = MP1,N,3                                                 
        STEMP = SX(I)                                                   
        SX(I) = SY(I)                                                   
        SY(I) = STEMP                                                   
        STEMP = SX(I + 1)                                               
        SX(I + 1) = SY(I + 1)                                           
        SY(I + 1) = STEMP                                               
        STEMP = SX(I + 2)                                               
        SX(I + 2) = SY(I + 2)                                           
        SY(I + 2) = STEMP                                               
   50 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      INTEGER FUNCTION ISAMAX(N,SX,INCX)                                
!                                                                       
!     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.            
!     JACK DONGARRA, LINPACK, 3/11/78.                                  
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 SX(1),SMAX                                                 
      INTEGER I,INCX,IX,N                                               
!                                                                       
      ISAMAX = 0                                                        
      IF( N .LT. 1 ) RETURN                                             
      ISAMAX = 1                                                        
      IF(N.EQ.1)RETURN                                                  
      IF(INCX.EQ.1)GO TO 20                                             
!                                                                       
!        CODE FOR INCREMENT NOT EQUAL TO 1                              
!                                                                       
      IX = 1                                                            
      SMAX = ABS(SX(1))                                                 
      IX = IX + INCX                                                    
      DO 10 I = 2,N                                                     
         IF(ABS(SX(IX)).LE.SMAX) GO TO 5                                
         ISAMAX = I                                                     
         SMAX = ABS(SX(IX))                                             
    5    IX = IX + INCX                                                 
   10 CONTINUE                                                          
      RETURN                                                            
!                                                                       
!        CODE FOR INCREMENT EQUAL TO 1                                  
!                                                                       
   20 SMAX = ABS(SX(1))                                                 
      DO 30 I = 2,N                                                     
         IF(ABS(SX(I)).LE.SMAX) GO TO 30                                
         ISAMAX = I                                                     
         SMAX = ABS(SX(I))                                              
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
