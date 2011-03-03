      SUBROUTINE GBASS(RBAS,GBAS)                                       
!     **                                                              **
!     **  CALCULATE RECIPROCAL LATTICE VECTORS FROM REAL SPACE        **
!     **  LATTICE VECTORS OR VICE VERSA (without 2 Pi factor !!!      **
!     **                                                              **
      INTEGER i,j
      DOUBLE PRECISION PI, DET
      DOUBLE PRECISION RBAS(3,3),GBAS(3,3) ,help(3,3)                                 
      PI=4.D0*DATAN(1.D0)                                               
      GBAS(1,1)=RBAS(2,2)*RBAS(3,3)-RBAS(3,2)*RBAS(2,3)                 
      GBAS(2,1)=RBAS(3,2)*RBAS(1,3)-RBAS(1,2)*RBAS(3,3)                 
      GBAS(3,1)=RBAS(1,2)*RBAS(2,3)-RBAS(2,2)*RBAS(1,3)                 
      GBAS(1,2)=RBAS(2,3)*RBAS(3,1)-RBAS(3,3)*RBAS(2,1)                 
      GBAS(2,2)=RBAS(3,3)*RBAS(1,1)-RBAS(1,3)*RBAS(3,1)                 
      GBAS(3,2)=RBAS(1,3)*RBAS(2,1)-RBAS(2,3)*RBAS(1,1)                 
      GBAS(1,3)=RBAS(2,1)*RBAS(3,2)-RBAS(3,1)*RBAS(2,2)                 
      GBAS(2,3)=RBAS(3,1)*RBAS(1,2)-RBAS(1,1)*RBAS(3,2)                 
      GBAS(3,3)=RBAS(1,1)*RBAS(2,2)-RBAS(2,1)*RBAS(1,2)                 
      DET=0.D0                                                          
      DO 100 I=1,3                                                      
      DET=DET+GBAS(I,1)*RBAS(I,1)                                       
100   CONTINUE                                                          
      DO 110 I=1,3                                                      
      DO 110 J=1,3     
 110     help(i,j)=gbas(i,j)/det                                                 
!      GBAS(I,J)=GBAS(I,J)*2.D0*PI/DET                                   
      DO 111 I=1,3                                                      
      DO 111 J=1,3     
      GBAS(I,J)=help(j,i)                                   
111   CONTINUE                                                          
      RETURN                                                            
      END                                                               
