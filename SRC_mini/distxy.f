      SUBROUTINE distxy(xx1,xx2,xx3,xx4,xx5,xx6,dists1)
!                                                                       
!                                                                       
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'param.inc'
!        Common blocks
!
!
!        LATTIC  - lattice type
!                  'P' ...... primitive latice (cubic, tetragonal,
!                             orthorhombic, monoclinic, triclin)
!                  'FC' ..... face centered
!                  'BC' ..... body centered
!                  'HEX' .... hexagonal
!                  'CXY' .... c-base centered (only orthorombic)
!                  'CYZ' .... a-base centered (only orthorombic)
!                  'CXZ' .... b-base centered (only orthorombic)
!                  'R' ...... rhombohedral
!
      CHARACTER*4        LATTIC
      COMMON  /CHAR/     LATTIC
      SAVE    /CHAR/
!
!        ORTH - .TRUE. for orthogonal lattice
!               .FALSE. otherwise
!
      LOGICAL            ORTHO
      COMMON  /ORTH/     ORTHO
      SAVE    /ORTH/
!
!        ALAT(1:3)   - lattice constants (a,b,c)
!        ALPHA(1:3)  - angles of the unit-cell's unit-vectors
!        IATNR(i)    - atom index of (inequivalent) atom i
!                      also indicates cubic and non-cubic symmetries
!                      IATNR(i) .GT. 0 ... cubic symmetry
!                      IATNR(i) .LT. 0 ... non-cubic symmetry
!        MULT(i)     - number of equivalent atoms for inequiv. atom i
!        PIA(1:3)    - reciprocal lattice constants (2*pi/a, 2*pi/b,
!                      2*pi/c)
!        POS(1:3,i)  - position vector of atom i (including equivalent
!                      atoms) in the unit-cell
!        RMT(i)      - muffin tin radius of atom i
!        V(i)        - relative muffin tin spherevolume for atom i
!        VI          - inverse volume of the direct lattice unit-cell
!
      DOUBLE PRECISION   VI
      DOUBLE PRECISION   ALAT(3), ALPHA(3), PIA(3)
      COMMON  /STRUK/    ALAT, ALPHA, PIA, VI
      SAVE    /STRUK/
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
      DIMENSION DIF(3),XX(3),PP(3),P(3),DISTS(NNN),NR(NNN),PNN(3,NNN)   
      DIMENSION NNAT(NNN),help(3),pos2(3,2)                                     
      DIMENSION BRnn(3,3)                                        
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!                                                                       
!     INDEX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE    
!     NONEQUIVALENT ATOMS                                               
!                                                                       
!                                                                       
!.....SET UP LATTICE, AND IF REQUIRED ROTATION MATRICES                 
      CALL DIRLAT (NAT,alpha,brnn)                                
!                                                                       
      pos2(1,1)=xx1/alat(1)
      pos2(2,1)=xx2/alat(2)
      pos2(3,1)=xx3/alat(3)
      pos2(1,2)=xx4/alat(1)
      pos2(2,2)=xx5/alat(2)
      pos2(3,2)=xx6/alat(3)
      INDEX=0                                                           
      DO 200 JATOM=1,2                                                
      INDEX=INDEX+1                                                     
      DO 150 J=1,3                                                      
  150 XX(J)=POS2(J,INDEX)                                                
      NC=0                                                              
          DO 180 I1=-2,2                                                
          DO 180 I2=-2,2                                                
          DO 180 I3=-2,2                                                
          IF(ortho) THEN                                   
            P(1)=I1*BRnn(1,1)+I2*BRnn(2,1)+I3*BRnn(3,1)                    
            P(2)=I1*BRnn(1,2)+I2*BRnn(2,2)+I3*BRnn(3,2)                    
            P(3)=I1*BRnn(1,3)+I2*BRnn(2,3)+I3*BRnn(3,3)                    
            ELSE                                                        
            P(1)=I1                                                     
            P(2)=I2                                                     
            P(3)=I3                                                     
          ENDIF                                                         
          K=0                                                           
      DO 120 JAT=1,2                                                  
      K=K+1                                                             
      DIST=0.d0                                                           
      DO 100 L=1,3                                                      
      PP(L)=POS2(L,K)+P(L)
  100 DIF(L)=XX(L)-PP(L)                                                
      IF (.not.ortho) THEN                                       
        help(1)=dif(1)  
        help(2)=dif(2)  
        help(3)=dif(3)  
      if(lattic(1:1).eq.'R') then
        dif(1)=help(1)*BRnn(1,1)+help(2)*BRnn(2,1)+help(3)*BRnn(3,1)
        dif(2)=help(1)*BRnn(1,2)+help(2)*BRnn(2,2)+help(3)*BRnn(3,2)
        dif(3)=help(1)*BRnn(1,3)+help(2)*BRnn(2,3)+help(3)*BRnn(3,3)           
      else
        dif(1)=(help(1)*BRnn(1,1)*ALAT(1)+help(2)*BRnn(2,1)*ALAT(2) &
               +help(3)*BRnn(3,1)*ALAT(3))/ALAT(1)           
        dif(2)=(help(1)*BRnn(1,2)*ALAT(1)+help(2)*BRnn(2,2)*ALAT(2) &
               +help(3)*BRnn(3,2)*ALAT(3))/ALAT(2)          
        dif(3)=(help(1)*BRnn(1,3)*ALAT(1)+help(2)*BRnn(2,3)*ALAT(2) &
               +help(3)*BRnn(3,3)*ALAT(3))/ALAT(3)           
      endif
      ENDIF                                                             
      DO 103 L=1,3                                                      
  103 DIST=DIST+DIF(L)*DIF(L)*ALAT(L)*ALAT(L)                                 
      DIST=SQRT(DIST)                                                   
      IF(DIST.GT.10.01d0) GO TO 120                                        
      IF(DIST.LT..000001) GO TO 120                                        
      NC=NC+1         
      if(nc.gt.nnn) goto 900
      DISTS(NC)=DIST                                                    
      NNAT(NC)=JAT                                                      
      DO 105 L=1,3                                                      
  105 PNN(L,NC)=PP(L)                                                   
!     WRITE(6,2)JAT,NAME(JAT),PP(1),PP(2),PP(3),DIST                    
!   2 FORMAT(' TO ATOM:',I2,2X,A10,' AT',3F8.4,                         
!    * ' IS ',F10.5,' A.U.')                                            
  120 CONTINUE                                                          
  180 CONTINUE                                                          
      CALL ORD2(DISTS,NR,NC)                                            
      N1=1                                                              
      N2=NR(N1)                                                         
      N3=NNAT(N2)                                                       
      dists1=DISTS(N1)
  200 CONTINUE                                                          
!                                                                       
                                                                        
      RETURN                                                            

!        Error messages
!
  900 CALL OUTERR('NN','nnn too small')
      STOP 'NN - Error'
 9000 FORMAT('RMT(',I2,')=',F7.5,' AND RMT(',I2,')=',F7.5)                
 9010 FORMAT('SUMS TO',F8.5,' GT NNN-DIST=',F8.5)                            
!
!        End of 'NN'
!
      END                                                               
!                                                                       
