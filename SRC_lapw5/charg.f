      SUBROUTINE CHARG(CHG,V,NDAT,POS,IATNR,Z,nat)
      use atomgrid                          
      use atpos                          
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      DIMENSION V(3),POS(3,*),VT(3),IATNR(NAT),Z(NAT)                      
      CHG=0.d0                                                            
      DO 30 IAT=1,NDAT                                                  
      JATOM=IABS(IATNR(IAT))                                            
      IF(Z(JATOM).LT.0.5) GOTO 30                                       
      DO 40 N=1,NPOS                                                    
      DO 35 II=1,3                                                      
   35 VT(II)=V(II)-(ATP(II,N)+POS(II,IAT))                              
      R=VNORM(VT)                                                       
      IF(R.LT.1.d-5) R=1.d-5                                            
      IR=1+LOG(R/RNOT(JATOM))/DX(JATOM)                                        
      IF(IR.GE.NRMAX(JATOM)) GOTO 40                                    
      IF(IR.LT.1) IR=1                                                  
      R1=Rn(IR,jatom)                                                         
      R2=Rn(IR+1,jatom)                                                       
      DR=R2-R1                                                          
      DDR=R-R1                                                          
      C1=RHO1(IR,JATOM)                                                 
      C2=RHO1(IR+1,JATOM)                                               
      DC=C2-C1                                                          
      CHG1=C1+(DDR/DR)*DC                                               
      CHG=CHG+CHG1    
   40 CONTINUE                                                          
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
