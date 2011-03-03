      SUBROUTINE RADIAL(CLM,RHO,R,IR,JATOM,rm,jri,nat)                      
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      DIMENSION CLM(NRAD)                                                
      dimension RM(NRAD,NAT),JRI(NAT)                                     
        ir1=ir   
        if(ir1.lt.2) ir1=2
        if(ir1.gt.jri(jatom)-2) ir1=jri(jatom)-2                           
      CALL INTERP(Rm(Ir1-1,jatom),clm(Ir1-1),4,R,SI,DSI,.FALSE.)      
      RHO=si/(R*R)                                                     
      RETURN                                                            
      END                                                               
