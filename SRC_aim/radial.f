      SUBROUTINE RADIAL(CLM,RHO,R,IR,JATOM)                             
      use rad
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      DIMENSION CLM(NRAD)                                                
!      COMMON /RAD/ RM(NRAD,NATO),JRI(NATO)                                 
        ir1=ir   
        if(ir1.lt.2) ir1=2
        if(ir1.gt.jri(jatom)-2) ir1=jri(jatom)-2                           
      CALL INTERP(Rm(Ir1-1,jatom),clm(Ir1-1),4,R,SI,DSI,.FALSE.)      
      RHO=si/(R*R)                                                     
      RETURN                                                            
      END                                                               
