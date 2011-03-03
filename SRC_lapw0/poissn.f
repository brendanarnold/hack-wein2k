!@PROCESS DC(WORK,DENSIT)                                               
      SUBROUTINE POISSN (LM,LMMAX)
!                                                                       
!     SUBROUTINE POISSON CALCULATES THE CHARGE DENSITY OF A GIVEN       
!     POTENTIAL BY SOLVING THE POISSON EQUATION.                        
!     THE DERIVATION OF THE RADIAL POTENTIAL COEFFICIENTS IS DONE       
!     BY A LINEAR DIFFERENCE APPROXIMATION  TO THE FUNCTION VLM,        
!               VLM=VLM(I) :  1 < I < JRI(JATOM)                       
!     WHERE THE I-VALUES ARE EQUALLY SPACED. TO TRANSFORM THESE         
!     DERIVATIVES D(VLM)/DI, D2(VLM)/DI2 TO THE REQUIRED FORM           
!     D(VLM)/DR, D2(VLM)/DR2 THE FOLLOWING FORMULAS ARE USED :          
!                                                                       
!        DV     DV     1.                                               
!        --  =  -- *  ------                                            
!        DR     DI    R*DX                                              
!                                                                       
!        D2V    D2V    DV                                               
!        --- =  --- -  -- * DX                                          
!        DR2    DI2    DI                                               
!              -----------------                                        
!                R*R*DX*DX                                              
!                                                                       
!     ^ VLM                                                             
!     :                                                                 
!     :                                                                 
!     :        .  .                 .  .                                
!     :     .      .             .     .  .                             
!     :    .                   .       :  .    .                        
!     :  .          .        .         :  :         .                   
!     : .              .  .            :  :                .            
!     --------------------------------------------------------> I       
!                                      I1 I2                            
!              VLM = VLM( R(I) )                                        
!              R(I)= R0*EXP( DX*(I-1) )                                 
!                                                                       
!     THE RESULT  IS  OF THE FORM CLM(R)*R*R*SQRT(4*PI) FOR L=0,        
!     AND CLM(R)*R*R FOR L>0. IT CAN BE COMPARED WITH THE ORIGINAL      
!     CHARGE DENSITY.                                                   
!     TAKE CARE OF POSSIBLY WANTED COMBINATION OF VLM                   
!                                                                       
!     LAST UPDATE :  20. JAN 87                                         
!     LAST UPDATE :  22. JAN 87 (IMPORT OF FILE OLDCLM)                 
!     LAST UPDATE :  26. JAN 87 (SUBROUTINE MODE)                       
!     LAST UPDATE :   9. MAR 87 (NDIF PARAMETER)                        
!                                                                       
      use densit
      use struct
      use work, only : yka,v
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      DIMENSION        LMMAX(NAT),LM(2,NCOM+3,NAT)                    
      DIMENSION        R(NRAD,NAT)        
      DIMENSION        QPOISS(NRAD)                                     
      DIMENSION        DVDI(NRAD),DER1(NRAD),DER2(NRAD)                 
      DIMENSION        SIGMA(NCOM,NAT),PSVAR(NCOM,NAT),X(2),Y(2)      
!-----------------------------------------------------------------------
      PI=ACOS(-1.)                                                      
!.....GENERATE RADIUS-MESH                                              
      DO 15 JATOM=1,NAT                                                 
         DO 17 JR=1,JRI(JATOM)                                         
            R(JR,JATOM)=R0(JATOM)*EXP( DX(JATOM)*(JR-1) )               
 17      CONTINUE                                                       
 15   CONTINUE                                                          
!                                                                       
      ICOLUM=5                                                          
!      WRITE(9,*) '    ORIGINAL AND CALCULATED CHARGE DENSITY'           
!      WRITE(9,*) '    NORM: CLM*R*R*SQRT(4*PI) FOR L=0, ELSE CLM*R*R'   
!      WRITE(9,1970)                                                     
!      WRITE(9,1970)                                                     
      DO 40 JATOM = 1,NAT                                               
!         WRITE(9,1970)                                                  
!         WRITE(9,2001) LMMAX(JATOM)                                     
         DO 50 LM1= 1,LMMAX(JATOM)                                      
!            WRITE(9,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)               
!...........CALCULATE 1.DERIVATIVE                                      
            DO 60 J=2,JRI(JATOM)-1                                     
               Y(1)=V(J-1,LM1,JATOM)                                    
               Y(2)=V(J+1,LM1,JATOM)                                    
               RAD=R(J,JATOM)                                           
               DVDI(J)=(Y(2)-Y(1))/2.                                   
               DER1(J)=DVDI(J)/RAD/DX(JATOM)                            
 60         CONTINUE                                                    
               Y(1)=V(1,LM1,JATOM)                                      
               Y(2)=V(2,LM1,JATOM)                                      
               RAD=R(1,JATOM)                                           
               DVDI(1)=Y(2)-Y(1)                                        
               DER1(1)=DVDI(1)/RAD/DX(JATOM)                            
               Y(1)=V(JRI(JATOM)-1,LM1,JATOM)                          
               Y(2)=V(JRI(JATOM),LM1,JATOM)                            
               RAD=R(JRI(JATOM),JATOM)                                 
               DVDI(JRI(JATOM))=Y(2)-Y(1)                              
               DER1(JRI(JATOM))=DVDI(JRI(JATOM))/RAD/DX(JATOM)        
!...........CALCULATE 2.DERIVATIVE                                      
            DO 70 J=2,JRI(JATOM)-1                                     
               Y(1)=DVDI(J-1)                                           
               Y(2)=DVDI(J+1)                                           
               RAD=R(J,JATOM)                                           
               D2VDI2=(Y(2)-Y(1))/2.                                    
               DER2(J)=(D2VDI2-DVDI(J)*DX(JATOM))                        &
                       /(RAD*RAD*DX(JATOM)*DX(JATOM))                   
 70         CONTINUE                                                    
               Y(1)=DVDI(1)                                             
               Y(2)=DVDI(2)                                             
               RAD=R(1,JATOM)                                           
               D2VDI2=Y(2)-Y(1)                                         
               DER2(1)=(D2VDI2-DVDI(1)*DX(JATOM))                        &
                       /(RAD*RAD*DX(JATOM)*DX(JATOM))                   
               Y(1)=DVDI(JRI(JATOM)-1)                                 
               Y(2)=DVDI(JRI(JATOM))                                   
               RAD=R(JRI(JATOM),JATOM)                                 
               D2VDI2=Y(2)-Y(1)                                         
               DER2(JRI(JATOM))=(D2VDI2-DVDI(JRI(JATOM))*DX(JATOM))    &
                                 /(RAD*RAD*DX(JATOM)*DX(JATOM))         
!...........===> C A L C U L A T E   CLM'S <===                         
            DO 80 JR=1,JRI(JATOM)                                      
              RAD=R(JR,JATOM)                                           
              D1=DER1(JR)                                               
              D2=DER2(JR)                                               
              VC=V(JR,LM1,JATOM)                                        
              L =ABS( LM(1,LM1,JATOM) )                                 
              RHO=( 2./RAD*D1 + D2 - VC/RAD/RAD*L*(L+1) )/(-8.*PI)      
              QPOISS(JR)=RHO*RAD*RAD                                    
              IF(L.EQ.0) QPOISS(JR)=QPOISS(JR)*SQRT(4.*PI)              
 80         CONTINUE                                                    
!...........CALULATE MEAN SQUARE DEVIATIONS                             
            SIGMA(LM1,JATOM)=0.                                         
            PSVAR(LM1,JATOM)=0.                                         
            DO 85 JR=3,JRI(JATOM)-2                                    
              RAD=R(JR,JATOM)                                           
              DELTA=CLM(JR,LM1,JATOM)-QPOISS(JR)                        
              SIGMA(LM1,JATOM)=SIGMA(LM1,JATOM) + DELTA*DELTA           
              PSVAR(LM1,JATOM)=PSVAR(LM1,JATOM) + ABS( DELTA/            &
                               CLM(JR,LM1,JATOM) )                      
 85         CONTINUE                                                    
            SIGMA(LM1,JATOM)=SIGMA(LM1,JATOM)/(JRI(JATOM)-4)           
            IF(L.EQ.0) SIGMA(LM1,JATOM)=SIGMA(LM1,JATOM)/SQRT(4.*PI)    
            PSVAR(LM1,JATOM)=PSVAR(LM1,JATOM)/(JRI(JATOM)-4)           
            IF(L.EQ.0) PSVAR(LM1,JATOM)=PSVAR(LM1,JATOM)/SQRT(4.*PI)    
!...........WRITE OUT RESULTS                                           
!            JEND=JRI(JATOM)-ICOLUM+1                                   
!            DO 90 JR=1,JEND,ICOLUM                                      
!                WRITE(9,2020) ( CLM(J,LM1,JATOM), J=JR,JR+ICOLUM-1 )    
!                WRITE(9,2020) ( QPOISS(J),        J=JR,JR+ICOLUM-1 )    
!                WRITE(9,1970)                                           
! 90         CONTINUE                                                    
!                JREST=JRI(JATOM) - MOD(JRI(JATOM),ICOLUM) + 1         
!                WRITE(9,2020) (CLM(J,LM1,JATOM), J=JREST,JRI(JATOM))   
!                WRITE(9,2020) (QPOISS(J),        J=JREST,JRI(JATOM))   
!                WRITE(9,2031)                                           
 50      CONTINUE                                                       
!         WRITE(9,2030)                                                  
 40   CONTINUE                                                          
!                                                                       
!.....WRITE OUT  MEAN SQUARE DEVIATIONS                                 
      DO 100 JATOM=1,NAT                                                
!         WRITE(9,1000) JATOM                                            
         DO 110 LM1=1,LMMAX(JATOM)                                      
!            WRITE(9,1010) LM(1,LM1,JATOM),LM(2,LM1,JATOM),              
!     *                    SIGMA(LM1,JATOM),PSVAR(LM1,JATOM)             
 110     CONTINUE                                                       
 100  CONTINUE                                                          
      RETURN                                                            
!                                                                       
!                                                                       
 1000 FORMAT(3X,I2,'. ATOM :   MEAN SQUARE DEV., PSEUDO VAR.COEF.')     
 1010 FORMAT(16X,'CLM(R) FOR L=',I2,3X,'M=',I2,' :   SIGMA=',E12.5,      &
             /,43X,'PSVAR=',E12.5)                                      
 1970 FORMAT(1X,A70)                                                    
 1980 FORMAT(3X)                                                        
 1990 FORMAT(3X,'ATOMNUMBER  =',I2,5X,10A4)                             
 2000 FORMAT(16X,I2//)                                                  
 2001 FORMAT(3X,'NUMBER OF LM=',I2//)                                   
 2010 FORMAT(16X,I2,5X,I2/)                                             
 2011 FORMAT(3X,'VLM(R) FOR L=',I2,3X,'M=',I2/)                         
 2020 FORMAT(3X,5E13.7)                                                 
 2030 FORMAT(///)                                                       
 2031 FORMAT(/)                                                         
 2061 FORMAT(1X,'NUMBER OF PW',I6)                                     
 2062 FORMAT(16X,I3)                                                    
 2070 FORMAT(3X,3I5,2E15.7)                                             
!                                                                       
      END                                                               
