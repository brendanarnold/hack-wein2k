      SUBROUTINE POTFAC(VINF,VOUTF,L,M,IATNR)                           
!                                                                       
!-------------------------------------------------------------          
!---  POTFAC CALCULATES FACTORS                             --          
!---                                                        --          
!-------------------------------------------------------------          
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8   ::  c_kub(0:10,0:10)
      COMMON /norm_kub/ c_kub

      pi  =ACOS(-1.0D0)                                                 
!                                                                       
      VINF =1.D0                                                        
      IF(L.EQ.0) VINF =1./SQRT(4.D0*PI)                               
      VOUTF=1.D0                                                        
      IF(IATNR.LT.0) RETURN                                             
!                                                                       
         IF(l.EQ.0.AND.m.EQ.0) THEN
            voutf=1.d0
         ELSEIF (l.EQ.-3.AND.m.EQ.2) THEN  
            voutf=1.d0
         ELSE
            voutf=.5d0/(c_kub(abs(l),m)*c_kub(abs(l),m))
         ENDIF
!                                                                       
      RETURN                                                            
      END                                                               
