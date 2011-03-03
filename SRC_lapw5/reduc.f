      SUBROUTINE REDUC(V,ATP,NPOS,POS,IAT,rmt)                              
!                                                                       
!     REDUCES VECTOR V TO EQUIVALENT SMALLEST VECTOR                    
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V(3),ATP(3,*),VHELP(3),POS(3,*),VTEST(3)                
!        write(6,*) iat,v,pos(1,iat),pos(2,iat),pos(3,iat)
      DO  9 J=1,3                                                       
   9  VHELP(J)=V(J)                                                     
!      R=VNORM(V)                                                        
      R=99999.                                                        
      DO 10 I=1,NPOS                                                    
      DO 11 J=1,3                                                       
  11  VTEST(J)=V(J)-ATP(J,I)-POS(J,IAT)                                 
      RTEST=VNORM(VTEST)                                                
!         write(6,*) i,atp(1,i),atp(2,i),atp(3,i),rtest,r
      IF(RTEST.LT.R) THEN                                               
      R=RTEST                                                           
      DO 12 J=1,3                                                       
  12  VHELP(J)=VTEST(J)                                                 
      END IF                                                            
  10  CONTINUE                                                          
      DO 13 J=1,3                                                       
  13  V(J)=VHELP(J)                          
      if(r.gt.rmt) write(6,*) ' Reduction failed (wrong symmetry oper) r-reduc gt.rmt',r,rmt                           
      RETURN                                                            
      END                                                               
