      FUNCTION SININT(X)                                                
!CC SINE INTEGRAL FROM ZERO TO ARGUMENT X                               
      IMPLICIT REAL*8 (A-H,O-Z)
      IF (ABS(X).GE.16.d0) GOTO 1                                         
      Z = X/16.d0                                                         
      Y = 4.d0*Z**2-2.d0                                                    
      B =      +0.000000000000007d0                                     
      A = Y*B  -0.000000000000185d0                                     
      B = Y*A-B+0.000000000004185d0                                     
      A = Y*B-A-0.000000000084710d0                                     
      B = Y*A-B+0.000000001522370d0                                     
      A = Y*B-A-0.000000024100076d0                                     
      B = Y*A-B+0.000000332988589d0                                     
      A = Y*B-A-0.000003972908746d0                                     
      B = Y*A-B+0.000040420271419d0                                     
      A = Y*B-A-0.000345469155569d0                                     
      B = Y*A-B+0.002436221404749d0                                     
      A = Y*B-A-0.013867445589417d0                                     
      B = Y*A-B+0.062033679432003d0                                     
      A = Y*B-A-0.211263780976555d0                                     
      B = Y*A-B+0.53014884791652d0                                      
      A = Y*B-A-0.96832223698708d0                                      
      B = Y*A-B+1.38930877117188d0                                      
      A = Y*B-A-1.92656509115065d0                                      
      B = Y*A-B+2.77875638174266d0                                      
      A = Y*B-A-4.0639808449119d0                                       
      A = Y*A-B+8.1058529553612d0                                       
      SININT = Z*.5d0*(A-B)                                               
      RETURN                                                            
!                                                                       
    1 Z = 16.d0/X                                                         
      G = Z**2                                                          
      Y = 4.d0*G-2.d0                                                       
      B =      +0.000000000000002d0                                     
      A = Y*B  -0.000000000000014d0                                     
      B = Y*A-B+0.000000000000107d0                                     
      A = Y*B-A-0.000000000000964d0                                     
      B = Y*A-B+0.000000000010308d0                                     
      A = Y*B-A-0.000000000136096d0                                     
      B = Y*A-B+0.000000002356196d0                                     
      A = Y*B-A-0.000000058670317d0                                     
      B = Y*A-B+0.000002453755677d0                                     
      A = Y*B-A-0.000233756041393d0                                     
      A = Y*A-B+0.124527458057854d0                                     
      F = Z*.5d0*(A-B)                                                    
      B =      +0.000000000000002d0                                     
      A = Y*B  -0.000000000000012d0                                     
      B = Y*A-B+0.000000000000087d0                                     
      A = Y*B-A-0.000000000000717d0                                     
      B = Y*A-B+0.000000000006875d0                                     
      A = Y*B-A-0.000000000079604d0                                     
      B = Y*A-B+0.000000001169202d0                                     
      A = Y*B-A-0.000000023468225d0                                     
      B = Y*A-B+0.000000724995950d0                                     
      A = Y*B-A-0.000042644182622d0                                     
      A = Y*A-B+0.007725712193407d0                                     
      G = G*.5d0*(A-B)                                                    
      B = 1.57079632679489d0                                            
      IF (X.LT.0.d0) B = -B                                               
      SININT = B-F*COS(X)-G*SIN(X)                                      
      RETURN                                                            
      END                                                               
