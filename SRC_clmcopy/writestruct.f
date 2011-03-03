SUBROUTINE writestruct
  ! writes struct file stored in module struct
  ! 23/2-01 PB
  USE struct
  USE symetr
  USE reallocate

  REAL*8              :: test,ninety
  TEST=1.D-5
  NINETY=90.0D0
  
  write(21,1000) title                                               
  write(21,1010) lattic,nat,cform,irel    
  write(21,1020) aa,bb,cc,alpha(1),alpha(2),alpha(3)
!
  index=0                                                           
  DO jatom = 1,nat                                               
     index=index+1
     !
     write(21,1030) iatnr(jatom),( pos(j,index),j=1,3 ), &
          mult(jatom),isplit(jatom) 
     DO m=1,mult(jatom)-1                                     
        index=index+1                                            
        write(21,1031) iatnr(jatom),( pos(j,index),j=1,3)         
     ENDDO
     write(21,1050) aname(jatom),jri(jatom),r0(jatom),rmt(jatom), &
          zz(jatom)
     write(21,1051) ((rotloc(i,j,jatom),i=1,3),j=1,3)                
  ENDDO
  write(21,1151) iord
  DO j=1,iord
     write(21,1101) ( (iz(j1,j2,j),j1=1,3),tau(j2,j),j2=1,3 ),inum(j)
  ENDDO
  RETURN
1000 FORMAT(A80)                                                       
1010 FORMAT(A4,23X,I3,1x,a4,/,13X,A4,18X,A4)                                 
1020 FORMAT(6F10.6)                                          
1030 FORMAT('ATOM',I4,': X=',F10.8,' Y=',F10.8,' Z=',F10.8,/,10X,'MULT=',I2,10X,'ISPLIT=',I2)          
1031 FORMAT(4X,I4,': X=',F10.8,' Y=',F10.8,' Z=',F10.8)                          
1050 FORMAT(A10,' NPT=',I5,'  R0=',F10.8,' RMT=',F10.4,'   Z:',F10.5)                      
1051 FORMAT(20X,3F10.7)                                             
1151 FORMAT(I4,'      NUMBER OF SYMMETRY OPERATIONS')
1101 FORMAT(3(3I2,F11.8/),I8)
  
END SUBROUTINE writestruct
