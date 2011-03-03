      subroutine strwri(spos,mult,isplit,name,rotloc,iatnr,gt)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension spos(3,*),rotloc(3,3)
      character*79 name
      character*2 gt
 901  format(' ****** IATNR IN STRUCT_ST CHANGED TO A', &
       ' NEGATIVE NUMBER ******')
 902  format(' ****** IATNR IN STRUCT_ST CHANGED TO A', &
       ' POSITIVE NUMBER ******')
!
      if(gt.eq.'gt') then      
        if(iatnr.gt.0) then
              write(6,901)
              iatnr=-iatnr
        endif
      else       
        if(iatnr.lt.0) then
              write(6,902)
              iatnr=-iatnr
        endif
      endif
!
      write(21,1010) IATNR,( SPOS(J,1),J=1,3 )
      write(21,1012) mult,isplit
      DO 55 M=2,MULT                                     
 55   write(21,1011) IATNR,( SPOS(J,M),J=1,3)
      write(21,1510) NAME 
      write(21,1013) ((rotloc(i,j),j=1,3),i=1,3)          
 1013 format('LOCAL ROT MATRIX:   ',3f10.7,/,20x,3f10.7,/,20x,3f10.7)
 1510 FORMAT(A79)                                                       
 1010 FORMAT('ATOM',I4,': X=',F10.8,' Y=',F10.8,' Z=',F10.8)            
 1011 FORMAT(4X,I4,': X=',F10.8,' Y=',F10.8,' Z=',F10.8)                
 1012 FORMAT('          MULT=',i2,'          ISPLIT=',i2)
!
      end
